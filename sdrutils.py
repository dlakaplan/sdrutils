import pmt
from gnuradio.blocks import parse_file_metadata
import numpy as np,scipy,scipy.signal
from astropy.io import fits
import sys,os
from rtlsdr import RtlSdr  
from astropy import units as u, constants as c


from gnuradio import analog
from gnuradio import audio
from gnuradio import filter
from gnuradio import gr
from gnuradio import gr, blocks
from gnuradio import uhd
from gnuradio.fft import window
from gnuradio.filter import firdes
from optparse import OptionParser
import time
import datetime

##################################################
def gr2fits(filename, merge=True, verbose=False):

    try:
        handle=open(filename, 'rb')
    except:
        raise IOError('File %s does not exist' % filename)

    nheaders = 0
    nread = 0
    baseband=[]
    
    fitsout=fits.HDUList()
    
    while(True):
        """
        /opt/local/bin/gr_read_file_metadata
        note that there can be > 1 metadata blocks
        I think they can come every 1e6 items
        """
        # read out next header bytes
        hdr_start = handle.tell()
        header_str = handle.read(parse_file_metadata.HEADER_LENGTH)
        if(len(header_str) == 0):
            break

        # Convert from string to PMT (should be a dictionary)
        try:
            header = pmt.deserialize_str(header_str)
        except RuntimeError:
            raise IOError("Could not deserialize header: invalid or corrupt data file.\n")

        if verbose:
            print("HEADER {0}".format(nheaders))
        info = parse_file_metadata.parse_header(header, verbose)
        if (info["extra_len"] > 0):
            extra_str = handle.read(info["extra_len"])
            if(len(extra_str) == 0):
                break

            try:
                extra = pmt.deserialize_str(extra_str)
            except RuntimeError:
                sys.stderr.write("Could not deserialize extras: invalid or corrupt data file.\n")
                break

            if verbose:
                print("\nExtra Header:")
            extra_info = parse_file_metadata.parse_extra_dict(extra, info, verbose)

        nheaders += 1
        nread += parse_file_metadata.HEADER_LENGTH + info["extra_len"]
        handle.seek(nread, 0)
        h=extra_info
        if h['size']==8 and h['cplx']:
            dtype=scipy.complex64

        d=scipy.fromfile(handle, dtype=dtype, count=h['nitems'])
        t0=np.arange(2*len(d))/h['rx_rate']/2
        t=np.arange(len(d))/h['rx_rate']
        
        nread += info['nbytes']
        handle.seek(nread, 0)
        
        fitsout.append(fits.ImageHDU(data=np.c_[d.real,d.imag]))
        fitsout[-1].header['NITEMS']=(h['nitems'],'Number of complex samples')
        fitsout[-1].header['RATE']=(h['rx_rate'],'[Hz] sample rate')
        fitsout[-1].header['RX_FREQ']=(pmt.to_float(h['rx_freq'])/1e6,'[MHz] Radio frequency')
        fitsout[-1].header['RX_TIME']=(h['rx_time'],'[s] Time of start of block')
        

    if merge:
        totallen=0
        for i in xrange(0,len(fitsout)):
            totallen+=fitsout[i].header['NAXIS2']
        d=np.zeros((totallen,2),dtype=fitsout[1].data.dtype)
        nmax=0
        for i in xrange(0,len(fitsout)):
            d[nmax:nmax+fitsout[i].header['NAXIS2']]=fitsout[i].data
            nmax+=fitsout[i].header['NAXIS2']
        newfitsout=fits.HDUList()
        newfitsout.append(fits.PrimaryHDU(data=d))
        newfitsout[0].header=fitsout[1].header
        newfitsout[0].header['NITEMS']=totallen
        newfitsout[0].header['EXPTIME']=(d.shape[0]/newfitsout[0].header['RATE'],'[s] Duration of file')
        fitsout=newfitsout

    fitsout.verify('silentfix')
    if os.path.exists(filename + '.fits'):
        os.remove(filename + '.fits')
    fitsout.writeto(filename + '.fits')
    print('Wrote %s.fits' % filename)
    return fitsout

##################################################
class sdrdata():

    def __init__(self, file=None):

        self.rate=None
        self.file=None
        self.center_freq=None
        self.exptime=None
        self.data=None
        self.nitems=0

        if file is not None and file.endswith('.fits'):
            self.from_fits(file)

    ##############################
    def from_fits(self, fitsfile):
        self.file=fitsfile
        self.f=fits.open(self.file)
        # read in the data, convert to complex
        self.data=self.f[0].data[:,0]+1j*self.f[0].data[:,1]
        self.dtype=self.data.dtype
        self.exptime=self.f[0].header['EXPTIME']*u.s

        # this is the complex sampling rate
        self.rate=self.f[0].header['RATE']*u.Hz
        self.center_freq=self.f[0].header['RX_FREQ']*u.MHz
        self.nitems=len(self.data)
        
    ##############################
    def from_rtl(self, freq=93.3e6, rate=256000, samples=8192000):
        try:
            sdr=RtlSdr()

            # configure device
            # in Hz
            sdr.sample_rate = int(rate)
            sdr.center_freq = freq
            sdr.gain = 'auto'
    
            # Read samples
            samples = sdr.read_samples(samples)
            # Clean up the SDR device
            sdr.close()  
            del(sdr)

            self.data=np.array(samples).astype('complex64')
            self.dtype=self.data.dtype
            self.rate=rate*u.Hz
            self.center_freq=freq*u.Hz
            self.nitems=len(self.data)
            self.exptime=(self.nitems/self.rate)
        except:
            raise IOError("Could not read from RTLSDR")
        
    ##############################
    def tofm(self, outfile=None, offset=0):
        """
        https://witestlab.poly.edu/blog/capture-and-decode-fm-radio/
        """

        Fs=self.rate.value
        
        # the data could have been taken with a frequency offset
        # to avoid the DC spike
        # so mix down to shift that away
        
        # To mix the data down, generate a digital complex exponential 
        # (with the same length as x1) with phase -F_offset/Fs
        fc1 = np.exp(-1.0j*2.0*np.pi* ((offset/self.rate).decompose().value*np.arange(self.nitems)))
        # Now, just multiply x1 and the digital complex expontential
        self.x2 = self.data * fc1

        
        # An FM broadcast signal has  a bandwidth of 200 kHz
        f_bw = 200000  
        try:
            n_taps = 64  
            # Use Remez algorithm to design filter coefficients
            lpf = signal.remez(n_taps, [0, f_bw, f_bw+(Fs/2-f_bw)/4, Fs/2], [1,0], Hz=Fs)  
            self.x3 = signal.lfilter(lpf, 1.0, self.x2)

            dec_rate = int(Fs / f_bw)  
            self.x4 = self.x3[0::dec_rate]  
            # Calculate the new sampling rate
            self.Fs_y = Fs/dec_rate  
        except:
            self.x4=self.x2
            self.Fs_y=Fs

        ### Polar discriminator
        y5 = self.x4[1:] * np.conj(self.x4[:-1])  
        self.x5 = np.angle(y5)  
            
        # The de-emphasis filter
        # Given a signal 'x5' (in a numpy array) with sampling rate Fs_y
        d = self.Fs_y * 75e-6   # Calculate the # of samples to hit the -3dB point  
        x = np.exp(-1/d)   # Calculate the decay between each sample  
        b = [1-x]          # Create the filter coefficients  
        a = [1,-x]  
        self.x6 = scipy.signal.lfilter(b,a,self.x5)  

        
        # Find a decimation rate to achieve audio sampling rate between 44-48 kHz
        audio_freq = 44100.0  
        dec_audio = int(self.Fs_y/audio_freq)  
        self.Fs_audio = self.Fs_y / dec_audio
        
        self.x7 = scipy.signal.decimate(self.x6, dec_audio)
        # Scale audio to adjust volume
        self.x7 *= 10000 / np.max(np.abs(self.x7))  
        if outfile is not None:
            scipy.io.wavfile.write(outfile, self.Fs_audio, self.x7)
            print('Wrote %s' % outfile)
        return (self.x7,self.Fs_audio*u.Hz)
    
##################################################
class FX:

    def __init__(self, file1, file2,
                 df=0.1*u.kHz,
                 dt=0.1*u.s,
                 window=None,
                 shift=0):
        """
        f=FX(file1, file2,
        df=0.1*u.kHz,
        dt=0.1*u.s,
        window=None,
        shift=0)

        window should be a function that accepts the number of points in the FFT
        such as
        def w(n):
           return np.kaiser(n, 5)


        shift is in seconds
        if there is a phase slope of S radians/channel:

        shift=S/f.channel/2/np.pi
        
        """

        self.file1=file1
        self.file2=file2
        self.df=df
        self.dt=dt
        
        self.s1=sdrdata(self.file1)
        self.s2=sdrdata(self.file2)

        if self.s1.rate != self.s2.rate:
            raise ValueError('Sampling rate for Rx2 (%f Hz) does not match rate for Rx1 (%f Hz)' % (self.s1.rate.value,self.s2.rate.value))
        if self.s1.center_freq != self.s2.center_freq:
            raise ValueError('Frequency for Rx2 (%f MHz) does not match rate for Rx1 (%f MHz)' % (self.s1.center_freq.to(u.MHz).value,
                                                                                       self.s2.center_freq.to(u.MHz).value))
        

        # read in the data, convert to complex
        self.d1=self.s1.data
        self.d2=self.s2.data

        self.dtype=self.d1.dtype

        self.exptime=self.s1.exptime
        self.center_freq=self.s1.center_freq

        # this is the complex sampling rate
        self.rate=self.s1.rate
        
        # needs this long a fft
        self.nfft=int(round((self.rate/self.df).decompose()))
        # there are this many ffts
        self.nsamples=self.d1.shape[0]/self.nfft
        # so each chunk is this long
        self.chunk=(self.nfft/self.rate).decompose()
        # and we need to add this many together
        self.nint=int(round(self.dt/self.chunk))
        # and we have this many correlations
        self.ncorr=int((self.exptime/self.dt))

        # frequency of FFTs
        self.freq=np.fft.fftshift(np.fft.fftfreq(self.nfft,(1/self.rate).decompose().value))*u.Hz
        # and the actual frequency on the sky
        self.rf_freq=self.center_freq+self.freq

        # Hz per channel
        self.channel=np.diff(self.freq).mean()

        if window is None:
            self.window=np.kaiser(self.nfft,5).astype('complex')
        else:
            self.window=window(self.nfft).astype('complex')

        # apply a phase shift
        # to the second input
        phi=np.exp(2j*np.pi*(self.rf_freq*shift).decompose().value)
 
        # dynamic spectra
        self.DS1=np.zeros((self.nsamples,self.nfft),dtype=self.dtype)
        self.DS2=np.zeros((self.nsamples,self.nfft),dtype=self.dtype)

        # output correlation
        self.outcorr=np.zeros((self.ncorr,self.nfft),dtype=self.dtype)

        j=0
        for i in range(self.ncorr):
            corr=np.zeros((self.nfft),dtype=self.dtype)
            for k in range(self.nint):
                D1=self.d1[j*self.nfft:(j+1)*self.nfft]*self.window
                D2=self.d2[j*self.nfft:(j+1)*self.nfft]*self.window
        
                F1=np.fft.fftshift(np.fft.fft(D1))
                F2=np.fft.fftshift(np.fft.fft(D2))
                F2*=phi
                self.DS1[j]=F1
                self.DS2[j]=F2
        
                corr+=F1*np.conj(F2)
                j+=1

            self.outcorr[i]=corr
    
        self.phase=np.angle(self.outcorr)*u.rad
        self.amp=np.absolute(self.outcorr)

##################################################
class dual_recorder(gr.top_block):

    def __init__(self, samp_rate=256e3,
                 gain1=28,
                 gain2=28,
                 freq=89.7e6,
                 samples=1e6,
                 out=None,
                 verbose=True):
        gr.top_block.__init__(self)
        
        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate 
        self.gain2 = gain2 
        self.gain1 = gain1 
        self.freq = freq
        self.out = out
        self.verbose = verbose
        if samples is not None:
            self.samples=samples
        else:
            self.samples=1e6
        if self.verbose:
            print('Will record %d samples at %.2f MHz, sampling at %.1f kHz' % (self.samples,
                                                                                self.freq/1e6,
                                                                                self.samp_rate/1e3))

        ##################################################
        self.uhd_usrp_source_0 = uhd.usrp_source(
            ",".join(("master_clock_rate=30.72e6", "")),
            uhd.stream_args(
            cpu_format="fc32",
            channels=range(2),
            ),
            )
        
        self.uhd_usrp_source_0.set_samp_rate(self.samp_rate)
        self.uhd_usrp_source_0.set_time_now(uhd.time_spec(time.time()), uhd.ALL_MBOARDS)

        self.uhd_usrp_source_0.set_center_freq(self.freq, 0)
        self.uhd_usrp_source_0.set_gain(self.gain1, 0)
        self.uhd_usrp_source_0.set_antenna("TX/RX", 0)
        self.uhd_usrp_source_0.set_auto_dc_offset(True, 1)
        self.uhd_usrp_source_0.set_auto_iq_balance(True, 1)  
        self.uhd_usrp_source_0.set_center_freq(self.freq, 1)
        self.uhd_usrp_source_0.set_gain(self.gain2, 1)
        self.uhd_usrp_source_0.set_antenna("RX2", 1)
        self.uhd_usrp_source_0.set_auto_dc_offset(True, 1)
        self.uhd_usrp_source_0.set_auto_iq_balance(True, 1)
  


        self.filenamebase=datetime.datetime.now().strftime('%Y%m%dT%H%M%S')
        if self.out is not None:
            self.filenamebase=os.path.join(self.out,self.filenamebase)
        self.filename_1='%s_1.dat' % self.filenamebase
        self.filename_2='%s_2.dat' % self.filenamebase
        if self.verbose:
            print('Will write output to %s, %s' % (self.filename_1, self.filename_2))
        
        self.blocks_file_meta_sink_1 = blocks.file_meta_sink(gr.sizeof_gr_complex*1,
                                                             self.filename_1,
                                                             self.samp_rate, 1,
                                                             blocks.GR_FILE_FLOAT,
                                                             True, 1000000, "", False)

        self.blocks_file_meta_sink_1.set_unbuffered(False)
        self.blocks_file_meta_sink_2 = blocks.file_meta_sink(gr.sizeof_gr_complex*1,
                                                             self.filename_2,
                                                             self.samp_rate, 1,
                                                             blocks.GR_FILE_FLOAT,
                                                             True, 1000000, "", False)
        self.blocks_file_meta_sink_2.set_unbuffered(False)

        ##################################################
        # Connections
        ##################################################
        self._head1=blocks.head(gr.sizeof_gr_complex,
                                int(self.samples))
        self._head2=blocks.head(gr.sizeof_gr_complex,
                                int(self.samples))
        
        self.connect((self.uhd_usrp_source_0, 0), self._head1,
                     (self.blocks_file_meta_sink_1, 0))    
        self.connect((self.uhd_usrp_source_0, 1), self._head2,
                     (self.blocks_file_meta_sink_2, 0))    
##################################################
class single_recorder(gr.top_block):

    def __init__(self, samp_rate=256e3,
                 gain=28,
                 freq=89.7e6,
                 samples=1e6,
                 out=None,
                 verbose=True):
        gr.top_block.__init__(self)
        
        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate 
        self.gain = gain
        self.freq = freq
        self.out = out
        self.verbose = verbose
        if samples is not None:
            self.samples=samples
        else:
            self.samples=1e6
        if self.verbose:
            print('Will record %d samples at %.2f MHz, sampling at %.1f kHz' % (self.samples,
                                                                                self.freq/1e6,
                                                                                self.samp_rate/1e3))

        ##################################################
        self.uhd_usrp_source_0 = uhd.usrp_source(
            ",".join(("master_clock_rate=30.72e6", "")),
            uhd.stream_args(
            cpu_format="fc32",
            channels=range(1),
            ),
            )
        
        self.uhd_usrp_source_0.set_samp_rate(self.samp_rate)
        self.uhd_usrp_source_0.set_time_now(uhd.time_spec(time.time()), uhd.ALL_MBOARDS)

        self.uhd_usrp_source_0.set_center_freq(self.freq, 0)
        self.uhd_usrp_source_0.set_gain(self.gain, 0)
        self.uhd_usrp_source_0.set_antenna("TX/RX", 0)
        self.uhd_usrp_source_0.set_auto_dc_offset(True, 1)
        self.uhd_usrp_source_0.set_auto_iq_balance(True, 1)  
        #self.uhd_usrp_source_0.set_center_freq(self.freq, 1)
        #self.uhd_usrp_source_0.set_gain(self.gain2, 1)
        #self.uhd_usrp_source_0.set_antenna("RX2", 1)
        #self.uhd_usrp_source_0.set_auto_dc_offset(True, 1)
        #self.uhd_usrp_source_0.set_auto_iq_balance(True, 1)
  


        self.filenamebase=datetime.datetime.now().strftime('%Y%m%dT%H%M%S')
        if self.out is not None:
            self.filenamebase=os.path.join(self.out,self.filenamebase)
        self.filename_1='%s.dat' % self.filenamebase
        #self.filename_2='%s_2.dat' % self.filenamebase
        if self.verbose:
            print('Will write output to %s' % (self.filename_1))
        
        self.blocks_file_meta_sink_1 = blocks.file_meta_sink(gr.sizeof_gr_complex*1,
                                                             self.filename_1,
                                                             self.samp_rate, 1,
                                                             blocks.GR_FILE_FLOAT,
                                                             True, 1000000, "", False)

        self.blocks_file_meta_sink_1.set_unbuffered(False)
        #self.blocks_file_meta_sink_2 = blocks.file_meta_sink(gr.sizeof_gr_complex*1,
        #                                                     self.filename_2,
        #                                                     self.samp_rate, 1,
        #                                                     blocks.GR_FILE_FLOAT,
        #                                                     True, 1000000, "", False)
        #self.blocks_file_meta_sink_2.set_unbuffered(False)

        ##################################################
        # Connections
        ##################################################
        self._head1=blocks.head(gr.sizeof_gr_complex,
                                int(self.samples))
        #self._head2=blocks.head(gr.sizeof_gr_complex,
        #                        int(self.samples))
        
        self.connect((self.uhd_usrp_source_0, 0), self._head1,
                     (self.blocks_file_meta_sink_1, 0))    
        #self.connect((self.uhd_usrp_source_0, 1), self._head2,
        #             (self.blocks_file_meta_sink_2, 0))    
