#!/usr/bin/env python2

from gnuradio import eng_notation
from gnuradio.eng_option import eng_option
from optparse import OptionParser
import sys,os
import sdrutils



        
##################################################
if __name__ == '__main__':

    parser=OptionParser(option_class=eng_option, usage="%prog: [options]")
    parser.add_option("-f", "--freq", type="eng_float", default=89.7e6,
                      help="Center frequency [default=%default]")
    parser.add_option("--gain", type="eng_float", default=28,
                      help="Channel gain in dB [default=%default]")
    #parser.add_option("--gain2", type="eng_float", default=28,
    #                  help="Second channel gain in dB [default=%default]")
    parser.add_option("-s", "--sampling", type="eng_float", default=256e3,
                      help="Sampling rate [default=%default]")
    parser.add_option("-n", "--nsamples", type="eng_float", default=None,
                      help="Number of samples to collect [default=+inf]")
    parser.add_option("-t", "--time", type="float", default=None,
                      help="Time to collect samples [default=+inf]")
    parser.add_option("-o", "--out", type="str", default='./',
                      help="Output directory [default=curdir]")
    parser.add_option('--noconvert',action='store_false',default=True,
                      dest='convert',
                      help='Do not convert raw output to FITS [default=convert]')
    parser.add_option('-c','--clean',
                      action="store_true",dest="clean",default=False,
                      help="Remove temporary files?")
    #parser.add_option('-m','--merge',
    #                  action="store_true",dest="merge",default=False,
    #                  help="Merge separate output into single file?")
    parser.add_option('-v','--verbose',
                      action="store_true",dest="verbose",default=False,
                      help="Increase verbosity of output")

    (options, args) = parser.parse_args()

    if options.time is not None and options.nsamples is None:
        options.nsamples=int(options.time*options.sampling)

    elif options.nsamples is not None:
        options.time=options.nsamples/options.sampling

    if options.nsamples is None:
        if options.verbose:
            print('Will record indefinitely...')
    else:
        if options.verbose:
            print('Will record %d samples (%f s)' % (options.nsamples,options.time))

    rec = sdrutils.single_recorder(samp_rate=options.sampling,
                                   freq=options.freq,
                                   gain=options.gain,
                                   #gain2=options.gain2,
                                   samples=options.nsamples,
                                   out=options.out,
                                   verbose=options.verbose)
    try:
        rec.run()
    except KeyboardInterrupt:
        pass
    if options.convert:
        filename_1=rec.filename_1
        filenamebase=rec.filenamebase
        # make sure the buffered output is flushed
        del rec
        out1=sdrutils.gr2fits(filename_1,verbose=options.verbose)
        #out2=sdrutils.gr2fits(filename_2,verbose=options.verbose)
        if options.clean:
            if options.verbose:
                print('Deleting %s...' % (filename_1))
            os.remove(filename_1)
            #os.remove(filename_2)
        if options.verbose:
            print('Found %d samples in %s' % (out1[0].header['NITEMS'],
                                              filename_1))
            #print('Found %d samples in %s' % (out2[0].header['NITEMS'],
            #                                  filename_2))

    sys.exit(0)
    
