# sdrutils
Utilities for SDR

contents of `sdrutils`:
* `gr2fits`: convert gnuradio binary output to FITS
* `sdrdata`: store/process SDR data, from FITS or direct from RTLSDR
* `FX`: FX correlator
* `dual_recorder`: use 2 channels on USRP to record
* `single_recorder`: use 1 channel on USRP to record


executables:
* `gr2fits.py`
* `record1.py`: `python ./record1.py -v -t 5 -f 93.3M -c -o ../`
* `record2.py`: `python ./record1.py -v -t 5 -f 93.3M -c -o ../`
