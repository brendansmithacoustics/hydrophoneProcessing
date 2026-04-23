This repository contains MATLAB code to compute absolute sound level measurements from Ocean Networks Canada hydrophone data.

There are two .m files:
- calibrateONChydrophoneSpectra.m is a function which computes calibrated spectra when supplied with audio data (in units of counts) and calibration data. This returns a SpectData structure.
- processHydrophoneAudio.m is a master sript which points to example audio and sensitivity data, calls calibrateONChydrophoneSpectra.m, and plots the result.

The example data is from the Main Endeavour Field.