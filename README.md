# JWST-tools

These notebooks provide help for searching and analyzing JWST data.

## GetInfoPrograms
Create a csv file with relevant information on all observed MIRI programs -or a part- and download data from MAST.

## StoreMAST
Store JWST data into folders based on a logical architecture.

## GetJWSTTelemetry
Download and plot JWST telemetry data from a FITS file or for a specific observation date using mnemonics.

## DitherPattern
Plot the MIRI dither patterns.

## JWSTExplorers
Interactive JWST/MIRI explorers (demo).
- ExploreCube.py: MIRI/MRS cube explorer. Show the spectrum for each spaxel or the avarage pixel in a specific region.

How to run: python ExploreCube.py mrs_s3d_file
- ExploreFlags.py: Show the FITS image and translate the DQ value for a specific pixel.

How to run: python ExploreFlags.py jwst_2d_file
- ExploreRamps.py: Plot the ramps (groups and integrations) of the uncalibrated (if _uncal data available) and calibrated (if _ramp or _jump data available), show the flags for each group and each pixel and the fitted line (if _fitopt or _rateints data available). Files should be in the corresponding /stage0 and /stage1 folders.

How to run: python path/to/ExploreRamps.py path/to/fitsfile.fits

## MRSWCSCorrection
Improve JWST MIRI/MRS WCS using source detection on simultaneous imaging.

Created by: B. Trahin