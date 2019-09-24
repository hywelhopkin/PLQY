# PLQY
Routine for PLQY calculations

This file is based on the previous Matlab routine used to analyse the PLQY in the G48 setup.

The measurement procedure is based on the paper by De Mello et al. (https://doi.org/10.1002/adma.19970090308)

## Installation
First install a python distribution, for example Anaconda (download here: https://www.anaconda.com/distribution/)

Copy the main file and the calibration files into your working folder

## Usage
Usage of the program is described running `PLQE_bw.py --help`:

```
usage: PLQE_bw.py [-h] -sn SHORT_NAME [-d DIRECTORY] [-cfg CONFIG]
                  [-lr LASER_RANGE LASER_RANGE] [-plr PL_RANGE PL_RANGE]
                  [-st SHORT_TIME] [-ln LONG_NAME] [-lt LONG_TIME] [-c COMMON]
                  [-cb COMMON_BCKG] [-ce COMMON_EMPTY] [-clb COMMON_LONG_BCKG]
                  [-cle COMMON_LONG_EMPTY]

Calculation of PL quantum efficiency.

optional arguments:
  -h, --help            show this help message and exit
  -sn SHORT_NAME, --short_name SHORT_NAME
                        Prefix of the file names (e.g 'filename_') to be
                        followed by 'in.txt', 'out.txt', 'bckg.txt' or
                        'empty.txt'

optional arguments:
  -d DIRECTORY, --directory DIRECTORY
                        folder in which the raw files are stored (e.g.
                        'folder_name/sub_folder/'
  -cfg CONFIG, --config CONFIG
                        Fiber and spectrometer configurations. Choices: 'Maya
                        red', 'Maya steel', 'QE red 25', 'QE red 200', 'QE
                        steel 25' or 'QE steel 200', Default = 'Maya red'
  -lr LASER_RANGE LASER_RANGE, --laser_range LASER_RANGE LASER_RANGE
                        Laser wavelength range in nm, Default = [397, 407]
  -plr PL_RANGE PL_RANGE, --pl_range PL_RANGE PL_RANGE
                        PL detection range in nm, Default = [550, 850]
  -st SHORT_TIME, --short_time SHORT_TIME
                        Integration time for short measurement in ms
  -ln LONG_NAME, --long_name LONG_NAME
                        Prefix of the long file name (e.g 'filename_long_') to
                        be followed by 'in.txt', 'out_.txt', 'bckg.txt' or
                        'empty.txt'
  -lt LONG_TIME, --long_time LONG_TIME
                        Integration time for long measurement in ms
  -c COMMON, --common COMMON
                        Indicates that common background and empty files are
                        used. Default = True
  -cb COMMON_BCKG, --common_bckg COMMON_BCKG
                        Name of the common background file. Default =
                        'bckg.txt'
  -ce COMMON_EMPTY, --common_empty COMMON_EMPTY
                        Name of the common empty file. Default = 'empty.txt'
  -clb COMMON_LONG_BCKG, --common_long_bckg COMMON_LONG_BCKG
                        Name of the common long background file. Default =
                        'long_bckg.txt'
  -cle COMMON_LONG_EMPTY, --common_long_empty COMMON_LONG_EMPTY
                        Name of the common long empty file. Default =
                        'long_empty.txt'
  ```
  
  ## Example
  Run: 
  `python PLQE_bw.py  -d ../../Desktop/PL/pl2/ -st 15 -lt 100 -plr 700 850 -sn S10_sh_ -ln S10_lg_`
