#!/usr/bin/env pythonw
from __future__ import print_function, division
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import argparse
from astropy import log
from os import path
from glob import glob
from subprocess import check_call
import shutil
from astropy.table import Table
from astropy.io import fits

from nicer.values import *
from nicer.plotutils import plot_light_curve

############################################################################################
## Code based on method by S. Bogdanov
## Code still in development. 
##
## TO DO:
##  - fix the os.system() or check_call()
##  - add consistency checks:
##     1) cut is within (min,max) of RATE values in light curve
##     2) if --lcfile and -- timebin are given, check that time bin is consistent
##  - give total time after filtering
##  - consider making the FITS file manipulation (ftcalc, ftcreate, etc.) directly in python
##  - consider implementing call from psrpipe.py
############################################################################################

def runcmd(cmd):
    # CMD should be a list of strings since it is not processed by a shell
    log.info('CMD: '+" ".join(cmd))
    os.system(" ".join(cmd))
    ## Some ftools calls don't work properly with check_call...not sure why!
    ## so I am using os.system instead of check_call
    #check_call(cmd,env=os.environ)

def getgti(evf):
    # Read the GTIs from the  event FITS file
    gtitable = Table.read(evf,hdu=2)
    return gtitable


################################################
# Checking the presence of HEASOFT
try:
    check_call('nicerversion',env=os.environ)
except:
    print("You need to initialize FTOOLS/HEASOFT first (e.g., type 'heainit')!", file=sys.stderr)
    exit()

################################################
# Checking the presence of gti header and columns in data/
gticolumns = path.join(datadir,'gti_columns.txt')
gtiheader = path.join(datadir,'gti_header.txt')

if not os.path.isfile(gtiheader) or not os.path.isfile(gticolumns):
    log.error('The files gti_header.txt or gti_columns.txt are missing.  Check the {} directory'.format(os.path.abspath(datadir)))
    exit()

################################################
desc = """
Count rate cut on event file, using ftools (following method by S. Bogdanov).  Automatic if count rate cut is provided, ortherwise interactive (calling sci_plot) 
"""
parser = argparse.ArgumentParser(description = desc)
parser.add_argument("evfiles", help="event file", default = None)
parser.add_argument("--lcfile", help="Light curve file (optional)", type=str, default = None)
parser.add_argument("--cut", help="Count rate cut in cts/sec (optional)", type=float, default=None)
parser.add_argument("--timebin", help="Time bin width in sec (default = 4 sec)", type=float, default=4.0)
parser.add_argument("--plotfilt", help="Ploting filtered lightcurve at the end", default=False, action='store_true')

args = parser.parse_args()
pipedir = os.getcwd()


################################################
##  STEP 0 - open event file and get GTI
eventfile = path.join(pipedir,args.evfiles)
etable = Table.read(eventfile,hdu=1)
eventgti = getgti(eventfile)
log.info('Changing name of TIME column of event file to MET (this is just for the nicer.plotutils.plot_light_curve call)')
etable.columns['TIME'].name = 'MET'


################################################
##  STEP 1 -- making light curve
if not args.lcfile:
    log.info('No light curve file provided. Making light curve with timebin {0} sec'.format(args.timebin))    
    lcfile = path.splitext(eventfile)[0] + ".lcurve"
    cmd = ["extractor", eventfile, "eventsout=none", "imgfile=none",
        "phafile=none", "fitsbinlc={0}".format(lcfile),
        "binlc={0}".format(args.timebin), "regionfile=none", "timefile=none",
        "xcolf=RAWX", "ycolf=RAWY", "tcol=TIME", "ecol=PI", "gti=GTI"]
    runcmd(cmd)
else:
    lcfile = path.join(pipedir,args.lcfile)
    log.info('Using light curve file provided: {0}'.format(lcfile)) 


################################################
## STEP 2 - Setting count rate cut from args or from interactive
if  args.cut:
    log.info("The count rate cut will be performed at {0} cts/sec".format(args.cut))
    CRcut = args.cut
else:        
    log.info("No cut provided. I will now show you the light curve to choose: Please close the display and enter your CRcut")
    plt.subplot(1,1,1)
    meanrate, a = plot_light_curve(etable, False, eventgti, binsize=args.timebin)
    plt.title('Light Curve')
    plt.ylabel('Count rate (c/s)')
    plt.xlabel('Time Elapsed (s)')
    plt.grid()
    plt.show()
    plt.clf()
    CRcut = input('Choose your count rate cut: ')
    log.info("The count rate cut will be performed at {0} cts/sec".format(CRcut))

    
################################################
## STEP 3 - Making Cut with lcfile
lcfile_cut = path.splitext(lcfile)[0] + "_cut.lcurve"
cmd = ["ftcopy", "{0}[1][RATE<{1}]".format(lcfile,CRcut), lcfile_cut, "clobber=yes"]
## Somehow, this line does not work work with os.system().  This is all a mystery to me!
log.info('CMD: '+" ".join(cmd))
check_call(cmd,env=os.environ)


################################################
## STEP 4 - calculate start and end times of remaining bins
log.info("Calculating the start and end times of the remaining bins")
cmd = ['ftcalc', lcfile_cut, lcfile_cut, 'TSTART', '\"TIME-(0.5*{0})+#TIMEZERO\"'.format(args.timebin), "clobber=yes"]
runcmd(cmd)
cmd = ['ftcalc', lcfile_cut, lcfile_cut, 'TEND', '\"TIME+(0.5*{0})+#TIMEZERO\"'.format(args.timebin), "clobber=yes"]
runcmd(cmd)


################################################
## STEP 5 - dumping the TSTART and TEND into text file
log.info("Writing the calculated TSTART and TEND columns into a text file, necessary for ftcreate (in next step)")
cmd = ['ftlist', '{0}[1]'.format(lcfile_cut), 'columns=TSTART,TEND', 'rownum=no', 'colheader=no', 'opt=t', '>', 'gti_data.txt']
runcmd(cmd)


################################################
##  STEP 6 - Making the GTI file from the text file
log.info("Making the GTI file gti.fits from the GTI data textfile")
cmd = ['ftcreate', '{}'.format(gticolumns), 'gti_data.txt', 'gti.fits', 'headfile={}'.format(gtiheader), 'extname="GTI"', 'clobber=yes']
runcmd(cmd)


################################################
##  STEP 7 - Extracting the new event file using the new GTI file created
log.info("Making the filtered event file using niextract-event and gti.fits")
outevtfile = path.splitext(eventfile)[0] + "_filt.evt"
cmd = ['niextract-events', '{0}'.format(eventfile), '{0}'.format(outevtfile), 'timefile="gti.fits[GTI]"', 'clobber=yes']
runcmd(cmd)

if args.plotfilt:
    log.info("Showing the filtered light curve")
    filtetable = Table.read(outevtfile,hdu=1)
    filteventgti = getgti(outevtfile)
    filtetable.columns['TIME'].name = 'MET'
    plt.subplot(1,1,1)
    meanrate, a = plot_light_curve(filtetable, False, filteventgti, binsize=args.timebin)
    plt.title('Light Curve')
    plt.xlabel('Time Elapsed (s)')
    plt.grid()
    plt.show()
    plt.clf()

    
################################################
log.info('DONE')
        
