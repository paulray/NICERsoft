#!/usr/bin/env python
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
from nicer.plotutils import find_hot_detectors
desc = """
Merging NICER ObsID or event files.

Output will be written in merged directory. This script assumes that all filtering has been done.

Output files include:
* merged event file with it science plot (nicerql.py)
* merged lightcurve produced by extrator
* merged spectrum PHA produced by extrator


The PULSE_PHASE column will be ignored for the nicerql.py plot, since it requires a merged orbit file -- not implemented yet.
"""

parser = argparse.ArgumentParser(description = desc)
parser.add_argument("infile", help="list of event files or text file with list", nargs='+')
parser.add_argument("outroot", help="root of the output filenames", type=str)
parser.add_argument("outdir", help="Name of output directories", type=str)
parser.add_argument("--par", help="Par file to use for phases (in nicerql.py)")
parser.add_argument("--orb", help="text file containing the paths to all orbits files")
parser.add_argument("--cut", help="perform count rate cut at the end", action='store_true')
parser.add_argument("--clobber", help="replace existing merged files", action='store_true')
args = parser.parse_args()


os.environ['HEADASNOQUERY'] = ' '
os.environ['HEADASPROMPT'] = '/dev/null'

# Checking the presence of HEASOFT
try:
    check_call('nicerversion',env=os.environ)
except:
    print("You need to initialize FTOOLS/HEASOFT first (e.g., type 'heainit')!", file=sys.stderr)
    exit()

# For some reason if the script is called via #!/usr/bin/env python
# it does not inherit LD_LIBRARY_PATH so ftools don't run.
#print(os.environ['LD_LIBRARY_PATH'])
#print(os.environ['HEADAS'])
#os.environ['LD_LIBRARY_PATH'] = path.join(os.environ['HEADAS'],'lib')
#print(os.environ['LD_LIBRARY_PATH'])

all_evfiles = []

def runcmd(cmd):
    # CMD should be a list of strings since it is not processed by a shell
    log.info('CMD: '+" ".join(cmd))
    log.info(cmd)
    check_call(cmd,env=os.environ)

## Check if outdir contains 'None', 'NONE', or 'none' (causes bug in ni-extractevents)

names = ['none', 'None', 'NONE']
if any(st in args.outdir for st in names):
    log.error("Due to a current bug in ni-extractevents, outdir cannot contain 'none', 'None', or 'NONE'.  Existing...")
    exit()
if any(st in args.outroot for st in names):
    log.error("Due to a current bug in ni-extractevents, outroot cannot contain 'none', 'None', or 'NONE'.  Existing...")
    exit()
    
# Make directory for working files if it does not exist
pipedir = "{0}".format(args.outdir)
if not os.path.exists(pipedir):
    log.info("Creating merged directory: {}".format(args.outdir))
    os.makedirs(pipedir)

## The list of event files all_evfiles is created in the loop above
if len(args.infile)==1:
    if args.infile[0].startswith('@'):
        inputfile = args.infile[0].split('@')[1]
        log.info('Reading input file list: {}'.format(inputfile))
        all_evfiles = open(inputfile,'r').readlines()
        if len(all_evfiles)==1:
            log.error('Why would you try to merge a single event file. Exiting...')
            exit()

        ## if pa
    else:
        log.error('Why would you try to merge a single event file. Exiting...')
        exit()
else:
    all_evfiles = args.infile

log.info('Event Files to merge: {0}'.format("\n      " + "\n      ".join(all_evfiles)))
log.info('NOTE: the input event files will be merged in that order')

if args.par and args.orb:
    if args.orb.startswith('@'):    
        inputorb = args.orb.split('@')[1]
    else:
        inputorb = args.orb
    all_orbfiles = open(inputorb,'r').readlines()
    if len(all_evfiles)!=len(all_orbfiles):
        log.error('Different numbers of event files and orbit files! Exiting...')
        exit()
    else:
        log.info('Orbit Files to use: {0}'.format("\n      " + "\n      ".join(all_orbfiles)))
else:
    if (args.par and not args.orb) or (args.orb and not args.par):
        log.error('You need both --par and --orb for the phase calculations! Exiting...')
        exit()
        
# Build input file for niextract-events
evlistname=path.join(pipedir,'evfiles.txt')
fout = open(evlistname,'w')
for en in all_evfiles:
    print('{0}'.format(en),file=fout)
fout.close()

# Build selection expression for niextract-events
log.info('Merging of input event files into single event file with niextract-event')
outname = path.join(pipedir,"{}.evt".format(args.outroot))
cmd = ["niextract-events", "filename=@{0}".format(evlistname),
       "eventsout={0}".format(outname)]
if args.clobber is not None:
    cmd.append("clobber=yes")
runcmd(cmd)

# Call to cr_cut.py
if args.cut:
    cmd = ["cr_cut.py", "{}.evt".format(outname), "--lcfile {}".format(lcfile), "--plotfilt"]
    runcmd(cmd)

# Make final merged clean plot
cmd = ["nicerql.py", "--save", "--merge",
       "--sci", outname, "--lcbinsize", "4.0","--allspec","--alllc",
       "--basename", path.splitext(outname)[0]]
if args.par is not None:
    cmd.append("--par")
    cmd.append("{0}".format(args.par))
    cmd.append("--orb")
    cmd.append("@{0}".format(args.orb))
runcmd(cmd)

# Make Phaseogram 
if args.par is not None:
    plotfile = path.join(pipedir,"{}_phaseogram.png".format(args.outroot))
    cmd = ["photonphase", "--fix", "--orb", "@{0}".format(args.orb) ,
           "--plot", "--plotfile", plotfile, "--addphase", outname, args.par]
    runcmd(cmd)

# Extract simple PHA file and light curve
phafile = path.join(pipedir,"{}.pha".format(args.outroot)) 
lcfile = path.join(pipedir,"{}.lc".format(args.outroot)) 
cmd = ["extractor", outname, "eventsout=none", "imgfile=none",
       "phafile={0}".format(phafile), "fitsbinlc={0}".format(lcfile),
       "binlc=1.0", "regionfile=none", "timefile=none",
       "xcolf=RAWX", "ycolf=RAWY", "tcol=TIME", "ecol=PI", "xcolh=RAWX",
       "ycolh=RAWY", "gti=GTI"]
runcmd(cmd)


