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

from nicer.values import *
desc = """
Pipeline process NICER data.

Output will be written in current working directory in directories that end in '_pipe'.

Several diagnostic plots are produced, and the following data processing steps are run:
* Select good times according to the following:
  * (ANG_DIST.lt.0.015).and.(ELV>30.0)
  * (MODE.eq.1).and.(SUBMODE_AZ.eq.2).and.(SUBMODE_EL.eq.2)
  * SAT_LAT, SAT_LONG not in SAA or polar horn regions specified by region files
  * If --dark is specified then also filter on SUNSHINE.eq.0
Optionally, you can filter on overshoot rate or rate of bad ratio events

* Select events according to:
  * EVENT_FLAGS=bx1x000 (SLOW-only or SLOW+FAST events)
  * PI in the specified range (default is 0.25-12.0 keV)
  * Remove events from any DET_ID specified by the --mask parameter

 The final output is 'cleanfilt.evt' and extracted PHA and lightcurve FITS
 files produced by extrator. The event file will have a PULSE_PHASE column
 computed using PINT if the --par file is specified on the command line.
"""
parser = argparse.ArgumentParser(description = desc)
parser.add_argument("indirs", help="Input directories to process", nargs='+')
parser.add_argument("--emin", help="Minimum energy to include (keV, default=0.25)", type=float, default=0.25)
parser.add_argument("--emax", help="Maximum energy to include (kev, default=12.0)", type=float, default=12.0)
parser.add_argument("--mask",help="Mask these IDS", nargs = '*', type=int, default=None)
parser.add_argument("--nofiltpolar",help="Disable filtering polar horn regions from data",default=False,action='store_true')
parser.add_argument("--cormin",help="Set minimum cutoff rigidity (COR_SAX) for nimaketime filtering (typical value = 4)",default=None)
parser.add_argument("--maxovershoot",help="Select data where overshoot rate is below this limit (default: no filter)",
    type=float,default=-1)
parser.add_argument("--badcut",help="Select data where bad ratio event rate is below this limit (default: no filter)",
    type=float,default=-1)
parser.add_argument("--obsid", help="Use this as OBSID for directory and filenames",
    default=None)
parser.add_argument("--dark", help="Apply SUNSHINE=0 filter to get only data in Earth shadow", action='store_true')
parser.add_argument("--par", help="Par file to use for phases")
parser.add_argument("--outdir", help="Add name to output directories (by default: directories end in '_pipe')", default='pipe')
parser.add_argument("--merge", help="Merge all ObsIDs provided into single event list, lightcurve and spectrum (outputdir called 'merged')", action='store_true')
args = parser.parse_args()

os.environ['HEADASNOQUERY'] = ' '
os.environ['HEADASPROMPT'] = '/dev/null'

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

for obsdir in args.indirs:

    # Set up a basename and make a work directory
    if args.obsid is not None:
        basename = args.obsid
    else:
        # Trim trailing / if needed
        if obsdir.endswith("/"):
            obsdir = obsdir[:-1]
        basename = path.basename(obsdir)

    # Make directory for working files and output
    pipedir = "{0}_{1}".format(basename,args.outdir)
    if not os.path.exists(pipedir):
        os.makedirs(pipedir)

    log.info('Making initial QL plots')
    cmd = ["nicerql.py", "--save", "--filtall", "--lcbinsize", "4.0",
           "--lclog", "--useftools", "--extraphkshootrate", "--writebkf",
           "--emin", "{0}".format(args.emin), "--emax", "{0}".format(args.emax),
           "--sci", "--eng", "--bkg", "--map", "--obsdir", obsdir,
           "--basename", path.join(pipedir,basename)+'_prefilt']
    if args.mask is not None:
        cmd.append("--mask")
        for detid in args.mask:
            cmd.append("{0}".format(detid))
    if args.par is not None:
        cmd.append("--par")
        cmd.append("{0}".format(args.par))
    if (args.maxovershoot>0) or (args.badcut>0):
        cmd.append("--writebkf")
    runcmd(cmd)


    # Get event filenames (could be just one)
    evfiles = glob(path.join(obsdir,'xti/event_cl/ni*mpu7_cl.evt'))
    evfiles.sort()
    log.info('Cleaned Event Files: {0}'.format("\n" + "    \n".join(evfiles)))

    # Get ufa file (unfiltered events)
    ufafiles = glob(path.join(obsdir,'xti/event_cl/ni*mpu7_ufa.evt'))
    ufafiles.sort()
    log.info('Unfiltered Event Files: {0}'.format("\n" + "    \n".join(evfiles)))

    # Get orbit file
    orbfile = glob(path.join(obsdir,'auxil/ni*.orb'))[0]
    log.info('Orbit File: {0}'.format(orbfile))
    # Copy orbit file to results dir for pulsar analysis
    shutil.copy(orbfile,pipedir)

    # Get filter file
    mkfile = glob(path.join(obsdir,'auxil/ni*.mkf'))[0]
    log.info('MKF File: {0}'.format(mkfile))
    # Copy orbit file to results dir for pulsar analysis
    shutil.copy(mkfile,pipedir)

    #  Get ATT hk files
    attfile = glob(path.join(obsdir,'auxil/ni*.att'))[0]
    log.info('ATT HK File: {0}'.format(attfile))

    #  Get BKF file for filtering based on background indicators
    bkffile = path.join(pipedir,basename)+'_prefilt.bkf'
    log.info('BKF File: {0}'.format(bkffile))

    #  Get MPU hk files
    hkfiles = glob(path.join(obsdir,'xti/hk/ni*.hk'))
    hkfiles.sort()
    log.info('MPU HK Files: {0}'.format("\n" + "    \n".join(hkfiles)))

    # Create any additional GTIs beyond what nimaketime does...
    extragtis="NONE"
    if not args.nofiltpolar:
        saafile = path.join(datadir,'saa.reg')
        mkf_expr = 'regfilter("{0}",SAT_LON,SAT_LAT)'.format(saafile)
        phfile = path.join(datadir,'polarhorns.reg')
        mkf_expr += '.and.regfilter("{0}",SAT_LON,SAT_LAT)'.format(phfile)
        gtiname2 = path.join(pipedir,'extra.gti')
        cmd = ["pset", "maketime", "expr={0}".format(mkf_expr)]
        runcmd(cmd)
        cmd = ["maketime", "infile={0}".format(mkfile), "outfile={0}".format(gtiname2),
            "compact=no", "time=TIME",  "prefr=0", "postfr=0", "clobber=yes"]
        runcmd(cmd)
        if len(Table.read(gtiname2,hdu=1))==0:
            log.error('No good time left after filtering!')
            continue
        extragtis = gtiname2

    gtiname3 = None
    # Create GTI from overshoot file using overshoot rate
    if args.maxovershoot > 0:
        gtiname3 = path.join(pipedir,'bkf.gti')
        bkf_expr = 'EV_OVERSHOOT.lt.{0}'.format(args.maxovershoot)
        cmd = ["maketime", bkffile, gtiname3, 'expr={0}'.format(bkf_expr),
            "compact=no", "time=TIME", "prefr=0", "postfr=0", "clobber=yes"]
        runcmd(cmd)
        if len(Table.read(gtiname3,hdu=1))==0:
            log.error('No good time left after filtering!')
            continue

    # Create GTI from overshoot file using bad event lightcurve
    if args.badcut > 0:
        gtiname3 = path.join(pipedir,'bkf.gti')
        bkf_expr = 'BAD_RATIO.lt.{0}'.format(args.badcut)
        cmd = ["maketime", bkffile, gtiname3, 'expr={0}'.format(bkf_expr),
            "compact=no", "time=TIME", "prefr=0", "postfr=0", "clobber=yes"]
        runcmd(cmd)
        if len(Table.read(gtiname3,hdu=1))==0:
            log.error('No good time left after filtering!')
            continue
    # If either of the bkf filters were used, include that GTI
    # in the extragtis passed to nimaketime
    if gtiname3 is not None:
        if extragtis == "NONE":
            extragtis = gitname3
        else:
            extragtis = extragtis + ',{0}'.format(gtiname3)

    # Make final merged GTI using nimaketime
    gtiname_merged = path.join(pipedir,"tot.gti")
    extra_expr="NONE"
    if args.dark:
        extra_expr = "(SUNSHINE.eq.0)"
    cor_string="-"
    if args.cormin is not None:
        cor_string = "{0}-".format(args.cormin)
    cmd = ["nimaketime",  "infile={0}".format(mkfile),
        'outfile={0}'.format(gtiname_merged), 'nicersaafilt=YES',
        'saafilt=NO', 'trackfilt=YES', 'ang_dist=0.015', 'elv=30',
        'br_earth=40', 'cor_range={0}'.format(cor_string), 'min_fpm=38',
        'ingtis={0}'.format(extragtis), "clobber=yes",
        'expr={0}'.format(extra_expr),
        'outexprfile={0}'.format(path.join(pipedir,"psrpipe_expr.txt"))]
    runcmd(cmd)


    ###  Extract filtered events and apply merged GTI
    filteredname = path.join(pipedir,"cleanfilt.evt")
    intermediatename = path.join(pipedir,"intermediate.evt")

    # Build input file for niextract-events
    evlistname=path.join(pipedir,'evfiles.txt')
    fout = file(evlistname,'w')
    for en in evfiles:
        print('{0}'.format(en),file=fout)
    fout.close()

    # Build selection expression for niextract-events
    # Select events with PI in the selected range, require SLOW trigger (FAST optional)
    # and filter all undershoot, overshoot, and force trigger events
    evfilt_expr = 'PI={0}:{1},EVENT_FLAGS=bx1x000'.format(
        int(args.emin*KEV_TO_PI), int(args.emax*KEV_TO_PI))

    cmd = ["niextract-events", "filename=@{0}[{1}]".format(evlistname,evfilt_expr),
        "eventsout={0}".format(intermediatename), "timefile={0}".format(gtiname_merged),
        "gti=GTI", "clobber=yes"]
    runcmd(cmd)

    evfilt_expr = '(EVENT_FLAGS==bx1x000)'
    if args.mask is not None:
        for detid in args.mask:
            evfilt_expr += ".and.(DET_ID!={0})".format(detid)
    cmd = ["ftcopy", "{0}[{1}]".format(intermediatename,evfilt_expr), filteredname,
        "clobber=yes", "history=yes"]
    runcmd(cmd)
    # Remove intermediate file
    os.remove(intermediatename)

    # Make final clean plot
    cmd = ["nicerql.py", "--save",
           "--orb", path.join(pipedir,path.basename(orbfile)),
           "--sci", filteredname, "--lcbinsize", "4.0",
           "--basename", path.join(pipedir,basename)+"_cleanfilt"]
    if args.par is not None:
        cmd.append("--par")
        cmd.append("{0}".format(args.par))
    runcmd(cmd)

    # Add phases
    if args.par is not None:
        plotfile = path.join(pipedir,"phaseogram.png")
        cmd = ["photonphase", "--fix", "--orb", orbfile, "--plot", "--plotfile",
            plotfile, "--planets", "--addphase", filteredname, args.par]
        runcmd(cmd)

    # Extract simple PHA file and light curve
    phafile = path.splitext(filteredname)[0] + ".pha"
    lcfile = path.splitext(filteredname)[0] + ".lc"
    cmd = ["extractor", filteredname, "eventsout=none", "imgfile=none",
        "phafile={0}".format(phafile), "fitsbinlc={0}".format(lcfile),
        "binlc=1.0", "regionfile=none", "timefile=none",
        "xcolf=RAWX", "ycolf=RAWY", "tcol=TIME", "ecol=PI", "xcolh=RAWX",
        "ycolh=RAWY", "gti=GTI"]
    runcmd(cmd)

    # if --merge option, Add clean evt file to list of files to merge
    if args.merge:
        all_evfiles.append(filteredname)


        
# Merging all ObsIDs
if args.merge and (len(all_evfiles)>1) :
    
    # Make directory for working files and output
    pipedir = "merged_{0}".format(args.outdir)
    if not os.path.exists(pipedir):
        os.makedirs(pipedir)
    log.info('Merging all ObsIDs into single event file with niextract-event')

    ## The list of event files all_evfiles is created in the loop above
    all_evfiles.sort()
    log.info('Cleaned Event Files to Merge: {0}'.format("\n" + "    \n".join(all_evfiles)))

    # Build input file for niextract-events
    evlistname=path.join(pipedir,'evfiles.txt')
    fout = file(evlistname,'w')
    for en in all_evfiles:
        print('{0}'.format(en),file=fout)
    fout.close()

    # Build selection expression for niextract-events
    outname = path.join(pipedir,"merged.evt")
    cmd = ["niextract-events", "filename=@{0}".format(evlistname),
        "eventsout={0}".format(outname), "clobber=yes"]
    runcmd(cmd)

    # Make final merged clean plot
    cmd = ["nicerql.py", "--save",
           "--sci", outname, "--lcbinsize", "4.0",
           "--basename", path.splitext(outname)[0]]
    if args.par is not None:
        log.info('The use of par files requires a merged orbit file -- not implemented yet')
        #cmd.append("--par")
        #cmd.append("{0}".format(args.par))
    runcmd(cmd)

    # Extract simple PHA file and light curve
    phafile = path.splitext(outname)[0] + ".pha"
    lcfile = path.splitext(outname)[0] + ".lc"
    cmd = ["extractor", outname, "eventsout=none", "imgfile=none",
        "phafile={0}".format(phafile), "fitsbinlc={0}".format(lcfile),
        "binlc=1.0", "regionfile=none", "timefile=none",
        "xcolf=RAWX", "ycolf=RAWY", "tcol=TIME", "ecol=PI", "xcolh=RAWX",
        "ycolh=RAWY", "gti=GTI"]
    runcmd(cmd)
