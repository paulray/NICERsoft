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

parser = argparse.ArgumentParser(description = "Process NICER data.  Output will be written in current working directory.")
parser.add_argument("indirs", help="Input directories to process", nargs='+')
parser.add_argument("--emin", help="Minimum energy to include (keV, default=0.4)", type=float, default=0.4)
parser.add_argument("--emax", help="Maximum energy to include (kev, default=8.0)", type=float, default=8.0)
parser.add_argument("--mask",help="Mask these IDS", nargs = '*', type=int, default=None)
parser.add_argument("--filtpolar",help="Filter polar horn regions from data",action='store_true')
parser.add_argument("--maxovershoot",help="Select data where overshoot rate is below this limit (default: no filter)",
    type=float,default=-1)
parser.add_argument("--badcut",help="Select data where bad ratio event rate is below this limit (default: no filter)",
    type=float,default=-1)
parser.add_argument("--obsid", help="Use this as OBSID for directory and filenames",
    default=None)
parser.add_argument("--ultraclean", help="Apply ultraclean background filters", action='store_true')
parser.add_argument("--par", help="Par file to use for phases")
args = parser.parse_args()

os.environ['HEADASNOQUERY'] = ' '
os.environ['HEADASPROMPT'] = '/dev/null'

# For some reason if the script is called via #!/usr/bin/env python
# it does not inherit LD_LIBRARY_PATH so ftools don't run.
#print(os.environ['LD_LIBRARY_PATH'])
#print(os.environ['HEADAS'])
#os.environ['LD_LIBRARY_PATH'] = path.join(os.environ['HEADAS'],'lib')
#print(os.environ['LD_LIBRARY_PATH'])

def runcmd(cmd):
    # CMD should be a list of strings since it is not processed by a shell
    log.info('CMD: '+" ".join(cmd))
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
    pipedir = "{0}_pipe".format(basename)
    if not os.path.exists(pipedir):
        os.makedirs(pipedir)

    log.info('Making initial QL plots')
    cmd = ["master_plotter.py", "--save", "--filtall",
           "--writebkf", "--eventshootrate",
           "--guessobj", "--lclog", "--useftools",
           "--emin", "{0}".format(args.emin), "--emax", "{0}".format(args.emax),
           "--sci", "--eng", "--bkg", "--obsdir", obsdir,
           "--basename", path.join(pipedir,basename)+'_prefilt']
    if args.par is not None:
        cmd.append("--par")
        cmd.append("{0}".format(args.par))
    runcmd(cmd)


    # Get event filenames (could be just one)
    evfiles = glob(path.join(obsdir,'xti/event_cl/ni*mpu?_cl.evt'))
    evfiles.sort()
    log.info('Event Files: {0}'.format("\n" + "    \n".join(evfiles)))

    # Get orbit file
    orbfile = glob(path.join(obsdir,'auxil/ni*.orb'))[0]
    log.info('Orbit File: {0}'.format(orbfile))
    # Copy orbit file to results dir for pulsar analysis
    shutil.copy(orbfile,pipedir)

    # Get filter file
    mkfile = glob(path.join(obsdir,'auxil/ni*.mkf'))[0]
    log.info('MKF File: {0}'.format(mkfile))

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

    # Create GTI from .mkf file
    # Not filtering on SAA for now since it is computed from an old map
    if args.ultraclean:
        mkf_expr='(ANG_DIST.lt.0.01).and.(ELV>30.0).and.(SUNSHINE.eq.0)'
    else:
        mkf_expr='(ANG_DIST.lt.0.01).and.(ELV>30.0)'
    saafile = path.join(datadir,'saa.reg')
    mkf_expr += '.and.regfilter({0},SAT_LON,SAT_LAT)'.format(saafile)
    if args.filtpolar:
        phfile = path.join(datadir,'polarhorns.reg')
        mkf_expr += '.and.regfilter({0},SAT_LON,SAT_LAT)'.format(phfile)
    gtiname1 = path.join(pipedir,'mkf.gti')
    cmd = ["maketime", mkfile, gtiname1, 'expr={0}'.format(mkf_expr),
        "compact=no", "time=TIME",  "prefr=0", "postfr=0", "clobber=yes"]
    runcmd(cmd)
    if len(Table.read(gtiname1,hdu=1))==0:
        log.error('No good time left after filtering!')
        continue

    # Create GTI from attitude data
    gtiname2 = path.join(pipedir,'att.gti')
    att_expr = '(MODE.eq.1).and.(SUBMODE_AZ.eq.2).and.(SUBMODE_EL.eq.2)'
    cmd = ["maketime", attfile, gtiname2, 'expr={0}'.format(att_expr),
        "compact=no", "time=TIME", "prefr=0", "postfr=0", "clobber=yes"]
    runcmd(cmd)
    if len(Table.read(gtiname2,hdu=1))==0:
        log.error('No good time left after filtering!')
        continue

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

    gtiname_merged = path.join(pipedir,"tot.gti")
    try:
        os.remove(gtiname_merged)
    except:
        pass

    if args.maxovershoot > 0 or args.badcut>0:
        cmd = ["mgtime", "{0},{1},{2}".format(gtiname1,gtiname2,gtiname3), gtiname_merged, "AND"]
    else:
        cmd = ["mgtime", "{0},{1}".format(gtiname1,gtiname2), gtiname_merged, "AND"]
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
    evfilt_expr = 'PI={0}:{1},EVENT_FLAGS=bx1x000'.format(
        int(args.emin*KEV_TO_PI), int(args.emax*KEV_TO_PI))

    cmd = ["niextract-events", "filename=@{0}[{1}]".format(evlistname,evfilt_expr),
        "eventsout={0}".format(intermediatename), "timefile={0}".format(gtiname_merged),
        "gti=gti", "clobber=yes"]
    runcmd(cmd)


    ### Final step filters and masked detector and does ratio filter
    maxratio = 1.4
    if args.ultraclean:
        maxratio=1.14
    evfilt_expr = '((float)PHA/(float)PHA_FAST < {0})'.format(maxratio)
    if args.mask is not None:
        for detid in args.mask:
            evfilt_expr += ".and.(DET_ID!={0})".format(detid)
    cmd = ["ftcopy", "{0}[{1}]".format(intermediatename,evfilt_expr), filteredname,
        "clobber=yes", "history=yes"]
    runcmd(cmd)
    # Remove intermediate file
    os.remove(intermediatename)

    # Add phases and plot, if requested
    cmd = ["master_plotter.py", "--save", "--filtratio", "{0}".format(maxratio),
           "--emin", "{0}".format(args.emin), "--emax", "{0}".format(args.emax),
           "--orb", path.join(pipedir,path.basename(orbfile)),
           "--sci", filteredname,
           "--basename", path.join(pipedir,basename)]
    if args.par is not None:
        cmd.append("--par")
        cmd.append("{0}".format(args.par))
    runcmd(cmd)
