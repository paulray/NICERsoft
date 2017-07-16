#!/usr/bin/env python
from __future__ import print_function, division
import os
import matplotlib.pyplot as plt
import numpy as np
import argparse
from astropy import log
from os import path
from glob import glob
from subprocess import check_call
import shutil

from nicer.values import *

parser = argparse.ArgumentParser(description = "Process NICER pulsar data.  Output will be written in current working directory.")
parser.add_argument("indirs", help="Input directories to process", nargs='+')
parser.add_argument("--emin", help="Minimum energy to include (keV)", type=float, default=0.3)
parser.add_argument("--emax", help="Maximum energy to include (kev)", type=float, default=8.0)
parser.add_argument("--obsid", help="Use this as OBSID for directory and filenames",
    default=None)
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

    #  Get MPU hk files
    hkfiles = glob(path.join(obsdir,'xti/hk/ni*.hk'))
    hkfiles.sort()
    log.info('MPU HK Files: {0}'.format("\n" + "    \n".join(hkfiles)))

    # Create GTI from .mkf file
    mkf_expr='(SAA.eq.0).and.(ANG_DIST.lt.0.02)'
    gtiname = path.join(pipedir,'good.gti')
    cmd = ["maketime", mkfile, gtiname, 'expr={0}'.format(mkf_expr),
        "compact=no", "time=TIME", "clobber=yes"]
    runcmd(cmd)

    # Build input file for ftmerge
    evlistname=path.join(pipedir,'evfiles.txt')
    fout = file(evlistname,'w')
    evfilt_expr = '(PI>{0}).and.(PI<{1}).and.(EVENT_FLAGS==bx1x000)'.format(
        args.emin*KEV_TO_PI, args.emax*KEV_TO_PI)
    for en in evfiles:
        print('{0}[{1}]'.format(en,evfilt_expr),file=fout)
    fout.close()

    # Run ftmerge
    mergedname = path.join(pipedir,"merged.evt")
    cmd = ["ftmerge", "@{0}".format(evlistname), "outfile={0}".format(mergedname),
        "clobber=yes"]
    runcmd(cmd)

    # Now apply the good GTI to remove SAA and slew time ranges
    filteredname = path.join(pipedir,"clean.evt")
    cmd = ["fltime", mergedname, gtiname, filteredname, "copyall=yes", "clobber=yes"]
    runcmd(cmd)

    # Add phases and plot, if requested
    cmd = ["python", "master_plotter.py", "--filtall", "--emin", "{0}".format(args.emin),
        "--emax", "{0}".format(args.emax), "--obsdir", "{0}".format(args.pipedir),
        "--sci", "--eng", "--bkg"]
    if args.par is not None:
        cmd.append(["--par", "{0}".args.par])
