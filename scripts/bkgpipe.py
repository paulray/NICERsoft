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

parser = argparse.ArgumentParser(description = "Process NICER background data.  Output will be written in current working directory.")
parser.add_argument("indirs", help="Input directories to process", nargs='+')
parser.add_argument("--emin", help="Minimum energy to include (keV, default=0.4)",
    type=float, default=0.5)
parser.add_argument("--emax", help="Maximum energy to include (kev, default=8.0)",
    type=float, default=8.0)
parser.add_argument("--mask",help="Mask these IDS", nargs = '*', type=int, default=None)
parser.add_argument("--obsid", help="Use this as OBSID for directory and filenames",
    default=None)
parser.add_argument("--ultraclean", help="Apply ultraclean background filters", action='store_true')
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
    pipedir = "{0}_bkg".format(basename)
    if not os.path.exists(pipedir):
        os.makedirs(pipedir)

    log.info('Making initial QL plots')
    cmd = ["master_plotter.py", "--save", "--filtall",
           "--writebkf", 
           "--guessobj", "--useftools",
           "--emin", "{0}".format(args.emin), "--emax", "{0}".format(args.emax),
           "--sci", "--eng", "--map", "--bkg", "--obsdir", obsdir,
           "--basename", path.join(pipedir,basename)+'_prefilt']
    # I don't think this mask is applied to the  .bkf data generated, but it should be
    if args.mask is not None:
        cmd.append(["--mask"].append(["{0}".format(m) for m in args.mask]))
    runcmd(cmd)


    # Get filter file
    mkfile = glob(path.join(obsdir,'auxil/ni*.mkf'))[0]
    log.info('MKF File: {0}'.format(mkfile))

    #  Get ATT hk files
    attfile = glob(path.join(obsdir,'auxil/ni*.att'))[0]
    log.info('ATT HK File: {0}'.format(attfile))

    #  Get BKF file for filtering based on background indicators
    bkffile = path.join(pipedir,basename)+'_prefilt.bkf'
    log.info('BKF File: {0}'.format(bkffile))

    # Create GTI from .mkf file
    if args.ultraclean:
        mkf_expr='(SAA.eq.0).and.(ANG_DIST.lt.0.01).and.(ELV>30.0).and.(SUNSHINE.eq.0)'
    else:
        mkf_expr='(SAA.eq.0).and.(ANG_DIST.lt.0.01).and.(ELV>30.0)'
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

    gtiname_merged = path.join(pipedir,"tot.gti")
    try:
        os.remove(gtiname_merged)
    except:
        pass

    cmd = ["mgtime", "{0},{1}".format(gtiname1,gtiname2), gtiname_merged, "AND"]
    runcmd(cmd)

    ### Final step filters and masked detector and does ratio filter
    filteredname = path.join(pipedir,"filtered.bkf")
    expr = 'gtifilter("{0}")'.format(gtiname_merged)
    cmd = ["ftcopy", '{0}[{1}]'.format(bkffile,expr), filteredname,
        "clobber=yes", "history=yes"]
    runcmd(cmd)
