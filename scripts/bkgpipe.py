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
from nicer.NicerFileSet import *

parser = argparse.ArgumentParser(description = "Process NICER background data.  Output will be written in current working directory.")
parser.add_argument("indirs", help="Input directories to process", nargs='+')
parser.add_argument("--emin", help="Minimum energy to include (keV, default=2.0)",
    type=float, default=0.5)
parser.add_argument("--emax", help="Maximum energy to include (kev, default=8.0)",
    type=float, default=8.0)
parser.add_argument("--filtpolar",help="Filter polar horn regions from data",default=False,action='store_true')
parser.add_argument("--cormin",help="Set minimum cutoff rigidity (COR_SAX) for nimaketime filtering (typical value = 4)",default=None)
parser.add_argument("--dark", help="Apply SUNSHINE=0 filter to get only data in Earth shadow", action='store_true')
parser.add_argument("--mask",help="Mask these IDS", nargs = '*', type=int, default=None)
parser.add_argument("--obsid", help="Use this as OBSID for directory and filenames",
    default=None)
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

    args.obsdir = obsdir
    args.orb = None
    args.sps = None
    args.useftools = True
    args.applygti = None
    args.extraphkshootrate = True
    args.eventshootrate = False
    args.writebkf = True
    args.object = None
    args.gtirows = None
    args.basename = path.join(pipedir,basename)+'_prefilt'
    args.filtall = True
    try:
        nfs = NicerFileSet(args)
        nfs.writebkffile()
    except:
        log.error('NicerFileSet failed for some reason on {0}'.format(obsdir))
        continue

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

    # Create any additional GTIs beyond what nimaketime does...
    extragtis="NONE"
    if args.filtpolar:
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

    # nimaketime gives an output GTI file with 1 row and START==STOP in the
    # case where no good time is selected.  This differs from the normal
    # maketime, which produces a GTI file with no rows in that case

    ### Final step filters and masked detector and does ratio filter
    filteredname = path.join(pipedir,"filtered.bkf")
    expr = 'gtifilter("{0}")'.format(gtiname_merged)
    cmd = ["ftcopy", '{0}[{1}]'.format(bkffile,expr), filteredname,
        "clobber=yes", "history=yes"]
    runcmd(cmd)
