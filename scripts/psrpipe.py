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

from nicer.values import *

parser = argparse.ArgumentParser(description = "Process NICER pulsar data.  Output will be written in current working directory.")
parser.add_argument("indirs", help="Input directories to process", nargs='+')
parser.add_argument("--obsid", help="Use this as OBSID for directory and filenames",
    default=None)
args = parser.parse_args()

for obsdir in indirs:

    # Set up a basename and make a work directory
    if args.obsid is not None:
        basename = args.obsid
    else:
        # Trim trailing / if needed
        if args.obsid.endswith("/"):
            args.obsid = args.obsid[:-1]
        basename = path.basename(args.obsid).strip()

    # Make directory for working files and output
    pipedir = "{0}_pipe".format(basename)
    if not os.path.exists(pipedir):
        os.makedirs(pipedir)

    # Get event filenames (could be just one)
    evfiles = glob(path.join(obsdir,'xti/event_cl/ni*mpu?_cl.evt'))
    evfiles.sort()

    # Get orbit file
    orbfile = glob(path.join(obsdir,'auxil/ni*.orb'))[0]

    # Get filter file
    mkfile = glob(path.join(obsdir,'auxil/ni*.mkf'))[0]

    #  Get MPU hk files
    hkfiles = glob(path.join(obsdir,'xti/hk/ni*.hk'))
    hkfiles.sort()

    # Create GTI from .mkf file
    cmd = "maketime {0} good.gti expr="(SAA.eq.0).and.(FOV_FLAG.eq.0).and.(ANG_DIST.lt.0.02)" compact=no time=TIME clobber=yes".format(mkfile)
    check_call(cmd,shell=True)
