#!/usr/bin/env python
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


def runcmd(cmd):
    # CMD should be a list of strings since it is not processed by a shell
    log.info("CMD: " + " ".join(cmd))
    os.system(" ".join(cmd))
    ## Some ftools calls don't work properly with check_call...not sure why!
    ## so I am using os.system instead of check_call
    # check_call(cmd,env=os.environ)


################################################
# Checking the presence of HEASOFT
try:
    check_call("nicerversion", env=os.environ)
except:
    print(
        "You need to initialize FTOOLS/HEASOFT first (e.g., type 'heainit')!",
        file=sys.stderr,
    )
    exit()

################################################
# Checking the presence of gti header and columns in data/
gticolumns = path.join(datadir, "gti_columns.txt")
gtiheader = path.join(datadir, "gti_header.txt")

if not os.path.isfile(gtiheader) or not os.path.isfile(gticolumns):
    log.error(
        "The files gti_header.txt or gti_columns.txt are missing.  Check the {} directory".format(
            os.path.abspath(datadir)
        )
    )
    exit()


desc = """
Create a simple GTI file from a pair of NICER METs. This is handy as an input file to niextract-events timefile=xxx.gti
"""
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("startmet", help="Starting MET for GTI", type=float)
parser.add_argument("stopmet", help="Ending MET for GTI", type=float)
parser.add_argument(
    "--gtiname",
    help="Name of output GTI FITS file (default gti.fits)",
    default="gti.fits",
)

args = parser.parse_args()

################################################
## STEP 5 - dumping the TSTART and TEND into text file
## Made the following alteration on July 15 2019 (Mason Ng):
## Wrapped "repr" around the arguments, in addition to explicitly specifying double-precision for startmet and stopmet.
## This is because fp.write was somehow truncating the number of digits to millisecond precision, instead of sub-ms!
## Also note that in Python3, repr is renamed to reprlib.

import tempfile

fp = tempfile.NamedTemporaryFile(mode="w")
fp.write("{0} {1}\n".format(repr(args.startmet), repr(args.stopmet)))
fp.flush()

################################################
##  STEP 6 - Making the GTI file from the text file
log.info("Making the GTI file gti.fits from the GTI data textfile")
cmd = [
    "ftcreate",
    "{}".format(gticolumns),
    fp.name,
    args.gtiname,
    "headfile={}".format(gtiheader),
    'extname="GTI"',
    "clobber=yes",
]
runcmd(cmd)

fp.close()
