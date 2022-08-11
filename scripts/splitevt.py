#!/usr/bin/env python
import numpy as np
import astropy.units as u
import os, sys
from astropy import log
import astropy.io.fits as pyfits
import argparse
import tempfile
from nicer.values import *
from subprocess import check_call
import random
import string

# MJDREF from the GTI header
# MJDREF = 56658.0 + 7.775925925925930E-04

# Checking the presence of HEASOFT
try:
    check_call("nicerversion", env=os.environ)
except:
    print(
        "You need to initialize FTOOLS/HEASOFT first (e.g., type 'heainit')!",
        file=sys.stderr,
    )
    sys.exit(1)


def mjd2met(m):
    return (m - MJDREF) * 86400.0


def runcmd(cmd):
    # CMD should be a list of strings since it is not processed by a shell
    log.info("CMD: " + " ".join(cmd))
    # log.info(cmd)
    check_call(cmd, env=os.environ)


def randomString(stringLength=10):
    """Generate a random string of fixed length"""
    letters = string.ascii_lowercase
    return "".join(random.choice(letters) for i in range(stringLength))


def write_chunk(evname, met0, met1, basename, index):
    # Build GTI file
    gticolumns = path.join(datadir, "gti_columns.txt")
    if tel.startswith("SWIFT"):
        gtiheader = path.join(datadir, "gti_header_swift.txt")
    else:
        gtiheader = path.join(datadir, "gti_header.txt")
    fp = tempfile.NamedTemporaryFile(mode="w+")
    fp.write("{0} {1}\n".format(met0, met1))
    fp.flush()
    gtiname = "temp{0}.gti".format(randomString())
    log.debug(gtiname)
    cmd = [
        "ftcreate",
        "{}".format(gticolumns),
        fp.name,
        gtiname,
        "headfile={}".format(gtiheader),
        'extname="GTI"',
        "clobber=yes",
    ]
    runcmd(cmd)
    fp.close()
    # Call nextract-events to extract the file
    cmd = [
        "niextract-events",
        "filename={0}".format(evname),
        "eventsout={0}{1:04d}.evt".format(basename, index),
        "timefile={0}".format(gtiname),
        "gti=GTI",
        "fpmsel=no",
        "clobber=yes",
    ]
    runcmd(cmd)
    os.remove(gtiname)


desc = """Read NICER evt file (often a merged one from many ObsIDs) and split into shorter chunks of a specified length."""

parser = argparse.ArgumentParser(description=desc)
parser.add_argument("eventname", help="FITS file to read events from")
parser.add_argument("--outbase", help="Base name of output file", default="split")
parser.add_argument(
    "--minspan",
    help="Minimum span of each chunk (s) to not be discarded.",
    default=0.0,
    type=float,
)
parser.add_argument(
    "--maxspan", help="Maximum span of each chunk (s).", default=86400.0, type=float
)
parser.add_argument(
    "--maxgap",
    help="Maximum gap to allow in a chunk (s, -1 for don't enforce)",
    default=-1.0,
    type=float,
)
parser.add_argument(
    "--minexp",
    help="Minimum exposure for each chunk to not be discarded (s).",
    default=0.0,
    type=float,
)
# parser.add_argument("--mkf",help="MKF file to use for FPM information",default=None)

## Parse arguments
args = parser.parse_args()

# If no maxgap specified, make it equal to the full span (disabling it)
if args.maxgap < 0.0:
    args.maxgap = args.maxspan

# Load GTIs
f = pyfits.open(args.eventname)

tel = f["EVENTS"].header["TELESCOP"]
log.info("TEL {}".format(tel))
GTINAM = "GTI"
MJDREF = float(f[GTINAM].header["MJDREFI"]) + f[GTINAM].header["MJDREFF"]
try:
    TIMEZERO = float(f[GTINAM].header["TIMEZERO"])
except KeyError:
    TIMEZERO = 0.0

gti_t0 = f[GTINAM].data.field("start")
gti_t0 = MJDREF + (gti_t0 + TIMEZERO) / 86400.0

gti_t1 = f[GTINAM].data.field("stop")
gti_t1 = MJDREF + (gti_t1 + TIMEZERO) / 86400.0

gti_dt = (gti_t1 - gti_t0) * 86400

f.close()

# Loop over the GTIs
fileidx = 0
i = 0
t0 = gti_t0[i]
t1 = t0
span = 0.0
exp = 0.0
while i < len(gti_t0):
    print(i, gti_t0[i], gti_t1[i], gti_dt[i])
    span = (gti_t1[i] - t0) * 86400
    gap = (gti_t0[i] - t1) * 86400
    # log.info("     {0}  {1:.6f} {2:.6f} span {3} gap {4}".format(i,gti_t0[i],gti_t1[i],span,gap))
    if span < args.maxspan and gap < args.maxgap:
        # If still within MAX, just move end time to end of this GTI
        t1 = gti_t1[i]
        exp += gti_dt[i]
        i += 1
    else:
        span = (t1 - t0) * 86400
        # log.info("DEBUG: {0} {1} {2} {3}".format(exp,args.minexp,span,args.minspan))
        if (exp > args.minexp) and (span > args.minspan):
            # This GTI pushes us beyond MAX, so write chunk and start new
            log.info(
                "Writing chunk {0} - {1} (span {2}, exp {3})".format(
                    t0, t1, (t1 - t0) * 86400, exp
                )
            )
            log.debug(
                "METs {} - {}, MJDREF {}".format(mjd2met(t0), mjd2met(t1), MJDREF)
            )
            write_chunk(args.eventname, mjd2met(t0), mjd2met(t1), args.outbase, fileidx)
            fileidx += 1
        else:
            log.info(
                "DISCARDING chunk {0} - {1} (span {2}, exp {3})".format(
                    t0, t1, (t1 - t0) * 86400, exp
                )
            )
        t0 = gti_t0[i]
        t1 = gti_t0[i]
        exp = 0.0
# If there is one segment left, do it.
if (exp > args.minexp) and (span > args.minspan):
    log.info(
        "Writing final chunk {0} - {1} (span {2}, exp {3})".format(
            t0, t1, (t1 - t0) * 86400, exp
        )
    )
    write_chunk(args.eventname, mjd2met(t0), mjd2met(t1), args.outbase, fileidx)
    fileidx += 1
else:
    log.info(
        "DISCARDING final chunk {0} - {1} (span {2}, exp {3})".format(
            t0, t1, (t1 - t0) * 86400, exp
        )
    )
