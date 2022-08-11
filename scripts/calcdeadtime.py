#!/usr/bin/env python
from glob import glob
import astropy.io.fits as pyfits
import sys, os
from os import path, remove
from astropy import log
from astropy.table import Table
from subprocess import check_call
import argparse
import re
import numpy as np

# from nicer.values import *
# Array of DET_IDs that are used
IDS = np.array(
    [
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        30,
        31,
        32,
        33,
        34,
        35,
        36,
        37,
        40,
        41,
        42,
        43,
        44,
        45,
        46,
        47,
        50,
        51,
        52,
        53,
        54,
        55,
        56,
        57,
        60,
        61,
        62,
        63,
        64,
        65,
        66,
        67,
    ]
)

import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(
    description="Compute deadtime correction to an EXPOSURE defined by a GTI extension, for a single OBSID."
)
parser.add_argument("obsdir", help="Directory containing the raw data for this OBSID")
parser.add_argument(
    "gtifile",
    help="FITS file containing a GTI extension to be used. Can be an event file, PHA file or any FITS file with a 'GTI' extension.",
)
parser.add_argument(
    "--mask", help="Mask particular FPMs", nargs="+", type=int, default=[]
)
parser.add_argument("--plot", help="Plot deadtime per FPM", action="store_true")
args = parser.parse_args()

# The GTI file is assumed to apply to all FPMs. This is normally the case since the user
# is operating on a merged event file whose GTI is the AND of all the individual MPU GTIs
# then they may make additional GTI selections that are more restrictive than that.
# So, we can go over each MPU file and apply the GTI before counting up the deadtime.

# Get the names of all the individual MPU files
gstr = path.join(args.obsdir, "xti/event_uf/ni*mpu?_uf.evt*")
log.debug("Glob string {}".format(gstr))
ufiles = glob(gstr)
ufiles.sort()
log.info(
    "Reading unfiltered events from these files :\n\t{}".format("\n\t".join(ufiles))
)
if len(ufiles) != 7:
    log.error("Did not find 7 MPU files!")
fpm_deadtime = np.zeros(len(IDS))
t_mpu = -1
log.info("Mask {}".format(args.mask))
for i, det_id in enumerate(IDS):
    if det_id in args.mask:
        continue
    mpu = det_id // 10
    log.debug("{} DET_ID {} MPU {} File {}".format(i, det_id, mpu, ufiles[mpu]))
    # Only read the raw MPU file once per MPU since all the FPMs for this MPU are in this file
    if mpu != t_mpu:
        cmd = (
            "niextract-events {0} eventsout={1} timefile='{2}[GTI]' clobber=yes".format(
                ufiles[mpu], "tmp.evt", args.gtifile
            )
        )
        st = check_call(cmd, shell=True)
        if st != 0:
            log.error("niextract-events failed!")
        t = Table.read("tmp.evt", hdu=1)
        t_mpu = mpu
    dets = t["DET_ID"]
    if not np.any(dets == det_id):
        fpm_deadtime[i] = 0.0
    else:
        fpm_deadtime[i] = (t["DEADTIME"][dets == det_id]).sum()

gtitable = Table.read("{}".format(args.gtifile), hdu="GTI")
exp = (gtitable["STOP"] - gtitable["START"]).sum()
log.debug("exp {}".format(exp))
active = np.where(fpm_deadtime > 0)[0]
if not np.any(fpm_deadtime > 0):
    deadtime = 0.0
    mindead = 0.0
    maxdead = 0.0
    stddead = 0.0
else:
    deadtime = fpm_deadtime[active].mean()
    mindead = fpm_deadtime[active].min()
    maxdead = fpm_deadtime[active].max()
    stddead = fpm_deadtime[active].std()
if args.plot:
    if exp > 0:
        plt.plot(IDS, 100 * fpm_deadtime / exp, "s")
    plt.xlabel("DET_ID")
    plt.ylabel("Deadtime %")
    plt.title(t.meta["OBS_ID"])
    # plt.savefig("deadtimeplots/{0}_deadtimes.png".format(t.meta["OBS_ID"]))
    plt.show()
if exp == 0.0:
    percent_frac = 0.0
else:
    percent_frac = 100.0 * deadtime / exp
print(
    "\nFile {} Exposure {:12.5f}, Mean Deadtime {:12.5f} ({:.3f} %) -> Livetime {:12.5f}".format(
        args.gtifile, exp, deadtime, percent_frac, exp - deadtime
    )
)
print(
    "Deadtime Statistics for {} FPM: Min {:12.5f} Max {:12.5f} Std {:12.5f}".format(
        len(active), mindead, maxdead, stddead
    )
)
