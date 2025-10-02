#!/usr/bin/env python

import numpy as np
import argparse
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from subprocess import check_call
import os, sys
from nicer.values import *
from loguru import logger as log

log.remove()
log.add(sys.stdout, format="{level}: {message}")


def runcmd(cmd):
    # CMD should be a list of strings since it is not processed by a shell
    log.debug("CMD: " + " ".join(cmd))
    # log.info(cmd)
    check_call(cmd, env=os.environ)


parser = argparse.ArgumentParser(
    description="Plot per-detector exposure from an .evt or mkf file"
)
parser.add_argument("infile", help="Input file to process")

args = parser.parse_args()

outfname = "expo_plot.txt"
cmd = ["niexposum", f"infile={args.infile}", f"outfile={outfname}"]
runcmd(cmd)

detids, exp = np.loadtxt(outfname, unpack=True)
detids = detids.astype(int)

if args.infile.endswith(".evt"):
    total_exp = pyfits.getval(args.infile, "EXPOSURE", ext=1)
    log.info(f"Total exposure in events file, for scaling = {total_exp}")
else:
    log.info(
        "Not using EVENT file, so scaling to detector with maximum exposure instead of total EXPOSURE."
    )
    total_exp = exp.max()

# for d, e in zip(detid, exp):
#    print(d, e / total_exp)

dets_off = set(detids[exp == 0])

dead_dets_not_off = DEAD_DETS - dets_off
if len(dead_dets_not_off) > 0:
    log.warning(f"Found dead detectors unexpectedly on: {dead_dets_not_off}")

live_dets_off = dets_off - DEAD_DETS
if len(live_dets_off) > 0:
    log.info(f"Live detectors with no exposure: {live_dets_off}")

log.info(f"Minimum live detector on fraction: {exp[exp>0].min()/total_exp}")
log.info(f"Maximum live detector on fraction: {exp.max()/total_exp}")

fig, ax = plt.subplots()
ax.bar(detids, exp / total_exp)
ax.set_title("Fractional Exposure by DET_ID")
ax.set_ylabel("DET_ID")
ax.set_xlabel("Fraction Active")
plt.show()
