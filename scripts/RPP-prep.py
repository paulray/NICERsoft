#!/usr/bin/env python
import argparse
import sys, os
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import time

# This script is to be run from dir where we want *pipe dirs to be written

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""Make files needed for RPP catalog spectral analysis.
    NOTE: This script is to be called from within the dir where we want *_pipe output dirs to be written.
    """,
)

parser.add_argument(
    "parfile", help="Full path to the par file (possibly in /mnt/nicer/pars).", type=str
)
parser.add_argument(
    "datadir",
    help="Full path to the NICER data dir (usually a subdir in /mnt/nicer/data).",
    type=str,
)
parser.add_argument(
    "--nopipe",
    help="Don't run psrpipe. Use this is all the _pipe dirs are already there",
    default=False,
    action="store_true",
)
# parser.add_argument(
#    "--manual",
#    help="Prepare for manual processing with Xspec (default: prepare for automated processing with pyxspec).",
#    default=False,
#    action="store_true",
# )
args = parser.parse_args()

par = args.parfile
ddir = args.datadir

# -------------------------------------

if not args.nopipe:
    cmd = (
        "psrpipe.py --nomap --merge --emin 0.22 --emax 15.01 --tidy --cormin 1.5 --kpmax 5 --mingti 100 --maxovershoot 1.5 --maxundershoot 600 --medianundershoot=100 --par "
        + par
        + " --ephem DE421 "
        + ddir
        + "/[0123456789]*"
    )
    print(cmd)
    os.system(cmd)

# -------------------------------------

cmd = "ls -1 *pipe/cleanfilt.evt > files.txt"
print(cmd)
os.system(cmd)

cmd = "remove_empty_evtfiles.py files.txt cleanfiles.txt"
print(cmd)
os.system(cmd)

cmd = "niextract-events @cleanfiles.txt merged.evt"
print(cmd)
os.system(cmd)

cmd = "niextlc merged.evt merged_2-10keV.lc timebin=32 pirange=200:1000 lcthresh=0.9 clobber=yes"
print(cmd)
os.system(cmd)

# -------------------------------------

lcfile = "merged_2-10keV.lc"
lctable = Table.read(lcfile, hdu=1)
t = lctable["TIME"]
rate = lctable["RATE"]
p95 = np.percentile(rate, 95)
p99 = np.percentile(rate, 99)
av = np.mean(rate)
med = np.median(rate)
std = np.std(rate)
plt.plot(t, rate, ".")
plt.axhline(y=av, color="black", label="Mean")
plt.axhline(y=med, color="green", label="Median")
plt.axhline(y=av + std, color="red", label="Mean+1/2/3sigma")
plt.axhline(y=av + 2 * std, color="red")
plt.axhline(y=av + 3 * std, color="red")
plt.axhline(y=p95, color="orange", label="95/99 percentile")
plt.axhline(y=p99, color="orange")
plt.legend()
plt.xlabel("Time interval (32 s)")
plt.ylabel("Counts/s")
plt.title("ENTER DESIRED CPS CUT AFTER CLOSING THIS PLOT")
plt.show()
cps_cut = input(
    "Desired cps cut (hit Enter for default=cps of 99th percentile of time intervals): "
).strip() or str(p99)
print("Count rate (cps) cut is: " + str(cps_cut))

cmd = (
    "cr_cut.py merged.evt --outname merged_cut.evt --cut "
    + cps_cut
    + " --filterbinsize=32.0 --lcfile merged_2-10keV.lc"
)
print(cmd)
os.system(cmd)

# -------------------------------------

cmd = "ls -1 " + ddir + "/*/auxil/*.mkf > rawmkffiles.txt"
print(cmd)
os.system(cmd)

cmd = "ftmerge @rawmkffiles.txt columns=TIME,QUATERNION,MPU_UNDERONLY_COUNT,MPU_OVERONLY_COUNT,FPM_ON,SAT_LAT,SAT_LON,FPM_UNDERONLY_COUNT,MEDIAN_UNDERONLY_COUNT,FPM_OVERONLY_COUNT,FPM_RATIO_REJ_COUNT,KP,COR_SAX,SOLAR_PHI copyall=NO clobber=YES mergedraw.mkf"
print(cmd)
os.system(cmd)

logfile = "nicerl3-spect-" + time.strftime("%Y%m%d-%H%M%S") + ".log"


cmd = (
    "nicerl3-spect suffix=RPP bkgmodeltype=SCORPEON bkgformat=script noticerange=22:1500 indir=. cldir=. clfile=merged_cut.evt chatter=3 ufafile=none clobber=YES mkfile=mergedraw.mkf >> "
    + logfile
    + " 2>&1"
)
print(cmd)
os.system(cmd)

cmd = (
    "nicerl3-spect suffix=RPP bkgmodeltype=SCORPEON bkgformat=script noticerange=22:1500 indir=. cldir=. clfile=merged_cut.evt chatter=3 ufafile=none outlang=python clobber=YES mkfile=mergedraw.mkf >> "
    + logfile
    + " 2>&1"
)
print(cmd)
os.system(cmd)
