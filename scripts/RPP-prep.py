#!/usr/bin/env python
import argparse
import sys, os
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import time
from subprocess import *
from glob import glob

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
    "--srcname",
    help="Source name (used in output .yml file).",
    default=None,
)
parser.add_argument(
    "--no_psrpipe",
    help="Skip psrpipe.py, if *pipe dirs are already made.",
    default=False,
    action="store_true",
)
parser.add_argument(
    "--rephase",
    help="Skip psrpipe.py, rephase .evt files in existing *pipe dirs with the provided par file.",
    default=False,
    action="store_true",
)
parser.add_argument(
    "--par_range",
    help="Use only data within the time range covered by the provided par file.",
    default=False,
    action="store_true",
)
parser.add_argument(
    "--cps_cut_percentile",
    help="Data bins above this percentile in terms of counts/s will be discarded.",
    type=float,
    default=95
)
parser.add_argument(
    "--emin",
    help="Emin (keV) of the optimal energy range, if known (if Emin and Emax are not both specified, nioptcuts.py will be used to determine them).",
    type=float,
    default=None
)
parser.add_argument(
    "--emax",
    help="Emax (keV) of the optimal energy range, if known.",
    type=float,
    default=None
)
# NEW: expose Bayesian Blocks p0 as a CLI option (passed through to RPP-profile.py)
parser.add_argument(
    "--p0",
    help="Bayesian Blocks false-positive rate (passed to RPP-profile.py). Typical values: 0.005–0.3.",
    type=float,
    default=0.005
)

args = parser.parse_args()

par = args.parfile
ddir = args.datadir

# -------------------------------------

if not args.no_psrpipe and not args.rephase:
    cmd = (
        "psrpipe.py --nomap --emin 0.22 --emax 15.01 --tidy --cormin 1.5 --kpmax 5 --mingti 100 --maxovershoot 1.5 --maxundershoot 600 --medianundershoot=100 --par "
        + par
        + " --ephem DE421 "
        + ddir
        + "/[0123456789]*"
    )
    print(cmd)
    os.system(cmd)

# -------------------------------------
    
if args.rephase:
    pipe_dirs = glob('[0-9]*_pipe')
    
    for p in pipe_dirs:
        eventfile = p+'/cleanfilt.evt'
        # Accept both plain .orb and gzipped .orb.gz (use wildcard .orb*)
        orbfile = glob(p+'/*.orb*')[0]

        cmd = 'photonphase --orbfile '+orbfile+' --addphase '+eventfile+' '+args.parfile
        print(cmd)
        os.system(cmd)
    
# -------------------------------------

if not args.par_range:
    cmd = "ls -1 [0-9]*_pipe/cleanfilt.evt > files.txt"
    print(cmd)
    os.system(cmd)
else:
    # Get the START/FINISH from the par file:
    cmd = 'grep START '+args.parfile
    pr = Popen(cmd,shell=True,stdout=PIPE)
    pr.wait()
    output = pr.communicate()[0]
    par_start = float(output.split()[-1])

    cmd = 'grep FINISH '+args.parfile
    pr = Popen(cmd,shell=True,stdout=PIPE)
    pr.wait()
    output = pr.communicate()[0]
    par_end = float(output.split()[-1])

    print('Par file start: ',par_start,' end: ',par_end)

    # Check which cleanfilt.evt files are within the par time range and include only them in the list to be merged.
    evt_files = glob('[0-9]*_pipe/cleanfilt.evt')
    evt_files.sort()
    
    fd = open('files.txt','w')
    
    for evt in evt_files:
        hdulist = fits.open(evt,mode='denywrite')
        mjdi = hdulist[0].header['MJDREFI']
        mjdf = hdulist[0].header['MJDREFF']
        
        met_start = hdulist[0].header['TSTART']
        evt_start = mjdi + mjdf + met_start/(24.*3600.)
        hdulist.close()

        if evt_start-par_start > -1 and par_end-evt_start > -1:
            fd.write(evt+'\n')
        else:
            print(evt+' (MJD: '+str(evt_start)+') outside par range, will not use it.')
            
    fd.close()

# -------------------------------------

# Collect both .orb and .orb.gz into the list (wildcard .orb*)
cmd = 'ls -1 [0-9]*_pipe/*.orb* > orbfiles.txt'
print(cmd)
os.system(cmd)

cmd = "remove_empty_evtfiles.py files.txt cleanfiles.txt"
print(cmd)
os.system(cmd)

cmd = "niextract-events @cleanfiles.txt merged.evt"
print(cmd)
os.system(cmd)

#cmd = "niextlc merged.evt merged_2-10keV.lc timebin=32 pirange=200:1000 lcthresh=0.9 clobber=yes"
cmd = "nicerl3-lc pirange=200:1000 timebin=32.0 indir=.  mkfile=merged.evt lcfile=merged_2-10keV.lc bkgmodeltype=NONE detnormtype=arr52 clfile=merged.evt ufafile=none clobber=yes lcthresh=0.5"
print(cmd)
os.system(cmd)

# Count rate cut -------------------------------------

lcfile = "merged_2-10keV.lc"
lctable = Table.read(lcfile, hdu=1)
t = lctable["TIME"]
rate = lctable["RATE"]
p95 = np.percentile(rate, 95)
p99 = np.percentile(rate, 99)
av = np.mean(rate)
med = np.median(rate)
std = np.std(rate)
cps_cut = np.percentile(rate,args.cps_cut_percentile)
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
#plt.title("ENTER DESIRED CPS CUT AFTER CLOSING THIS PLOT")
plt.title('CPS av: %.2f  med: %.2f  95th%%: %.2f  99th%%: %.2f' % (av,med,p95,p99))
plt.suptitle('Chosen CPS cut percentile: %.2f value: %.2f counts/s' % (args.cps_cut_percentile, cps_cut))
#plt.show()
plt.savefig(args.srcname+'_crcut.png')

print('CPS av: %.2f  med: %.2f  95th%%: %.2f  99th%%: %.2f' % (av,med,p95,p99))
print('Chosen CPS cut percentile: %.2f value: %.2f counts/s' % (args.cps_cut_percentile, cps_cut))

cmd = (
    "cr_cut.py merged.evt --outname merged_cut.evt --cut "
    + str(cps_cut)
    + " --filterbinsize=32.0 --lcfile merged_2-10keV.lc"
)

print(cmd)
os.system(cmd)

# -------------------------------------

cmd = (
    'niphaseogram.py --subtract_min --outfile '
    + args.srcname
    + '_phaseogram.png merged_cut.evt'
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

# Prep for pulsed analysis -------------------------------------

# Find optimal energy range
if args.emin is None or args.emax is None:
    cmd = 'nioptcuts.py --noplot merged_cut.evt'
    print(cmd)
    print('(This may take a long time - a few hours if there is lots of data.)')
    pr = Popen(cmd,shell=True,stdout=PIPE,stderr=sys.stderr)
    pr.wait()
    output = pr.communicate()[0]
    output = output.split()
    emin_kev = output[-6].decode()
    emax_kev = output[-2].decode()
else:
    emin_kev = args.emin
    emax_kev = args.emax
    
print('Emin,Emax (keV): ',emin_kev,emax_kev)

# Build the RPP-profile command and inject user-selected p0
# NOTE: --p0 controls the Bayesian Blocks false-positive rate in off-pulse detection.
if args.srcname is None:
    cmd = 'RPP-profile.py --outbase prof --optemin '+str(emin_kev)+' --optemax '+str(emax_kev)+' --p0 '+str(args.p0)+' merged_cut.evt'
else:
    cmd = 'RPP-profile.py --srcname '+args.srcname+' --outbase '+args.srcname+'_prof'+' --optemin '+str(emin_kev)+' --optemax '+str(emax_kev)+' --p0 '+str(args.p0)+' merged_cut.evt'
print(cmd)
os.system(cmd)

# This makes the file load_pulsedspec.py for further automated processing (and load_pulsedspec.xcm that can be used for further manual processing)
if args.srcname is None:
    cmd = 'RPP-pulsedspec.py merged_cut.evt merged_cut_profinfo.yml merged_cutmpu7RPP.rmf merged_cutmpu7RPP.arf'
else:
    cmd = 'RPP-pulsedspec.py merged_cut.evt '+args.srcname+'_profinfo.yml merged_cutmpu7RPP.rmf merged_cutmpu7RPP.arf'
print(cmd)
os.system(cmd)