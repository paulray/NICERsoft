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
from astropy.io import fits
import tempfile

from nicer.values import *
from nicer.plotutils import find_hot_detectors
desc = """
Pipeline process NICER data.

Output will be written in current working directory in directories that end in '_pipe'.

Several diagnostic plots are produced, and the following data processing steps are run:
* Select good times according to the following:
  * (ANG_DIST.lt.0.015).and.(ELV>30.0)
  * (MODE.eq.1).and.(SUBMODE_AZ.eq.2).and.(SUBMODE_EL.eq.2)
  * SAT_LAT, SAT_LONG not in SAA or polar horn regions specified by region files
  * If --dark is specified then also filter on SUNSHINE.eq.0
Optionally, you can filter on overshoot rate or rate of bad ratio events

* Select events according to:
  * EVENT_FLAGS=bx1x000 (SLOW-only or SLOW+FAST events)
  * PI in the specified range (default is 0.25-12.0 keV)
  * Remove events from any DET_ID specified by the --mask parameter

 The final output is 'cleanfilt.evt' and extracted PHA and lightcurve FITS
 files produced by extrator. The event file will have a PULSE_PHASE column
 computed using PINT if the --par file is specified on the command line.
"""
parser = argparse.ArgumentParser(description = desc)
parser.add_argument("indirs", help="Input directories to process", nargs='+')
parser.add_argument("--emin", help="Minimum energy to include (keV, default=0.25)", type=float, default=0.25)
parser.add_argument("--emax", help="Maximum energy to include (kev, default=12.0)", type=float, default=12.0)
parser.add_argument("--mask",help="Mask these IDS", nargs = '*', type=int, default=None)
parser.add_argument("--filtpolar",help="Turn on  filtering polar horn regions from data",default=False,action='store_true')
parser.add_argument("--cormin",help="Set minimum cutoff rigidity (COR_SAX) for nimaketime filtering (default=no COR filtering, typical value = 4)",default=None)
parser.add_argument("--minfpm",help="Set minimum of FPMs active for nimaketime filtering (default=38)",default=38)
parser.add_argument("--uocut",help="Apply Teru's undershoot/overshoot parameter space cut",default=False,action='store_true')
parser.add_argument("--maxovershoot",help="Select data where overshoot rate is below this limit (default: no filter)",
    type=float,default=-1)
parser.add_argument("--badcut",help="Select data where bad ratio event rate is below this limit (default: no filter)",
    type=float,default=-1)
parser.add_argument("--angdist",help="Set threshold for ANG_DIST in call to nimaketime (degrees, default=0.015)", type=float, default=0.015)
parser.add_argument("--obsid", help="Use this as OBSID for directory and filenames",
    default=None)
parser.add_argument("--dark", help="Apply SUNSHINE=0 filter to get only data in Earth shadow", action='store_true')
parser.add_argument("--par", help="Par file to use for phases")
parser.add_argument("--ephem", help="Ephem to use with photonphase", default="DE421")
parser.add_argument("--outdir", help="Add name to output directories (by default: directories end in '_pipe')", default='pipe')
parser.add_argument("--merge", help="Merge all ObsIDs provided into single event list, lightcurve and spectrum (outputdir called 'merged')", action='store_true')
parser.add_argument("--crcut", help="perform count rate cut on merged event file (only if --merge)", action='store_true')
args = parser.parse_args()

os.environ['HEADASNOQUERY'] = ' '
os.environ['HEADASPROMPT'] = '/dev/null'

# Checking the presence of HEASOFT
try:
    check_call('nicerversion',env=os.environ)
except:
    print("You need to initialize FTOOLS/HEASOFT first (e.g., type 'heainit')!", file=sys.stderr)
    exit()

# Create a temporary dir for pfiles and set the environment 
tempdir = tempfile.mkdtemp()
os.mkdir(tempdir+'/pfiles')
headas = os.environ['HEADAS']
os.environ["PFILES"] = tempdir+'/pfiles;'+headas+'/syspfiles'
#print(os.environ["PFILES"])
# Called before each exit() and at end of program to clean up:
#shutil.rmtree(tempdir)

# For some reason if the script is called via #!/usr/bin/env python
# it does not inherit LD_LIBRARY_PATH so ftools don't run.
#print(os.environ['LD_LIBRARY_PATH'])
#print(os.environ['HEADAS'])
#os.environ['LD_LIBRARY_PATH'] = path.join(os.environ['HEADAS'],'lib')
#print(os.environ['LD_LIBRARY_PATH'])

all_evfiles = []
all_orbfiles = []

def runcmd(cmd):
    # CMD should be a list of strings since it is not processed by a shell
    log.info('CMD: '+" ".join(cmd))
    log.info(cmd)
    check_call(cmd,env=os.environ)

# Check if outdir contains 'None', 'NONE', or 'none' (causes bug in ni-extractevents)
if args.outdir:
    names = ['none', 'None', 'NONE']
    if any(st in args.outdir for st in names):
        log.error("Due to a current bug in ni-extractevents, --outdir cannot contain 'none', 'None', or 'NONE'.  Existing...")
        shutil.rmtree(tempdir)
        exit()

# Check if ObsIDs are listed or to be read from file
if len(args.indirs)==1:
    if args.indirs[0].startswith('@'):
        inputdirs = args.indirs[0].split('@')[1]
        log.info('Reading input ObsID list: {}'.format(inputdirs))
        all_obsids = np.loadtxt(inputdirs,dtype=str)
    else:
        all_obsids = args.indirs
else:
    all_obsids = args.indirs

# Start processing all ObsIDs
for obsdir in all_obsids:

    # Set up a basename and make a work directory
    if args.obsid is not None:
        basename = args.obsid
    else:
        # Trim trailing / if needed
        if obsdir.endswith("/"):
            obsdir = obsdir[:-1]
        basename = path.basename(obsdir)

    # Make directory for working files and output
    pipedir = "{0}_{1}".format(basename,args.outdir)
    if not os.path.exists(pipedir):
        os.makedirs(pipedir)

    log.info('Making initial QL plots')

    # Get event filenames (could be just one)
    evfiles = glob(path.join(obsdir,'xti/event_cl/ni*mpu7_cl.evt'))
    evfiles.sort()
    log.info('Cleaned Event Files: {0}'.format("\n" + "    \n".join(evfiles)))

    # Get ufa file (unfiltered events)
    ufafiles = glob(path.join(obsdir,'xti/event_cl/ni*mpu7_ufa.evt'))
    ufafiles.sort()
    log.info('Unfiltered Event Files: {0}'.format("\n" + "    \n".join(evfiles)))

    ufaevents = 0
    for ufafile in ufafiles:
        hdulist = fits.open(ufafile, memmap=True)
        nevents = hdulist[1].header['NAXIS2']
        hdulist.close()
        ufaevents = ufaevents + nevents

    if ufaevents < 10000000:
        cmd = ["nicerql.py", "--save", "--filtall", "--lcbinsize", "4.0","--allspec","--alllc",
               "--lclog", "--useftools", "--extraphkshootrate", "--writebkf",
               "--eventshootrate",
               "--emin", "{0}".format(args.emin), "--emax", "{0}".format(args.emax),
               "--sci", "--eng", "--bkg", "--map", "--obsdir", obsdir,
               "--basename", path.join(pipedir,basename)+'_prefilt']
    else:
        cmd = ["nicerql.py", "--save", "--filtall", "--lcbinsize", "4.0","--allspec","--alllc",
               "--lclog", "--useftools", "--extraphkshootrate", "--writebkf",
               "--emin", "{0}".format(args.emin), "--emax", "{0}".format(args.emax),
               "--sci", "--eng", "--bkg", "--map", "--obsdir", obsdir,
               "--basename", path.join(pipedir,basename)+'_prefilt']
    if args.mask is not None:
        cmd.append("--mask")
        for detid in args.mask:
            cmd.append("{0}".format(detid))
#    if args.par is not None:
#        cmd.append("--par")
#        cmd.append("{0}".format(args.par))
#    if (args.maxovershoot>0) or (args.badcut>0):
#        cmd.append("--writebkf")
    runcmd(cmd)

    # Get orbit file
    orbfile = glob(path.join(obsdir,'auxil/ni*.orb'))[0]
    log.info('Orbit File: {0}'.format(orbfile))
    # Copy orbit file to results dir for pulsar analysis
    shutil.copy(orbfile,pipedir)

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

    #  Get MPU hk files
    hkfiles = glob(path.join(obsdir,'xti/hk/ni*.hk'))
    hkfiles.sort()
    log.info('MPU HK Files: {0}'.format("\n" + "    \n".join(hkfiles)))

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
    # Create GTI from overshoot file using overshoot rate
    if args.maxovershoot > 0:
        gtiname3 = path.join(pipedir,'bkf.gti')
        bkf_expr = 'EV_OVERSHOOT.lt.{0}'.format(args.maxovershoot)
        cmd = ["maketime", bkffile, gtiname3, 'expr={0}'.format(bkf_expr),
            "compact=no", "time=TIME", "prefr=0", "postfr=0", "clobber=yes"]
        runcmd(cmd)
        if len(Table.read(gtiname3,hdu=1))==0:
            log.error('No good time left after filtering!')
            continue

    # Create GTI from overshoot file using bad event lightcurve
    if args.badcut > 0:
        gtiname3 = path.join(pipedir,'bkf.gti')
        bkf_expr = 'BAD_RATIO.lt.{0}'.format(args.badcut)
        cmd = ["maketime", bkffile, gtiname3, 'expr={0}'.format(bkf_expr),
            "compact=no", "time=TIME", "prefr=0", "postfr=0", "clobber=yes"]
        runcmd(cmd)
        if len(Table.read(gtiname3,hdu=1))==0:
            log.error('No good time left after filtering!')
            continue
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
    if args.uocut:
        if extra_expr == "NONE":
            extra_expr = "((float)TOT_OVER_COUNT<50.0+5.5*((float)TOT_UNDER_COUNT/1000.0)**2)"
        else:
            extra_expr += ".and.((float)TOT_OVER_COUNT<50.0+5.5*((float)TOT_UNDER_COUNT/1000.0)**2)"
    cor_string="-"
    if args.cormin is not None:
        cor_string = "{0}-".format(args.cormin)
    cmd = ["nimaketime",  "infile={0}".format(mkfile),
        'outfile={0}'.format(gtiname_merged), 'nicersaafilt=YES',
        'saafilt=NO', 'trackfilt=YES', 'ang_dist={0:.3f}'.format(args.angdist), 'elv=30',
        'br_earth=40', 'cor_range={0}'.format(cor_string), 'min_fpm={0}'.format(args.minfpm),
        'ingtis={0}'.format(extragtis), "clobber=yes",
        'expr={0}'.format(extra_expr),
        'outexprfile={0}'.format(path.join(pipedir,"psrpipe_expr.txt"))]
    runcmd(cmd)

    # nimaketime gives an output GTI file with 1 row and START==STOP in the
    # case where no good time is selected.  This differs from the normal
    # maketime, which produces a GTI file with no rows in that case

    ###  Extract filtered events and apply merged GTI
    filteredname = path.join(pipedir,"cleanfilt.evt")
    intermediatename = path.join(pipedir,"intermediate.evt")

    # Build input file for niextract-events
    evlistname=path.join(pipedir,'evfiles.txt')
    fout = open(evlistname,'w')
    for en in evfiles:
        print('{0}'.format(en),file=fout)
    fout.close()

    # Build selection expression for niextract-events
    # Select events with PI in the selected range, require SLOW trigger (FAST optional)
    # and filter all undershoot, overshoot, and force trigger events
    evfilt_expr = 'PI={0}:{1},EVENT_FLAGS=bx1x000'.format(
        int(args.emin*KEV_TO_PI), int(args.emax*KEV_TO_PI))

    cmd = ["niextract-events", "filename=@{0}[{1}]".format(evlistname,evfilt_expr),
        "eventsout={0}".format(intermediatename), "timefile={0}".format(gtiname_merged),
        "gti=GTI", "clobber=yes"]
    runcmd(cmd)

    # Here we can automatically mask detectors by parsing the intermediate file
    bad_dets = None
    if args.mask is not None and args.mask[0] < 0:
        etable = Table.read(intermediatename,hdu=1)
        # Apply TIMEZERO if needed
        if 'TIMEZERO' in etable.meta:
            log.info('Applying TIMEZERO of {0} to etable'.format(etable.meta['TIMEZERO']))
            etable['TIME'] += etable.meta['TIMEZERO']
            etable.meta['TIMEZERO'] = 0.0
        log.info('Auto-masking detectors')
        bad_dets = find_hot_detectors(etable)
        if bad_dets is not None:
            log.info('Found hot detectors {0}!!'.format(bad_dets))
        # Make intermediate eng plot to show bad detectors
        cmd = ["nicerql.py", "--save",
               "--eng", intermediatename, "--lcbinsize", "4.0","--allspec","--alllc",
               "--basename", path.join(pipedir,basename)+"_intermediate"]
        runcmd(cmd)

    # Now filter any bad detectors
    evfilt_expr = '(EVENT_FLAGS==bx1x000)'
    # Filter any explicitly specified masked detectors
    if args.mask is not None:
        for detid in args.mask:
            evfilt_expr += ".and.(DET_ID!={0})".format(detid)
    # Filter any automatically identified hot detectors
    if bad_dets is not None:
        fout = open(path.join(pipedir,"bad_dets.txt"),"w")
        print("{0}".format(bad_dets),file=fout)
        fout.close()
        for detid in bad_dets:
            evfilt_expr += ".and.(DET_ID!={0})".format(detid)

    cmd = ["ftcopy", "{0}[{1}]".format(intermediatename,evfilt_expr), filteredname,
        "clobber=yes", "history=yes"]
    runcmd(cmd)
    # Remove intermediate file
    #os.remove(intermediatename)

    # Make final clean plot
    cmd = ["nicerql.py", "--save",
           "--orb", path.join(pipedir,path.basename(orbfile)),
           "--sci", "--eng", filteredname, "--lcbinsize", "4.0","--allspec","--alllc",
           "--basename", path.join(pipedir,basename)+"_cleanfilt"]
    if args.par is not None:
        cmd.append("--par")
        cmd.append("{0}".format(args.par))
    runcmd(cmd)

    # Add phases
    if args.par is not None:
        plotfile = path.join(pipedir,"phaseogram.png")
        log.info('Applying phases to {0}'.format(filteredname))
        if fits.getval(filteredname,"TIMEZERO",ext=1) < 0.0:
            log.info('Event file has TIMZERO < 0, so not applying --fix in photonphase!')
            cmd = ["photonphase", "--ephem", args.ephem, "--orb", orbfile, "--plot", "--plotfile",
                plotfile, "--addphase", filteredname, args.par]
        else:
            cmd = ["photonphase", "--ephem", args.ephem, "--fix", "--orb", orbfile, "--plot", "--plotfile",
                plotfile, "--addphase", filteredname, args.par]
        runcmd(cmd)

    # Extract simple PHA file and light curve
    phafile = path.splitext(filteredname)[0] + ".pha"
    lcfile = path.splitext(filteredname)[0] + ".lc"
    cmd = ["extractor", filteredname, "eventsout=none", "imgfile=none",
        "phafile={0}".format(phafile), "fitsbinlc={0}".format(lcfile),
        "binlc=1.0", "regionfile=none", "timefile=none",
        "xcolf=RAWX", "ycolf=RAWY", "tcol=TIME", "ecol=PI", "xcolh=RAWX",
        "ycolh=RAWY", "gti=GTI"]
    runcmd(cmd)

    # if --merge option, Add clean evt file to list of files to merge
    if args.merge:
        all_evfiles.append(filteredname)
        all_orbfiles.append(orbfile)


# Merging all ObsIDs
if args.merge and (len(all_evfiles)>1) :
    ## save list of orbit files
    orbfiles_list = path.join(os.getcwd(),"list_orbfiles.txt")
    np.savetxt(orbfiles_list,all_orbfiles,fmt=['%s'])

    ## Call merge.py
    cmd = ["merge.py"]
    for evt in all_evfiles:
        cmd.append(evt)
    cmd.append("merged")
    cmd.append("merged_{0}".format(args.outdir))
    cmd.append("--clobber")
    if args.par is not None:
        cmd.append("--par")
        cmd.append("{0}".format(args.par))
        cmd.append("--orb")
        cmd.append("{0}".format(orbfiles_list))
    if args.crcut:
        cmd.append("--cut")
    runcmd(cmd)

else:
    if args.crcut:
        log.warning("Count rate cuts are only performed on merged event files (add the option --merge)")

shutil.rmtree(tempdir)
