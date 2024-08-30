#!/usr/bin/env python
import os, sys

# import pint.logging
from loguru import logger as log

# pint.logging.setup(level=pint.logging.script_level)

import numpy as np
import argparse
from os import path
from glob import glob
from subprocess import check_call, CalledProcessError
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
* Select good times according to customizable criteria, including:
  * (ANG_DIST.lt.0.015).and.(ELV>20.0)
  * (MODE.eq.1).and.(SUBMODE_AZ.eq.2).and.(SUBMODE_EL.eq.2)
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
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("indirs", help="Input directories to process", nargs="+")
parser.add_argument(
    "--emin",
    help="Minimum energy to include (keV, default=0.25)",
    type=float,
    default=0.25,
)
parser.add_argument(
    "--emax",
    help="Maximum energy to include (kev, default=15.0)",
    type=float,
    default=15.0,
)
parser.add_argument("--mask", help="Mask these IDS", nargs="*", type=int, default=None)
parser.add_argument(
    "--filtpolar",
    help="Turn on  filtering polar horn regions from data",
    default=False,
    action="store_true",
)
parser.add_argument(
    "--cormin",
    help="Set minimum cutoff rigidity (COR_SAX) for nimaketime filtering (default=no COR filtering, typical value = 1.5 or 2)",
    default=None,
)
parser.add_argument(
    "--kpmax",
    help="Set maximum KP value for nimaketime filtering (default=no KP filtering, typical value = 5)",
    default=None,
)
parser.add_argument(
    "--minfpm",
    help="Set minimum of FPMs active for nimaketime filtering (default=38)",
    default=38,
)
parser.add_argument(
    "--mingti",
    help="Set minimum duration (in seconds) of a GTI to be kept (default=16)",
    default=16,
)
parser.add_argument(
    "--maxovershoot",
    help="Select data where overshoot rate (per FPM) is below this limit (default: 1.5)",
    type=float,
    default=1.5,
)
parser.add_argument(
    "--badcut",
    help="Select data where bad ratio event rate is below this limit (default: no filter)",
    type=float,
    default=-1,
)
parser.add_argument(
    "--angdist",
    help="Set threshold for ANG_DIST in call to nimaketime (degrees, default=0.015)",
    type=float,
    default=0.015,
)
parser.add_argument(
    "--obsid", help="Use this as OBSID for directory and filenames", default=None
)
parser.add_argument(
    "--dark",
    help="Apply SUNSHINE=0 filter to get only data in Earth shadow",
    action="store_true",
)
parser.add_argument(
    "--maxundershoot",
    help="Select data where undershoot rate is below this limit (default: 200)",
    type=float,
    default=200,
)
parser.add_argument(
    "--medianundershoot",
    help="Select data where MEDIAN undershoot rate is below this limit",
    type=float,
    default=None,
)
parser.add_argument(
    "--nounderfilt",
    help="Don't filter good times based on UNDERONLY rate",
    action="store_true",
)
parser.add_argument(
    "--minsun",
    help="Set minimum sun angle (SUN_ANGLE) for nimaketime filtering (default=no SUN_ANGLE filtering, typical values = 60, 70, 80, 90 deg). Note: Allows dark time at any Sun angle!",
    default=None,
)
parser.add_argument(
    "--day",
    help="Apply SUNSHINE=1 filter to get only data in ISS-day",
    action="store_true",
)
parser.add_argument("--par", help="Par file to use for phases")
parser.add_argument("--ephem", help="Ephem to use with photonphase", default="DE421")
parser.add_argument(
    "--outdir",
    help="Add name to output directories (by default: directories end in '_pipe')",
    default="pipe",
)
parser.add_argument(
    "--merge",
    help="Merge all ObsIDs provided into single event list, lightcurve and spectrum (outputdir called 'merged')",
    action="store_true",
)
parser.add_argument(
    "--crcut",
    help="perform count rate cut on merged event file (only if --merge)",
    action="store_true",
)
parser.add_argument(
    "--lcbinsize",
    help="Lightcurve bin size (sec, default=4.0)",
    type=float,
    default=4.0,
)
parser.add_argument(
    "--filterbinsize",
    help="Bin size for Count rate and Overshoot rate filtering (sec, default=16.0)",
    type=float,
    default=16.0,
)
parser.add_argument(
    "--overonlyexpr",
    help="Turn on the expression to filter on overshoot rate based on COR_SAX",
    action="store_true",
)
parser.add_argument(
    "--tidy",
    help="Make fewer plots and remove some intermediate files, to run faster and save disk space",
    action="store_true",
)
parser.add_argument(
    "--nomap",
    help="Don't make maps (for use if cartopy is not working)",
    action="store_true",
)
parser.add_argument(
    "--copymkf",
    help="Copy ni*.mkf into the _pipe output directory, for convenience. Warning, mkf files are large!",
    action="store_true",
)
parser.add_argument(
    "--log-level",
    type=str,
    choices=("TRACE", "DEBUG", "INFO", "WARNING", "ERROR"),
    default="INFO",
    help="Logging level",
    dest="loglevel",
)

# parser.add_argument("--crabnorm", help="normalize the spectrum with the crab (only if --merge)", action='store_true')
args = parser.parse_args()
log.remove()
log.add(
    sys.stderr,
    level=args.loglevel,
    colorize=True,
    format="<level>{level: <8}</level> ({name: <30}): <level>{message}</level>",
    # filter=pint.logging.LogFilter(),
)


os.environ["HEADASNOQUERY"] = " "
os.environ["HEADASPROMPT"] = "/dev/null"

# Checking the presence of HEASOFT
try:
    check_call("nicerversion", env=os.environ)
except:
    print(
        "You need to initialize FTOOLS/HEASOFT first (e.g., type 'heainit')!",
        file=sys.stderr,
    )
    exit()

# Create a temporary dir for pfiles and set the environment
tempdir = tempfile.mkdtemp()
os.mkdir(tempdir + "/pfiles")
headas = os.environ["HEADAS"]
os.environ["PFILES"] = tempdir + "/pfiles;" + headas + "/syspfiles"
# print(os.environ["PFILES"])
# Called before each exit() and at end of program to clean up:
# shutil.rmtree(tempdir)

# For some reason if the script is called via #!/usr/bin/env python
# it does not inherit LD_LIBRARY_PATH so ftools don't run.
# print(os.environ['LD_LIBRARY_PATH'])
# print(os.environ['HEADAS'])
# os.environ['LD_LIBRARY_PATH'] = path.join(os.environ['HEADAS'],'lib')
# print(os.environ['LD_LIBRARY_PATH'])

all_evfiles = []
all_orbfiles = []


def runcmd(cmd):
    # CMD should be a list of strings since it is not processed by a shell
    log.info("CMD: " + " ".join(cmd))
    # log.info(cmd)
    try:
        check_call(cmd, env=os.environ)
    except CalledProcessError:
        log.error("Command failed! CMD: " + " ".join(cmd))
        raise


# Check if outdir contains 'None', 'NONE', or 'none' (causes bug in ni-extractevents)
if args.outdir:
    names = ["none", "None", "NONE"]
    if any(st in args.outdir for st in names):
        log.error(
            "Due to a current bug in ni-extractevents, --outdir cannot contain 'none', 'None', or 'NONE'.  Exiting..."
        )
        shutil.rmtree(tempdir)
        exit()

# Check if ObsIDs are listed or to be read from file
if len(args.indirs) == 1:
    if args.indirs[0].startswith("@"):
        inputdirs = args.indirs[0].split("@")[1]
        log.info("Reading input ObsID list: {}".format(inputdirs))
        all_obsids = np.loadtxt(inputdirs, dtype=str)
    else:
        all_obsids = args.indirs
else:
    all_obsids = args.indirs

#
# Start processing all ObsIDs
#
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
    pipedir = "{0}_{1}".format(basename, args.outdir)
    if not os.path.exists(pipedir):
        os.makedirs(pipedir)

    log.info("Making initial QL plots")

    # Get filter file
    # Searching for *mkf file
    list_mkfiles = glob(path.join(obsdir, "auxil/ni*.mkf"))
    if len(list_mkfiles) == 0:
        # Searching for *mkf.gz file
        list_mkfiles = glob(path.join(obsdir, "auxil/ni*.mkf.gz"))
        log.info("Using the *.mkf.gz file found instead of a *.mkf file.")
        if len(list_mkfiles) == 0:
            log.error("No *mkf or *mkf.gz files found. Exiting.")
            exit()

    if len(list_mkfiles) == 1:
        mkfile = list_mkfiles[0]
    else:
        # When more than one mkf or mkf.gz file found, using the most recent one
        mkfile = max(list_mkfiles, key=os.path.getmtime)
        log.warning("Multiple MKF files found. Using the most recent one.")

    log.info("MKF File: {0}".format(mkfile))

    # Check if MKF file contains the new columns (try opening one of the new columns)
    try:
        mkftest = fits.open(mkfile)[1].data["FPM_OVERONLY_COUNT"]
    except:
        log.error(
            "New *mkf files needed in {}/auxil/. Please use niprefilter2.".format(
                obsdir
            )
        )
        exit()
    try:
        mkftest = fits.open(mkfile)[1].data["KP"]
        has_KP = True
    except:
        log.warning("No KP column in MKF file. Will not use any KP based filters!")
        has_KP = False

    # Copy MKF file to results dir for pulsar analysis
    if args.copymkf:
        shutil.copy(mkfile, pipedir)

    if args.dark and args.day:
        log.warning("Both --dark and --day are requested")
        args.dark = False
        args.day = False

    cmd = [
        "nicerql.py",
        "--save",
        "--filtall",
        "--lcbinsize",
        "{}".format(args.lcbinsize),
        "--lclog",
        "--useftools",
        "--filterbinsize",
        "{}".format(args.filterbinsize),
        "--emin",
        "{0}".format(args.emin),
        "--emax",
        "{0}".format(args.emax),
        "--sci",
        "--eng",
        "--bkg",
        "--obsdir",
        obsdir,
        "--basename",
        path.join(pipedir, basename) + "_prefilt",
    ]
    if not args.nomap:
        cmd.append("--map")

    if args.mask is not None:
        cmd.append("--mask")
        for detid in args.mask:
            cmd.append("{0}".format(detid))
    if args.badcut > 0:
        cmd.append("--writebkf")

    if not args.tidy:
        runcmd(cmd)

    # Get orbit file
    orbfile = glob(path.join(obsdir, "auxil/ni*.orb*"))[0]
    if len(glob(path.join(obsdir, "auxil/ni*.orb*"))) > 1:
        log.info("Multiple orbit files found :")
        for f in glob(path.join(obsdir, "auxil/ni*.orb*")):
            log.info("  -> {}".format(path.basename(f)))
    log.info("Orbit File: {0}".format(orbfile))
    # Copy orbit file to results dir for pulsar analysis
    shutil.copy(orbfile, pipedir)

    #  Get BKF file for filtering based on background indicators
    if args.badcut > 0:
        bkffile = path.join(pipedir, basename) + "_prefilt.bkf"
        log.info("BKF File: {0}".format(bkffile))
    else:
        bkffile = None

    # Get merged unfiltered event filename (should just be one)
    evfiles = glob(path.join(obsdir, "xti/event_cl/ni*mpu7_cl.evt"))
    evfiles.sort()
    log.info("Cleaned Event Files: {0}".format("\n" + "    \n".join(evfiles)))

    # Build input file for niextract-events
    # I don't think there should ever be more than one. Removing this code.
    # evlistname = path.join(pipedir, "evfiles.txt")
    # fout = open(evlistname, "w")
    # for en in evfiles:
    #     print("{0}".format(en), file=fout)
    # fout.close()
    # Just take first one!
    evfilename = evfiles[0]

    intermediatename = path.join(pipedir, "intermediate.evt")
    filteredname = path.join(pipedir, "cleanfilt.evt")

    # First filter any bad detectors.
    # Start with launch detectors
    detfilt_expr = "launch"

    # Filter any explicitly specified masked detectors
    if args.mask is not None:
        for detid in args.mask:
            if detid >= 0:
                detfilt_expr += ",-{0}".format(detid)

    if detfilt_expr == "launch":
        cmd = ["cp", evfilename, intermediatename]
    else:
        cmd = [
            "nifpmsel",
            "infile={}".format(evfilename),
            "outfile={}".format(intermediatename),
            "detlist={}".format(detfilt_expr),
            "mkfile={}".format(mkfile),
            "clobber=yes",
            "history=yes",
        ]
    runcmd(cmd)

    # Create any additional GTIs beyond what nimaketime does...
    extragtis = "NONE"
    if args.filtpolar:
        saafile = path.join(datadir, "saa.reg")
        mkf_expr = 'regfilter("{0}",SAT_LON,SAT_LAT)'.format(saafile)
        phfile = path.join(datadir, "polarhorns.reg")
        mkf_expr += '.and.regfilter("{0}",SAT_LON,SAT_LAT)'.format(phfile)
        gtiname2 = path.join(pipedir, "extra.gti")
        cmd = ["pset", "maketime", "expr={0}".format(mkf_expr)]
        runcmd(cmd)
        cmd = [
            "maketime",
            "infile={0}".format(mkfile),
            "outfile={0}".format(gtiname2),
            "compact=no",
            "time=TIME",
            "prefr=0",
            "postfr=0",
            "clobber=yes",
        ]
        runcmd(cmd)
        if len(Table.read(gtiname2, hdu=1)) == 0:
            log.error("No good time left after filtering!")
            continue
        extragtis = gtiname2

    gtiname3 = None

    # Create GTI from overshoot file using bad event lightcurve
    if args.badcut > 0:
        gtiname3 = path.join(pipedir, "bkf.gti")
        bkf_expr = "BAD_RATIO.lt.{0}".format(args.badcut)
        cmd = [
            "maketime",
            bkffile,
            gtiname3,
            "expr={0}".format(bkf_expr),
            "compact=no",
            "time=TIME",
            "prefr=0",
            "postfr=0",
            "clobber=yes",
        ]
        runcmd(cmd)
        if len(Table.read(gtiname3, hdu=1)) == 0:
            log.error("No good time left after filtering!")
            continue

    # If either of the bkf filters were used, include that GTI
    # in the extragtis passed to nimaketime
    if gtiname3 is not None:
        if extragtis == "NONE":
            extragtis = gtiname3
        else:
            extragtis = extragtis + ",{0}".format(gtiname3)

    # Make final merged GTI using nimaketime
    gtiname_merged = path.join(pipedir, "tot.gti")
    gtiname_merged_and_eventgti = path.join(pipedir, "tot_and_eventgti.gti")
    gtiname_clipped = path.join(pipedir, "tot_clipped.gti")

    # Manage extra_expr for nimaketime (ST_VALID, DARK/DAY, and FPM_OVER_ONLY filters from --KEITH)
    list_extra_expr = ["ST_VALID.eq.1"]

    if args.dark:
        list_extra_expr.append("SUNSHINE.eq.0")
    if args.day:
        list_extra_expr.append("SUNSHINE.eq.1")

    if args.kpmax and has_KP:
        list_extra_expr.append("KP.lt.{0}".format(args.kpmax))

    if args.minsun:
        # Exclude data that is at a Sun angle less that some value, unless it is in eclipse
        list_extra_expr.append(
            "(SUN_ANGLE.gt.{0}.or.SUNSHINE.eq.0)".format(args.minsun)
        )

    if args.overonlyexpr:
        list_extra_expr.append("FPM_OVERONLY_COUNT<1.52*COR_SAX**(-0.633)")
    if args.medianundershoot:
        # Exclude data when the MEDIAN undershoot is above this level
        list_extra_expr.append(f"(MEDIAN_UNDERONLY_COUNT.lt.{args.medianundershoot})")

    extra_expr = "(" + " && ".join("%s" % expr for expr in list_extra_expr) + ")"

    cor_string = "-"
    if args.cormin is not None:
        cor_string = "{0}-".format(args.cormin)
    ## Default is now 20 and 30 since that is standard processing as of 2019 May.
    elvcut = 20.0
    brcut = 30.0
    ## Might want to change this.
    maxunder = args.maxundershoot
    if args.nounderfilt:
        maxunder = 650.0
    cmd = [
        "nimaketime",
        "infile={0}".format(mkfile),
        "outfile={0}".format(gtiname_merged),
        "nicersaafilt=YES",
        "saafilt=NO",
        "trackfilt=YES",
        "ang_dist={0:.3f}".format(args.angdist),
        "elv={0}".format(elvcut),
        "br_earth={0}".format(brcut),
        "cor_range={0}".format(cor_string),
        "min_fpm={0}".format(args.minfpm),
        "underonly_range=0-{0}".format(maxunder),
        "overonly_range=0-{0}".format(args.maxovershoot),
        "ingtis={0}".format(extragtis),
        "clobber=yes",
        "expr={0}".format(extra_expr),
        "outexprfile={0}".format(path.join(pipedir, "psrpipe_expr.txt")),
    ]
    runcmd(cmd)

    # Make a GTI file that is the AND of gtiname_merged and the intermediate file GTI
    # ftmgtime will overwrite a file if the file exists and if clobber=YES,
    # For now, however, delete the old file if it exists and don't use clobber=YES.
    if path.isfile(gtiname_merged_and_eventgti):
        os.remove(gtiname_merged_and_eventgti)
    cmd = [
        "ftmgtime",
        f"{intermediatename}[GTI],{gtiname_merged}",
        f"outgti={gtiname_merged_and_eventgti}",
        "merge=AND",
    ]
    runcmd(cmd)

    # nimaketime gives an output GTI file with 1 row and START==STOP in the
    # case where no good time is selected.  This differs from the normal
    # maketime, which produces a GTI file with no rows in that case

    # Now run ftadjustgti to discard short GTIs
    cmd = [
        "ftadjustgti",
        "infile={0}".format(gtiname_merged_and_eventgti),
        "outfile={0}".format(gtiname_clipped),
        "mingti={}".format(args.mingti),
        "copyall=NO",
        "clobber=YES",
    ]
    runcmd(cmd)

    # Build selection expression for niextract-events
    # Select events with PI in the selected range, require SLOW trigger (FAST optional)
    # and filter all undershoot, overshoot, and force trigger events
    evfilt_expr = "PI={0}:{1},EVENT_FLAGS=bx1x000".format(
        int(args.emin * KEV_TO_PI), int(args.emax * KEV_TO_PI)
    )

    cmd = [
        "niextract-events",
        "filename={0}[{1}]".format(intermediatename, evfilt_expr),
        "eventsout={0}".format(filteredname),
        "timefile={0}".format(gtiname_clipped),
        "gti=GTI",
        "clobber=yes",
    ]
    runcmd(cmd)

    # Remove intermediate file
    if args.tidy:
        os.remove(intermediatename)

    # Check that there are events left
    if len(Table.read(filteredname, hdu=1)) == 0:
        log.error("No events left !!!")
        continue

    # Make cleanfilt.mkf file from ObsID .mkf file and merged_GTI
    cleanfilt_mkf = path.join(pipedir, "cleanfilt.mkf")
    log.info("Applying the GTI filtering to the *mkf file")
    cmd = [
        "fltime",
        "infile={}[1]".format(mkfile),
        "gtifile={}".format(gtiname_merged),
        "outfile={}".format(cleanfilt_mkf),
        "clobber=yes",
    ]
    runcmd(cmd)

    # Make final clean plot
    cmd = [
        "nicerql.py",
        "--save",
        "--orb",
        orbfile,
        #        path.join(pipedir, path.basename(orbfile)),
        "--sci",
        "--eng",
        filteredname,
        # "--allspec",
        # "--alllc",
        "--lcbinsize",
        "{}".format(args.lcbinsize),
        "--filterbinsize",
        "{}".format(args.filterbinsize),
        "--mkf",
        cleanfilt_mkf,
        "--bkg",
        "--basename",
        path.join(pipedir, basename) + "_cleanfilt",
    ]
    if not args.nomap:
        cmd.append("--map")
    if args.par is not None:
        cmd.append("--par")
        cmd.append("{0}".format(args.par))
    runcmd(cmd)

    # Add phases
    if args.par is not None:
        plotfile = path.join(pipedir, "phaseogram.png")
        log.info("Applying phases to {0}".format(filteredname))
        log.info("Event file has TIMZERO < 0, so not applying --fix in photonphase!")
        cmd = [
            "photonphase",
            "--ephem",
            args.ephem,
            "--orb",
            orbfile,
            "--plot",
            "--plotfile",
            plotfile,
            "--addphase",
            filteredname,
            args.par,
        ]
        runcmd(cmd)

    # Extract simple PHA file and light curve
    phafile = path.splitext(filteredname)[0] + ".pha"
    lcfile = path.splitext(filteredname)[0] + ".lc"
    cmd = [
        "extractor",
        filteredname,
        "eventsout=none",
        "imgfile=none",
        "phafile={0}".format(phafile),
        "fitsbinlc={0}".format(lcfile),
        "binlc={}".format(args.lcbinsize),
        "regionfile=none",
        "timefile=none",
        "xcolf=RAWX",
        "ycolf=RAWY",
        "tcol=TIME",
        "ecol=PI",
        "xcolh=RAWX",
        "ycolh=RAWY",
        "gti=GTI",
    ]
    runcmd(cmd)

    # if --merge option, Add clean evt file to list of files to merge
    if args.merge:
        all_evfiles.append(filteredname)
        all_orbfiles.append(orbfile)


# Merging all ObsIDs
if args.merge and (len(all_evfiles) > 1):
    ## save list of orbit files
    orbfiles_list = path.join(os.getcwd(), "list_orbfiles.txt")
    np.savetxt(orbfiles_list, all_orbfiles, fmt=["%s"])

    ## Call merge.py
    cmd = ["merge.py"]
    for evt in all_evfiles:
        cmd.append(evt)
    cmd.append("merged")
    cmd.append("merged_{0}".format(args.outdir))
    cmd.append("--clobber")
    cmd.append("--lcbinsize")
    cmd.append("{}".format(args.lcbinsize))

    if args.par is not None:
        cmd.append("--par")
        cmd.append("{0}".format(args.par))
        cmd.append("--orb")
        cmd.append("{0}".format(orbfiles_list))
    if args.crcut:
        cmd.append("--cut")
        cmd.append("--filterbinsize")
        cmd.append("{}".format(args.filterbinsize))

    # if args.crabnorm:
    #     cmd.append("--crabnorm")
    # if args.dark:
    #     cmd.append("--dark")
    runcmd(cmd)

else:
    if args.crcut:
        log.warning(
            "Count rate cuts are only performed on merged event files (add the option --merge)"
        )

    # # Perform Crab Normalization (in the case of single ObsID bring processed)
    # --- NOT IMPLEMENTED --- #
    # if args.crabnorm:
    #     log.info("Performing normalization of the cleanfilt spectrum with crab residuals")
    #     normed_spec = path.join(phafile.strip('.pha'),"_crabnorm.pha")
    #     if args.dark:
    #         cmd = ["mathpha","{}/{}".format(phafile,CRAB_RES_NIGHT),"R","{}".format(normed_spec),"{} % POISS-0 0".format(phafile)]
    #         runcmd(cmd)
    #     else:
    #         cmd = ["mathpha","{}/{}".format(phafile,CRAB_RES_TOT),"R","{}".format(normed_spec),"{} % POISS-0 0".format(phafile)]
    #         runcmd(cmd)

shutil.rmtree(tempdir)
