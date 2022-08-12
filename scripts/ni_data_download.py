#!/usr/bin/env python

import os
import argparse
import numpy as np
import pandas as pd
from fuzzywuzzy import process
from astropy import log
from getpass import getpass
from datetime import datetime
import sys
import subprocess
import shutil
import fnmatch


# ----------------------------------------------------------------------
#   command line interface
# ----------------------------------------------------------------------

desc = """
NICER data downloader. Script written by P. Bult. Minor tweaks by S. Guillot, to make it more automatic, and removing multiproccessing (it was causing problems).
"""
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("sourcename", help="Give source name", default=None)
parser.add_argument("heasarc_user", help="Enter heasarc username", default=None)
parser.add_argument("heasarc_pwd", help="Enter heasarc password", default=None)
parser.add_argument(
    "--outdir",
    help="Directory to place all ObsIDs (relative to current)",
    type=str,
    default=None,
)
parser.add_argument("--clobber", help="Clobber existing files", action="store_true")
parser.add_argument(
    "--cl-only",
    help="Grab the clean event files only",
    action="store_true",
    dest="cl_only",
)
parser.add_argument(
    "--mkf-only",
    help="Grab the prefilter data files only",
    action="store_true",
    dest="mkf_only",
)
parser.add_argument(
    "--decryptkey",
    help="Add decryption key to run GPG after download",
    type=str,
    default=None,
)
parser.add_argument("--unzip", help="Gunzip after decrypting", action="store_true")
args = parser.parse_args()

# ----------------------------------------------------------------------
# Checking the presence of GPG or GPG2
# ----------------------------------------------------------------------
gpg_version = "gpg"
if args.decryptkey:
    try:
        cmd = ["{}".format(gpg_version), "--version"]
        subprocess.call(cmd)
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            log.warning("gpg not found...Trying gpg2")
            gpg_version = "gpg2"
            cmd = ["{}".format(gpg_version), "--version"]
            try:
                subprocess.call(cmd)
            except OSError as e:
                if e.errno == os.errno.ENOENT:
                    log.warning("gpg2 not found...No decryption will be performed")
                    exit()
                else:
                    log.error("Something went wrong while trying to run gpg2")
                    raise
        else:
            log.error("Something went wrong while trying to run gpg")
            raise
log.info("{} will be used".format(gpg_version))

# ----------------------------------------------------------------------
#   HEASARC scraper
# ----------------------------------------------------------------------


def print_nicer_segment(
    url="https://heasarc.gsfc.nasa.gov/docs/nicer/team_schedule/nicer_seg_team.html",
    username=None,
    password=None,
):
    """
    This prints out the segment detail table in text format

    usage: % print_nicer_segment(username = "nicer_user_name" password = "nicer password")

    outputs: prints the nicer segment table to the terminal

    :param url: location of the segment detail page
    :param username: nicer team username
    :param password: nicer team password
    :return:
    """
    from bs4 import BeautifulSoup
    import requests

    if (not username) or (not password):
        raise ValueError(
            "must supply username and password to access the NICER obs page"
        )
    req = requests.get(url, auth=(username, password))
    if req.status_code != 200:
        raise ValueError(
            "Problem accessing {0} with ({1}, {2}) \nReturn code: {3}".format(
                url, username, password, req.status_code
            )
        )

    soup = BeautifulSoup(req.text, "lxml")
    tabs = soup.find_all("table")[1]
    df = pd.read_html(str(tabs))
    return df[0]


# ----------------------------------------------------------------------
#   file handling
# ----------------------------------------------------------------------
def untar(obsid, outdir):
    cmd = "tar -xf {}.tar --directory {}".format(obsid, outdir)
    os.system(cmd)


def decrypt(passwd, obsid_dir, gpg_v):
    cmd = "find {} -name '*.gpg' -exec {} --batch --yes --passphrase '{}'  {{}} \; -delete".format(
        obsid_dir, gpg_v, passwd
    )
    try:
        os.remove("{}/ni-gpg-call.log".format(obsid_dir))
    except OSError:
        pass
    os.system(cmd + " 2>> {}/ni-gpg-call.log".format(obsid_dir))


def unzip(obsid_dir):
    cmd = "find {} -name '*.gz' -exec gunzip {{}} \;".format(obsid_dir)
    os.system(cmd)


def check_files(teststring, checkdir):
    filesfound = []
    for root, dirnames, filenames in os.walk(checkdir):
        for filename in fnmatch.filter(filenames, teststring):
            filesfound.append(os.path.join(root, filename))
    return filesfound


# ----------------------------------------------------------------------
#   main routine
# ----------------------------------------------------------------------

if args.unzip:
    if not args.decryptkey:
        log.warning("You need to set --decryptkey if you wish to unzip files. Exiting")
        print("")
        parser.print_help()
        exit()

# Set up work directory
if args.outdir:
    workingdir = os.path.join(os.getcwd(), args.outdir)
    if not os.path.exists(workingdir):
        log.info("Making data directory: {}".format(workingdir))
        os.makedirs(workingdir)
    else:
        log.info("Data will be placed in: {}".format(workingdir))

else:
    workingdir = os.getcwd()
    log.info("Data will be placed in current directory: {}".format(workingdir))


# Making heasarc query
log.info(" --| QUERYING HEASARC FOR {0} |-- ".format(args.sourcename))
username = args.heasarc_user  ###input  ("username: ")
password = args.heasarc_pwd  ###getpass("password: ")
df = print_nicer_segment(username=username, password=password)

# Filter on source name
source = df[df["Target Name"] == args.sourcename]

# Ensure data selection
if len(source) == 0:
    # Source not found, do fuzzy logic search
    log.info("Source not found. Using fuzzy logic search...")
    choices = np.unique(df["Target Name"])
    matches = process.extract(args.sourcename, choices, limit=3)

    # Ask for new input
    log.info("Did you mean:")
    for n, match in enumerate(matches):
        print("  [{}] : {}".format(n, match[0]))
    print("  [{}] : {}".format(3, "Quit"))

    answer = int(input("(0-3) : "))
    if answer == 3:
        log.info("You picked: 'Quit'. Quitting now...")
        exit()
    else:
        log.info("You picked: {}".format(matches[answer][0]))
        source = df[df["Target Name"] == matches[answer][0]]

# ObsID output
log.info("The following ObsIDs were found:")
log.info("  # ObsIDs: {}".format(len(source)))
log.info(" Good expo: {} s".format(np.sum(source["Good Expo[s]"])))
log.info(source[["Target Name", "Start TimeUTC", "Observation ID", "Good Expo[s]"]])

# Grab the relevant columns
cols = source[["Target Name", "Start TimeUTC", "Observation ID", "Good Expo[s]"]]

# Print ObsID logs to file
log.info("Writing the ObsID log file: ObsLogs.txt")
ObsIDsLogFile = open(os.path.join(workingdir, "ObsLogs.txt"), "w")
ObsIDsLogFile.write(" ".join(str(x) for x in sys.argv))
ObsIDsLogFile.write("\n\n # ObsIDs: {}".format(len(source)))
ObsIDsLogFile.write("\n Good expo: {} s".format(np.sum(source["Good Expo[s]"])))
ObsIDsLogFile.write("\n\n")
ObsIDsLogFile.write(cols.to_csv(sep="\t", index=False, header=False))
ObsIDsLogFile.close()

log.info("Downloading from heasarc and file handling begins...")

for n, [no, row] in enumerate(source.iterrows()):
    # Extract the relevant columns
    obsid = "{:010d}".format(row["Observation ID"])
    gexpo = row["Good Expo[s]"]
    try:
        dtime = datetime.strptime(row["Start TimeUTC"], "%Y-%m-%dT%H:%M:%S")
    except:
        log.warning("Time is NaN...Ignoring ObsID {}".format(obsid))
        continue

    # Set the archive
    base = "ftp://heasarc.gsfc.nasa.gov/nicer/data/.obs2validate/"  ## TEAM ARCHIVE
    # base = "ftp://heasarc.gsfc.nasa.gov/nicer/data/obs/"          ## PUBLIC ARCHIVE

    # Select the date folder
    base += "{}_{:02d}/".format(dtime.year, dtime.month)

    # Select the files to download
    if args.cl_only:
        target = "{0}/xti/event_cl/ni{0}_0mpu7_cl.evt.gz.gpg".format(obsid)
    elif args.mkf_only:
        target = "{0}/auxil/ni{0}.mkf.gz.gpg".format(obsid)
    else:
        target = obsid + ".tar"

    # Handle clobbering
    # First check if files are present in workingdir
    matchedfiles = check_files("*{}*".format(obsid), os.path.join(workingdir, obsid))

    if not matchedfiles or args.clobber:
        # Build the download command
        cmd = "curl --retry 10 {0}{1} --ftp-ssl -k --create-dirs -o {1}".format(
            base, target
        )
        log.info(
            "{:3d} / {:3d} :: Trying: curl obsid {}.tar".format(
                n + 1, len(source), obsid
            )
        )
        os.system(cmd)

        # File management
        log.info(
            "{:3d} / {:3d} :: Un-archiving obsid {}.tar".format(
                n + 1, len(source), obsid
            )
        )
        untar(obsid, workingdir)
        try:
            os.remove("{}.tar".format(obsid))
        except:
            log.warning(
                "{:3d} / {:3d} :: Cannot remove {}.tar. File not found".format(
                    n + 1, len(source), obsid
                )
            )
            continue

        if args.decryptkey:
            decryptOK = True
            log.info(
                "{:3d} / {:3d} :: De-crypting ObsID {}".format(
                    n + 1, len(source), obsid
                )
            )
            decrypt(args.decryptkey, os.path.join(workingdir, obsid), gpg_version)
            if ("decryption failed") in open(
                "{}/ni-gpg-call.log".format(os.path.join(workingdir, obsid))
            ).read():
                log.warning("A bad key was provided. Decrypting was not performed!")
                if args.unzip:
                    log.warning("A bad key was provided. Skipping unzip!")
                decryptOK = False
        else:
            decryptOK = False

        if args.unzip and decryptOK:
            log.info(
                "{:3d} / {:3d} :: unzipping all files of #{}".format(
                    n + 1, len(source), obsid
                )
            )
            unzip(os.path.join(workingdir, obsid))
    else:
        log.info(
            "{:3d} / {:3d} :: Obsid {} has already been downloaded".format(
                n + 1, len(source), obsid
            )
        )

        if args.decryptkey:
            decryptOK = True
            # Check there are files that have --not-- been decrypted
            matchedfiles = check_files("*gpg".format(obsid), workingdir)
            if matchedfiles:
                log.info(
                    "{:3d} / {:3d} :: De-crypting ObsID {}".format(
                        n + 1, len(source), obsid
                    )
                )
                decrypt(args.decryptkey, os.path.join(workingdir, obsid), gpg_version)
                if (
                    "decryption failed"
                    in open(
                        "{}/ni-gpg-call.log".format(os.path.join(workingdir, obsid))
                    ).read()
                ):
                    log.warning("A bad key was provided. Decrypting was not performed!")
                    if args.unzip:
                        log.warning("A bad key was provided. Skipping unzip!")
                    decryptOK = False
        else:
            decryptOK = False

        if args.unzip and decryptOK:
            matchedfiles = check_files("*gz".format(obsid), workingdir)
            if matchedfiles:
                log.info(
                    "{:3d} / {:3d} :: unzipping all files of #{}".format(
                        n + 1, len(source), obsid
                    )
                )
                unzip(os.path.join(workingdir, obsid))

log.info(" --| DONE DOWNLOADING DATA FOR {0}|-- ".format(args.sourcename))
print("")

if (args.decryptkey and not decryptOK) or not args.decryptkey:
    log.warning(
        "Data files not decrypted. The following command can decrypt all files:"
    )
    print(
        " find {} -name '*.gpg' -exec gpg --batch --yes --passphrase 'your_passphrase'  {{}} \; -delete".format(
            workingdir
        )
    )

if (args.unzip and not decryptOK) or not args.unzip:
    log.warning(
        "Data files not unzipped. The following command can unzip all (decrypted) files:"
    )
    print(" find {} -name '*.gz' -exec gunzip {{}} \;".format(workingdir))
