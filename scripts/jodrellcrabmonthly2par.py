#!/usr/bin/env python
# Code to write out Tempo2-format PAR files based on Jodrell Crab monthly ephemeris parameters
# Update the text file with the monthly ephemeris:
# % curl -O https://www.jb.man.ac.uk/pulsar/crab/crab2.txt
import numpy as np
import astropy.units as u
import argparse
from astropy.table import Table
from astropy import log
import sys
import os.path as path
import urllib.request
import urllib.error

desc = "Write out Tempo2 par files based on Jodrell montly ephemeris"
CRAB2_URL = "https://www.jb.man.ac.uk/pulsar/crab/crab2.txt"


def download_crab2(dest="crab2.txt", prompt=True):
    """Download crab2.txt from Jodrell Bank, optionally prompting the user."""
    if prompt:
        try:
            resp = input(
                f"{dest} does not exist. Download from {CRAB2_URL}? [Y/n]: "
            )
        except (EOFError, KeyboardInterrupt):
            print()
            return False
        if resp.strip().lower() not in ("", "y", "yes"):
            return False

    log.info(f"Downloading {CRAB2_URL} to {dest}...")
    try:
        urllib.request.urlretrieve(CRAB2_URL, dest)
        log.info("Download complete.")
        return True
    except urllib.error.URLError as e:
        log.error(f"Failed to download {CRAB2_URL}: {e}")
        return False
    except Exception as e:
        log.error(f"Unexpected error downloading {CRAB2_URL}: {e}")
        return False


def get_parser():
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("MJD", help="MJD to extract ephemeris for", type=float)
    parser.add_argument(
        "--auto-download",
        "-y",
        action="store_true",
        help="Automatically download crab2.txt without prompting if not found",
    )
    return parser


def main(argv=None):
    ## Parse arguments
    args = get_parser().parse_args(argv)

    if not path.exists("crab2.txt"):
        success = download_crab2("crab2.txt", prompt=not args.auto_download)
        if not success or not path.exists("crab2.txt"):
            log.error(
                f'crab2.txt does not exist!  Download with "curl -O {CRAB2_URL}"'
            )
            sys.exit(2)

    validmonths = [
        "JAN",
        "FEB",
        "MAR",
        "APR",
        "MAY",
        "JUN",
        "JUL",
        "JLY",
        "AUG",
        "SEP",
        "OCT",
        "NOV",
        "DEC",
    ]

    day = []
    mon = []
    yr = []
    mjd = []
    tJPL = []
    tacc = []
    f0 = []
    f1 = []
    dm = []
    dm1 = []
    for line in open("crab2.txt"):
        line = line.strip()
        if len(line) == 0:
            continue
        line = (line.replace("(", " ")).replace(")", " ")
        cols = line.split()
        # Make sure second column is a valid month, or skip line. This skips headers, etc.
        if not cols[1] in validmonths:
            continue
        # Skip entries before 55910 when DM1 column was added.
        if float(cols[3]) < 55910:
            continue
        # print(line)
        day.append(int(cols[0]))
        mon.append(cols[1])
        yr.append(int(cols[2]))
        mjd.append(float(cols[3]))
        tJPL.append(float(cols[4]))
        tacc.append(float(cols[5]))
        f0.append(float(cols[6]))
        f1.append(float(cols[8]) * 1.0e-15)
        dm.append(float(cols[10]))
        dm1.append(float(cols[11]))

    jtab = Table(
        [day, mon, yr, mjd, tJPL, tacc, f0, f1, dm, dm1],
        names=("day", "month", "year", "mjd", "tJPL", "tacc", "f0", "f1", "dm", "dm1"),
        dtype=("i", "S3", "i4", "f8", "f8", "f8", "f8", "f8", "f8", "f8"),
    )

    diff = np.abs(jtab["mjd"] - args.MJD)
    idx = diff.argmin()

    if diff[idx] > 18:
        log.warning(
            "Too large a difference ({} days) to closest ephemeris entry!".format(
                diff[idx]
            )
        )
        # sys.exit(1)

    f0 = jtab[idx]["f0"]
    f1 = jtab[idx]["f1"]
    p0 = 1.0 / f0
    p1 = -f1 / (f0 * f0)
    # This calculation of F2 comes from the C code at the end of the
    # explanatory notes for the Jodrell ephemeris
    f2 = 2.0 * p1 * p1 / (p0 * p0 * p0)

    skel = """
PSR B0531+21
RAJ 05:34:31.97232
DECJ +22:00:52.069
EPHEM DE200
F0 {0:.15f}
F1 {1:.15g}
C F2 computed as F2=2.0*P1*P1/(P0*P0*P0) according to JBO ephemeris
F2 {2:.15g}
PEPOCH {3:.15f}
TZRSITE @
TZRFRQ 0.0
TZRMJD {4:.15f}
START {5:.4f}
FINISH {6:.4f}
UNITS TDB
CLOCK TT(TAI)
DM {7:.5f}
DM1 {8:.5g}
DMEPOCH {9:.15f}
TRES {10:.3f}
"""

    # Note 5: The observed barycentric frequency and its first derivative
    # are quoted at the arrival time, using the DE200 ephemeris
    # So we use this value for both the period epoch and for TZRMJD
    # This is in TDB units
    mjd = jtab[idx]["mjd"]
    pepoch = mjd + jtab[idx]["tJPL"] / 86400.0

    print(
        skel.format(
            f0,
            f1,
            f2,
            pepoch,
            pepoch,
            mjd - 16,
            mjd + 16,
            jtab[idx]["dm"],
            jtab[idx]["dm1"],
            pepoch,
            jtab[idx]["tacc"],
        )
    )


if __name__ == "__main__":
    main()
