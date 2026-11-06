"""Helpers that run FTOOLS to merge and filter event files.

All of these need HEASoft to be installed and initialized.
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import argparse
from loguru import logger as log
from os import path
from glob import glob
from subprocess import check_call
import shutil
import tempfile
from astropy.table import Table

from nicer.values import *


def runcmd(cmd):
    """Run a command, logging it first.

    Parameters
    ----------
    cmd : list of str
        The command and its arguments. It is not passed through a shell.

    Raises
    ------
    OSError
        If the command could not be run (for example, if FTOOLS is not
        initialized).
    subprocess.CalledProcessError
        If the command returns a non-zero exit status.
    """
    # CMD should be a list of strings since it is not processed by a shell
    log.info("CMD: " + " ".join(cmd))
    try:
        check_call(cmd, env=os.environ)
    except OSError:
        log.error("Command failed! Did you install and initialize FTOOLS?")
        raise


def filtallandmerge_ftools(evfiles, workdir=None):
    """Merge event files with ``ftmerge``, keeping only standard science events.

    Events are selected with ``EVENT_FLAGS==bx1x000``.

    Parameters
    ----------
    evfiles : list of str
        Event files to merge.
    workdir : str, optional
        Directory in which to make the temporary working directory.

    Returns
    -------
    astropy.table.Table
        The selected events, with ``TIMEZERO`` added to ``TIME``.
    """

    tmpdir = tempfile.mkdtemp(dir=workdir)

    # Build input file for ftmerge
    evlistname = path.join(tmpdir, "evfiles.txt")
    fout = open(evlistname, "w")
    # evfilt_expr = '(PI>20).and.(EVENT_FLAGS==bx1x000)'
    evfilt_expr = "(EVENT_FLAGS==bx1x000)"
    for en in evfiles:
        print("{0}[{1}]".format(en, evfilt_expr), file=fout)
    fout.close()

    # Run ftmerge
    mergedname = path.join(tmpdir, "merged.evt")
    cmd = [
        "ftmerge",
        "@{0}".format(evlistname),
        "outfile={0}".format(mergedname),
        "clobber=yes",
    ]
    runcmd(cmd)

    # Read merged event FITS file into a Table
    etable = Table.read(mergedname, hdu=1)
    # Apply TIMEZERO if needed
    if "TIMEZERO" in etable.meta:
        log.info(
            "Applying TIMEZERO of {0} to etable in filtallandmerge".format(
                etable.meta["TIMEZERO"]
            )
        )
        etable["TIME"] += etable.meta["TIMEZERO"]

    # Clean up
    os.remove(evlistname)
    os.remove(mergedname)
    os.rmdir(tmpdir)

    return etable


def get_eventovershoots_ftools(evfiles, workdir=None):
    """Merge event files with ``ftmerge``, keeping only overshoot-only events.

    Events are selected with ``EVENT_FLAGS==bxxx010``.

    Parameters
    ----------
    evfiles : list of str
        Event files to merge.
    workdir : str, optional
        Directory in which to make the temporary working directory.

    Returns
    -------
    astropy.table.Table
        The selected events, with ``TIMEZERO`` added to ``TIME``.
    """

    tmpdir = tempfile.mkdtemp(dir=workdir)

    # Build input file for ftmerge
    evlistname = path.join(tmpdir, "evfiles.txt")
    fout = open(evlistname, "w")
    evfilt_expr = "(EVENT_FLAGS==bxxx010)"
    for en in evfiles:
        print("{0}[{1}]".format(en, evfilt_expr), file=fout)
    fout.close()

    # Run ftmerge
    mergedname = path.join(tmpdir, "overshoots.evt")
    cmd = [
        "ftmerge",
        "@{0}".format(evlistname),
        "outfile={0}".format(mergedname),
        "clobber=yes",
    ]
    runcmd(cmd)

    # Read merged event FITS file into a Table
    etable = Table.read(mergedname, hdu=1)
    if "TIMEZERO" in etable.meta:
        log.info(
            "Applying TIMEZERO of {0} to etable in get_eventovershoots".format(
                etable.meta["TIMEZERO"]
            )
        )
        etable["TIME"] += etable.meta["TIMEZERO"]

    # Clean up
    os.remove(evlistname)
    os.remove(mergedname)
    os.rmdir(tmpdir)

    return etable


def get_eventundershoots_ftools(evfiles, workdir=None):
    """Merge event files with ``ftmerge``, keeping only undershoot-only events.

    Events are selected with ``EVENT_FLAGS==bxxx001``.

    Parameters
    ----------
    evfiles : list of str
        Event files to merge.
    workdir : str, optional
        Directory in which to make the temporary working directory.

    Returns
    -------
    astropy.table.Table
        The selected events, with ``TIMEZERO`` added to ``TIME``.
    """

    tmpdir = tempfile.mkdtemp(dir=workdir)

    # Build input file for ftmerge
    evlistname = path.join(tmpdir, "evfiles.txt")
    fout = open(evlistname, "w")
    evfilt_expr = "(EVENT_FLAGS==bxxx001)"
    for en in evfiles:
        print("{0}[{1}]".format(en, evfilt_expr), file=fout)
    fout.close()

    # Run ftmerge
    mergedname = path.join(tmpdir, "undershoots.evt")
    cmd = [
        "ftmerge",
        "@{0}".format(evlistname),
        "outfile={0}".format(mergedname),
        "clobber=yes",
    ]
    runcmd(cmd)

    # Read merged event FITS file into a Table
    etable = Table.read(mergedname, hdu=1)
    if "TIMEZERO" in etable.meta:
        log.info(
            "Applying TIMEZERO of {0} to etable in get_eventundershoots".format(
                etable.meta["TIMEZERO"]
            )
        )
        etable["TIME"] += etable.meta["TIMEZERO"]

    # Clean up
    os.remove(evlistname)
    os.remove(mergedname)
    os.rmdir(tmpdir)

    return etable


def get_eventbothshoots_ftools(evfiles, workdir=None):
    """Merge event files with ``ftmerge``, keeping only events flagged as both overshoot and undershoot.

    Events are selected with ``EVENT_FLAGS==bxxx011``.

    Parameters
    ----------
    evfiles : list of str
        Event files to merge.
    workdir : str, optional
        Directory in which to make the temporary working directory.

    Returns
    -------
    astropy.table.Table
        The selected events, with ``TIMEZERO`` added to ``TIME``.
    """

    tmpdir = tempfile.mkdtemp(dir=workdir)

    # Build input file for ftmerge
    evlistname = path.join(tmpdir, "evfiles.txt")
    fout = open(evlistname, "w")
    evfilt_expr = "(EVENT_FLAGS==bxxx011)"
    for en in evfiles:
        print("{0}[{1}]".format(en, evfilt_expr), file=fout)
    fout.close()

    # Run ftmerge
    mergedname = path.join(tmpdir, "bothshoots.evt")
    cmd = [
        "ftmerge",
        "@{0}".format(evlistname),
        "outfile={0}".format(mergedname),
        "clobber=yes",
    ]
    runcmd(cmd)

    # Read merged event FITS file into a Table
    etable = Table.read(mergedname, hdu=1)
    if "TIMEZERO" in etable.meta:
        log.info(
            "Applying TIMEZERO of {0} to etable in get_eventbothshoots".format(
                etable.meta["TIMEZERO"]
            )
        )
        etable["TIME"] += etable.meta["TIMEZERO"]

    # Clean up
    os.remove(evlistname)
    os.remove(mergedname)
    os.rmdir(tmpdir)

    return etable


def get_badratioevents_ftools(evfiles, workdir=None):
    """Merge event files with ``ftmerge``, keeping only events above the slow/fast ratio "trumpet".

    Events are selected with ``EVENT_FLAGS==bx11000`` and
    ``PI/PI_FAST > 1.1 + 120/PI``.

    Parameters
    ----------
    evfiles : list of str
        Event files to merge.
    workdir : str, optional
        Directory in which to make the temporary working directory.

    Returns
    -------
    astropy.table.Table
        The selected events, with ``TIMEZERO`` added to ``TIME``.
    """

    tmpdir = tempfile.mkdtemp(dir=workdir)

    # Build input file for ftmerge
    evlistname = path.join(tmpdir, "evfiles.txt")
    fout = open(evlistname, "w")
    # evfilt_expr = '(EVENT_FLAGS==bx11000).and.((float)PI/(float)PI_FAST > (1.0 + 25.0/(float)PI + 4.0E-11*(float)PI**3))'
    evfilt_expr = "(EVENT_FLAGS==bx11000).and.((float)PI/(float)PI_FAST > (1.1 + 120.0/(float)PI + 0.0E-11*(float)PI**3))"
    for en in evfiles:
        print("{0}[{1}]".format(en, evfilt_expr), file=fout)
    fout.close()

    # Run ftmerge
    mergedname = path.join(tmpdir, "badratios.evt")
    cmd = [
        "ftmerge",
        "@{0}".format(evlistname),
        "outfile={0}".format(mergedname),
        "clobber=yes",
    ]
    runcmd(cmd)

    # Read merged event FITS file into a Table
    etable = Table.read(mergedname, hdu=1)
    if "TIMEZERO" in etable.meta:
        log.info(
            "Applying TIMEZERO of {0} to etable in get_badratioevents".format(
                etable.meta["TIMEZERO"]
            )
        )
        etable["TIME"] += etable.meta["TIMEZERO"]

    # Clean up
    os.remove(evlistname)
    os.remove(mergedname)
    os.rmdir(tmpdir)

    return etable
