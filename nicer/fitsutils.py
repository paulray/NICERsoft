from __future__ import print_function, division
import os
import matplotlib.pyplot as plt
import numpy as np
import argparse
from astropy import log
from os import path
from glob import glob
from subprocess import check_call
import shutil
import tempfile
from astropy.table import Table

from nicer.values import *

def runcmd(cmd):
    # CMD should be a list of strings since it is not processed by a shell
    log.info('CMD: '+" ".join(cmd))
    try:
        check_call(cmd,env=os.environ)
    except OSError:
        log.error('Command failed! Did you install and initialize FTOOLS?')
        raise


def filtallandmerge_ftools(evfiles,workdir=None):
    'Merges and filters a set of event files, returning an etable'

    tmpdir = tempfile.mkdtemp(dir=workdir)

    # Build input file for ftmerge
    evlistname=path.join(tmpdir,'evfiles.txt')
    fout = open(evlistname,'w')
    #evfilt_expr = '(PI>20).and.(EVENT_FLAGS==bx1x000)'
    evfilt_expr = '(EVENT_FLAGS==bx1x000)'
    for en in evfiles:
        print('{0}[{1}]'.format(en,evfilt_expr),file=fout)
    fout.close()

    # Run ftmerge
    mergedname = path.join(tmpdir,"merged.evt")
    cmd = ["ftmerge", "@{0}".format(evlistname), "outfile={0}".format(mergedname),
        "clobber=yes"]
    runcmd(cmd)

    # Read merged event FITS file into a Table
    etable = Table.read(mergedname,hdu=1)
    # Apply TIMEZERO if needed
    if 'TIMEZERO' in etable.meta:
        log.info('Applying TIMEZERO of {0} to etable in filtallandmerge'.format(etable.meta['TIMEZERO']))
        etable['TIME'] += etable.meta['TIMEZERO']

    # Clean up
    os.remove(evlistname)
    os.remove(mergedname)
    os.rmdir(tmpdir)

    return etable

def get_eventovershoots_ftools(evfiles,workdir=None):
    'Merges and filters a set of event files, returning an etable with only the OVERSHOOT-only events'

    tmpdir = tempfile.mkdtemp(dir=workdir)

    # Build input file for ftmerge
    evlistname=path.join(tmpdir,'evfiles.txt')
    fout = open(evlistname,'w')
    evfilt_expr = '(EVENT_FLAGS==bxxx010)'
    for en in evfiles:
        print('{0}[{1}]'.format(en,evfilt_expr),file=fout)
    fout.close()

    # Run ftmerge
    mergedname = path.join(tmpdir,"overshoots.evt")
    cmd = ["ftmerge", "@{0}".format(evlistname), "outfile={0}".format(mergedname),
        "clobber=yes"]
    runcmd(cmd)

    # Read merged event FITS file into a Table
    etable = Table.read(mergedname,hdu=1)
    if 'TIMEZERO' in etable.meta:
        log.info('Applying TIMEZERO of {0} to etable in get_eventovershoots'.format(etable.meta['TIMEZERO']))
        etable['TIME'] += etable.meta['TIMEZERO']

    # Clean up
    os.remove(evlistname)
    os.remove(mergedname)
    os.rmdir(tmpdir)

    return etable

def get_eventundershoots_ftools(evfiles,workdir=None):
    'Merges and filters a set of event files, returning an etable with only the UNDERSHOOT-only events'

    tmpdir = tempfile.mkdtemp(dir=workdir)

    # Build input file for ftmerge
    evlistname=path.join(tmpdir,'evfiles.txt')
    fout = open(evlistname,'w')
    evfilt_expr = '(EVENT_FLAGS==bxxx001)'
    for en in evfiles:
        print('{0}[{1}]'.format(en,evfilt_expr),file=fout)
    fout.close()

    # Run ftmerge
    mergedname = path.join(tmpdir,"undershoots.evt")
    cmd = ["ftmerge", "@{0}".format(evlistname), "outfile={0}".format(mergedname),
        "clobber=yes"]
    runcmd(cmd)

    # Read merged event FITS file into a Table
    etable = Table.read(mergedname,hdu=1)
    if 'TIMEZERO' in etable.meta:
        log.info('Applying TIMEZERO of {0} to etable in get_eventundershoots'.format(etable.meta['TIMEZERO']))
        etable['TIME'] += etable.meta['TIMEZERO']

    # Clean up
    os.remove(evlistname)
    os.remove(mergedname)
    os.rmdir(tmpdir)

    return etable

def get_eventbothshoots_ftools(evfiles,workdir=None):
    '''Merges and filters a set of event files, returning an etable with
    events that are both OVERSHOOT and UNDERSHOOT'''

    tmpdir = tempfile.mkdtemp(dir=workdir)

    # Build input file for ftmerge
    evlistname=path.join(tmpdir,'evfiles.txt')
    fout = open(evlistname,'w')
    evfilt_expr = '(EVENT_FLAGS==bxxx011)'
    for en in evfiles:
        print('{0}[{1}]'.format(en,evfilt_expr),file=fout)
    fout.close()

    # Run ftmerge
    mergedname = path.join(tmpdir,"bothshoots.evt")
    cmd = ["ftmerge", "@{0}".format(evlistname), "outfile={0}".format(mergedname),
        "clobber=yes"]
    runcmd(cmd)

    # Read merged event FITS file into a Table
    etable = Table.read(mergedname,hdu=1)
    if 'TIMEZERO' in etable.meta:
        log.info('Applying TIMEZERO of {0} to etable in get_eventbothshoots'.format(etable.meta['TIMEZERO']))
        etable['TIME'] += etable.meta['TIMEZERO']

    # Clean up
    os.remove(evlistname)
    os.remove(mergedname)
    os.rmdir(tmpdir)

    return etable

def get_badratioevents_ftools(evfiles,workdir=None):
    'Merges and filters a set of event files, returning an etable with only the events that are above the PI ratio trumpet'

    tmpdir = tempfile.mkdtemp(dir=workdir)

    # Build input file for ftmerge
    evlistname=path.join(tmpdir,'evfiles.txt')
    fout = open(evlistname,'w')
    evfilt_expr = '(EVENT_FLAGS==bx11000).and.((float)PI/(float)PI_FAST > (1.0 + 25.0/(float)PI + 4.0E-11*(float)PI**3))'
    for en in evfiles:
        print('{0}[{1}]'.format(en,evfilt_expr),file=fout)
    fout.close()

    # Run ftmerge
    mergedname = path.join(tmpdir,"badratios.evt")
    cmd = ["ftmerge", "@{0}".format(evlistname), "outfile={0}".format(mergedname),
        "clobber=yes"]
    runcmd(cmd)

    # Read merged event FITS file into a Table
    etable = Table.read(mergedname,hdu=1)
    if 'TIMEZERO' in etable.meta:
        log.info('Applying TIMEZERO of {0} to etable in get_badratioevents'.format(etable.meta['TIMEZERO']))
        etable['TIME'] += etable.meta['TIMEZERO']

    # Clean up
    os.remove(evlistname)
    os.remove(mergedname)
    os.rmdir(tmpdir)

    return etable
