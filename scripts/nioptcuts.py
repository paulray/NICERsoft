#!/usr/bin/env python
from __future__ import (print_function, division, unicode_literals, absolute_import)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
from astropy import log
from nicer.values import *
from astropy.table import Table, vstack
import astropy.units as u
from pint.eventstats import hm

def cached_hm(mask):
    nph = mask.sum()
    if nph == 0:
        return 0
    s = (cached_hm._cache[...,mask].sum(axis=2)**2).sum(axis=1)
    return ((2./nph)*np.cumsum(s)-4*np.arange(0,20)).max()

parser = argparse.ArgumentParser(description="Read event file with phase and optimize cuts")
parser.add_argument("evfiles", help="Name of event files to process", nargs='+')
parser.add_argument("--noplot", help="Suppress plotting.", action="store_true",default=False)
args = parser.parse_args()

tlist = []
for fn in args.evfiles:
    log.info('Reading file {0}'.format(fn))
    tlist.append(Table.read(fn,hdu=1))
log.info('Concatenating files')
if len(tlist) == 1:
    etable = tlist[0]
else:
    etable = vstack(tlist,metadata_conflicts='silent')
del tlist

phasesinitial = etable['PULSE_PHASE'].astype(np.float32)
hbest = hm(phasesinitial)
eminbest = 0.0
emaxbest = 100.0
print("Initial Htest = {0}".format(hbest))

# assemble cache
m = 20
cache = np.empty([m,2,len(phasesinitial)],dtype=np.float32)
for i in xrange(m):
    cache[i,0] = np.cos(phasesinitial*(2*np.pi*(i+1)))
    cache[i,1] = np.sin(phasesinitial*(2*np.pi*(i+1)))
cached_hm._cache = cache

emins = np.arange(0.25,4.0,0.025)
for emin in emins:
    emaxs = np.arange(emin+0.1,12.01,0.1)
    for emax in emaxs:
        mask = np.logical_and(etable['PI']*PI_TO_KEV>emin,
            etable['PI']*PI_TO_KEV<emax)
        h = cached_hm(mask)
        if h>=hbest:
            hbest=h
            eminbest=emin
            emaxbest=emax


print("Final Htest = {0}".format(hbest))
print("Best emin {0} emax {1}".format(eminbest,emaxbest))

if not args.noplot:
    fig,ax = plt.subplots()
    ax.hist(phasesinitial, bins = 32)
    idx = np.where(np.logical_and(etable['PI']*PI_TO_KEV>eminbest,
         etable['PI']*PI_TO_KEV<emaxbest))[0]
    phases = etable['PULSE_PHASE'][idx]
    ax.hist(phases,bins=32)
    ax.text(0.1, 0.1, 'H = {0:.2f}'.format(hbest), transform=ax.transAxes)
    ax.set_ylabel('Counts')
    ax.set_xlabel('Pulse Phase')
    ax.set_title('Pulse Profile')

    plt.show()
