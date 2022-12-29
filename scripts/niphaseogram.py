#!/usr/bin/env python
# Program: niphaseogram.py
# Author: Paul S. Ray <paul.ray@nrl.navy.mil>
# Description:
# Reads a FITS event file with  event times and a pulse phase column and makes a
# 2-D phaseogram plot
import sys, math, os
from optparse import OptionParser
import numpy as np
import pylab
import astropy.io.fits as pyfits
import scipy.stats
import datetime
import matplotlib.pyplot as plt
import astropy.units as u
from astropy import log
from nicer.values import *
from pint.fits_utils import read_fits_event_mjds
from astropy.time import Time

desc = (
    "Read a FITS event file containing a PULSE_PHASE column and make a 2-d phaseogram"
)
parser = OptionParser(usage=" %prog [options] [FT1_FILENAME]", description=desc)
parser.add_option(
    "-n",
    "--ntoa",
    type="int",
    default=60,
    help="Number of TOAs to produce between TSTART and TSTOP.",
)
parser.add_option(
    "-b", "--nbins", type="int", default=32, help="Number of bins in each profile."
)
parser.add_option(
    "-e", "--emin", type="float", default=0.25, help="Minimum energy to include."
)
parser.add_option(
    "-x", "--emax", type="float", default=12.0, help="Maximum energy to include."
)
parser.add_option(
    "-o",
    "--outfile",
    type="string",
    default=None,
    help="File name for plot file.  Type comes from extension.",
)
parser.add_option(
    "-r", "--radio", type="string", default=None, help="Radio profile to overplot."
)
parser.add_option(
    "-B",
    "--radiobins",
    type="int",
    default=128,
    help="Number of bins in each profile (default=128)",
)
parser.add_option(
    "-w",
    "--weights",
    type="string",
    default=None,
    help="FITS column to use as photon weight.",
)
parser.add_option(
    "-t", "--tmin", type="float", default=0.0, help="Minimum time to include (MJD)"
)
# parser.add_option("-t","--tmax",type="float",default=0.0,help="Maximum time to include (MJD)")
## Parse arguments
(options, args) = parser.parse_args()
if len(args) != 1:
    parser.error("event FILTS file argument is required.")

evname = args[0]

# Read FT1 file
hdulist = pyfits.open(evname)
evhdr = hdulist[1].header
evdat = hdulist[1].data

# Hack for XMM barycentered data
if evhdr["TELESCOP"].startswith("XMM"):
    log.info("XMM data, setting MET0 to 50814.0")
    MET0 = Time(50814.0, format="mjd", scale="tdb")

# Hack for NuSTAR barycentered data
if evhdr["TELESCOP"].startswith("NuSTAR"):
    log.info("NuSTAR data, setting MET0 to 55197.00076601852")
    MET0 = Time(55197.00076601852, format="mjd", scale="tdb")

if evhdr["TELESCOP"].startswith("SWIFT"):
    log.info("SWIFT data, setting MET0 to 51910.0007428703700000")
    MET0 = Time(51910.000742870370, format="mjd", scale="tdb")

# mets = evdat.field('TIME')
# mets += evhdr['TIMEZERO']
mjds = read_fits_event_mjds(hdulist[1])
# log.info('MJDs {0}'.format(mjds[:10]))
# evtimes = MET0 + mets*u.s
# mjds = evtimes.mjd
# log.info('Evtimes {0}'.format(evtimes[:10]))
try:
    phases = evdat.field("PULSE_PHASE")
except:
    phases = evdat.field("PHASE")

if options.weights is not None:
    weights = evdat.field(options.weights)

# FILTER ON ENERGY
try:
    en = evdat.field("PI") * PI_TO_KEV
    idx = np.where(np.logical_and((en >= options.emin), (en <= options.emax)))
    mjds = mjds[idx]
    phases = phases[idx]
    log.info("Energy cuts left %d out of %d events." % (len(mjds), len(en)))
except:
    log.warning("PI column not found. No energy cuts applied")

TSTART = float(evhdr["TSTART"])
TSTARTtime = MET0 + TSTART * u.s
TSTOP = float(evhdr["TSTOP"])
TSTOPtime = MET0 + TSTOP * u.s

# Compute MJDSTART and MJDSTOP in MJD
MJDSTART = TSTARTtime.mjd
MJDSTOP = TSTOPtime.mjd

if options.tmin != 0:
    MJDSTART = options.tmin

# Compute observation duration for each TOA
toadur = (MJDSTOP - MJDSTART) / options.ntoa
log.info("MJDSTART {0}, MJDSTOP {1}, toadur {2}".format(MJDSTART, MJDSTOP, toadur))

mjdstarts = MJDSTART + toadur * np.arange(options.ntoa, dtype=float)
mjdstops = mjdstarts + toadur

# Build array of bins for histogram call. Is n+1 long because it includes enpoint (1.0)
profbins = np.linspace(0.0, 1.0, options.nbins + 1, endpoint=True)

# Loop over blocks to process
a = None
for tstart, tstop in zip(mjdstarts, mjdstops):

    if options.tmin != 0 and tstart < options.tmin:
        continue

    idx = (mjds > tstart) & (mjds < tstop)

    if options.weights is not None:
        profile, edges = np.histogram(phases[idx], bins=profbins, weights=weights[idx])
    else:
        profile, edges = np.histogram(phases[idx], bins=profbins)

    if a is None:
        a = profile
        fullprof = profile
    else:
        a = np.append(a, profile)
        fullprof += profile

log.info("Total photons = {0}".format(fullprof.sum()))
log.info(
    "DC level = {0} counts per bin ({1} bins)".format(fullprof.min(), options.nbins)
)
try:
    log.info("Exposure = {0}".format(evhdr["EXPOSURE"]))
    log.info("Total count rate = {0}".format((fullprof.sum()) / evhdr["EXPOSURE"]))
    log.info(
        "Pulsed count rate (estimated from minimum profile bin) = {0}".format(
            (fullprof.sum() - options.nbins * fullprof.min()) / evhdr["EXPOSURE"]
        )
    )
except:
    pass

b = a.reshape(options.ntoa, options.nbins)

c = np.hstack([b, b])

pylab.figure(1, figsize=(6.0, 8.0))

# This axis is the 2-d phaseogram
ax1 = pylab.axes([0.15, 0.05, 0.75, 0.6])
ax1.imshow(
    c,
    interpolation="nearest",
    origin="lower",
    cmap=pylab.cm.binary,
    extent=(0, 2.0, MJDSTART, MJDSTOP),
    aspect=2.0 / (MJDSTOP - MJDSTART),
)
pylab.xlabel("Pulse Phase")
pylab.ylabel("Time (MJD)")

# This axis is the summed profile plot
ax2 = pylab.axes([0.15, 0.65, 0.75, 0.3])
bbins = np.concatenate((profbins, profbins[1:] + 1.0))
ax2.step(
    bbins,
    np.concatenate((fullprof, fullprof, np.array([fullprof[0]]))),
    where="post",
    color="k",
    lw=1.5,
)
py = np.concatenate((fullprof, fullprof))
pylab.ylabel("Photons")
pylab.setp(ax2.get_xticklabels(), visible=False)
pymin = py.min() - 0.1 * (py.max() - py.min())
pymax = py.max() + 0.1 * (py.max() - py.min())
# pylab.ylim(ymin=0.0)
pylab.ylim((pymin, pymax))
pylab.xlim((0.0, 2.0))
# ax2.set_xticks(np.arange(20.0)/10.0)
ax2.minorticks_on()
# pylab.xlabel('Pulse Phase')

# Add radio profile
if options.radio is not None:
    x, y = np.loadtxt(options.radio, unpack=True)
    # import psr_utils
    # x = np.arange(len(y))/len(y)
    # y = psr_utils.fft_rotate(y,0.2498*len(x))
    # y=y-y[20:50].mean()

    if len(y) > options.radiobins:
        # Average every N points to bring the profile down to number of radiobins
        y = y.reshape((options.radiobins, int(len(y) / options.radiobins)))
        y = np.mean(y, 1)
    y = y - y.min()
    y = y / y.max()
    y = y * (py.max() - py.min()) + py.min()

    pylab.plot(
        np.arange(2.0 * len(y)) / len(y),
        np.concatenate((y, y)),
        linewidth=1.5,
        color="r",
    )

# pylab.subplots_adjust(hspace=0)

if options.outfile is not None:
    log.info("writing outfile")
    pylab.savefig(options.outfile)
else:
    pylab.show()
