#!/usr/bin/env python
# Author: Paul S. Ray <paul.ray@nrl.navy.mil> / Deepto Chakrabarty
# Description:
# Compute fractional RMS amplitude of a pulsation
#
from __future__ import division, print_function
import sys, math, os
from optparse import OptionParser
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import astropy.units as u
from astropy import log
from nicer.values import *
from pint.fits_utils import read_fits_event_mjds
from astropy.time import Time

desc="Calc RMS fractional amplitude of a pulsation"
parser=OptionParser(usage=" %prog [options] [FILENAME]",
                                        description=desc)
parser.add_option("-b","--nbins",type="int",default=32,help="Number of bins in each profile.")
parser.add_option("--bkg",type="float",default=0.0,help="Background rate to subtract (in counts/second)")
parser.add_option("--random",default=False,help="Simulate white noise")
## Parse arguments
(options,args) = parser.parse_args()
if len(args) != 1:
        parser.error("event FILTS file argument is required.")

evname = args[0]

# Read FT1 file
hdulist = pyfits.open(evname)
evhdr=hdulist[1].header
evdat=hdulist[1].data

try:
	phases = evdat.field('PULSE_PHASE')
except:
	phases = evdat.field('PHASE')

TSTART = float(evhdr['TSTART'])
TSTARTtime = MET0 + TSTART*u.s
TSTOP = float(evhdr['TSTOP'])
TSTOPtime = MET0 + TSTOP*u.s

# Try out totally random phases to simulate white noise
if options.random:
    log.warning('Replacing phases with RANDOM phases')
    phases = np.random.rand(len(phases))

# Compute MJDSTART and MJDSTOP in MJD
MJDSTART = TSTARTtime.mjd
MJDSTOP = TSTOPtime.mjd

# Build array of bins for histogram call. Is n+1 long because it includes enpoint (1.0)
profbins = np.linspace(0.0,1.0,options.nbins+1,endpoint=True)

fullprof, edges = np.histogram(phases,bins=profbins)

log.info('Total photons = {0}'.format(fullprof.sum()))
log.info('DC level = {0} counts per bin ({1} bins)'.format(fullprof.min(),options.nbins))
try:
	log.info('Exposure = {0}'.format(evhdr['EXPOSURE']))
	log.info('Total count rate = {0}'.format((fullprof.sum())/evhdr['EXPOSURE']))
	log.info('Pulsed count rate (estimated from minimum profile bin) = {0}'.format((fullprof.sum()-options.nbins*fullprof.min())/evhdr['EXPOSURE']))
except:
	pass


fullprof = np.asarray(fullprof,dtype=np.float)
if options.bkg > 0.0:
    bkgperbin = options.bkg*evhdr['EXPOSURE']/options.nbins
    log.info('Subtracting {0} counts/bin of background'.format(bkgperbin))
    fullprof -= bkgperbin

# Compute Fvar, which is the fractional RMS variability amplitude
# Equation 10 of Vaughan et al (2003, MNRAS, 345, 1271)
fracrms = np.sqrt(fullprof.var()-fullprof.mean())/fullprof.mean()
print('Fractional RMS is {0}'.format(fracrms))
