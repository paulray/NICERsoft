#!/usr/bin/env python
# Program: photon_toa.py
# Authors: Paul S. Ray <paul.ray@nrl.navy.mil>
#          Matthew Kerr <matthew.kerr@gmail.com>
# Description:
# Reads a FITS file of photon event times (from NICER or another X-ray mission)
# and generates TOAs from the unbined times using a pulsar timing model
# and an analytic template. The TOAs can be output them in Tempo2 format.
from __future__ import division, print_function
import os, sys
import argparse
import numpy as np
from astropy import log

# import astropy.units as u
import astropy.io.fits as pyfits
from pint.eventstats import hmw, hm, h2sig
from pint.templates.lctemplate import LCTemplate, prim_io
from pint.templates import lcfitters
import pickle

desc = """Perform a template fit using a stored phase column and print source and background properties."""

parser = argparse.ArgumentParser(description=desc)
parser.add_argument("eventname", help="FITS file to read events from")
parser.add_argument("templatename", help="Name of file to read template from")
parser.add_argument(
    "--phasecol", help="FITS column to use for pulse phase.", default="PULSE_PHASE"
)

## Parse arguments
args = parser.parse_args()

# Load Template objects
try:
    template = pickle.load(open(args.templatename, "rb"), encoding="latin1")
except:
    primitives, norms = prim_io(args.templatename)
    template = LCTemplate(primitives, norms)

log.info("Template properties: \n{0}".format(str(template)))

# load up events and exposure information
f = pyfits.open(args.eventname)
phases = f[1].data.field(args.phasecol)
try:
    exposure = f[1].header["exposure"]
except:
    exposure = 0


h = float(hm(phases))
log.info("Htest : {0:.2f} ({1:.2f} sigma)".format(h, h2sig(h)))

# make sure template shape is fixed
lcf = lcfitters.LCFitter(template, phases)
for prim in lcf.template.primitives:
    prim.free[:] = False

# do iterative fit of position and background
for i in range(2):
    lcf.fit_position(unbinned=False)
    lcf.fit_background(unbinned=False)

# get exposure information
try:
    f = pyfits.open(args.eventname)
    f.close()
except:
    exposure = 0

# Use PINT's TOA writer to save the TOA
nsrc = lcf.template.norm() * len(lcf.phases)
nbkg = (1 - lcf.template.norm()) * len(lcf.phases)
log.info("Exposure = {0:.2f} s, Nsrc = {1}, Nbkg = {2}".format(exposure, nsrc, nbkg))
log.info(
    "Src rate = {0:.3f} c/s, Bkg rate = {1:.3f} c/s".format(
        nsrc / exposure, nbkg / exposure
    )
)
