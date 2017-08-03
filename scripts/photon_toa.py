#!/usr/bin/env python
# Program: photon_toa.py
# Authors: Paul S. Ray <paul.ray@nrl.navy.mil>
#          Matthew Kerr <matthew.kerr@gmail.com>
# Description:
# Reads a FITS file of photon event times (from NICER or another X-ray mission)
# and generates TOAs from the unbined times using a pulsar timing model
# and an analytic template. The TOAs can be output them in Tempo2 format.
from __future__ import division, print_function
import os,sys
import argparse
import numpy as np
from astropy import log
import astropy.units as u
import astropy.io.fits as pyfits
import pint.residuals
from pint.event_toas import load_NICER_TOAs
from pint.event_toas import load_RXTE_TOAs
from pint.event_toas import load_XMM_TOAs
from pint.plot_utils import phaseogram_binned
from pint.observatory.nicer_obs import NICERObs
from pint.observatory.rxte_obs import RXTEObs
import pint.toa, pint.models
from pint.eventstats import hmw, hm, h2sig
from astropy.time import Time, TimeDelta
#from os.path import exists
#from nicer.toabinner import TOABinner,UniformLogBinner
#from nicer.toagen import UnbinnedTOAGenerator
from pint.templates.lctemplate import LCTemplate,prim_io
from pint.templates import lcfitters

desc="""Generate TOAs from photon event data."""

parser=argparse.ArgumentParser(description=desc)
parser.add_argument("eventname",help="FITS file to read events from")
parser.add_argument("templatename",help="Name of file to read template from")
parser.add_argument("parname",help="Timing model file name")
parser.add_argument("--orbfile",help="Name of orbit file", default=None)
parser.add_argument("--planets",help="Use planetary Shapiro delay in calculations (default=False)", default=False, action="store_true")
parser.add_argument("--ephem",help="Planetary ephemeris to use (default=DE421)", default="DE421")
parser.add_argument("--plot",help="Show phaseogram plot.", action='store_true', default=False)
parser.add_argument("--plotfile",help="Output figure file name (default=None)", default=None)
#parser.add_argument("--ntoa","-n",type=int,default=10,help="Number of TOAs to produce between TSTART and TSTOP.")
# parser.add_argument("-b","--nbins",type=int,default=32,help="Number of bins in each profile.")
# parser.add_argument("-e","--emin",type=float,default=0.0,help="Minimum energy to include.")
# parser.add_argument("-x","--emax",type=float,default=300000.0,help="Maximum energy to include.")
# parser.add_argument("-p","--plot",default=None,help="Set to base name of files to generate likelihood surface plots for each TOA.")
# parser.add_argument("--wmin",type=float,default=0.0,help="Minimum weight value to include.")
# parser.add_argument("--wmax",type=float,default=1.0,help="Maximum weight value to include.")
# parser.add_argument("-t","--template",default=None,help="File name containing template (LCTemplate compatible)")
# parser.add_argument("-w","--weights",default=None,help="Specify column of FT1 file to use as weights for TOA computation")
# parser.add_argument("-u","--uniform_sigma",action="store_true",default=False,help="Instead of fixed time bins, make a TOA each time the H-test gets to a sufficient sigma.")
# parser.add_argument("--blind",action="store_true",default=False,help="Force blind search for TOAs rather than tracking.")
# parser.add_argument("-o","--output",default=None,help="File for output of .tim file.  Otherwise output to STDOUT.")

## Parse arguments
args = parser.parse_args()

# Load PINT model objects
modelin = pint.models.get_model(args.parname)
log.info(str(modelin))

# check for consistency between ephemeris and options
if (not args.planets) and (
    'SolarSystemShapiro' in modelin.components.keys()):
    log.error(
            "Ephemeris specifies Shapiro delay calculation.  Removing it.")
    components = modelin.components
    components.pop('SolarSystemShapiro')
    modelin.setup_components(components.values())

# Load Template objects
primitives,norms = prim_io(args.templatename)
template = LCTemplate(primitives,norms)
print(template)

# Load photons as PINT toas, and weights, if specified
# Here I might loop over the files specified
# Read event file header to figure out what instrument is is from
hdr = pyfits.getheader(args.eventname,ext=1)

log.info('Event file TELESCOPE = {0}, INSTRUMENT = {1}'.format(hdr['TELESCOP'],
    hdr['INSTRUME']))

if hdr['TELESCOP'] == 'NICER':
    # Instantiate NICERObs once so it gets added to the observatory registry
    if args.orbfile is not None:
        log.info('Setting up NICER observatory')
        NICERObs(name='NICER',FPorbname=args.orbfile,tt2tdb_mode='none')
    # Read event file and return list of TOA objects
    tl  = load_NICER_TOAs(args.eventname)
elif hdr['TELESCOP'] == 'XTE':
    # Instantiate RXTEObs once so it gets added to the observatory registry
    if args.orbfile is not None:
        # Determine what observatory type is.
        log.info('Setting up RXTE observatory')
        RXTEObs(name='RXTE',FPorbname=args.orbfile,tt2tdb_mode='none')
    # Read event file and return list of TOA objects
    tl  = load_RXTE_TOAs(args.eventname)
elif hdr['TELESCOP'].startswith('XMM'):
    # Not loading orbit file here, since that is not yet supported.
    tl  = load_XMM_TOAs(args.eventname)
else:
    log.error("FITS file not recognized, TELESCOPE = {0}, INSTRUMENT = {1}".format(
        hdr['TELESCOP'], hdr['INSTRUME']))
    sys.exit(1)

# Now convert to TOAs object and compute TDBs and posvels
ts = pint.toa.TOAs(toalist=tl)
ts.filename = args.eventname
ts.compute_TDBs()
ts.compute_posvels(ephem=args.ephem,planets=args.planets)

print(ts.get_summary())
mjds = ts.get_mjds()
print(mjds.min(),mjds.max())

# Compute model phase for each TOA
phss = modelin.phase(ts.table)[1].value # discard units
# ensure all postive
phases = np.where(phss < 0.0, phss + 1.0, phss)
mjds = ts.get_mjds()
h = float(hm(phases))
print("Htest : {0:.2f} ({1:.2f} sigma)".format(h,h2sig(h)))
if args.plot:
    phaseogram_binned(mjds,phases,bins=100,plotfile = args.plotfile)

# Given some subset of the event times, phases, and weights, compute
# the TOA based on a reference event near the middle of the span.
# Build the TOA as a PINT TOA() object
lcf = lcfitters.LCFitter(template,phases)
dphi,dphierr = lcf.fit_position()

# find MJD closest to center of observation and turn it into a TOA
argmid = np.searchsorted(mjds,0.5*(mjds.min()+mjds.max()))
tmid = ts.table['tdb'][argmid]
tplus = tmid + TimeDelta(1*u.s,scale='tdb')
toamid = pint.toa.TOA(tmid)
toaplus = pint.toa.TOA(tplus)
toas = pint.toa.TOAs(toalist=[toamid,toaplus])
toas.compute_TDBs()
toas.compute_posvels(ephem=args.ephem,planets=args.planets)
phsi,phsf = modelin.phase(toas.table)
fbary = (phsi[1]-phsi[0]) + (phsf[1]-phsf[0])
fbary._unit = u.Hz
tfinal = tmid + TimeDelta((-(dphi+phsf[0].value)/fbary),scale='tdb')

# get exposure information
try:
    f = pyfits.open(args.eventname)
    exposure = f[1].header['exposure']
    f.close()
except:
    exposure = 0

# Use PINT's TOA writer to save the TOA
nsrc = lcf.template.norm()*len(lcf.phases)
nbkg = (1-lcf.template.norm())*len(lcf.phases)
toafinal = pint.toa.TOA(tfinal,
        nsrc='%.2f'%nsrc,nbkg='%.2f'%nbkg,exposure='%.2f'%exposure)
toas = pint.toa.TOAs(toalist=[toafinal])
toas.table['error'][:] = dphierr/fbary*1e6
toas.write_TOA_file(sys.stdout,name='nicer',format='tempo2')


