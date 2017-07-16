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
from astropy.time import Time
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
phss = modelin.phase(ts.table)[1]
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

print(lcf)

print(dphi,dphierr)


# Use PINT's TOA writer to save the TOA



### JUNK BELOW HERE
sys.exit()
# poly = Polyco(polyconame,recalc_polycos=options.recalc_polycos)
# data = PhaseData(ft1name,poly,use_weights=(options.weights is not None),we_col_name=options.weights,wmin=options.wmin,emin=options.emin,emax=options.emax)
# if options.addphase: data.write_phase()

primitives,norms = prim_io(args.template)
template = LCTemplate(primitives,norms)
if args.weights is not None:
    logl = np.log(1+data.weights*(template(data.ph)-1)).sum()
else: logl = np.log(template(data.ph)).sum()
print('Total log likelihood:  %.2f'%(logl))

tg = UnbinnedTOAGenerator(data,poly,template,plot_stem=options.plot,good_ephemeris=(not options.blind))

if args.uniform_sigma:
    binner = UniformLogBinner(args.ntoa,data,template)
else:
    binner = TOABinner(args.ntoa,data)

toas,err_toas,tim_strings = tg.get_toas(binner)

if args.output is not None:
    f = open(args.output,'w')
    f.write('\n'.join(tim_strings))
    f.close()
    print('Wrote TOAs to %s'%(args.output))
else: print('\n'.join(tim_strings))
