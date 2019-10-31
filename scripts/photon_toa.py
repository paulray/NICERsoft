#!/usr/bin/env python
# Program: photon_toa.py
# Authors: Paul S. Ray <paul.ray@nrl.navy.mil>
#          Matthew Kerr <matthew.kerr@gmail.com>
#          Julia Deneva <julia.deneva@gmail.com>
# Description:
# Reads a FITS file of photon event times (from NICER or another X-ray mission)
# and generates TOAs from the unbined times using a pulsar timing model
# and an analytic template. The TOAs can be output them in Tempo2 format.
from __future__ import division, print_function
#from future import standard_library
#standard_library.install_aliases()
from builtins import str
from builtins import zip
from builtins import range
import os,sys
import argparse
import numpy as np
from astropy import log
import astropy.units as u
import astropy.io.fits as pyfits
import pint.residuals
from pint.event_toas import load_NICER_TOAs
from pint.event_toas import load_RXTE_TOAs
from pint.event_toas import load_NuSTAR_TOAs
from pint.event_toas import load_XMM_TOAs
from pint.plot_utils import phaseogram_binned
from pint.observatory.nicer_obs import NICERObs
from pint.observatory.rxte_obs import RXTEObs
import pint.toa, pint.models
from pint.eventstats import hmw, hm, h2sig
from astropy.time import Time, TimeDelta
from pint.templates.lctemplate import LCTemplate,prim_io
from pint.templates import lcfitters
from copy import deepcopy
import pickle
import io
from collections import deque
import astropy.constants as const
from pint.observatory import get_observatory
from pint.observatory.special_locations import SpacecraftObs, BarycenterObs

def local_load_NICER_TOAs(eventname):
    """ Local override to add MET field to each TOA object."""
    # TODO -- add this to PINT method ?
    tl  = load_NICER_TOAs(eventname)
    f = pyfits.open(eventname)
    mets = f['events'].data.field('time')
    f.close()
    for t,met in zip(tl,mets):
        t.met = met
    return tl
# The returned tl has topocentric TT MJD photon times; TIMESYS=TT, TIMEREF=LOCAL in the .evt file

desc="""Generate TOAs from photon event data."""

parser=argparse.ArgumentParser(description=desc)
parser.add_argument("eventname",help="FITS file to read events from")
parser.add_argument("templatename",help="Name of file to read template from")
parser.add_argument("parname",help="Timing model file name")
parser.add_argument("--orbfile",help="Name of orbit file", default=None)
parser.add_argument("--ephem",help="Planetary ephemeris to use (default=DE421)", default="DE421")
parser.add_argument("--plot",help="Show phaseogram plot.", action='store_true', default=False)
parser.add_argument("--plotfile",help="Output figure file name (default=None)", default=None)
#parser.add_argument("--fitbg",help="Fit an overall background level (e.g. for changing particle background level (default=False).",action='store_true',default=False)
parser.add_argument("--unbinned",help="Fit position with unbinned likelihood.  Don't use for large data sets. (default=False)",action='store_true',default=False)
#parser.add_argument("--fix",help="Adjust times to fix 1.0 second offset in NICER data (default=False)", action='store_true',default=False)
parser.add_argument("--tint",help="Integrate for tint seconds for each TOA, or until the total integration exceeds maxint.  The algorithm is based on GTI, so the integration will slightly exceed tint (default None; see maxint.)",default=None)
parser.add_argument("--maxint",help="Maximum time interval to accumulate exposure for a single TOA (default=2*86400s)",type=float, default=2*86400.)
parser.add_argument("--minexp",help="Minimum exposure (s) for which to include a TOA (default=0.0).",default=0.0,type=float)
parser.add_argument("--track",help="Assume model is close to good and only search near 0 phase (to avoid getting TOAs off by 0.5 in double peaked pulsars)",action="store_true",default=False)
parser.add_argument("--dice",help="Dice up long GTIs into chunks of length <= tint",action="store_true",default=False)
parser.add_argument("--use_bipm",help="Use BIPM clock corrections",action="store_true",default=False)
parser.add_argument("--use_gps",help="Use GPS to UTC clock corrections",action="store_true",default=False)
parser.add_argument("--topo",help="Make topocentric TOAs; include the spacecraft ECI position on the TOA line",action="store_true",default=False)
parser.add_argument("--outfile",help="Name of file to save TOAs to (default is STDOUT)",default=None)

## Parse arguments
args = parser.parse_args()

# Load PINT model objects
modelin = pint.models.get_model(args.parname)
log.info(str(modelin))

# check for consistency between ephemeris and options
if modelin.PLANET_SHAPIRO.quantity:
    planets=True
else:
    planets=False

# Load Template objects
try:
    template = pickle.load(file(args.templatename))
except:
    primitives,norms = prim_io(args.templatename)
    template = LCTemplate(primitives,norms)
#print(template)

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
        obs = NICERObs(name='NICER',FPorbname=args.orbfile)
    # Read event file and return list of TOA objects
    try:
        tl  = local_load_NICER_TOAs(args.eventname)
    except KeyError:
        log.error('Failed to load NICER TOAs. Make sure orbit file is specified on command line!')
        raise
elif hdr['TELESCOP'] == 'XTE':
    # Instantiate RXTEObs once so it gets added to the observatory registry
    if args.orbfile is not None:
        # Determine what observatory type is.
        log.info('Setting up RXTE observatory')
        obs = RXTEObs(name='RXTE',FPorbname=args.orbfile)
    # Read event file and return list of TOA objects
    tl  = load_RXTE_TOAs(args.eventname)
elif hdr['TELESCOP'].startswith('XMM'):
    # Not loading orbit file here, since that is not yet supported.
    tl  = load_XMM_TOAs(args.eventname)
elif hdr['TELESCOP'].startswith('NuSTAR'):
    # Not loading orbit file here, since that is not yet supported.
    tl  = load_NuSTAR_TOAs(args.eventname)
    f = pyfits.open(args.eventname)
    mets = f['events'].data.field('time')
    f.close()
    for t,met in zip(tl,mets):
        t.met = met
else:
    log.error("FITS file not recognized, TELESCOPE = {0}, INSTRUMENT = {1}".format(
        hdr['TELESCOP'], hdr['INSTRUME']))
    sys.exit(1)

if args.topo: #for writing UTC topo toas
    SpacecraftObs(name='spacecraft')

if not args.topo and (args.orbfile is None):
    obs = 'Barycenter'

if len(tl) <= 0:
    log.error('No TOAs found. Aborting.')
    sys.exit(1)

# Now convert to TOAs object and compute TDBs and (SSB) posvels
ts = pint.toa.get_TOAs_list(tl,ephem=args.ephem,planets=planets,include_bipm=False,include_gps=False)
ts.filename = args.eventname

#print(ts.get_summary())
mjds = ts.get_mjds() # TT topocentric MJDs as floats; only used to find the index of the photon time closest to the middle of the MJD range

# Compute model phase for each TOA; 
#(TODO--'tdbs' var. needs to be renamed since it's TT and not TDB for the topocentric photon times in the if block)
phss = modelin.phase(ts,abs_phase=True)[1].value # discard units
if args.topo: 
    tdbs = ts.table['mjd']
else:
    tdbs = ts.table['tdb']

# ensure all positive
phases = np.where(phss < 0.0, phss + 1.0, phss)

h = float(hm(phases))
print("Htest : {0:.2f} ({1:.2f} sigma)".format(h,h2sig(h)))
if args.plot:
    phaseogram_binned(mjds,phases,bins=100,plotfile = args.plotfile)

# get exposure information
try:
    f = pyfits.open(args.eventname)
    exposure = f[1].header['exposure']
    f.close()
except:
    exposure = 0


def estimate_toa(mjds,phases,tdbs,topo,obs):
    """ Return a pint TOA object for the provided times and phases."""

    # Given some subset of the event times, phases, and weights, compute
    # the TOA based on a reference event near the middle of the span.
    # Build the TOA as a PINT TOA() object
    lcf = lcfitters.LCFitter(deepcopy(template),phases)
# fitbg does not work!  Disabling.
#    if args.fitbg:
#        for i in xrange(2):
#            lcf.fit_position(unbinned=False)
#            lcf.fit_background(unbinned=False)
    dphi,dphierr = lcf.fit_position(unbinned=args.unbinned,track=args.track)
    log.info('Measured phase shift dphi={0}, dphierr={1}'.format(dphi,dphierr))

    # find MJD closest to center of observation and turn it into a TOA
    argmid = np.searchsorted(mjds,0.5*(mjds.min()+mjds.max()))
    tmid = tdbs[argmid] # Should we used tdbld?

    if topo: #tdbs are topocentric TT MJD times
        tplus = tmid + TimeDelta(1*u.s,scale='tt')
        toamid = pint.toa.TOA(tmid,obs=obs.name)
        toaplus = pint.toa.TOA(tplus,obs=obs.name)
    else: #tdbs are TDB MJD times
        tplus = tmid + TimeDelta(1*u.s,scale='tdb')
        toamid = pint.toa.TOA(tmid)
        toaplus = pint.toa.TOA(tplus)

    toas = pint.toa.get_TOAs_list([toamid,toaplus],include_gps=args.use_gps,
        include_bipm=args.use_bipm,ephem=args.ephem,planets=planets)

    phsi,phsf = modelin.phase(toas,abs_phase=True)
    if topo:
        sc = 'tt'
    else:
        sc = 'tdb'
    # Compute frequency = d(phase)/dt
    f = (phsi[1]-phsi[0]) + (phsf[1]-phsf[0])
    f._unit = u.Hz

    # First delta is to get time of phase 0.0 of initial model
    # Second term corrects for the measured phase offset to align with template
    tfinal = tmid + TimeDelta(-phsf[0].value/f,scale=sc) + TimeDelta(dphi/f,scale=sc)
    
    # Use PINT's TOA writer to save the TOA
    nsrc = lcf.template.norm()*len(lcf.phases)
    nbkg = (1-lcf.template.norm())*len(lcf.phases)

    if args.topo: #tfinal is a topocentric TT MJD
        telposvel = obs.posvel_gcrs(tfinal)
        x = telposvel.pos[0].to(u.km)
        y = telposvel.pos[1].to(u.km)
        z = telposvel.pos[2].to(u.km)
        vx = telposvel.vel[0].to(u.km/u.s)
        vy = telposvel.vel[1].to(u.km/u.s)
        vz = telposvel.vel[2].to(u.km/u.s)
        
        toafinal = pint.toa.TOA(tfinal.utc,scale='utc',obs='spacecraft',nsrc='%.2f'%nsrc,nbkg='%.2f'%nbkg,exposure='%.2f'%exposure,dphi='%.5f'%dphi,mjdTT='%.8f'%tfinal.tt.mjd,
                                telx='%.8f'%x.value,tely='%.8f'%y.value,telz='%.8f'%z.value,
                                vx='%.8f'%vx.value,vy='%.8f'%vy.value,vz='%.8f'%vz.value)

    else:
        toafinal = pint.toa.TOA(tfinal,
                                nsrc='%.2f'%nsrc,nbkg='%.2f'%nbkg,exposure='%.2f'%exposure,dphi='%.5f'%dphi)


    log.info("Src rate = {0} c/s, Bkg rate = {1} c/s".format(nsrc/exposure, nbkg/exposure))
    return toafinal,dphierr/f.value*1.e6

if args.tint is None:

    # do a single TOA for table
    toafinal,toafinal_err = estimate_toa(mjds,phases,tdbs,args.topo,obs)
    if 'OBS_ID' in hdr:
        # Add ObsID to the TOA flags
        toafinal.flags['obsid'] = hdr['OBS_ID']
        toafinal.flags['htest'] = "{0:.2f}".format(hm(phases))
    toafinal = [toafinal]
    toafinal_err = [toafinal_err]
else:
    # Load in GTIs
    f = pyfits.open(args.eventname)
    # Warning:: This is ignoring TIMEZERO!!!!
    gti_t0 = f['gti'].data.field('start')
    gti_t1 = f['gti'].data.field('stop')
    gti_dt = gti_t1-gti_t0
    mets = np.asarray([t.met for t in tl])

    tint = float(args.tint)

    if args.dice:
        # Break up larger GTIs into small chunks
        new_t0s = deque()
        new_t1s = deque()
        for t0,t1 in zip(gti_t0,gti_t1):
            dt = t1-t0
            if dt < tint:
                new_t0s.append(t0)
                new_t1s.append(t1)
            else:
                # break up GTI in such a way to avoid losing time (to tmin) and
                # to avoid having pieces longer than tint
                npiece = int(np.floor(dt/tint))+1
                new_edges = np.linspace(t0,t1,npiece+1)
                for it0,it1 in zip(new_edges[:-1],new_edges[1:]):
                    new_t0s.append(it0)
                    new_t1s.append(it1)
        gti_t0 = np.asarray(new_t0s)
        gti_t1 = np.asarray(new_t1s)
        gti_dt = gti_t1-gti_t0

    # the algorithm here is simple -- go through the GTI and add them up
    # until either the good time exceeds tint, or until the total time
    # interval exceeds maxint
    i0 = 0
    current = 0.0
    toas = deque()
    maxint = float(args.maxint)
    for i in range(len(gti_t0)):
        current += gti_dt[i]
        #print('iteration=%d, current=%f'%(i,current))
        if (current >= tint) or ((gti_t1[i]-gti_t0[i0])>maxint) or (i==len(gti_t0)-1):
            # make a TOA
            ph0,ph1 = np.searchsorted(mets,[gti_t0[i0],gti_t1[i]])
            m,p,t = mjds[ph0:ph1],phases[ph0:ph1],tdbs[ph0:ph1]
            #print('Generating TOA ph0={0}, ph1={1}, len(m)={2}, i0={3}, i={4}'.format(ph0,ph1,len(m),i0,i))
            #print('m[0]={0}, m[1]={1}'.format(m[0],m[-1]))
            if len(m) > 0:
                toas.append(estimate_toa(m,p,t,args.topo,obs))
                toas[-1][0].flags['htest'] = "{0:.2f}".format(hm(p))
                # fix exposure
                toas[-1][0].flags['exposure'] = current
            current = 0.0
            i0 = i+1
    toafinal,toafinal_err = list(zip(*toas))

if args.minexp > 0.0:
    x = [(t,e) for t,e in zip(toafinal,toafinal_err) if float(t.flags['exposure']) > args.minexp]
    if len(x) > 0:
        toafinal,toafinal_err= list(zip(*x))
    else:
        print('No TOAs passed exposure cut!')
        sys.exit(0)

for t in toafinal:
    t.flags["-t"] = hdr['TELESCOP']
toas = pint.toa.TOAs(toalist=toafinal)
toas.table['error'][:] = np.asarray(toafinal_err)
sio = io.StringIO()
toas.write_TOA_file(sio,name='nicer',format='tempo2')
output = sio.getvalue()

if args.topo:
    output = output.replace('spacecraft','STL_GEO')
else:
    output = output.replace('bat','@')

if args.outfile is not None:
    print(output,file=open(args.outfile,"w"))
else:
    print(output)
