#!/usr/bin/env python
import sys
# Hack to add this to pythonpath
#sys.path.append('/Users/paulr/src/NICERsoft')
import matplotlib.pyplot as plt
import numpy as np
import argparse
from astropy import log
from astropy.table import Table, vstack
import astropy.io.fits as pyfits
import astropy.units as u
from astropy.time import Time
from os import path
from nicer.values import *
from nicer.cartographer import *
from nicer.plotutils import *
from nicer.sci_plots import sci_plots
from nicer.eng_plots import eng_plots
from glob import glob
from nicer.bkg_plots import *
from nicer.fitsutils import *

parser = argparse.ArgumentParser(description = "Plot the NICER data nicely.")
parser.add_argument("infiles", help="Input files", nargs='*')
parser.add_argument("--obsdir",help = "Find alllllllll the files!")
parser.add_argument("--object", help="Override object name", default=None)
parser.add_argument("--guessobj", help="Try to guess object from directory name", action="store_true")
parser.add_argument("--mask",help="Mask these IDS", nargs = '*', type=int, default=None)
parser.add_argument("-s", "--save", help = "Save plots to file", action = "store_true")
parser.add_argument("--sci", help = "Makes some nice science plots", action = "store_true")
parser.add_argument("--eng", help = "Makes some nice engineering plots", action = "store_true")
parser.add_argument("--ratio", help = "Display extra figure with diagnostic data", action = 'store_true')
parser.add_argument("--filtswtrig", help = "Filter SW TRIG events", action = "store_true")
parser.add_argument("--filtovershoot", help = "Filter OVERSHOOT events", action = "store_true")
parser.add_argument("--filtundershoot", help = "Filter UNDERSHOOT events", action = "store_true")
parser.add_argument("--filtratio", help="Filter PHA/PHA_FAST ratio (argument is ratio to cut at)", type=float, default=0.0)
parser.add_argument("--filtall", help = "Filter SWTRIG, UNDERSHOOT and OVERSHOOT events", action = "store_true")
parser.add_argument("--emin", help="Minimum energy (keV) to keep", default=-1.0, type=float)
parser.add_argument("--emax", help="Minimum energy (keV) to keep", default=-1.0, type=float)
parser.add_argument("--tskip", help="Seconds to skip at beginning of data", default=0.0, type=float)
parser.add_argument("--lcbinsize", help="Light curve bin size (s)", default=1.0, type=float)
parser.add_argument("--pi", help="Force use of internal PHA to PI conversion", action='store_true')
parser.add_argument("--basename", help="Basename for output plots", default=None)
parser.add_argument("--lclog", help = "make light curve log axis", action = "store_true")
parser.add_argument("--foldfreq", help="Make pulse profile by folding at a fixed freq (Hz)",
    default=0.0,type=float)
parser.add_argument("--nyquist", help="Nyquist freq for power spectrum (Hz)",
    default=100.0,type=float)
parser.add_argument("--map", help= "Creates a map with some stuff on it", action = 'store_true')
parser.add_argument("--orb", help="Path to orbit FITS filed", default = None)
parser.add_argument("--par", help="Path to par file", default = None)
parser.add_argument("--sps", help="Path to SPS HK file (_apid0260.hk)",default=None)
parser.add_argument("--powspec",help = "Display power spectrum (replaces ratio plot)", action = 'store_true')
parser.add_argument("--pslog", help = "make power spectrum log axis", action = "store_true")
parser.add_argument("--writeps", help = "write out power spectrum", action = "store_true")
parser.add_argument("--writeovershoot",help="Write summed overshoot rates to FITS file", action='store_true')
parser.add_argument("--applygti",help="Read GTI from provided FITS file", default=None)
parser.add_argument("--filtou",help="Filter Over/Undershoot Events for both flags", action='store_true')
args = parser.parse_args()

args.hkfiles = []
mkfiles = []

#Get file names from the directory
if args.obsdir:
    # Get names of event files from obsdir
    if len(args.infiles) == 0:
        # Here we could grab the clean data (cl) or the unfiltered merged (ufa) data.
        # Clean will have no flagged events
        # args.infiles = glob(path.join(args.obsdir,'xti/event_cl/ni*mpu?_ufa.evt'))
        args.infiles = glob(path.join(args.obsdir,'xti/event_cl/ni*mpu?_cl.evt'))
        args.infiles.sort()
    if len(args.infiles) == 0:
        log.error("No event files found!")
        sys.exit(1)
    log.info('Found event files: {0}'.format("\n" + "    \n".join(args.infiles)))

    # Get name of orbit file from obsdir
    try:
        args.orb = glob(path.join(args.obsdir,'auxil/ni*.orb'))[0]
    except:
        log.error("Orbit file not found!")
    log.info('Found the orbit file: {0}'.format(args.orb))

    # Get name of SPS HK file (apid0260)
    if args.sps is None:
        try:
            args.sps = glob(path.join(args.obsdir,'auxil/ni*_apid0260.hk'))[0]
        except:
            args.sps = None

    # Get name of MPU housekeeping files
    hkfiles = glob(path.join(args.obsdir,'xti/hk/ni*.hk'))
    hkfiles.sort()
    log.info('Found the MPU housekeeping files: {0}'.format("\n"+"\t\n".join(hkfiles)))
    args.hkfiles = hkfiles

    mkfiles = glob(path.join(args.obsdir,'auxil/ni*.mkf'))[0]

#Get the concatenated and filtered data
etable = filtandmerge(args.infiles,workdir=None)

if args.filtall:
    args.filtswtrig=True
    args.filtovershoot=True
    args.filtundershoot=True
    args.filtratio=1.4

if not args.sci and not args.eng and not args.map and not args.ratio:
    log.info("No specific plot requested, making all")
    args.sci = True
    args.eng = True
    args.map = True
    args.ratio = True

# Read the GTIs from the first event FITS file
gtitable = Table.read(args.infiles[0],hdu=2)
log.info('Got the good times from GTI')
gtitable['DURATION'] = gtitable['STOP']-gtitable['START']
# Only keep GTIs longer than 16 seconds
idx = np.where(gtitable['DURATION']>16.0)[0]
gtitable = gtitable[idx]
print(gtitable)

#Making the MK Table
log.info('Getting MKTable')
if len(mkfiles) > 0:
    mktable = Table.read(mkfiles,hdu=1)
else:
    mktable = None

# Change TIME column name to MET to reflect what it really is
etable.columns['TIME'].name = 'MET'
# Update exposure to be sum of GTI durations
etable.meta['EXPOSURE'] = gtitable['DURATION'].sum()

# Sort table by MET
etable.sort('MET')
log.info("Event MET Range : {0} to {1}".format(etable['MET'].min(),
    etable['MET'].max(), etable['MET'].max()-etable['MET'].min()))
log.info("TSTART {0}  TSTOP {1} (Span {2} seconds)".format(etable.meta['TSTART'],
    etable.meta['TSTOP'], etable.meta['TSTOP']-etable.meta['TSTART'] ))
log.info("DATE Range {0} to {1}".format(etable.meta['DATE-OBS'],
    etable.meta['DATE-END']))
if args.object is not None:
    etable.meta['OBJECT'] = args.object
'''
# Hack to trim first chunk of data
if args.tskip > 0.0:
    t0 = gtitable['START'][0]
    etable = etable[etable['MET']>t0+args.tskip]
    # Correct exposure (approximately)
    etable.meta['TSTART'] += args.tskip
    if gtitable['START'][0]+args.tskip < gtitable['STOP'][0]:
        gtitable['START'][0] += args.tskip
    else:
        log.error('Trying to skip more than first GTI segment!  **NOT IMPLEMENTED**')
        sys.exit(1)

'''
# If there are no PI columns, add them with approximate calibration
if args.pi or not ('PI' in etable.colnames):
    log.info('Adding PI')
    calfile = path.join(datadir,'gaincal_linear.txt')
    pi = calc_pi(etable,calfile)
    etable['PI'] = pi

#filtering out chosen IDS
if args.mask != None:
    log.info('Masking IDS')
    for id in args.mask:
        etable = etable[np.where(etable['DET_ID'] != id)]

if args.applygti is not None:
    g = Table.read(args.applygti)
    log.info('Applying external GTI from {0}'.format(args.applygti))
    g['DURATION'] = g['STOP']-g['START']
    # Only keep GTIs longer than 16 seconds
    g = g[np.where(g['DURATION']>16.0)]
    log.info('Applying external GTI')
    print g
    etable = apply_gti(etable,g)
    # Replacing this GTI does not work. It needs to be ANDed with the existing GTI
    #filttable.meta['EXPOSURE'] = g['DURATION'].sum()
    #gtitable = g

log.info('Exposure : {0:.2f}'.format(etable.meta['EXPOSURE']))

# Add Time column with astropy Time for ease of use and for PINT TOAs
if args.par is not None:
    log.info('Adding time column')
    # This should really be done the FITS way using MJDREF etc...
    # For now, just using MET0
    etime = filttable.columns['MET'] + MET0
    filttable['T'] = etime

# Set up the light curve bins, so we can have them for building
# light curves of various quantities, like overshoot rate and ratio filtered events
# Hmmm. The lc_elapsed_bins and lc_met_bins are never used, but CUMTIME is
startmet = gtitable['START'][0]
stopmet = gtitable['STOP'][0]
duration = stopmet-startmet
# Add 1 bin to make sure last bin covers last events
lc_elapsed_bins = np.arange(0.0,duration+args.lcbinsize,args.lcbinsize)
lc_met_bins = startmet+lc_elapsed_bins
cumtimes = [ 0.0 ]
cumtime = lc_elapsed_bins[-1]+args.lcbinsize
for i in range(1,len(gtitable['START'])):
    startmet = gtitable['START'][i]
    stopmet = gtitable['STOP'][i]
    duration = stopmet-startmet
    myelapsedbins = np.arange(0.0,duration+args.lcbinsize,args.lcbinsize)
    lc_elapsed_bins = np.append(lc_elapsed_bins,
        cumtime+myelapsedbins)
    lc_met_bins = np.append(lc_met_bins,np.arange(startmet,stopmet+args.lcbinsize,args.lcbinsize))
    mylcduration = myelapsedbins[-1]+args.lcbinsize
    cumtimes.append(cumtime)
    cumtime += mylcduration
gtitable['CUMTIME'] = np.array(cumtimes)

# Define basename for writing output files
bn = path.basename(args.infiles[0]).split('_')[0]
log.info('OBS_ID {0}'.format(etable.meta['OBS_ID']))
if etable.meta['OBS_ID'].startswith('000000'):
    log.info('Overwriting OBS_ID with {0}'.format(bn))
    etable.meta['OBS_ID'] = bn
    etable.meta['OBS_ID'] = bn

if args.basename is None:
    basename = '{0}'.format(bn)
    args.basename = basename
else:
    basename = args.basename

# Overwrite bad OBJECT name, if requested (early FITS all have OBJECT=Crab)
if args.guessobj and args.obsdir:
    # Trim trailing slash, if needed
    if args.obsdir[-1] == '/':
        args.obsdir = args.obsdir[:-1]
    objname = path.basename(args.obsdir)[11:]
    log.info('Guessing Object name {0}'.format(objname))
    etable.meta['OBJECT'] = objname
    etable.meta['OBJECT'] = objname

# getting the overshoot and undershoot rate from HK files.  Times are hkmet
log.info('Getting overshoot and undershoot rates')
if len(args.hkfiles) > 0:
    log.info('Reading '+hkfiles[0])
    hdulist = pyfits.open(hkfiles[0])
    td = hdulist[1].data
    hkmet = td['TIME']
    log.info("HK MET Range {0} to {1} (Span = {2:.1f} seconds)".format(hkmet.min(),
        hkmet.max(),hkmet.max()-hkmet.min()))
    overshootrate = td['MPU_OVER_COUNT'].sum(axis=1)
    undershootrate = td['MPU_UNDER_COUNT'].sum(axis=1)
    reset_rates = td['MPU_UNDER_COUNT'].sum(axis=0)
    for fn in hkfiles[1:]:
	log.info('Reading '+fn)
	hdulist = pyfits.open(fn)
	mytd = hdulist[1].data
	mymet = mytd['TIME']
	myovershootrate = mytd['MPU_OVER_COUNT'].sum(axis=1)
	myundershootrate = mytd['MPU_UNDER_COUNT'].sum(axis=1)
	myreset = mytd['MPU_UNDER_COUNT'].sum(axis=0)
        if not np.all(mymet == hkmet):
	    log.error('TIME axes are not compatible')
	    sys.exit(1)
	overshootrate += myovershootrate
	undershootrate += myundershootrate
	reset_rates= np.append(reset_rates,myreset)
    del hdulist

    if args.filtou:
        b1 = np.where(etable['EVENT_FLAGS'][:,FLAG_UNDERSHOOT] == True)
        b2 = np.where(etable['EVENT_FLAGS'][:,FLAG_OVERSHOOT] == True)
        b3 = np.intersect1d(b1,b2)
        time = etable['MET'][b3]
        counts,a = np.histogram(time,np.append(hkmet,(hkmet[-1]+hkmet[1]-hkmet[0])))
        overshootrate = overshootrate - counts
        del b1, b2, b3, time, counts, a
    # Write overshoot and undershoot rates to file for filtering
    if args.writeovershoot:
        tcol = pyfits.Column(name='TIME',unit='S',array=hkmet,format='D')
        ocol = pyfits.Column(name='OVERSHOOT',array=overshootrate,format='D')
        ucol = pyfits.Column(name='UNDERSHOOT',array=undershootrate,format='D')
        ovhdu = pyfits.BinTableHDU.from_columns([tcol,ocol,ucol], name='HKP')
        ovhdu.writeto("{0}.ovs".format(basename),overwrite=True,checksum=True)
else:
    hkmet = None
    overshootrate=None
    undershootrate = None
    nresets = reset_rate(etable, IDS)
    reset_rates = nresets/etable.meta['EXPOSURE']

filttable = etable

# Background plots are diagnostics for background rates and filtering
if args.ratio:
    figure4 = ratio_plots(etable, overshootrate, gtitable, args, hkmet, undershootrate, mktable)
    figure4.set_size_inches(16,12)
    if args.save:
        log.info('Writing ratio plot {0}'.format(basename))
        figure4.savefig('{0}_bkg.png'.format(basename), dpi = 100)

# Engineering plots are reset rates, count rates by detector, and deadtime
if args.eng:
    figure1 = eng_plots(etable, args, reset_rates, filttable)
    figure1.set_size_inches(16,12)
    if args.save:
    	log.info('Writing eng plot {0}'.format(basename))
    	if args.filtall:
        	figure1.savefig('{0}_eng_clean_{1:.1f}-{2:.1f}keV.png'.format(basename,args.emin,args.emax), dpi = 100)
    	else:
        	figure1.savefig('{0}_eng.png'.format(basename), dpi = 100)

#DELETING ETABLE HERE. USE FILTTABLE FROM NOW ON

# Science plot is light curve, spectrum, pulse profile, and PHA ratio plot (or poweer spectrum)
if args.sci:
    # Make science plots using filtered events
    figure2 = sci_plots(filttable, gtitable, args)
    figure2.set_size_inches(16,12)
    if args.save:
    	log.info('Writing sci plot {0}'.format(basename))
    	if args.filtall:
        	figure2.savefig('{0}_sci_clean_{1:.1f}-{2:.1f}keV.png'.format(basename,args.emin,args.emax), dpi = 100)
    	else:
        	figure2.savefig('{0}_sci.png'.format(basename), dpi = 100)

# Map plot is overshoot and undershoot rates on maps
if args.map:
    log.info("I'M THE MAP I'M THE MAP I'M THE MAAAAP")
    figure3 = cartography(hkmet, overshootrate, args, undershootrate, filttable)

    if args.save:
        log.info('Writing MAP {0}'.format(basename))
        figure3.savefig('{0}_map.png'.format(basename), dpi = 100)

# Show all plots at the end, if not saving
if not args.save:
    log.info('Showing plots...')
    plt.show()
