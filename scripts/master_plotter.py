#!/usr/bin/env python
from __future__ import (print_function, division, unicode_literals, absolute_import)
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
import sys
from nicer.bkg_plots import *
from nicer.fitsutils import *
from InteractiveLC import *
from nicer.NicerFileSet import *

parser = argparse.ArgumentParser(description = "Plot the NICER data nicely.")
parser.add_argument("infiles", help="Input files", nargs='*', default = None)
parser.add_argument("--obsdir",help = "Find alllllllll the files!", default = None)
parser.add_argument("--object", help="Override object name", default=None)
parser.add_argument("--guessobj", help="Try to guess object from directory name", action="store_true")
parser.add_argument("--useftools", help="Use FTOOLS for filter and merge", action="store_true")
parser.add_argument("--mask",help="Mask these IDS", nargs = '*', type=int, default=None)
parser.add_argument("-s", "--save", help = "Save plots to file", action = "store_true")
parser.add_argument("--sci", help = "Makes some nice science plots", action = "store_true")
parser.add_argument("--eng", help = "Makes some nice engineering plots", action = "store_true")
parser.add_argument("--bkg", help = "Display background diagnostic plots", action = 'store_true')
parser.add_argument("--filtswtrig", help = "Filter SW TRIG events", action = "store_true")
parser.add_argument("--filtovershoot", help = "Filter OVERSHOOT events", action = "store_true")
parser.add_argument("--filtundershoot", help = "Filter UNDERSHOOT events", action = "store_true")
parser.add_argument("--filtratio", help="Filter PHA/PHA_FAST ratio (argument is ratio to cut at)", type=float, default=1.4)
parser.add_argument("--filtall", help = "Filter SWTRIG, UNDERSHOOT and OVERSHOOT events", action = "store_true")
parser.add_argument("--emin", help="Minimum energy (keV) to keep", default=-1.0, type=float)
parser.add_argument("--emax", help="Minimum energy (keV) to keep", default=-1.0, type=float)
parser.add_argument("--tskip", help="Seconds to skip at beginning of data", default=0.0, type=float)
parser.add_argument("--lcbinsize", help="Light curve bin size (s)", default=1.0, type=float)
parser.add_argument("--pi", help="Force use of internal PHA to PI conversion", action='store_true')
parser.add_argument("--basename", help="Basename for output plots", default=None)
parser.add_argument("--lclog", help = "make light curve log axis", action = "store_true")
parser.add_argument("--foldfreq", help="Make pulse profile by folding at a fixed freq (Hz)",default=0.0,type=float)
parser.add_argument("--nyquist", help="Nyquist freq for power spectrum (Hz)",default=100.0,type=float)
parser.add_argument("--map", help= "Creates a map with overshoots and undershoots", action = 'store_true')
parser.add_argument("--orb", help="Path to orbit FITS filed", default = None)
parser.add_argument("--par", help="Path to par file", default = None)
parser.add_argument("--sps", help="Path to SPS HK file (_apid0260.hk)",default=None)
parser.add_argument("--powspec",help = "Display power spectrum (replaces ratio plot)", action = 'store_true')
parser.add_argument("--pslog", help = "make power spectrum log axis", action = "store_true")
parser.add_argument("--writeps", help = "write out power spectrum", action = "store_true")
parser.add_argument("--writeovershoot",help="Write summed overshoot rates to FITS file", action='store_true')
parser.add_argument("--applygti",help="Read GTI from provided FITS file", default=None)
parser.add_argument("--eventshootrate",help="Gets over/undershoot rates from the events", action='store_true')
parser.add_argument("--interactive", help= "TEST FOR INTERACTIVE LC", action = 'store_true')
parser.add_argument("--readovs", help = "Filters events with overshoot > input number", nargs = '*', type = float, default = None)
parser.add_argument("--andrea",help = "Help out Andrea!!!", action = 'store_true')
args = parser.parse_args()

#------------------------------Getting the data and concatenating------------------------------
if np.logical_or(args.obsdir is not None, args.infiles is not None):
    if args.obsdir is not None:
        #Create the data structure
        data = NicerFileSet(args)
        etable = data.etable
        gtitable = data.gtitable
        #Some Definitions
        mktable = data.mktable
        hkmet = data.hkmet
        basename = data.basename

    else:
        #Creating the data table for each separate file
        if args.useftools:
            etable = filtandmerge(args.infiles,workdir=None)
        else:
            log.info('Reading files')
            tlist = []
            for fn in args.infiles:
                log.info('Reading file {0}'.format(fn))
                tlist.append(Table.read(fn,hdu=1))
                if len(tlist[0]) > 3000000:
                    log.error('There is too much data to handle. Not processing...')
                    sys.exit(3)
            log.info('Concatenating files')

            if len(tlist) == 1:
                etable = tlist[0]
            else:
                etable = vstack(tlist,metadata_conflicts='silent')
            del tlist

        # Read the GTIs from the first event FITS file
        gtitable = Table.read(args.infiles[0],hdu=2)
        log.info('Got the good times from GTI')
        gtitable['DURATION'] = gtitable['STOP']- gtitable['START']
        # Only keep GTIs longer than 16 seconds
        idx = np.where(gtitable['DURATION']>16.0)[0]
        gtitable = gtitable[idx]
        print(gtitable)

        # Change TIME column name to MET to reflect what it really is
        etable.columns['TIME'].name = 'MET'

        # Update exposure to be sum of GTI durations
        etable.meta['EXPOSURE'] = gtitable['DURATION'].sum()

        log.info('Got the good times from GTI')
        # Sort table by MET
        etable.sort('MET')
        log.info("Event MET Range : {0} to {1}".format(etable['MET'].min(),etable['MET'].max(), etable['MET'].max()-etable['MET'].min()))
        log.info("TSTART {0}  TSTOP {1} (Span {2} seconds)".format(etable.meta['TSTART'],etable.meta['TSTOP'], etable.meta['TSTOP']-etable.meta['TSTART'] ))
        log.info("DATE Range {0} to {1}".format(etable.meta['DATE-OBS'],etable.meta['DATE-END']))

        if args.object is not None:
            etable.meta['OBJECT'] = args.object

        bn = path.basename(args.infiles[0]).split('_')[0]
        log.info('OBS_ID {0}'.format(etable.meta['OBS_ID']))
        if etable.meta['OBS_ID'].startswith('000000'):
            log.info('Overwriting OBS_ID with {0}'.format(bn))
            etable.meta['OBS_ID'] = bn
            etable.meta['OBS_ID'] = bn
        if args.basename is None:
            basename = '{0}'.format(bn)
        else:
            basename = args.basename

        reset_rates = None
else:
    log.warning('You have not specified any files, please input the path to the files you want to see. Exiting.')
    sys.exit()


#---------------------Options for data filtering / Plotting -------------------
if not args.sci and not args.eng and not args.map and not args.bkg and not args.interactive:
    log.warning("No specific plot requested, making all")
    args.sci = True
    args.eng = True
    args.map = True
    args.bkg = True

if args.filtall:
    args.filtswtrig=True
    args.filtovershoot=True
    args.filtundershoot=True
    args.filtratio=1.4

if args.andrea:
    args.sci = True
    args.eng = True
    args.bkg = True
#--------------------Editing / Filtering the event data Options-----------------
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

#filtering out chosen IDS
if args.mask is not None:
    log.info('Masking IDS')
    for id in args.mask:
        etable = etable[np.where(etable['DET_ID'] != id)]
        # If there are no PI columns, add them with approximate calibration
if args.pi or not ('PI' in etable.colnames):
    log.info('Adding PI')
    calfile = path.join(datadir,'gaincal_linear.txt')
    pi = calc_pi(etable,calfile)
    etable['PI'] = pi

if args.applygti is not None:
    g = Table.read(args.applygti)
    log.info('Applying external GTI from {0}'.format(args.applygti))
    g['DURATION'] = g['STOP']-g['START']
    # Only keep GTIs longer than 16 seconds
    g = g[np.where(g['DURATION']>16.0)]
    print(g)
    etable = apply_gti(etable,g)
    # Replacing this GTI does not work. It needs to be ANDed with the existing GTI
    etable.meta['EXPOSURE'] = g['DURATION'].sum()
    gtitable = g
log.info('Exposure : {0:.2f}'.format(etable.meta['EXPOSURE']))


# Set up the light curve bins, so we can have them for building
# light curves of various quantities, like overshoot rate and ratio filtered events
# Hmmm. The lc_elapsed_bins and lc_met_bins are never used, but CUMTIME is
if gtitable is not None:
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
        lc_elapsed_bins = np.append(lc_elapsed_bins,cumtime+myelapsedbins)
        lc_met_bins = np.append(lc_met_bins,np.arange(startmet,stopmet+args.lcbinsize,args.lcbinsize))
        mylcduration = myelapsedbins[-1]+args.lcbinsize
        cumtimes.append(cumtime)
        cumtime += mylcduration
    gtitable['CUMTIME'] = np.array(cumtimes)

# Overwrite bad OBJECT name, if requested (early FITS all have OBJECT=Crab)
if args.guessobj and args.obsdir:
    # Trim trailing slash, if needed
    if args.obsdir[-1] == '/':
        args.obsdir = args.obsdir[:-1]
    objname = path.basename(args.obsdir)[11:]
    log.info('Guessing Object name {0}'.format(objname))
    etable.meta['OBJECT'] = objname
    etable.meta['OBJECT'] = objname

#Getting over/undershoot rate from event data.
if args.eventshootrate:
    eventovershoot = data.eventovershoot
    eventundershoot = data.eventundershoot
    bothshoots = data.bothshoots
else:
    bothshoots = None
    eventundershoot = None
    eventovershoot = None

if args.obsdir is not None:
    hkovershoots = data.hkovershoots
    hkundershoots = data.hkundershoots
    bothshoots = data.eventbothshoots
    reset_rates = data.reset_rates

# Write overshoot and undershoot rates to file for filtering
if args.writeovershoot:

    temptable = etable[np.logical_and(etable['EVENT_FLAGS'][:,FLAG_SLOW],etable['EVENT_FLAGS'][:,FLAG_FAST])]
    ratio = np.array(temptable['PHA'],dtype=np.float)/np.array(temptable['PHA_FAST'],dtype=np.float)
    badtable = temptable[np.where(ratio > args.filtratio)[0]]
    # This is no good.  Needs to be binned on hkmet bins.
    r, badlightcurve = plot_light_curve(badtable, args.lclog, gtitable,
        binsize=16.0, noplot=True)

    badlightcurve = np.array(badlightcurve)
    data.writeovsfile(badlightcurve)
    del r, temptable, ratio

if np.logical_and(args.readovs is not None, args.writeovershoot == True):
    ovsfile = "{0}.ovs".format(basename)
    ovstable = Table.read(ovsfile,hdu=1)
    print(ovstable)



#---------------------------------------------Filting all the data as necessary!---------------------------------------------------------------
log.info('Filtering...')
filt_str = 'Filter: {0:.2f} < E < {1:.2f} keV'.format(args.emin,args.emax)
if args.emin >= 0:
    b4 = etable['PI'] > args.emin/PI_TO_KEV
else:
    b4 = np.ones_like(etable['PI'],dtype=np.bool)
if args.emax >= 0:
    b4 = np.logical_and(b4, etable['PI']< args.emax/PI_TO_KEV)

if args.filtswtrig:
    b1 = etable['EVENT_FLAGS'][:,FLAG_SWTRIG] == False
    filt_str += ", not SWTRIG"
else:
    b1 = np.ones_like(etable['PI'],dtype=np.bool)
if args.filtundershoot:
    b2 = etable['EVENT_FLAGS'][:,FLAG_UNDERSHOOT] == False
    filt_str += ", not UNDERSHOOT"
else:
    b2 = np.ones_like(etable['PI'],dtype=np.bool)
if args.filtovershoot:
    b3 = etable['EVENT_FLAGS'][:,FLAG_OVERSHOOT] == False
    filt_str += ", not OVERSHOOT"
else:
    b3 = np.ones_like(etable['PI'],dtype=np.bool)

idx = np.where(b1 & b2 & b3 & b4)[0]
del b1, b2, b3, b4
filttable = etable[idx]
filttable.meta['FILT_STR'] = filt_str
etable.meta['FILT_STR'] = filt_str


# Add Time column with astropy Time for ease of use and for PINT TOAs
if args.par is not None:
    log.info('Adding time column')
    # This should really be done the FITS way using MJDREF etc...
    # For now, just using MET0
    etime = filttable.columns['MET'] + MET0
    filttable['T'] = etime


#------------------------------------------------------PLOTTING HAPPENS BELOW HERE ------------------------------------------------------
# Background plots are diagnostics for background rates and filtering
if args.bkg:
    if hkmet is None:
        log.error("Can't make background plots without MPU HKP files")
    else:
        if eventovershoot is not None:

            figure4 = bkg_plots(etable, eventovershoot, gtitable, args, hkmet, eventundershoot, mktable, bothshoots)
        else:
            figure4 = bkg_plots(etable, hkovershoots, gtitable, args, hkmet, hkundershoots, mktable, bothshoots)
        figure4.set_size_inches(16,12)
        if args.save:
            log.info('Writing bkg plot {0}'.format(basename))
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
    if eventovershoot is not None:
        figure3 = cartography(hkmet, eventovershoot, args, eventundershoot, filttable, mktable)
    else:
        figure3 = cartography(hkmet, hkovershoots, args, hkundershoots, filttable, mktable)

    if args.save:
        log.info('Writing MAP {0}'.format(basename))
        figure3.savefig('{0}_map.png'.format(basename), dpi = 100)

if args.interactive:
    log.info("Interaction is coming")
    figure4 = plot.figure()
    ILC = InteractiveLC(etable, args.lclog, gtitable, figure4, basename, binsize=1.0)
    ILC.getgoodtimes()
    ILC.writegti()

# Show all plots at the end, if not saving
if not args.save:
    log.info('Showing plots...')
    plt.show()
