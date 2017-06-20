#!/usr/bin/env python
import sys
# Hack to add this to pythonpath
sys.path.append('/Users/paulr/src/NICERsoft')
import matplotlib.pyplot as plt
import numpy as np
import argparse
from astropy import log
from astropy.table import Table, vstack
import astropy.io.fits as pyfits
import astropy.units as u
from astropy.time import Time
#CHANGE THIS BEFORE COMMITTING
from nicer.values import *

from functionality import *
from sci_plots import sci_plots
from eng_plots import eng_plots


parser = argparse.ArgumentParser(description = "Plot the NICER data nicely.")
parser.add_argument("infiles", help="Input files", nargs='+')
parser.add_argument("-s", "--save", help = "Save plots to file", action = "store_true")
parser.add_argument("--sci", help = "Makes some nice science plots", action = "store_true")
parser.add_argument("--eng", help = "Makes some nice engineering plots", action = "store_true")
parser.add_argument("--filtswtrig", help = "Filter SW TRIG events", action = "store_true")
parser.add_argument("--filtovershoot", help = "Filter OVERSHOOT events", action = "store_true")
parser.add_argument("--filtundershoot", help = "Filter UNDERSHOOT events", action = "store_true")
parser.add_argument("--filtall", help = "Filter SWTRIG, UNDERSHOOT and OVERSHOOT events", action = "store_true")
parser.add_argument("--emin", help="Minimum energy (keV) to keep", default=-1.0, type=float)
parser.add_argument("--emax", help="Minimum energy (keV) to keep", default=-1.0, type=float)
parser.add_argument("--tskip", help="Seconds to skip at beginning of data", default=0.0, type=float)
parser.add_argument("--lcbinsize", help="Light curve bin size (s)", default=0.5, type=float)
parser.add_argument("--pi", help="Force use of internal PHA to PI conversion", action='store_true')
parser.add_argument("--basename", help="Basename for output plots", default=None)
parser.add_argument("--lclog", help = "make light curve log axis", action = "store_true")
parser.add_argument("--foldfreq", help="Make pulse profile by folding at a fixed freq (Hz)", 
    default=0.0,type=float)
args = parser.parse_args()

if args.filtall:
    args.filtswtrig=True
    args.filtovershoot=True
    args.filtundershoot=True
# Load files and build events table
log.info('Reading files')
tlist = []
for fn in args.infiles:
    log.info('Reading file {0}'.format(fn))
    tlist.append(Table.read(fn,hdu=1))

log.info('Concatenating files')
etable = vstack(tlist,metadata_conflicts='silent')
del tlist
# Change TIME column name to MET to reflect what it really is
etable.columns['TIME'].name = 'MET'

# Sort table by MET
etable.sort('MET')
log.info("MET Range : {0} to {1}".format(etable['MET'].min(), etable['MET'].max()))

# Hack to trim first chunk of data
if args.tskip > 0.0:
    t0 = etable['MET'].min()
    etable = etable[etable['MET']>t0+args.tskip]
    # Correct exposure (approximately)
    etable.meta['EXPOSURE'] -= args.tskip

# Add Time column with astropy Time for ease of use
log.info('Adding time column')
# This should really be done the FITS way using MJDREF etc...
# For now, just using MET0
etime = etable.columns['MET'] + MET0
etable['T'] = etime

# If there are no PI columns, add them with approximate calibration
if args.pi or not ('PI' in etable.colnames):
    log.info('Adding PI')
    calfile = 'data/gaincal_linear.txt'
    pi = calc_pi(etable,calfile)
    etable['PI'] = pi

# Note: To access event flags, use etable['EVENT_FLAGS'][:,B], where B is
# the bit number for the flag (e.g. FLAG_UNDERSHOOT)

exposure = etable.meta['EXPOSURE']
log.info('Exposure : {0:.2f}'.format(exposure))
log.info('Filtering...')

filt_str = 'Filter: {0:.2f} < E < {1:.2f} keV'.format(args.emin,args.emax)
if args.emin >= 0:
    b4 = etable['PI'] > args.emin/PI_TO_KEV
else:
    b4 = np.ones_like(etable['PI'],dtype=np.bool)
if args.emax >= 0:
    b4 = np.logical_and(b4, etable['PI']< args.emax/PI_TO_KEV)

# Apply filters for good events
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

log.info("Filtering cut {0} events to {1} ({2:.2f}%)".format(len(etable),
    len(filttable), 100*len(filttable)/len(etable)))

if args.basename is None:
    basename = '{0}-{1}'.format(etable.meta['OBJECT'],etable.meta['OBS_ID'])
else:
    basename = args.basename

#Making the plots
if args.eng:
    figure1 = eng_plots(filttable)
    figure1.set_size_inches(16,12)
    if args.save:
        log.info('Writing eng plot {0}'.format(basename))
        figure1.savefig('{0}_eng.png'.format(basename), dpi = 100)
    else:
        plt.show()

if args.sci:
    # Make science plots using filtered events
    figure2 = sci_plots(filttable, args.lclog, args.lcbinsize, args.foldfreq)
    figure2.set_size_inches(11,8.5)
    if args.save:
    	log.info('Writing sci plot {0}'.format(basename))
    	figure2.savefig('{0}_sci.png'.format(basename), dpi = 100)
    else:
    	plt.show()
