import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys
from astropy import log
from astropy.table import Table, vstack
import astropy.io.fits as pyfits
import astropy.units as u
from astropy.time import Time
from nicer.values import *

from functionality import *
from sci_plots import sci_plots
from eng_plots import eng_plots

parser = argparse.ArgumentParser(description = "Plot the NICER data nicely.")
parser.add_argument("infiles", help="Input files", nargs='+')
parser.add_argument("-s", "--save", help = "Save plots to file", action = "store_true")
parser.add_argument("--sci", help = "Makes some nice science plots", action = "store_true")
parser.add_argument("--eng", help = "Makes some nice engineering plots", action = "store_true")
parser.add_argument("--emin", help="Minimum energy (keV) to keep", default=0.2, type=float)
parser.add_argument("--emax", help="Minimum energy (keV) to keep", default=12.0, type=float)
args = parser.parse_args()

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
# Add Time column with astropy Time for ease of use
log.info('Adding time column')
etime = etable.columns['MET'] + MET0
etable['T'] = etime

# If there are no PI columns, add them with approximate calibration
if not ('PI' in etable.colnames):
    log.info('Adding PI')
    pi = PI(etable)

# Note: To access event flags, use etable['EVENT_FLAGS'][:,B], where B is
# the bit number for the flag (e.g. FLAG_UNDERSHOOT)

#data1, event_flags, info = smush(args.infiles)
exposure = etable.meta['EXPOSURE'] * u.s

log.info('Computing reset rates')
nresets = reset_rate(etable, IDS)
reset_rates = nresets/exposure

log.info('Filtering...')
# Apply filters for good events
b1 = etable['EVENT_FLAGS'][:,FLAG_SWTRIG] == False
b2 = etable['EVENT_FLAGS'][:,FLAG_UNDERSHOOT] == False
b3 = etable['EVENT_FLAGS'][:,FLAG_OVERSHOOT] == False
b4 = np.logical_and(etable['PI'] > args.emin*PI_TO_KEV, etable['PI']< args.emax*PI_TO_KEV)
idx = np.where(b1 & b2 & b3 )[0]
del b1, b2, b3, b4
filttable = etable[idx]

log.info("Filtering cut {0} events to {1} ({2:.2f}%)".format(len(etable),
    len(filttable), 100*len(filttable)/len(etable)))

#Making the plots
if args.eng:
    figure1 = eng_plots(etable, reset_rates)
    figure1.set_size_inches(16,12)
    if args.save:
        figure1.savefig(str(etable.meta["OBJECT"])+ str(etable.meta["DATE-OBS"]) +  "_ENG_.png", dpi = 100)
    else:
        plt.show()

if args.sci:
    figure2 = sci_plots(etable)
    figure2.set_size_inches(11,8.5)
    if args.save:
        figure2.savefig(str(etable.meta["OBJECT"])+ str(etable.meta["DATE-OBS"]) +  "_SCI_.png", dpi = 100)
    else:
        plt.show()
