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
parser.add_argument("-c", "--sci", help = "Makes some nice science plots", action = "store_true")
parser.add_argument("-e", "--eng", help = "Makes some nice engineering plots", action = "store_true")
parser.add_argument("-p", "--pi", help = "Use approximate calibration to compute PI", action = "store_true")
args = parser.parse_args()

log.info('Reading files')
tlist = []
hdr = pyfits.getheader(args.infiles[0])
for fn in args.infiles:
    tlist.append(Table.read(fn,hdu=1))

etable = vstack(tlist,metadata_conflicts='silent')
# Change TIME column name to MET to reflect what it really is
etable.columns['TIME'].name = 'MET'
etime = etable.columns['MET'] + MET0
# Add Time column with astropy Time for computations
etable['T'] = etime

# Note: To access event flags, use etable['EVENT_FLAGS'][:,B], where B is
# the bit number for the flag (e.g. FLAG_UNDERSHOOT)

#data1, event_flags, info = smush(args.infiles)
exposure = etable.meta['EXPOSURE'] * u.s

log.info('Computing reset rates')
nresets = reset_rate(etable, IDS)
reset_rates = nresets/exposure

#Making the plots
if args.eng:
    figure1 = eng_plots(etable, reset_rates)
    figure1.set_size_inches(16,12)
    if args.save:
        figure1.savefig(str(etable.meta["OBJECT"])+ str(etable.meta["DATE-OBS"]) +  "_ENG_.png", dpi = 100)
    else:
        plt.show()

if args.sci:
    IDS, num_events, stdev, colors = hist_use(etable)
    figure2 = sci_plots(etable)
    figure2.set_size_inches(16,12)
    if args.save:
        figure2.savefig(str(etable.meta["OBJECT"])+ str(etable.meta["DATE-OBS"]) +  "_SCI_.png", dpi = 100)
    else:
        plt.show()
