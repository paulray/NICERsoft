import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from astropy import log

from plotutils import *

def bkg_plots(etable, overshootrate, gtitable, args, hkmet, undershootrate, mktable, bothrate):
    figure = plt.figure(figsize = (8.5, 11), facecolor = 'white')
    bkg_grid = gridspec.GridSpec(25,4)

    #Lightcurve of Rejected events
    log.info('Building Rejected Event Light curve')
    plt.subplot(bkg_grid[1:5,:4])
    etable = etable[np.logical_and(etable['EVENT_FLAGS'][:,FLAG_SLOW],etable['EVENT_FLAGS'][:,FLAG_FAST])]
    ratio = np.array(etable['PHA'],dtype=np.float)/np.array(etable['PHA_FAST'],dtype=np.float)
    badtable = etable[np.where(ratio > 1.4)[0]]
    r, lc = plot_light_curve(badtable, args.lclog, gtitable, binsize=16.0)
    plt.legend(handles = [lc], loc = 2)
    #Overshoot rate plot -- use --lclog to make it a log y axis
    log.info('Building overshoot plot')
    plt.subplot(bkg_grid[5:9,:4])
    if overshootrate is not None:
        plot_overshoot(etable, overshootrate, gtitable, args, hkmet, bothrate)
        plot_SAA(mktable, gtitable, overshootrate)

    #Plot of undershoot rate
    log.info('Building undershoot plot')
    plt.subplot(bkg_grid[9:13,:4])
    if undershootrate is not None:
        plot_undershoot(etable, undershootrate, gtitable, args, hkmet, mktable)

    #Plot of Sun / Moon
    log.info('Building Sun / Moon / Earth angle Plot')
    plt.subplot(bkg_grid[13:17,:4])
    plot_angles(mktable, gtitable)

    #Pointing plot
    plt.subplot(bkg_grid[17:21,:4])
    plot_pointing(mktable, gtitable)

    #Plot of LAT / Long
    plt.subplot(bkg_grid[21:25,:4])
    plot_latlon(mktable, gtitable)


    figure.suptitle('Object: {0} at {1}'.format(etable.meta['OBJECT'],etable.meta['DATE-OBS'].replace('T', ' at ')),
                    fontsize=18)
    #plt.subplots_adjust(left = .07, right = .99, bottom = .05, top = .9, wspace = .95, hspace = .95)
    plt.tight_layout()
    return figure
