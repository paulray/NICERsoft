import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from astropy import log

from functionality import *

def sci_plots(etable, lclog):
    #GRID SET UP
    figure2 = plt.figure(figsize = (11, 8.5), facecolor = 'white')
    sci_grid = gridspec.GridSpec(5,7)
    plt.style.use('grayscale')



    #Light Curve
    log.info('Building light curve')
    ax_count = plt.subplot(sci_grid[3:5,:5])
    ax_rate = ax_count.twinx()
    meanrate = plot_light_curve(etable,ax_rate, ax_count, lclog, binsize=1.0)

    #Fast / Slow (Slow x, Fast y)
    log.info('Building fast/slow')
    plt.subplot(sci_grid[1:3,2:5])
    plot_slowfast(etable)

    #Energy Spectrum
    log.info('Building energy spectrum')
    plt.subplot(sci_grid[1:3,:2])
    plot_energy_spec(etable)

    #Power Spectrum
    log.info('Building power spectrum')
    fourier = plt.subplot(sci_grid[3:5,5:7])
    #plot_fft_of_power(etable)

    #PULSE PROFILE
    log.info('Building pulse profile')
    pulse = plt.subplot(sci_grid[1:3,5:7])
##    pulse = pulse_profile(data1, pulse, event_flags)

    #Making the plot all nice and stuff
    plt.subplots_adjust(left = .07, right = .99, bottom = .05, top = .9, wspace = .8, hspace = .8)

    figure2.suptitle('ObsID {0}: {1} at {2}'.format(etable.meta['OBS_ID'],
            etable.meta['OBJECT'],etable.meta['DATE-OBS']),
            fontsize=18)

    # Add text info here:
    plt.figtext(.07, .9, 'Mean count rate {0:.3f} c/s'.format(meanrate), fontsize = 10)
    plt.figtext(0.07, 0.8, etable.meta['FILT_STR'], fontsize=10)

    return figure2
