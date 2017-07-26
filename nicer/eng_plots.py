from __future__ import (print_function, division, unicode_literals, absolute_import)
import matplotlib.pyplot as plot
import numpy as np
import pylab
import matplotlib.gridspec as gridspec
import argparse
from nicer.plotutils import *
from nicer.values import *

def eng_plots(etable, args, reset_rates, filttable):
    #GRID SET UP
    figure1 = plot.figure(figsize = (11, 8.5), facecolor = 'white')
    sci_grid = gridspec.GridSpec(5,6)
    exposure = etable.meta['EXPOSURE']

    #Total Count Histogram
    ax_counts = plot.subplot(sci_grid[2:5, :3])
    ax_rate = ax_counts.twinx()
    num_events = plot_total_count_hist(filttable, ax_rate, ax_counts)

    #Detector Array Chart Thing
    ax_map = plot.subplot(sci_grid[2:5, 3:6])
    plot_detector_chart(filttable, num_events, ax_map)

    #Deadtime
    dead = plot.subplot(sci_grid[:2, 4:6])
    dead_plot = plot_deadtime(filttable)

    #RESET RATE PER DETECTOR
    log.info('Computing reset rates')
    plot.subplot(sci_grid[:2, 2:4])
    if reset_rates is not None:
        plot_resetrate(IDS, reset_rates)

    #Making the plot all nice and stuff
    plot.subplots_adjust(left = .07, right = .99, bottom = .05, top = .9, wspace = .7, hspace = .8)

    figure1.suptitle('Object: {0} at {1}'.format(etable.meta['OBJECT'],etable.meta['DATE-OBS'].replace('T', ' at ')),
        fontsize=18)


    #plot.figtext(0.02, 0.9, etable.meta['FILT_STR'], fontsize=10)
    #average = "Mean events per detector is {0:.2f}".format(num_events.mean())
    #ext1 = plot.figtext(.02, .9, average, fontsize = 12.5)
    plot.figtext(.02, .8,
      "Mean reset rate is {0:.2f}/s".format(reset_rates.mean()),
      fontsize = 12.5)
    if args.mask:
        plot.figtext(.07, .81, 'IDS {0} are masked'.format(args.mask), fontsize=10)

    return figure1
