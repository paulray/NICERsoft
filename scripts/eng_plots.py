import matplotlib.pyplot as plot
import numpy as np
import pylab
import matplotlib.gridspec as gridspec
import argparse
from functionality import *
from nicer.values import *

def eng_plots(etable):
    #GRID SET UP
    figure1 = plot.figure(figsize = (11.0, 8.5))
    sci_grid = gridspec.GridSpec(5,6)

    exposure = etable.meta['EXPOSURE']

    #Total Count Histogram
    ax_counts = plot.subplot(sci_grid[2:5, :3])
    ax_rate = ax_counts.twinx()
    num_events = plot_total_count_hist(etable, ax_rate, ax_counts)

    #Detector Array Chart Thing
    ax_map = plot.subplot(sci_grid[2:5, 3:6])
    plot_detector_chart(etable, num_events, ax_map)

    #Deadtime
    dead = plot.subplot(sci_grid[:2, 4:6])
    dead_plot = plot_deadtime(etable)

    #RESET RATE PER DETECTOR
    log.info('Computing reset rates')
    nresets = reset_rate(etable, IDS)
    reset_rates = nresets

    plot.subplot(sci_grid[:2, 2:4])
    plot_resetrate(IDS, reset_rates)

    #Making the plot all nice and stuff
    plot.subplots_adjust(left = .07, right = .99, bottom = .05, top = .9, wspace = .7, hspace = .8)

    figure1.suptitle('Object: {0} at {1}'.format(etable.meta['OBJECT'],etable.meta['DATE-OBS']),
        fontsize=18)

    average = "Mean events per detector is {0:.2f}".format(num_events.mean())

    text1 = plot.figtext(.02, .9, average, fontsize = 12.5)
    text3 = plot.figtext(.02, .8,
      "Mean reset rate is {0:.2f}/s".format(reset_rates.mean()),
      fontsize = 12.5)

    return figure1
