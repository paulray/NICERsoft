import matplotlib.pyplot as plot
import numpy as np
import pylab
import matplotlib.gridspec as gridspec
import argparse
from functionality import *
from nicer.values import *

def eng_plots(etable, reset_rates):
    #GRID SET UP
    figure1 = plot.figure(figsize = (8.0, 5.0))
    sci_grid = gridspec.GridSpec(5,6)

    #Total Count Histogram
    total_counts = plot.subplot(sci_grid[2:5, :3]) #totcount hist
    countrate = total_counts.twinx()
    #tc, IDS, num_events, rate = plot_total_count_hist(data1, countrate, total_counts)
    num_events = 10

    #Detector Array Chart Thing
    detector_map = plot.subplot(sci_grid[2:5, 3:6])
    #detector_map = plot_detector_chart(data1, IDS, num_events, sci_grid, detector_map)

    #Deadtime
    dead = plot.subplot(sci_grid[:2, 4:6])
    #dead_plot = plot_deadtime(data1, sci_grid, dead)

    #RESET RATE PER DETECTOR
    reset = plot.subplot(sci_grid[:2, 2:4])
    reset_plot = plot_resetrate(IDS, reset_rates)

    #Making the plot all nice and stuff
    plot.subplots_adjust(left = .07, right = .99, bottom = .05, top = .9, wspace = .7, hspace = .8)

    average = "The average total number of events per detector is " + str(round(np.mean(num_events),2)) + "."
    stddev = "Standard deviation for total event count is " + str(round(np.std(num_events),2)) + "."
    reset = "Average detector reset rate is "  + " Hz per detector."

    figure1.suptitle('Target ID: ' + str(etable.meta['OBJECT']) + ', Observed on ' + str(etable.meta['DATE-OBS'][0:9]) + ' at ' + str(etable.meta['DATE-OBS'][11:19]), fontsize = 25)

    text1 = plot.figtext(.02, .9, average, fontsize = 12.5)
    text2 = plot.figtext(.02, .85, stddev, fontsize = 12.5)
    text3 = plot.figtext(.02, .8,
      "Average detector reset rate is {0:.2f} Hz per detector.".format(reset_rates.mean()),
      fontsize = 12.5)

    return figure1
