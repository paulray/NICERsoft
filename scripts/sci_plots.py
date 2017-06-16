import matplotlib.pyplot as plot
import numpy as np
import matplotlib.gridspec as gridspec
from astropy import log

from functionality import *

def sci_plots(data1, event_flags, info, num_events, avg_rate, PI_flag, ID_rates):
    #GRID SET UP
    figure2 = plot.figure(figsize = (8.0, 5.0), facecolor = 'white')
    sci_grid = gridspec.GridSpec(5,7)

    plot.style.use('grayscale')
    #Light Curve
    log.info('Building light curve')
    light_curve_plot = plot.subplot(sci_grid[1:5,:2])
    light_curve_plot = plot_light_curve(data1, sci_grid, light_curve_plot, event_flags)
    
    #Fast / Slow (Slow x, Fast y)
    log.info('Building fast/slow')
    fastslow_ratio = plot.subplot(sci_grid[1:3,2:5])
    fast, slow, total = slow_fast(data1)
    fastslow_ratio = plot_slowx_fasty(data1, sci_grid, fastslow_ratio, fast, slow, total)

    #Energy Spectrum
    
    log.info('Building spectrum')
    power_spec = plot.subplot(sci_grid[3:5,2:5])
    power_spec = plot_power_spec(data1, sci_grid, power_spec, event_flags, PI_flag)

    #Power SPpectrum
    log.info('Building power spectrum')
    fourier = plot.subplot(sci_grid[3:5,5:7])
    time = no_more_resets(data1[0], event_flags)
    #power_spec = plot_fft_of_power(time, data1)

    #PULSE PROFILE
    log.info('Building pulse profile')
    pulse = plot.subplot(sci_grid[1:3,5:7])
##    pulse = pulse_profile(data1, pulse, event_flags)
    
    #Making the plot all nice and stuff
    plot.subplots_adjust(left = .07, right = .99, bottom = .05, top = .9, wspace = .7, hspace = .8)

    #list of times without reset events in it
    log.info('Running light_curve')
    sums, count_rate = light_curve(data1, event_flags)
    average = "Mean detector count rate is " + str(round(count_rate,4)) + "."
    stddev = "Standard deviation for event number is " + str(round(np.std(ID_rates[1]),2)) + "."

    figure2.suptitle('Target ID: ' + str(info['OBJECT']) + ', Observed on ' + str(info['DATE-OBS'][0:9]) + ' at ' + str(info['DATE-OBS'][11:19]), fontsize = 25)
    
    text1 = plot.figtext(.07, .9, average, fontsize = 13)
    text2 = plot.figtext(.07, .85, stddev, fontsize = 13)

    return figure2

