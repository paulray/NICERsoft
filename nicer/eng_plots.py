from __future__ import (print_function, division, unicode_literals, absolute_import)
import matplotlib.pyplot as plot
import numpy as np
import pylab
import matplotlib.gridspec as gridspec
import argparse
from nicer.plotutils import *
from nicer.values import *

#def eng_plots(etable, args, reset_rates, filttable, gtitable):
def eng_plots(etable, args, mktable, filttable, gtitable):
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
    if mktable is not None:
        log.info('Computing reset rates')
        plot.subplot(sci_grid[:2, 2:4])
        reset_rates = calc_nresets(mktable, IDS)/ etable.meta['EXPOSURE']
        if reset_rates is not None:
            plot_resetrate(IDS, reset_rates)
    else:
        reset_rates = None 

    # Hot detector
    num_events, colors = hist_use(filttable)
    max_id = IDS[np.where(num_events == num_events.max())[0]][0]
    log.info('max_id {0}'.format(max_id))
    # Plot spectrum and lightcurve of hottest detector
    plot.subplot(sci_grid[:1,:2])
    idx = np.where(filttable['DET_ID']==max_id)[0]
    plot_energy_spec(filttable[idx],binscale=4.0)
    plot.title('PI Spectrum of DET_ID {0}'.format(max_id))
    plot.subplot(sci_grid[1:2,:2])
    meanrate, a = plot_light_curve(filttable[idx], args.lclog, gtitable, binsize=args.lcbinsize)
    plot.title('Light Curve of DET_ID {0}'.format(max_id))
    plot.xlabel('Time (s)')


    #Making the plot all nice and stuff
    plot.subplots_adjust(left = .07, right = .99, bottom = .05, top = .9, wspace = .7, hspace = .8)

    figure1.suptitle('Object: {0} at {1}'.format(etable.meta['OBJECT'],etable.meta['DATE-OBS'].replace('T', ' at ')),
        fontsize=18)


    #plot.figtext(0.02, 0.9, etable.meta['FILT_STR'], fontsize=10)
    #average = "Mean events per detector is {0:.2f}".format(num_events.mean())
    #ext1 = plot.figtext(.02, .9, average, fontsize = 12.5)
    if reset_rates is not None:
        plot.figtext(.42, .87,
            "Mean reset rate is {0:.2f}/s".format(reset_rates.mean()),
            fontsize = 12.5)
    if args.mask:
        plot.figtext(.55, .05, 'IDS {0} are masked'.format(args.mask), fontsize=10)

    return figure1



def plot_all_spectra(etable, args, filttable, gtitable):
    #GRID SET UP
    fig_all_spec = plot.figure(figsize = (11, 8.5), facecolor = 'white')    
    ncols = 8
    nbDET = len(IDS)
    spec_grid = gridspec.GridSpec(nbDET // ncols , ncols)

    # Plot spectrum of all detectors
    for i, detid in enumerate(IDS):
        row = (i // ncols)
        col = i % ncols
        idx = np.where(filttable['DET_ID']==detid)[0]
        plot.subplot(spec_grid[row, col])
        if col==0:
            plot_energy_spec(filttable[idx],binscale=4.0,plot_pos='left')
        if row == ((nbDET // ncols)-1):
            plot_energy_spec(filttable[idx],binscale=4.0,plot_pos='bottom')
        if (col==0) and (row == ((nbDET // ncols)-1) ):
            plot_energy_spec(filttable[idx],binscale=4.0,plot_pos='corner')
        if (col!=0) and (row != ((nbDET // ncols)-1) ):
            plot_energy_spec(filttable[idx],binscale=4.0,plot_pos='center')
        plot.title('DET_ID {0}'.format(detid))
        
    plot.subplots_adjust(left = .04, right = .99, bottom = .05, top = .94, wspace = .2, hspace = .2)
    fig_all_spec.suptitle('Object: {0} at {1}'.format(etable.meta['OBJECT'],etable.meta['DATE-OBS'].replace('T', ' at ')),fontsize=18)
    
    return fig_all_spec


def plot_all_lc(etable, args, filttable, gtitable):
    #GRID SET UP
    fig_all_lc = plot.figure(figsize = (11, 8.5), facecolor = 'white')    
    ncols = 8
    nbDET = len(IDS)
    lc_grid = gridspec.GridSpec(nbDET // ncols , ncols)
    lc_grid.update(hspace=0.5)

    # Plot lightcurve of all detectors
    for i, detid in enumerate(IDS):
        row = (i // ncols)
        col = i % ncols
        idx = np.where(filttable['DET_ID']==detid)[0]
        plot.subplot(lc_grid[row, col])
        if col==0:
            meanrate, a = plot_light_curve(filttable[idx], args.lclog, gtitable, binsize=5*args.lcbinsize,plot_pos='left')

        if row == ((nbDET // ncols)-1):
            meanrate, a = plot_light_curve(filttable[idx], args.lclog, gtitable, binsize=5*args.lcbinsize,plot_pos='bottom')
            plot.xlabel('Time (s)')
            
        if (col==0) and (row == ((nbDET // ncols)-1) ):
            meanrate, a = plot_light_curve(filttable[idx], args.lclog, gtitable, binsize=5*args.lcbinsize,plot_pos='corner')
            plot.xlabel('Time (s)')
            
        if (col!=0) and (row != ((nbDET // ncols)-1) ):
            meanrate, a = plot_light_curve(filttable[idx], args.lclog, gtitable, binsize=5*args.lcbinsize,plot_pos='center')

        plot.title('DET_ID {0}'.format(detid))

    plot.subplots_adjust(left = .04, right = .99, bottom = .05, top = .94, wspace = .2, hspace = .0)
    fig_all_lc.suptitle('Object: {0} at {1}'.format(etable.meta['OBJECT'],etable.meta['DATE-OBS'].replace('T', ' at ')),fontsize=18)
    
    return fig_all_lc
