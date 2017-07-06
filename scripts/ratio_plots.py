import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from astropy import log

from functionality import *

def ratio_plots(etable, overshootrate, gtitable, elapsedbins, metbins, args, hkmet):
	figure = plt.figure(figsize = (11, 8.5), facecolor = 'white')
   	ratio_grid = gridspec.GridSpec(4,4)

	#Fast / Slow (Slow x, Fast y)
	log.info('Building fast/slow subplot')
	plt.subplot(ratio_grid[2:4,2:4])
	log.info('Building actual slow fast data')
	#plot_slowfast(etable)

	#Overshoot rate plot
	plt.subplot(ratio_grid[2:3,:2])
	plot_overshoot(etable, overshootrate, gtitable, elapsedbins, args, hkmet)

        #Plot of Pointing
	plt.subplot(ratio_grid[3:4,:2])

	#Plot of ISS LatLong vs Time
	plt.subplot(ratio_grid[1:2,2:4])

	#Plot of events rejected
	log.info('Building overshoot plot')
	plt.subplot(ratio_grid[1:2,:2])

	
	
	figure.suptitle('Object: {0} at {1}'.format(etable.meta['OBJECT'],etable.meta['DATE-OBS']),
        fontsize=18)
	plt.subplots_adjust(left = .07, right = .99, bottom = .05, top = .9, wspace = .4, hspace = .4)

	return figure
