import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from astropy import log

from functionality import *

def ratio_plots(etable, overshootrate, gtitable, args, hkmet, undershootrate, mktable):
	figure = plt.figure(figsize = (8.5, 11), facecolor = 'white')
   	ratio_grid = gridspec.GridSpec(12,4)

	#Lightcurve of Rejected events
	log.info('Building Rejected Event Light curve')
	plt.subplot(ratio_grid[0:2,:4])


	#Overshoot rate plot -- use --lclog to make it a log y axis
	log.info('Building overshoot plot')
	plt.subplot(ratio_grid[2:4,:4])
	plot_overshoot(etable, overshootrate, gtitable, args, hkmet)

	#Plot of undershoot rate
	log.info('Building undershoot plot')
	plt.subplot(ratio_grid[4:6,:4])
	undershoot, etime=plot_undershoot(etable, undershootrate, gtitable, args, hkmet)
	plot_sunshine(args, undershoot, etime, gtitable, hkmet)
        
	#Plot of Sun / Moon
	log.info('Building Sun / Moon / Earth angle Plot')
	plt.subplot(ratio_grid[6:8,:4])	
	plot_angles(mktable, gtitable)	

	#Plot of ISS LatLong vs Time
	plt.subplot(ratio_grid[8:10,:4])

	#Plot of events rejected
	plt.subplot(ratio_grid[10:12,:4])



	figure.suptitle('Object: {0} at {1}'.format(etable.meta['OBJECT'],etable.meta['DATE-OBS']),
        fontsize=18)
	plt.subplots_adjust(left = .07, right = .99, bottom = .05, top = .9, wspace = .9, hspace = .95)

	return figure

