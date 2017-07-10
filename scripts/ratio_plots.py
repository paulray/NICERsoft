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
	ratio = np.array(etable['PHA'],dtype=np.float)/np.array(etable['PHA_FAST'],dtype=np.float)
	badtable = etable[np.where(ratio > 1.4)]
        plot_light_curve(badtable, args.lclog, overshootrate, gtitable, binsize=1.0)

	#Overshoot rate plot -- use --lclog to make it a log y axis
	log.info('Building overshoot plot')
	plt.subplot(ratio_grid[2:4,:4])
	plot_overshoot(etable, overshootrate, gtitable, args, hkmet)
        plot_SAA(mktable, gtitable, overshootrate)

	#Plot of undershoot rate
	log.info('Building undershoot plot')
	plt.subplot(ratio_grid[4:6,:4])
	plot_undershoot(etable, undershootrate, gtitable, args, hkmet, mktable)
        
	#Plot of Sun / Moon
	log.info('Building Sun / Moon / Earth angle Plot')
	plt.subplot(ratio_grid[6:8,:4])	
	plot_angles(mktable, gtitable)	

	#Pointing plot
	plt.subplot(ratio_grid[8:10,:4])
	plot_pointing(mktable, gtitable)

	#Plot of LAT / Long
	plt.subplot(ratio_grid[10:12,:4])
	plot_latlon(mktable, gtitable)


	figure.suptitle('Object: {0} at {1}'.format(etable.meta['OBJECT'],etable.meta['DATE-OBS']),
        fontsize=18)
	#plt.subplots_adjust(left = .07, right = .99, bottom = .05, top = .9, wspace = .95, hspace = .95)
	plt.tight_layout()
	return figure

