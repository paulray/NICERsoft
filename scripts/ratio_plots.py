import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from astropy import log

from functionality import *

def ratio_plots(etable):
	figure = plt.figure(figsize = (11, 8.5), facecolor = 'white')
   	sci_grid = gridspec.GridSpec(5,7)

	#Fast / Slow (Slow x, Fast y)
	log.info('Building fast/slow subplot')
	plt.subplot(sci_grid[:3,2:5])
	log.info('Building actual slow fast data')
	plot_slowfast(etable)
	    
	figure.suptitle('Object: {0} at {1}'.format(etable.meta['OBJECT'],etable.meta['DATE-OBS']),
        fontsize=18)
	return figure
