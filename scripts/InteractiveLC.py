import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from astropy import log

from nicer.plotutils import *

class InteractiveLC(object):

    def init(self):
	self.fig = 1
	self.etable = None
	self.lclog = None
	self.gtitable = None
	self.binsize = 1.0

    def welcome(self):
        print 'Some code is running, and then you can click on teh plots!'
	print 'Click on the plot to choose an initial point'
	print 'Drag WITHOUT RELEASING to chose a second point'

    def __init__(self, etable, lclog, gtitable,  fig, binsize):
        self.init()
        #self.__dict__.update(**kwargs)
        #self.phases = phases
        #self.primitives = []
        #self.norms = []
        #self.dom = np.linspace(0,1,100)
        #self.welcome()
	self.fig = fig
        self.ax  = plot.gca()
        self.connect()
        meanrate, lc = plot_light_curve(etable, lclog, gtitable, binsize)
	plot.title('Light Curve')
	plot.xlabel('Time Elapsed (s)')
	plot.show()

    def connect(self):
        self.cidpress  = self.fig.canvas.mpl_connect('button_press_event',self.on_press)
        self.cidrelese = self.fig.canvas.mpl_connect('button_release_event',self.on_release)
	
    def on_press(self,event):
        self.x0 = event.xdata
        self.y0 = event.ydata

    def on_release(self,event):
        self.x1 = event.xdata
        self.y1 = event.ydata
	x = [self.x0, self.x1]
	y = [self.y0, self.y0]
	plot.plot(x,y,color = 'r', marker = '+')
	print('The times that you chose are {0} s'.format(x))
	plot.show()

    def writegti(self):
        scol = pyfits.Column(name='START',unit='S',array=np.array(self.x0),format='D')
        ecol = pyfits.Column(name='STOP',array=np.array(self.x1),format='D')
        dcol = pyfits.Column(name='DURATION',array=np.array(self.x1-self.x0),format='D')
        ovhdu = pyfits.BinTableHDU.from_columns([tcol,ocol,ucol], name='NEWGTI')
        ovhdu.writeto("{0}.ovs".format(basename),overwrite=True,checksum=True)
'''	
x = np.arange(0,20,1)
y = x * 4
z = InteractiveLC(x, y)
'''
