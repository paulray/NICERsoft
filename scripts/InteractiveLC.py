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
        print 'Some code is running, and then you can click on the plots!'
        print 'Click on the plot to choose an initial point'
        print 'Drag WITHOUT RELEASING to choose a second point'

    def __init__(self, etable, lclog, gtitable,  fig, name, binsize):
        self.init()
        self.welcome()
        self.fig = fig
        self.ax  = plot.gca()
        self.starts = np.array([]) #start times
        self.stops = np.array([]) #end times
        self.dcol = np.array([])
        self.scol = np.array([])
        self.ecol = np.array([])
        self.connect()
        self.gtitable = gtitable
        self.name = name
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
        self.getclicks()
        plot.show()

    def getclicks(self):
        self.starts =  np.append(self.starts, self.x0)
        self.stops = np.append(self.stops,self.x1)
        self.starts = np.sort(self.starts)
        self.stops = np.sort(self.stops)
        print(self.starts)
        print(self.stops)

    def getgoodtimes(self):
        if len(self.starts) > 0: #If there were clicks
            for idx in xrange(0,len(self.gtitable['CUMTIME'])): #For every GTI
                goodlist = np.array([])
                for sidx in xrange(0, len(self.starts)):#For every click
                    if self.gticheck(self.starts[sidx], idx):#Check if any of the clicks were in this GTI
                        goodlist = np.append(goodlist,sidx)
                    
                if len(goodlist) > 0:#If there were clicks
                    #The first interval will start here
                    start = self.gtitable['START'][idx]
                    for index in xrange(0,len(goodlist)):#for every click in this GTI
                        if np.logical_and(index == 0,index == (len(goodlist)-1)) : #if it's the first AND ONLY one in goodlist
                            end = self.starts[index] + self.gtitable['START'][idx]
                            self.addtolist(start, end) 
                            print('FIRST AND ONLY')
                            break

                        elif index == 0:
                            end = self.starts[index] + self.gtitable['START'][idx]
                            self.addtolist(start, end) 
                            print('first')
                        elif index == (len(goodlist)-1): #if it's the last one in goodlist
                            start = self.stops[index-1] + self.gtitable['START'][idx]
                            end = self.starts[index]+ self.gtitable['START'][idx]
                            self.addtolist(start,end)
                            print('last')

                        else:#if it's in the middle
                            start = self.stops[index] + self.gtitable['START'][idx]
                            end = self.starts[index]+ self.gtitable['START'][idx]
                            self.addtolist(start,end)
                            print('middle')
                    #Tacks the end on to the interval set
                    print('tacking on the end')
                    start = self.stops[index] + self.gtitable['START'][idx]
                    end = self.gtitable['STOP'][idx]
                    self.addtolist(start,end)
                    goodlist = np.array([])

                    del start, end

                else:#If there were no clicks 
                    print('no clicks here')
                    start = self.gtitable['START'][idx]
                    end = self.gtitable['STOP'][idx]
                    self.addtolist(start,end)

        self.scol = np.sort(self.scol)
        self.ecol = np.sort(self.ecol)

    def addtolist(self, start, end):
        self.scol = np.append(self.scol,start)
        self.ecol = np.append(self.ecol,end)

    def gticheck(self, clickstart, idx):
        if idx == (len(self.gtitable['CUMTIME'])-1):
            if np.logical_and(clickstart < (self.gtitable['CUMTIME'][idx] + self.gtitable['DURATION'][idx]), clickstart > (self.gtitable['CUMTIME'][idx])):
                ans = True
            else:
                ans = False
        else:
            if np.logical_and(clickstart < self.gtitable['CUMTIME'][idx+1], clickstart > self.gtitable['CUMTIME'][idx]):
                ans = True
            else:
                ans = False 
        return ans
        
    def writegti(self):
        scol = pyfits.Column(name='START',unit='S',array=self.scol,format='D')
        ecol = pyfits.Column(name='STOP',array=self.ecol,format='D')
        ovhdu = pyfits.BinTableHDU.from_columns([scol,ecol], name='NEWGTI')
        ovhdu.writeto("{0}.ovs".format(self.name),overwrite=True,checksum=True)
        

