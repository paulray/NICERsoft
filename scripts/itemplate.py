#!/usr/bin/env python
"""
Provide a method for interactively fitting an analytic template to photon phases

Authors:
         Matthew Kerr <matthew.kerr@nrl.navy.mil>
         Paul S. Ray <paul.ray@nrl.navy.mil>
"""
import numpy as np
import pylab as pl
import os
import astropy.io.fits as pyfits
from pint.templates.lcprimitives import LCGaussian,LCKernelDensity,LCEmpiricalFourier
from pint.templates.lcfitters import LCFitter
from pint.templates.lctemplate import LCTemplate
from optparse import OptionParser

def light_curve(phases,weights=None,nbins=25,ec='blue',ls='solid',label=None,axes=None,fignum=1,nmc=100,template=None):
    if axes is None:
        import pylab as pl; pl.figure(fignum); axes = pl.gca()
    bins = np.linspace(0,1,nbins+1)
    bcs = (bins[1:]+bins[:-1])/2
    nph = len(phases)

    cod = axes.hist(phases,bins=bins,weights=weights,normed=True,histtype='step',ec=ec)[0]

    if weights is None:
        err = (cod*float(nbins)/len(phases))**0.5
    else:
        err = np.empty([nbins,nmc])
        for i in xrange(nmc):
            rweights = weights[np.argsort(np.random.rand(nph))]
            err[:,i] = np.histogram(phases,bins=bins,weights=rweights,normed=True)[0]
        err = np.std(err,axis=1)
    axes.errorbar(bcs,cod,yerr=err,color=ec,capsize=0,ls=' ',marker=' ')

    if (weights is not None):
        bg_level = 1-(weights**2).sum()/weights.sum()
        axes.axhline(bg_level,color=ec)
    else: bg_level = 0

    if template is not None:
        dom = np.linspace(0,1,101)
        axes.plot(dom,template(dom)*(1-bg_level)+bg_level,color='red')
        #axes.plot(dom,template(dom,suppress_bg=(weights is not None))*(1-bg_level)+bg_level,color='red')


def get_phases(ft1file,get_weights=False,weightcol='WEIGHT',tmax=999999999):
    f = pyfits.open(ft1file)
    phases = np.asarray(f['EVENTS'].data.field('PULSE_PHASE'),dtype=float)
    mask = f['events'].data.field("TIME") < tmax
    print len(mask),mask.sum()
    phases = phases[mask]
    if get_weights:
        weights = np.asarray(f['EVENTS'].data.field(weightcol),dtype=float)
        weights = weights[mask]
    else: weights = None
    f.close()
    return phases,weights

class InteractiveFitter(object):

    def init(self):
        self.nbins = 50
        self.weights = None
        self.fignum = 1
        self.errors = False

    def welcome(self):
        print 'Welcome to the interactive unbinned template fitter!'
        print 'Displaying the profile... now, we will specify where to put Gaussians.'
        print 'For each peak, drag a horizontal line'
        print '         AT THE HIGHEST POINT ON THE PEAK'
        print '         from HALF-MAX to HALF-MAX'
        print 'After each drag, the plot will refresh with the current template.'
        print 'After all Gaussians are specified, close the plot, and fitting will start.'
        print '(Note -- if using interactively, you will start the fit with do_fit; but do close the plot!'

    def __init__(self,phases,**kwargs):
        self.init()
        self.__dict__.update(**kwargs)
        self.phases = phases
        self.primitives = []
        self.norms = []
        self.dom = np.linspace(0,1,100)
        self.welcome()
        pl.close(self.fignum)
        self.fig = pl.figure(self.fignum)
        self.ax  = pl.gca()
        self.connect()
        light_curve(self.phases,weights=self.weights,nbins=self.nbins,axes=self.ax)
        pl.show()

    def do_fit(self):
        print 'Fitting the template with unbinned likelihood...'
        template = LCTemplate(self.primitives,norms=self.norms)
        fitter   = LCFitter(template,self.phases,weights=self.weights)
        fitter.fit(estimate_errors=self.errors)
        print 'Fitting finished!'
        print fitter
        print 'Overlaying fitted template...'
        self.fig = pl.figure(self.fignum)
        self.ax = pl.gca()
        light_curve(self.phases,weights=self.weights,nbins=self.nbins,axes=self.ax,template=template)
        pl.show()
        self.fitter = fitter

    def connect(self):
        self.cidpress  = self.fig.canvas.mpl_connect('button_press_event',self.on_press)
        self.cidrelese = self.fig.canvas.mpl_connect('button_release_event',self.on_release)

    def on_press(self,event):
        self.x0 = event.xdata
        self.y0 = event.ydata

    def on_release(self,event):
        x1 = event.xdata
        y1 = event.ydata

        fwhm  = x1 - self.x0
        peak  = (y1 + self.y0)/2.
        phase = (x1 + self.x0)/2.

        # just Gaussian for now
        sigma = fwhm/(8 * np.log(2))**0.5
        ampl  = peak * sigma * (2*np.pi)**0.5

        self.primitives.append(LCGaussian(p=[sigma,phase]))
        self.norms.append(ampl)
        template = LCTemplate(self.primitives,norms=self.norms)
        self.ax.clear()
        light_curve(self.phases,weights=self.weights,nbins=self.nbins,axes=self.ax,template=template)
        pl.draw()

    def write_template(self,outfile):
        if not hasattr(self,'fitter'):
            print 'Must do fit first!'; return
        self.fitter.write_template(outfile)


if __name__ == '__main__':

    desc="Read an FT1 file containing PULSE_PHASE info and interactively fit a template."""
    parser=OptionParser(usage=" %prog [options] [FT1_FILENAME]", description=desc)
    parser.add_option('-n','--nbins',type='int',default=50,help="Number of bins to use in phase histogram.")
    parser.add_option('-w','--weights',action='store_true',default=False,help='Use weighted light curve')
    parser.add_option('-c','--weightcol',type='string',default='WEIGHT',help='Column in FITS file that holds the weight')
    parser.add_option('-p','--prof',type='string',default=None,help='Output name for profile')
    parser.add_option('-m','--min_weight',type='float',default=1e-2,help='Minimum weight to include in fit.')
    parser.add_option('-T','--tmax',type='float',default=999999999,help='Maximum time to include in fit.')
    parser.add_option('-e','--errors',action='store_true',default=False,help='Compute errors on components.')

    ## Parse arguments
    (options,args) = parser.parse_args()
    if len(args) < 1:
        raise ValueError('Must provide an input FITS file!')

    phases,weights = get_phases(args[0],get_weights=options.weights,weightcol=options.weightcol)

    if options.weights:
        phases = phases[weights > options.min_weight]
        print '%d of %d photons survive weight cut'%(len(phases),len(weights))
        weights = weights[weights > options.min_weight]

    print 'Welcome to the interactive unbinned template fitter!'
    print 'What type of template would you like to fit?'
    print 'gauss=Gaussian [default], kd=Kernel Density, ef [NHARM]=Empirical Fourier'
    line = raw_input()
    if line.startswith('kd'):
        dom = np.linspace(0.0,1.0,100)
        prim = LCKernelDensity(phases=phases)
        template = LCTemplate([prim],norms=None)
        pl.hist(phases,options.nbins,normed=True,histtype='step',edgecolor='k')
        pl.plot(dom,template(dom),color='red')
        pl.title('Kernel Density Template Fit')
        pl.show()
        s = raw_input('Enter a filename here to save template for future use.  Just hit ENTER to skip the step.\n')
        if len(s) > 0:
            prim.to_file(s)

    elif line.startswith('ef'):
        dom = np.linspace(0.0,1.0,100)
        if len(line.split()) > 1:
            nharm = int(line.split()[1])
        else:
            nharm = 16
        lcf = LCEmpiricalFourier(phases=phases,nharm=nharm)
        template = LCTemplate([lcf],norms=None)
        pl.hist(phases,options.nbins,normed=True,histtype='step',edgecolor='k')
        pl.plot(dom,template(dom),color='red')
        pl.title('Empirical Fourier Template with %d harmonics' % (nharm,))
        pl.show()
        s = raw_input('Enter a filename here to save template for future use.  Just hit ENTER to skip the step.\n')
        if len(s) > 0:
            lcf.to_file(s)
    else:
        intf = InteractiveFitter(phases,nbins=options.nbins,weights=weights,errors=options.errors)
        intf.do_fit()

        if options.prof is not None:
            # check that specified directory exists
            out = options.prof
            dname = os.path.dirname(out)
            if len(dname) > 0 and not os.path.exists(dname):
                raise IOError('Specified directory %s does not exist!'%(os.path.dirname(out)))
        else:
            out = ''
            out = raw_input('Enter filename for gaussian profile output file, or just hit ENTER to exit...:  ')
        if len(out) > 0:
            print 'Writing Gaussian-style template to %s...'%(out)
            intf.write_template(out)

    print 'Goodbye!'
