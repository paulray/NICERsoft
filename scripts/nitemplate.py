#!/usr/bin/env python
"""
Provide a method for interactively fitting an analytic template to photon phases

Authors:
         Matthew Kerr <matthew.kerr@nrl.navy.mil>
         Paul S. Ray <paul.ray@nrl.navy.mil>
"""
from __future__ import division
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import range
from builtins import object
import numpy as np
import pylab as pl
import os
import astropy.io.fits as pyfits
from pint.templates.lcprimitives import LCGaussian,LCKernelDensity,LCEmpiricalFourier
from pint.templates.lcfitters import LCFitter
from pint.templates.lctemplate import LCTemplate
from optparse import OptionParser
import pickle

def light_curve(phases,weights=None,nbins=25,ec='blue',ls='solid',label=None,axes=None,fignum=1,nmc=100,template=None):
    if axes is None:
        pl.figure(fignum); axes = pl.gca()
    bins = np.linspace(0,1,nbins+1)
    bcs = (bins[1:]+bins[:-1])/2.0
    nph = len(phases)

    cod = axes.hist(phases,bins=bins,weights=weights,normed=True,histtype='step',ec=ec)[0]

    if weights is None:
        err = (cod*float(nbins)/len(phases))**0.5
    else:
        err = np.empty([nbins,nmc])
        for i in range(nmc):
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

    axes.axis([0,1,axes.axis()[2],axes.axis()[3]])
    axes.set_ylabel('Profile Amplitude')
    axes.set_xlabel('Pulse Phase')

    if template is not None:
        axes_resid = axes.twinx()
        model = (template(bcs)*(1-bg_level)+bg_level)
        resids = cod/model-1
        axes_resid.errorbar(bcs,resids,yerr=err/model,
                color='green',ls=' ',marker='o',alpha=0.5)
        axes_resid.axhline(0)
        axes_resid.set_ylabel('Fractional Residuals')
        axes_resid.axis([0,1,axes_resid.axis()[2],axes_resid.axis()[3]])



def get_phases(ft1file,get_weights=False,weightcol='WEIGHT',tmax=999999999):
    f = pyfits.open(ft1file)
    phases = np.asarray(f['EVENTS'].data.field('PULSE_PHASE'),dtype=float)
    mask = f['events'].data.field("TIME") < tmax
    print(len(mask),mask.sum())
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
        self.unbinned = False

    def welcome(self):
        print('Welcome to the interactive unbinned template fitter!')
        print('Displaying the profile... now, we will specify where to put Gaussians.')
        print('For each peak, drag a horizontal line')
        print('         AT THE HIGHEST POINT ON THE PEAK')
        print('         from HALF-MAX to HALF-MAX')
        print('After (left button) each drag, the plot will refresh with the current template.')
        print('To perform a fit, right click.')
        print('To undo a gaussian addition, press u.')
        print('To exit fitting, press c.')

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

    def do_fit(self,doshow=False):
        ubstr = 'unbinned' if self.unbinned else 'binned'
        print('Fitting the template with %s likelihood...'%ubstr)
        template = LCTemplate(self.primitives,norms=self.norms)
        fitter   = LCFitter(template,self.phases,weights=self.weights)
        fitter.fit(estimate_errors=self.errors,use_gradient=True,
                unbinned=self.unbinned)
        print(fitter)
        self.fig = pl.figure(self.fignum)
        self.ax = pl.gca()
        light_curve(self.phases,weights=self.weights,nbins=self.nbins,template=template)
        if doshow:
            pl.show()
        self.fitter = fitter

    def connect(self):
        self.cidpress  = self.fig.canvas.mpl_connect('button_press_event',self.on_press)
        self.cidrelese = self.fig.canvas.mpl_connect('button_release_event',self.on_release)
        self.keypress  = self.fig.canvas.mpl_connect('key_press_event',self.on_key)

    def on_key(self,event):
        if event.key=='u':
            print('Undoing last primitive.')
            self.primitives = self.primitives[:-1]
            self.norms = self.norms[:-1]
            pl.clf()
            template = LCTemplate(self.primitives,norms=self.norms)
            light_curve(self.phases,weights=self.weights,
                    nbins=self.nbins,template=template)
            pl.draw()
        elif event.key=='c':
            print('Closing figure and saving profile.')
            pl.close()

    def on_press(self,event):
        if event.button > 1:
            if len(self.primitives)==0:
                print('Must have at least one component to do fit.')
                return
            # do a fit
            pl.clf()
            self.do_fit(doshow=False)
            self.norms = list(self.fitter.template.norms())
            pl.draw()

        self.x0 = event.xdata
        self.y0 = event.ydata

    def on_release(self,event):
        if event.button > 1:
            return

        x1 = event.xdata
        y1 = event.ydata

        fwhm  = x1 - self.x0
        peak  = (y1 + self.y0)/2.
        phase = (x1 + self.x0)/2.

        # just Gaussian for now
        sigma = fwhm/(8 * np.log(2))**0.5
        # attempt a rough correction for existing template
        if len(self.primitives) > 0:
            template = LCTemplate(self.primitives,norms=self.norms)
            if peak > template(phase):
                peak -= template(phase)
        ampl  = peak * sigma * (2*np.pi)**0.5

        self.primitives.append(LCGaussian(p=[sigma,phase]))
        self.norms.append(ampl)
        norms = np.asarray(self.norms)
        if norms.sum() > 1:
            norms *= 1./norms.sum()
        template = LCTemplate(self.primitives,norms=norms)
        pl.clf()
        light_curve(self.phases,weights=self.weights,nbins=self.nbins,template=template)
        pl.draw()

    def write_template(self,outfile):
        if not hasattr(self,'fitter'):
            print('Must do fit first!'); return
        self.fitter.write_template(outfile)

    def write_profile(self,outfile,nbin,integral=True,suppress_bg=True):
        if not hasattr(self,'fitter'):
            print('Must do fit first!'); return
        self.fitter.template.write_profile(outfile,nbin,integral=integral,
                suppress_bg=suppress_bg)


if __name__ == '__main__':

    desc="Read an FT1 file containing PULSE_PHASE info and interactively fit a template."""
    parser=OptionParser(usage=" %prog [options] [FT1_FILENAME]", description=desc)
    parser.add_option('--nhistbins',type='int',default=50,help="Number of bins to use in phase histogram.")
    parser.add_option('--nprofbins',type='int',default=256,help="Number of bins to use output tabular profile.")
    parser.add_option('-u','--unbinned',action='store_true',default=False,help="Perform fit with unbinned likelihood.")
    parser.add_option('-a','--align',action='store_true',default=False,help="Align template such that peak falls at phase 0.")
    parser.add_option('-w','--weights',action='store_true',default=False,help='Use weighted light curve')
    parser.add_option('-c','--weightcol',type='string',default='WEIGHT',help='Column in FITS file that holds the weight')
    parser.add_option('-p','--prof',type='string',default=None,help='Output name for products (default itemplate.*')
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
        print('%d of %d photons survive weight cut'%(len(phases),len(weights)))
        weights = weights[weights > options.min_weight]

    """
    print 'What type of template would you like to fit?'
    print 'gauss=Gaussian [default], kd=Kernel Density, ef [NHARM]=Empirical Fourier'
    line = raw_input()
    if line.startswith('kd'):
        dom = np.linspace(0.0,1.0,100)
        prim = LCKernelDensity(phases=phases)
        template = LCTemplate([prim],norms=None)
        pl.hist(phases,options.nhistbins,normed=True,histtype='step',edgecolor='k')
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
        pl.hist(phases,options.nhistbins,normed=True,histtype='step',edgecolor='k')
        pl.plot(dom,template(dom),color='red')
        pl.title('Empirical Fourier Template with %d harmonics' % (nharm,))
        pl.show()
        s = raw_input('Enter a filename here to save template for future use.  Just hit ENTER to skip the step.\n')
        if len(s) > 0:
            lcf.to_file(s)
    """
    if False:
        pass
    else:
        intf = InteractiveFitter(phases,nbins=options.nhistbins,weights=weights,errors=options.errors,unbinned=options.unbinned)
        intf.do_fit(doshow=True)

        if options.align:
            intf.fitter.template.align_peak()

        if options.prof is not None:
            # check that specified directory exists
            out = options.prof
            dname = os.path.dirname(out)
            if len(dname) > 0 and not os.path.exists(dname):
                raise IOError('Specified directory %s does not exist!'%(os.path.dirname(out)))
            prof = options.prof
        else:
            prof = 'itemplate'
        intf.write_template(prof+'.gauss')
        intf.write_profile(prof+'.prof',options.nprofbins)
        pickle.dump(intf.fitter.template,open(prof+'.pickle','wb'))

    print('Goodbye!')
