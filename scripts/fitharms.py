#!/usr/bin/env python
from __future__ import print_function,division

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from pint.templates import lctemplate,lcprimitives,lcfitters
from pint.eventstats import z2m,sf_z2m, hm, sf_hm, sig2sigma
import sys
from astropy import log
import scipy.stats

def compute_fourier(phases,nh=10,pow_phase=False):
    '''Compute Fourier amplitudes from an array of pulse phases
    phases should be [0,1.0)
    nh is the number of harmonics (1 = fundamental only)
    Returns: cos and sin component arrays, unless pow_phase is True
    then returns Fourier power (Leahy normalized) and phase arrays 
    DC bin is not computed or returned
    '''
    phis = 2.0*np.pi*phases  # Convert phases to radians
    n = len(phis)
    c = np.asarray([(np.cos(k*phis)).sum() for k in range(1,nh+1)])/n
    s = np.asarray([(np.sin(k*phis)).sum() for k in range(1,nh+1)])/n
    c *= 2.0
    s *= 2.0

    if pow_phase:
        # CHECK!  There could be errors here!
        # These should be Leahy normalized powers
        fourier_pow = (n/2)*(c**2+s**2)
        fourier_phases = np.arctan2(s,c)
        return n,fourier_pow,fourier_phases
    else:
        return n,c,s

def evaluate_fourier(n,c,s,nbins,k=None):
    # This should be updated to do a little integral over each bin. 
    # Currently evaluates the model at the center of each bin
    model = np.zeros(nbins)+n/nbins
    theta = 2.0*np.pi*np.arange(nbins,dtype=np.float)/nbins
    theta += theta[1]/2.0
    if k is not None:
        model += (n/nbins)*(c[k]*np.cos((k+1)*theta) + s[k]*np.sin((k+1)*theta))
    else:
        for k in range(len(c)):
            model += (n/nbins)*(c[k]*np.cos((k+1)*theta) + s[k]*np.sin((k+1)*theta)) 

    return model

def evaluate_chi2(hist,model):
    # Question here is whether error should be sqrt(data) or sqrt(model)
    return ((hist-model)**2/model).sum()

def compute_phist(phases,nbins=200):
    h, edges = np.histogram(phases,bins=np.linspace(0.0,1.0,nbins+1,endpoint=True))
    return edges[:-1], h
    
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "Fit a set of pulse phases to harmonics")
    parser.add_argument("evname", help="Input event file (must have PULSE_PHASE column)") 
    parser.add_argument("--white",help = "Replace phases with white random numbers, for testing", action="store_true")
    parser.add_argument("--txt",help = "Assume input file is .txt instead of FITS", action="store_true")
    parser.add_argument("--showcomps",help = "Show individual components of harmonic fit on plot", action="store_true")
    parser.add_argument("--output",help = "Save figures with basename", default=None)
    parser.add_argument("--numharm",help="Max harmonic to use in analysis (1=Fundamental only)",default=4,type=int)
    parser.add_argument("--numbins",help="Number of bins for histograms",default=200,type=int)
    parser.add_argument("--emin",help="Minimum energy to include (keV)",default=0.25,type=float)
    parser.add_argument("--emax",help="Maximum energy to include (keV)",default=12.0,type=float)
    args = parser.parse_args()

    if args.txt:
        exposure = None
        ph,en = np.loadtxt(args.evname,unpack=True,usecols=(1,2),skiprows=3)
        log.info("Read {0} phases from .txt file".format(len(ph)))
    else:
        f = fits.open(args.evname)
        en = f['events'].data.field('pi')
        ph = f['events'].data.field('pulse_phase')
        log.info("Read {0} phases from FITS file".format(len(ph)))
        exposure = float(f['events'].header['EXPOSURE'])
        log.info("Exposure = {0} s".format(exposure))

    if args.white:
        # Random phases uniform over [0,1)
        ph = np.random.random_sample(len(en))
        log.info("Replaced with {0} random phases".format(len(en)))

 
    matplotlib.rcParams['font.family'] = "serif"
    matplotlib.rcParams.update({'font.size': 13})
    matplotlib.rc('axes', linewidth=1.5)
        
    # Filter on energy
    idx = np.where(np.logical_and(en > int(args.emin*100), en < int(args.emax*100) ))[0]
    ph = ph[idx]
    en = en[idx]
    
    # Hack to manually split out a segment
    #q = 3  # Use 0, 1, 2, 3
    #qn = len(ph)//4
    #ph = ph[q*qn:(q+1)*qn]
    #en = en[q*qn:(q+1)*qn]
    
    nbins = args.numbins
    bins,phist = compute_phist(ph,nbins=nbins)
    fig,axs = plt.subplots(nrows=2,ncols=1)
    plt.subplots_adjust(left=0.15, bottom=0.1, right=0.97, top=0.94,hspace=0.001)
    ax=axs[0]
    ax.tick_params(direction='in', length=6, width=2, colors='k',top=True, right=True, labelbottom=False)
#    ax.text(.5,.8,'PSR J0030+0451', horizontalalignment='center', transform=ax.transAxes)
#    ax.text(.5,.8,'PSR J0437-4715', horizontalalignment='center', transform=ax.transAxes)
#    ax.text(.2,.8,'PSR J1231-1411', horizontalalignment='center', transform=ax.transAxes)
#    ax.text(.8,.8,'PSR J2124-3358', horizontalalignment='center', transform=ax.transAxes)
    ax.step(np.concatenate((bins,np.ones(1))),np.concatenate((phist,phist[-1:])),color='k',where='post')
    ax.set_xlim(0.0,1.0)
    ax.set_ylabel('Counts per bin')


    n,c,s = compute_fourier(ph,nh=args.numharm)
    model = evaluate_fourier(n,c,s,nbins)
    ax.plot(bins+bins[1]/2.0,model,color='r',lw=2)
    if args.showcomps:
        for k in range(len(c)):
            ax.plot(np.linspace(0.0,1.0,nbins),evaluate_fourier(n,c,s,nbins,k=k),ls='--')

    fn,fpow,fphase = compute_fourier(ph,nh=args.numharm,pow_phase=True)
    i=1
    log.info("Harm  LeahyPower   Phase(deg)")
    for fp, fph in zip(fpow,fphase):
        log.info("{0:2d} {1:12.3f} {2:9.3f} deg".format(i,fp,np.rad2deg(fph)))
        i+=1

    pcounts = (model-model.min()).sum()
    if exposure:
        log.info("Pulsed counts = {0:.3f}, count rate = {1:.3f} c/s".format(pcounts, pcounts/exposure))
        log.info("Total rate = {0:.3f} c/s, Unpulsed rate = {1:.3f} c/s".format(n/exposure, n/exposure-pcounts/exposure))

    ax = axs[1]
    ax.tick_params(direction='in', length=6, width=2, colors='k',top=True, right=True)
    ax.errorbar(np.linspace(0.0,1.0,nbins),phist-model,yerr=np.sqrt(phist),fmt='.',ecolor='k')
    chisq = evaluate_chi2(phist,model)
    nparams = 1 + 2*args.numharm # 1 for DC + 2 for each sinusoidal component
    ax.set_xlim(0.0,1.0)
    ax.set_xlabel('Pulse Phase')
    ax.set_ylabel('Residuals (counts)')
    ax.tick_params(direction='in', length=6, width=2, colors='k',top=True)
    ndof = len(phist)-nparams
    axs[0].set_title("NumHarm = {0}, Chisq = {1:.2f}, DOF = {2}".format(args.numharm,chisq,ndof))
    ax.grid(1)
#    ax.set_label("{0} Harmonic Fit to Profile".format(args.numharm))
    plt.tight_layout()
    
    if args.output:
        fig.savefig("{0}_harmfit.pdf".format(args.output))

    # Plot distribution of residuals to compare to a gaussian
    fig,ax = plt.subplots()
    ax.tick_params(direction='in', length=6, width=2, colors='k',top=True, right=True)
    chi = (phist-model)/np.sqrt(model)
    #x, y = np.histogram(chi,bins=np.linspace(-2.0,2.0,0.1))
    x = np.linspace(-3.0,3.0,32,endpoint=True)
    ax.hist(chi,bins=x,density=True)
    ax.set_title('Histogram of residuals')
    ax.plot(x,scipy.stats.norm.pdf(x))
    plt.tight_layout()
    
    #  Plot histogram of phase differences to see if they are Poisson
    fig,ax = plt.subplots()
    ax.tick_params(direction='in', length=6, width=2, colors='k',top=True, right=True)
    ph.sort()
    pdiffs = (ph[1:]-ph[:-1])*1.0
    x = np.linspace(0.0,50.0e-6,200,endpoint=True)
    histn, histbins, histpatches  = ax.hist(pdiffs,bins=x,density=True,log=True)
    ax.set_title('Histogram of phase differences')
    ax.set_xlabel('Phase diff')
    ax.plot(x,np.exp(-len(pdiffs)*(x*1.0))*n)
    plt.tight_layout()

# Compute number of significant harmonics
# First by plotting Leahy powers

    fig,axs = plt.subplots(nrows=2,ncols=1)
    ax = axs[0]
    ax.tick_params(direction='in', length=6, width=2, colors='k',top=True, right=True)
    n,pow,phases = compute_fourier(ph,nh=nbins//2,pow_phase=True)
    ax.semilogy(np.arange(len(pow))+1,pow,marker='o')
    # Leahy power of 5.99 corresponds to 2 sigma, I think
    ax.axhline(5.99,color='r')
    ax.axhline(2.0,color='b',ls='--')
    #ax.xaxis.set_ticks(np.arange(1,len(pow)+1))
    #ax.set_xlabel('Harmonic Number')
    ax.set_ylabel('Leahy Power')
    ax.set_title("Power Spectrum")
    plt.tight_layout()
    
    ax = axs[1]
    ax.tick_params(direction='in', length=6, width=2, colors='k',top=True, right=True)
    ax.plot(np.arange(len(pow))+1,pow,marker='o')
    ax.axhline(5.99,color='r')
    ax.axhline(2.0,color='b',ls='--')
    #ax.xaxis.set_ticks(np.arange(1,len(pow)+1))
    ax.set_ylim(0.0,10.0)
    ax.text(1.0,7.0,'Mean power {0:.3f}'.format(pow.mean()))
    ax.set_xlabel('Harmonic Number')
    ax.set_ylabel('Leahy Power')
    if args.output:
        fig.savefig("{0}_leahy.pdf".format(args.output))
    plt.tight_layout()

# Then by computing chisq as a function of number of harmonics in model

    chisq = []
    ndof = []
    maxharms = np.arange(1,min(33,nbins//4+1))
    n,c,s = compute_fourier(ph,nh=maxharms[-1])
    for maxharm in maxharms:
        model = evaluate_fourier(n,c[:maxharm],s[:maxharm],nbins)
        chisq.append(evaluate_chi2(phist,model))
        nparams = 1 + 2*maxharm # 1 for DC + 2 for each sinusoidal component
        ndof.append(len(phist)-nparams)
    chisq = np.asarray(chisq)
    ndof = np.asarray(ndof)
    fig,ax = plt.subplots()
    ax.tick_params(direction='in', length=6, width=2, colors='k',top=True, right=True)
    ax.plot(maxharms,chisq/ndof,'o',ls='-')
    ax.set_ylim(0.5,3.0)
    ax.axhline(1.0,color='r',ls='--')
    ax.set_xlabel('Number of Harmonics')
    ax.set_ylabel('Chisq')
    ax.set_title("Chisq/DOF vs. Number of Harmonics")
    #ax.xaxis.set_ticks(maxharms)
    #ax.semilogy(maxharms,ndof)
    plt.tight_layout()
    
    if args.output:
        fig.savefig("{0}_chisq.pdf".format(args.output))

# Then look at amplitudes and phases as a function of energy cuts

# Look at color oscillations
# Select photons above and below some energy cut and look at the ratio
    ensplit = 55
    softidx = np.where(en<ensplit)[0]
    hardidx = np.where(en>=ensplit)[0]
    colorbins = 32
    softbins, softn = compute_phist(ph[softidx],nbins=colorbins)
    hardbins, hardn = compute_phist(ph[hardidx],nbins=colorbins)
    softn = np.asarray(softn,dtype=np.float)
    hardn = np.asarray(hardn,dtype=np.float)
    fig,ax = plt.subplots()
    color = hardn/softn
    # Propagate Poisson errors to get error in ratio
    cerr = color*np.sqrt(1.0/softn + 1.0/hardn)
    #ax.step(np.concatenate((softbins,np.ones(1))),np.concatenate((color,color[-1:])),color='C0',where='post')
    ax.errorbar(softbins+0.5*softbins[1],color,yerr=cerr,color='k',fmt='.')
    ax.set_xlim(0.0,1.0)
    ax.set_xlabel('Pulse Phase')
    ax.set_ylabel('Spectral Color')
    

    plt.show()
    
