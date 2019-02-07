#!/usr/bin/env python
from __future__ import print_function,division

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from pint.templates import lctemplate,lcprimitives,lcfitters
from pint.eventstats import z2m,sf_z2m, hm, sf_hm, sig2sigma
import sys
from astropy import log

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

def evaluate_fourier(n,c,s,nbins):
    model = np.zeros(nbins)+n/nbins
    theta = np.linspace(0.0,2.0*np.pi,nbins)+np.pi/nbins
    for k in range(len(c)):
        model += (n/nbins)*(c[k]*np.cos((k+1)*theta) + s[k]*np.sin((k+1)*theta)) 

    return model

def compute_phist(phases,nbins=200):
    h, edges = np.histogram(ph,bins=np.linspace(0.0,1.0,nbins+1,endpoint=True))
    return edges[:-1], h
    
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "Fit a set of pulse phases to harmonics")
    parser.add_argument("evname", help="Input event file (must have PULSE_PHASE column)") 
    parser.add_argument("--white",help = "Replace phases with white random numbers, for testing", action="store_true")
    parser.add_argument("--numharm",help="Max harmonic to use in analysis (1=Fundamental only)",default=4,type=int)
    parser.add_argument("--numbins",help="Number of bins for histograms",default=200,type=int)
    args = parser.parse_args()

    f = fits.open(args.evname)
    en = f['events'].data.field('pi')
    if args.white:
        # Random phases uniform over [0,1)
        ph = np.random.random_sample(len(en))
        log.info("Generated {0} random phases".format(len(en)))
    else:
        ph = f['events'].data.field('pulse_phase')
        log.info("Read {0} phases".format(len(ph)))

    nbins = args.numbins
    bins,phist = compute_phist(ph,nbins=nbins)
    fig,axs = plt.subplots(nrows=2,ncols=1)
    ax=axs[0]
    ax.step(bins,phist,where='post')
    ax.set_xlim(0.0,1.0)
    ax.set_ylabel('Counts')

    n,c,s = compute_fourier(ph,nh=args.numharm)
    model = evaluate_fourier(n,c,s,nbins)
    ax.plot(np.linspace(0.0,1.0,nbins),model)

    ax = axs[1]
    ax.errorbar(np.linspace(0.0,1.0,nbins),phist-model,yerr=np.sqrt(phist),fmt='.')
    chisq = ((phist-model)**2/phist).sum()
    nparams = 1 + 2*args.numharm # 1 for DC + 2 for each sinusoidal component
    ax.set_xlabel('Pulse Phase')
    ndof = len(phist)-nparams
    axs[0].set_title("NumHarm = {0}, Chisq = {1:.2f}, DOF = {2}".format(args.numharm,chisq,ndof))
    ax.grid(1)
    #ax.set_label("{0} Harmonic Fit to Profile".format(args.numharm))

# Compute number of significant harmonics
# First by plotting Leahy powers

    fig,axs = plt.subplots(nrows=2,ncols=1)
    ax = axs[0]
    n,pow,phases = compute_fourier(ph,nh=50,pow_phase=True)
    ax.semilogy(np.arange(len(pow))+1,pow,marker='o')
    # Leahy power of 5.99 corresponds to 2 sigma, I think
    ax.axhline(5.99,color='r')
    ax.axhline(2.0,color='b',ls='--')
    #ax.xaxis.set_ticks(np.arange(1,len(pow)+1))
    #ax.set_xlabel('Harmonic Number')
    ax.set_ylabel('Leahy Normalized Power')
    ax.set_title("Power Spectrum")
    
    ax = axs[1]
    ax.plot(np.arange(len(pow))+1,pow,marker='o')
    ax.axhline(5.99,color='r')
    ax.axhline(2.0,color='b',ls='--')
    #ax.xaxis.set_ticks(np.arange(1,len(pow)+1))
    ax.set_ylim(0.0,8.0)
    ax.text(1.0,7.0,'Mean power {0:.3f}'.format(pow.mean()))
    ax.set_xlabel('Harmonic Number')
    ax.set_ylabel('Leahy Normalized Power')


# Then by computing chisq as a function of number of harmonics in model

    chisq = []
    ndof = []
    maxharms = np.arange(1,20)
    for maxharm in maxharms:
        n,c,s = compute_fourier(ph,nh=maxharm)
        model = evaluate_fourier(n,c,s,nbins)
        chisq.append(((phist-model)**2/phist).sum())
        nparams = 1 + 2*maxharm # 1 for DC + 2 for each sinusoidal component
        ndof.append(len(phist)-nparams)
    chisq = np.asarray(chisq)
    ndof = np.asarray(ndof)
    fig,ax = plt.subplots()
    ax.plot(maxharms,chisq/ndof)
    ax.set_ylim(0.5,3.0)
    ax.axhline(1.0,color='r',ls='--')
    ax.set_xlabel('Number of Harmonics')
    ax.set_ylabel('Chisq')
    ax.set_title("Chisq/DOF vs. Number of Harmonics")
    ax.xaxis.set_ticks(maxharms)
    #ax.semilogy(maxharms,ndof)


# Then look at amplitudes and phases as a function of energy cuts

    plt.show()
