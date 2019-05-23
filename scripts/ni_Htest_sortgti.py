#!/usr/bin/env python
# Version: 1.0
# Author: M. Kerr (updated by S. Guillot)

from __future__ import division, print_function
import sys, math, os
import numpy as np
import scipy.stats
import datetime
import matplotlib.pyplot as plt
import glob
from collections import deque
import astropy.units as u
from astropy import log
from astropy.io import fits
from astropy.time import Time
from nicer.values import *
from pint.fits_utils import read_fits_event_mjds
from pint.eventstats import h2sig,hm
import argparse

desc= """
Read one or more event files 
to sort GTI by background rate 
and evaluate H-test
"""

plt.rc('font', size=14)          # controls default text sizes
plt.rc('axes', labelsize=13)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=13)    # fontsize of the tick labels
plt.rc('ytick', labelsize=13)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rc('axes', linewidth=1.5)
plt.rc('xtick.major', size=4, width=1.5)
plt.rc('ytick.major', size=4, width=1.5)


parser = argparse.ArgumentParser(description = desc)
#parser.add_argument("evt", help="Input event files using glob.glob (e.g.'10010101*_Proc/cleanfilt_par.evt')")
parser.add_argument("infile", help="file or text file with list of event file", nargs='+')
parser.add_argument("outfile", help="name for output files")
parser.add_argument("--emin", help="Minimum energy to include (keV, default=0.25)", type=float, default=0.25)
parser.add_argument("--emax", help="Maximum energy to include (keV, default=2.00)", type=float, default=2.00)
parser.add_argument("--gridsearch", help="Search over energies to find max H-test", action="store_true",default=False)
parser.add_argument("--coarsegridsearch", help="Search over energies to find max H-test", action="store_true",default=False)
parser.add_argument("--savefile", help="Saving optimized event file", action="store_true",default=False)
parser.add_argument("--nbins", help="Number of bins for plotting pulse profile (default=16)", type=int, default=16)
parser.add_argument("--name", help="Pulsar name for output figure", type=str, default='')
args = parser.parse_args()

def load_files(fnames):
    # make a comprehensive list of gtis/rates
    
    gtis = deque()
    times = deque()
    pis = deque()
    phases = deque()
    for fname in fnames:
        f = fits.open(fname)
        t0s = f['gti'].data.field('start')
        t1s = f['gti'].data.field('stop')
        for t0,t1 in zip(t0s,t1s):
            gtis.append([t0,t1])
        times.append(f['events'].data.field('time'))
        pis.append(f['events'].data.field('pi'))
        try:
            phases.append(f['events'].data.field('pulse_phase'))
        except:
            pass

    times = np.concatenate(times)
    pis = np.concatenate(pis)
    if len(phases) > 0:
        phases = np.concatenate(phases)
    else:
        phases = None
    t0s,t1s = np.asarray(list(gtis)).transpose()
    return times,phases,pis,t0s,t1s

def dice_gtis(data,tmax=100):
    """ Break larger GTIs into small pieces to handle rate variations."""
    times,phases,pis,t0s,t1s = data
    new_t0s = deque()
    new_t1s = deque()
    for t0,t1 in zip(t0s,t1s):
        dt = t1-t0
        if dt < tmax:
            new_t0s.append(t0)
            new_t1s.append(t1)
        else:
            # break up GTI in such a way to avoid losing time (to tmin) and
            # to avoid having pieces longer than tmax
            npiece = int(np.floor(dt/tmax))+1
            new_edges = np.linspace(t0,t1,npiece+1)
            for it0,it1 in zip(new_edges[:-1],new_edges[1:]):
                new_t0s.append(it0)
                new_t1s.append(it1)
    return times,phases,pis,np.asarray(new_t0s),np.asarray(new_t1s)


def ensemble_htest(phases,slices,m=20,c=4):
    """ Calculate H-test statistic for subsets of a set of phases.
        Cache intermediate products to avoid O(N^2) complexity!
    """
    
    phases = np.asarray(phases)*(2*np.pi) # in radians and copy
    cache = np.empty((2*m,len(phases)))
    for i in range(m):
        cache[2*i] = np.cos((i+1)*phases)
        cache[2*i+1] = np.sin((i+1)*phases)

    rvals = np.empty(len(slices))
    penalty = c*np.arange(0,m)
    for isl,sl in enumerate(slices):
        x = cache[:,sl]
        nph = x.shape[1]
        if nph == 0:
            rvals[isl] = 0.01
            continue
        t = np.sum(cache[:,sl],axis=1)**2
        t = np.cumsum(t[::2] + t[1::2])*(2./nph)
        rvals[isl] = np.max(t-penalty)
    return rvals

def make_sn(data,mask=None,rate=0.1,min_gti=5):
    """ data -- output of load_local
        mask -- optional mask to select events (e.g. on PI)
        rate -- assumed rate for S/N calculation in ct/sec
        min_gti -- minimum GTI length in seconds
    """

    times,phases,pis,t0s,t1s = data
    if mask is not None:
        times = times[mask]
        phases = phases[mask]
        pis = pis[mask]
    # determine which gti each event belongs to
    gti_idx = np.searchsorted(t1s,times)
    # count events in each gti
    gti_cts = np.bincount(gti_idx,minlength=len(t1s))
    gti_len = t1s-t0s

    mask = gti_len > min_gti
    rates = (gti_cts / gti_len)[mask]
    a = np.argsort(rates)
    gti_len_s = gti_len[mask][a]
    gti_cts_s = gti_cts[mask][a]
    gti_rts_s = gti_cts_s/gti_len_s

    sn = rate*np.cumsum(gti_len_s)/np.sqrt(np.cumsum(gti_cts_s)+rate*gti_len_s)
    rate = 0
    sn0 = np.cumsum(gti_len_s)/np.sqrt(np.cumsum(gti_cts_s))

    if phases is not None:
        counter = 0
        ph_gti = deque()
        for i,ct in enumerate(gti_cts):
            ph_gti.append(phases[counter:counter+ct])
            counter += ct
        # apply mask
        ph_gti = [ph_gti[i] for i,imask in enumerate(mask) if imask]
        # apply sorting
        ph_gti = [ph_gti[i] for i in a]
        # make a set of slices
        nph = np.cumsum([len(phg) for phg in ph_gti])
        slices = [slice(0,n) for n in nph]
        assert(len(slices)==len(ph_gti))
        # calculate H test
        hs = ensemble_htest(np.concatenate(ph_gti),slices)

    else:
        hs = None
        ph_gti = None

    return sn,sn0,hs,ph_gti,gti_rts_s,gti_len_s

## UNUSED...
def get_optimal_cuts(data,pred_rate = 0.017):
    """ Determine rates from analytic prediction for both 0.25 and 0.40 keV
        cuts, and use these to estimate an optimal H test.
    """
    pi_mask = (data[2]>25) & (data[2]<100)
    sn,sn0,hs,ph_gti,gti_rts_s,gti_len_s = make_sn(data,mask=pi_mask,rate=pred_rate)
    amax = np.argmax(sn)

    pi_mask = (data[2]>40) & (data[2]<100)
    pred_rate *= 3./5
    sn_40,sn0_40,hs_40,ph_gti_40,gti_rts_s_40,gti_len_s_40 = make_sn(data,mask=pi_mask,rate=pred_rate)
    hsig_40 = [h2sig(h) for h in hs]
    exposure = np.cumsum(gti_len_s)
    amax_40 = np.argmax(sn_40)

    # get intersection of 0.25 keV data sets and 0.4 keV data sets
    rate_025_cut = gti_rts_s[amax]
    rate_040_cut = gti_rts_s_40[amax_40]
    mask_025 = gti_rts_s<=rate_025_cut
    mask_040 = ~mask_025 & (gti_rts_s_40<rate_040_cut)
    ph1 = np.concatenate(np.asarray(ph_gti)[mask_025])
    try:
        ph2 = np.concatenate(np.asarray(ph_gti_40)[mask_040])
    except ValueError:
        ph2 = np.asarray([])

    print('Found %d photons satisfying 0.25 keV cut; exposure = %.1fs'%(len(ph1),np.sum(np.asarray(gti_len_s)[mask_025])))
    print('Found %d photons satisfying 0.40 keV cut; exposure = %.1fs'%(len(ph2),np.sum(np.asarray(gti_len_s_40)[mask_040])))
    print('H-test 1: %.2f'%hm(ph1))
    if len(ph2)>0:
        print('H-test 2: %.2f'%hm(ph2))
    else:
        print('H-test 2: no photons!')
    print('H-test joint: %.2f'%hm(np.append(ph1,ph2)))



if len(args.infile)==1:
    if args.infile[0].startswith('@'):
        inputfile = args.infile[0].split('@')[1]
        log.info('Reading input ObsID list: {}'.format(inputfile))
        all_files = np.loadtxt(inputfile,dtype=str)
    else:
        all_files = args.infile
else:
    all_files = args.infile

data = load_files(all_files)
data_diced = dice_gtis(data)
#get_optimal_cuts(data_diced)

if args.gridsearch:
    all_emin = np.arange(0.24,1.0,0.01)
elif args.coarsegridsearch:
    all_emin = np.arange(0.3,2.0,0.1)
else:
    all_emin = np.array([args.emin])

hbest = 0.0
eminbest = 0.0
emaxbest = 100.0

eminlist = []
emaxlist = []
hgrid = []
    
for emin in all_emin:
    
    if args.gridsearch:
        all_emax = np.arange(1.0,3.0,0.01)
    elif args.coarsegridsearch:
        all_emax = np.arange(emin+0.1,7.0,0.2)
    else:
        all_emax = np.array([args.emax])

    for emax in all_emax:

        print("Energy range: {:0.2f}-{:0.2f} keV".format(emin,emax))
        
        pi_mask = (data[2]>emin*KEV_TO_PI) & (data[2]<emax*KEV_TO_PI)
        pred_rate = 0.05/10.0 # 2241
        sn,sn0,hs,ph_gti,gti_rts_s,gti_len_s = make_sn(data_diced,mask=pi_mask,rate=pred_rate)
        
        hsig = [h2sig(h) for h in hs]
        exposure = np.cumsum(gti_len_s)
        
        plt.figure(5); plt.clf()
        # scale exposure to the expected S/N
        amax = np.argmax(sn)
        exposure_scale = sn[amax]/exposure[amax]**0.5

        Hmax = np.argmax(hsig)
        
        if not args.gridsearch and not args.coarsegridsearch:
            #plt.plot(gti_rts_s,exposure**0.5*exposure_scale,label='scaled exposure')
            #plt.plot(gti_rts_s,sn,label='predicted S/N')
            plt.plot(gti_rts_s,hsig,label='H-test significance')
            plt.axvline(gti_rts_s[amax],color='k',ls='--')
            plt.axvline(gti_rts_s[Hmax],color='r',ls='--')
            plt.xlabel('Background Rate (ct/s)')
            plt.ylabel('Significance (sigma)')
            plt.title('{} - [{},{}]'.format(args.name,emin,emax))
            plt.legend(loc='lower right')
            plt.savefig('{}_sig.png'.format(args.outfile))
            
            plt.clf()
            nbins=args.nbins
            select_ph = np.concatenate(ph_gti[:Hmax]).ravel()
            profbins = np.linspace(0.0,1.0,nbins+1,endpoint=True)
            profile, edges = np.histogram(select_ph,bins=profbins)
            bbins = np.concatenate((profbins, profbins[1:]+1.0, profbins[1:]+2.0))
            fullprof = np.concatenate((profile,profile,profile,np.array([profile[0]])))
            plt.errorbar(bbins-(0.5/nbins),fullprof,
                         yerr=fullprof**0.5,
                         marker ='',
                         drawstyle='steps-mid',
                         linewidth=1.5,
                         color='k'
            )
            #plt.subplots_adjust(left=0.15, right=0.93)  #, bottom=0.1)
            plt.xlim((0.0,2.0))
            plt.ylabel('Photons')
            plt.xlabel('Phase')
            plt.title(args.name)
            plt.savefig('{}_profile.png'.format(args.outfile))
            plt.clf()
            
        else:
            eminlist.append(emin)
            emaxlist.append(emax)
            hgrid.append(hsig[Hmax])
            if hsig[Hmax]>=hbest:
                hbest=hsig[Hmax]
                eminbest=emin
                emaxbest=emax

if args.gridsearch or args.coarsegridsearch:

    pi_mask = (data[2]>eminbest*KEV_TO_PI) & (data[2]<emaxbest*KEV_TO_PI)
    pred_rate = 0.05/10.0 # 2241
    sn,sn0,hs,ph_gti,gti_rts_s,gti_len_s = make_sn(data_diced,mask=pi_mask,rate=pred_rate)
    hsig = [h2sig(h) for h in hs]
    exposure = np.cumsum(gti_len_s)
    
    plt.figure(5); plt.clf()
    # scale exposure to the expected S/N
    amax = np.argmax(sn)
    exposure_scale = sn[amax]/exposure[amax]**0.5
    
    Hmax = np.argmax(hsig)

    plt.plot(gti_rts_s,hsig,label='H-test significance')
    plt.axvline(gti_rts_s[amax],color='k',ls='--')
    plt.axvline(gti_rts_s[Hmax],color='r',ls='--')
    plt.xlabel('Background Rate (ct/s)')
    plt.ylabel('Significance (sigma)')
    plt.title('{} - [{},{}]'.format(args.name,eminbest,emaxbest))
    plt.legend(loc='lower right')
    plt.savefig('{}_sig_bestrange.png'.format(args.outfile))

    plt.clf()
    plt.scatter(eminlist,emaxlist, c=hgrid, s=10, edgecolor='')
    cbar = plt.colorbar()
    cbar.set_label('H-test')
    plt.xlabel('Low Energy Cut (keV)')
    plt.ylabel('High Energy Cut (keV)')
    plt.savefig('{}_grid.png'.format(args.outfile))
    plt.clf()
    
    nbins=args.nbins
    select_ph = np.concatenate(ph_gti[:Hmax]).ravel()
    profbins = np.linspace(0.0,1.0,nbins+1,endpoint=True)
    profile, edges = np.histogram(select_ph,bins=profbins)
    bbins = np.concatenate((profbins, profbins[1:]+1.0, profbins[1:]+2.0))
    fullprof = np.concatenate((profile,profile,profile,np.array([profile[0]])))
    plt.errorbar(bbins-(0.5/nbins),fullprof,
                 yerr=fullprof**0.5,
                 marker ='',
                 drawstyle='steps-mid',
                 linewidth=1.5,
                 color='k'
    )
    
    plt.xlim((0.0,2.0))
    plt.ylabel('Photons')
    plt.xlabel('Phase')
    plt.title(args.name)
    plt.savefig('{}_profile.png'.format(args.outfile))

    print("Maximum significance: {:0.3f} sigma".format(hsig[Hmax]))
    print("Maximum significance: {:0.3f} sigma".format(hbest))
    print("   obtained in {:0.3f} ksec".format(exposure[Hmax]))
    print("   for {} events".format(len(select_ph)))
    
else:
    
    print("Maximum significance: {:0.3f} sigma".format(hsig[Hmax]))
    print("   obtained in {:0.3f} ksec".format(exposure[Hmax]))
    print("   for {} events".format(len(select_ph)))
