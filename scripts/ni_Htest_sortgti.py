#!/usr/bin/env python
# Version: 1.0
# Author: M. Kerr (updated by S. Guillot)

from __future__ import division, print_function
import argparse
from collections import deque
import glob
import os
import sys
from subprocess import check_call

import astropy.units as u
from astropy import log
from astropy.io import fits
from astropy.time import Time
import numpy as np
from pint.fits_utils import read_fits_event_mjds
from pint.eventstats import h2sig,hm,sig2sigma
import scipy.stats
from scipy.stats import chi2

#from nicer.values import *
from nicer.values import datadir,KEV_TO_PI

desc= """
Read one or more event files 
to sort GTI by background rate 
and evaluate H-test
"""

def local_h2sig(h):
    h = np.atleast_1d(h)
    rvals = np.zeros_like(h)
    for ix,x in enumerate(h):
        if x > 0:
            rvals[ix] = h2sig(x)
    return rvals

parser = argparse.ArgumentParser(description = desc)
parser.add_argument("infile", help="file or text file with list of event file", nargs='+')
parser.add_argument("outfile", help="name for output files")
parser.add_argument("--emin", help="Minimum energy to include (keV, default=0.25)", type=float, default=0.25)
parser.add_argument("--emax", help="Maximum energy to include (keV, default=2.00)", type=float, default=2.00)
parser.add_argument("--maxemin", help="Maximum emin to use in fine grid search (keV, default=1.00)", type=float, default=1.00)
parser.add_argument("--minemax", help="Minimum emax to use in fine grid search (keV, default=1.00)", type=float, default=1.00)
parser.add_argument("--gridsearch", help="Search over energies to find max H-test", action="store_true",default=False)
parser.add_argument("--coarsegridsearch", help="Search over energies to find max H-test", action="store_true",default=False)
parser.add_argument("--minbw", help="Minimum fractional bandwidth used during energy grid searching.  E.g., --minbw=0.5 would allow a 1.0 to 1.5 (50%%) keV energy range, but not a 2.0 to 2.2 (10%%) range.", type=float,default=None)
parser.add_argument("--minexp", help="Minimum exposure allowed for a candidate cut, expressed as a fraction of the total.  E.g., --minexp=0.2 would allow a cut that throws away 80%% of the GTI.", type=float,default=None)
parser.add_argument("--mingti", help="Minimum GTI length to allow -- short GTIs don't give a reliable count rate. (seconds, default=10.0)", type=float, default=10.0)
parser.add_argument("--nopulsetest", help="Only use the predicted S/N to determine the GTIs to use.", action="store_true",default=False)
parser.add_argument("--verbosity", help="Verbosity (0=quiet,1=default,2=verbose,3=very verbose).", type=int, default=1)
parser.add_argument("--savefile", help="Saving optimized event file", action="store_true",default=False)
parser.add_argument("--writegti", help="Write out the GTI corresponding to the event selection.", action="store_true",default=False)
parser.add_argument("--writeevents", help="Write out the corresponding event file", action="store_true",default=False)
parser.add_argument("--remote", help="Disable interactive plotting backend", action="store_true",default=False)
parser.add_argument("--usez", help="Use Z^2_2 test instead of H test.", action="store_true",default=False)
parser.add_argument("--nbins", help="Number of bins for plotting pulse profile (default=16)", type=int, default=16)
parser.add_argument("--name", help="Pulsar name for output figure", type=str, default='')
args = parser.parse_args()

import matplotlib
if (args.remote):
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rc('font', size=14)          # controls default text sizes
plt.rc('axes', labelsize=13)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=13)    # fontsize of the tick labels
plt.rc('ytick', labelsize=13)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rc('axes', linewidth=1.5)
plt.rc('xtick.major', size=4, width=1.5)
plt.rc('ytick.major', size=4, width=1.5)

def merge_gtis(data):
    # gti starts and stops
    gti_start = data[-2]
    gti_stop = data[-1]
    # quick and dirty loop -- array
    out_gtis = deque()
    out_gtis.append([gti_start[0],gti_stop[0]])
    for start,stop in zip(gti_start[1:],gti_stop[1:]):
        # if start is same value is last stop, just update
        if start == out_gtis[-1][1]:
            out_gtis[-1][1] = stop
        else:
            out_gtis.append([start,stop])
    t0s,t1s = np.asarray(out_gtis).transpose()
    return data[0],data[1],data[2],t0s,t1s

def runcmd(cmd):
    # CMD should be a list of strings since it is not processed by a shell
    log.info('CMD: '+" ".join(cmd))
    os.system(" ".join(cmd))
    ## Some ftools calls don't work properly with check_call...not sure why!
    ## so I am using os.system instead of check_call
    #check_call(cmd,env=os.environ)

def load_files(fnames):
    """ Load in time stamps, PIs, GTIs, etc. from all files."""
    gtis = deque()
    times = deque()
    pis = deque()
    phases = deque()
    for fname in fnames:
        f = fits.open(fname)
        try:
            tzero = f['gti']._header['TIMEZERO']
        except KeyError:
            tzero = 0
        t0s = f['gti'].data.field('start') + tzero
        t1s = f['gti'].data.field('stop') + tzero
        for t0,t1 in zip(t0s,t1s):
            gtis.append([t0,t1])
        try:
            tzero = f['events']._header['TIMEZERO']
        except KeyError:
            tzero = 0
        times.append(f['events'].data.field('time') + tzero )
        pis.append(f['events'].data.field('pi'))
        try:
            phases.append(f['events'].data.field('pulse_phase'))
        except:
            pass
        f.close()

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

def write_gtis(gti_start, gti_stop, outfile, merge_gti=False):
    # write out GTIs -- sort them by time rather than bkg
    a = np.argsort(gti_start)
    gti_start = gti_start[a]
    gti_stop = gti_stop[a]
    # merge adjacent GTI -- quick and dirty loop
    out_gtis = deque()
    out_gtis.append([gti_start[0],gti_stop[0]])
    for start,stop in zip(gti_start[1:],gti_stop[1:]):
        # if start is same value is last stop, just update
        if merge_gti and (start == out_gtis[-1][1]):
            out_gtis[-1][1] = stop
        else:
            out_gtis.append([start,stop])
    out_gtis = np.asarray(out_gtis)
    np.savetxt("{}_OptimalGTI.txt".format(outfile),out_gtis)

    # Checking the presence of HEASOFT
    try:
        check_call('nicerversion',env=os.environ)
    except:
        print("You need to initialize FTOOLS/HEASOFT first (e.g., type 'heainit')!", file=sys.stderr)
        return

    # Checking the presence of gti header and columns in data/
    gticolumns = os.path.join(datadir,'gti_columns.txt')
    gtiheader = os.path.join(datadir,'gti_header.txt')

    # make sure TIMEZERO == 0
    lines = open(gtiheader).readlines()
    for line in lines:
        toks = line.split('=')
        if (len(toks) > 0) and (toks[0].strip()=='TIMEZERO'):
            if float(toks[1].strip().split()[0]) != 0:
                print('WARNING! TIMEZERO in output GTI not consistent.')
            break

    if not os.path.isfile(gtiheader) or not os.path.isfile(gticolumns):
        log.error('The files gti_header.txt or gti_columns.txt are missing.  Check the {} directory'.format(os.path.abspath(datadir)))
        exit()
        
    ## Making the GTI file from the text file
    log.info("Making the GTI file gti.fits from the GTI data textfile")
    cmd = ['ftcreate', '{}'.format(gticolumns), '{}_OptimalGTI.txt'.format(outfile), '{}_OptimalGTI.fits'.format(outfile), 'headfile={}'.format(gtiheader), 'extname="GTI"', 'clobber=yes']
    runcmd(cmd)

    ## Extracting the new event file using the new GTI file created
    if args.writeevents:
        if len(args.infile)==1:
            eventfile = args.infile[0]
            outevtfile = "{}_OptimalEvents.fits".format(outfile)
            cmd = ['niextract-events', '{0}'.format(eventfile), '{0}'.format(outevtfile), 'timefile={}_OptimalGTI.fits'.format(outfile), 'clobber=yes']
            runcmd(cmd)
        else:
            log.warning("Cannot create events file. niextract-events needs a single file or a list of events files (@list.txt)")

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

def ensemble_ztest(phases,slices,m=2):
    """ Calculate H-test statistic for subsets of a set of phases.
        Cache intermediate products to avoid O(N^2) complexity!
    """
    
    phases = np.asarray(phases)*(2*np.pi) # in radians and copy
    cache = np.empty((2*m,len(phases)))
    for i in range(m):
        cache[2*i] = np.cos((i+1)*phases)
        cache[2*i+1] = np.sin((i+1)*phases)

    rvals = np.empty(len(slices))
    for isl,sl in enumerate(slices):
        x = cache[:,sl]
        nph = x.shape[1]
        if nph == 0:
            rvals[isl] = 0.01
            continue
        t = np.sum(cache[:,sl],axis=1)**2
        t = np.sum(t[::2] + t[1::2])*(2./nph)
        rvals[isl] = t
    return rvals

def make_sn(data,mask=None,rate=0.1,min_gti=10,usez=False,snonly=False):
    """ data -- output of load_local
        mask -- optional mask to select events (e.g. on PI)
        rate -- assumed rate for S/N calculation in ct/sec
        min_gti -- minimum GTI length in seconds
        usez -- use Z^2 test instead of H test
        snonly -- skip computation of pulsed statistic, only do S/N
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
    gti_t0_s = t0s[mask][a]
    gti_t1_s = t1s[mask][a]
    gti_len_s = gti_len[mask][a]
    gti_cts_s = gti_cts[mask][a]
    gti_rts_s = gti_cts_s/gti_len_s

    sn = rate*np.cumsum(gti_len_s)/np.sqrt(np.cumsum(gti_cts_s)+rate*gti_len_s)
    rate = 0
    zero_mask = gti_cts_s > 0
    sn0 = np.empty(len(gti_len_s))
    sn0[zero_mask] = np.cumsum(gti_len_s[zero_mask])/np.sqrt(np.cumsum(gti_cts_s[zero_mask]))
    sn0[~zero_mask] = np.inf

    counter = 0
    pi_gti = deque()
    for i,ct in enumerate(gti_cts):
        pi_gti.append(pis[counter:counter+ct])
        counter += ct
    # apply mask
    pi_gti = [pi_gti[i] for i,imask in enumerate(mask) if imask]
    # apply sorting
    pi_gti = [pi_gti[i] for i in a]

    if (not snonly) and (phases is not None):
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
        if usez:
            hs = ensemble_ztest(np.concatenate(ph_gti),slices)
        else:
            hs = ensemble_htest(np.concatenate(ph_gti),slices)

    else:
        hs = None
        ph_gti = None

    return sn,sn0,hs,ph_gti,list(pi_gti),gti_rts_s,gti_len_s,gti_t0_s,gti_t1_s




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
#print('There are %d GTIs.'%(len(data[-1])))
#data = merge_gtis(data)
#print('There are now %d GTIs.'%(len(data[-1])))
data_diced = dice_gtis(data)
#import sys; sys.exit(0)

if args.writeevents:
    args.writegti=True

if args.gridsearch:
    all_emin = np.arange(max(0.24,args.emin),args.maxemin+0.005,0.01)
elif args.coarsegridsearch:
    all_emin = np.arange(max(0.24,args.emin),2.0,0.1)
else:
    all_emin = np.array([args.emin])

hbest = 0.0
eminbest = 0.0
emaxbest = 100.0
dosearch = args.gridsearch or args.coarsegridsearch

eminlist = []
emaxlist = []
hgrid = []
    
for emin in all_emin:
    
    if args.gridsearch:
        delta_e = 0.02
        all_emax = np.arange(args.minemax,min(3.0,args.emax)+0.005,delta_e)
    elif args.coarsegridsearch:
        delta_e = 0.1
        hilimit = min(7.0,args.emax)
        all_emax = np.arange(emin+0.1,hilimit+0.01,delta_e)
    else:
        delta_e = 0
        all_emax = np.array([args.emax])

    if len(all_emax) == 0:
        break

    if (args.verbosity >= 1):
        print("emin={:0.2f}, emax ranging from {:0.2f}-{:0.2f} by {:0.2f} keV".format(emin,all_emax[0],all_emax[-1],delta_e))
        

    for emax in all_emax:

        if (args.verbosity >= 3):
            print("    emin={:0.2f}, emax={:0.2f}".format(emin,emax))

        # test for energy bandwidth
        if dosearch and (args.minbw is not None):
            if emax/emin-1 < args.minbw:
                if (args.verbosity >= 2):
                    print('    excluding emin={:0.2f}, emax={:0.2f} because smaller than specified minbw'.format(emin,emax))
                continue

        pi_mask = (data[2]>emin*KEV_TO_PI) & (data[2]<emax*KEV_TO_PI)
        pred_rate = 0.05/10.0 # 2241
        sn,sn0,hs,ph_gti,pi_gti,gti_rts_s,gti_len_s,gti_t0_s,gti_t1_s = \
            make_sn(data_diced,mask=pi_mask,rate=pred_rate,usez=args.usez,
                    min_gti=args.mingti)

        exposure = np.cumsum(gti_len_s)
        #idx = np.searchsorted(exposure,np.arange(10)*1./10*exposure[-1])
        #print(gti_rts_s[idx])

        if args.minexp is not None:
            # depress S/N for values that do not satisfy eposure cuts
            mask = exposure/exposure[-1] < args.minexp
            sn[mask] = 0
            sn0[mask] = 0
            hs[mask] = 0

        if args.usez:
            hsig = sig2sigma(chi2.sf(hs,4))
        else:
            hsig = local_h2sig(hs)
        
        # scale exposure to the expected S/N
        amax = np.argmax(sn)
        exposure_scale = sn[amax]/exposure[amax]**0.5

        if args.nopulsetest:
            Hmax = amax
        else:
            Hmax = np.argmax(hsig)
        
        if not dosearch:
            # Make output plots after single iteration.
            plt.figure(5); plt.clf()
            if args.usez:
                plt.plot(gti_rts_s,hsig,label='Z-test significance')
            else:
                plt.plot(gti_rts_s,hsig,label='H-test significance')
            plt.axvline(gti_rts_s[amax],color='k',ls='--',label='No  H-test (sig={:0.3f})'.format(hsig[amax]))
            if not args.nopulsetest:
                plt.axvline(gti_rts_s[Hmax],color='r',ls='--',label='Max H-test (sig={:0.3f})'.format(hsig[Hmax]))
            plt.xlabel('Background Rate (ct/s)')
            plt.ylabel('Significance (sigma)')
            plt.title('{} - [{:0.2f},{:0.2f}]'.format(args.name,emin,emax))
            plt.legend(loc='lower right')
            plt.savefig('{}_sig.png'.format(args.outfile))
            
            plt.clf()
            nbins=args.nbins
            select_ph = np.concatenate(ph_gti[:Hmax+1]).ravel()
            profbins = np.linspace(0.0,1.0,nbins+1,endpoint=True)
            profile, edges = np.histogram(select_ph,bins=profbins)
            bbins = np.concatenate((profbins, profbins[1:]+1.0, profbins[1:]+2.0))
            fullprof = np.concatenate((profile,profile,profile,np.array([profile[0]])))
            plt.errorbar(bbins-(0.5/nbins),fullprof,
                         yerr=fullprof**0.5,
                         marker ='',
                         drawstyle='steps-mid',
                         linewidth=1.5,
                         color='k',
                         label= '{:0.2f}-{:0.2f} keV'.format(emin,emax)
            )
            #plt.subplots_adjust(left=0.15, right=0.93)  #, bottom=0.1)
            plt.tight_layout()
            plt.xlim((0.0,2.0))
            plt.ylabel('Photons')
            plt.xlabel('Phase')
            plt.title(args.name)
            plt.legend(loc='upper left')
            plt.savefig('{}_profile.png'.format(args.outfile))
            plt.clf()

            if args.writegti:
                write_gtis(gti_t0_s[:Hmax+1],gti_t1_s[:Hmax+1],args.outfile)
            
        else:
            # store data for future comparison
            eminlist.append(emin)
            emaxlist.append(emax)
            hgrid.append(hsig[Hmax])
            if hsig[Hmax]>=hbest:
                hbest=hsig[Hmax]
                eminbest=emin
                emaxbest=emax

if dosearch:

    # recreate data -- really need to encapsulate this!

    pi_mask = (data[2]>eminbest*KEV_TO_PI) & (data[2]<emaxbest*KEV_TO_PI)
    pred_rate = 0.05/10.0 # 2241
    sn,sn0,hs,ph_gti,pi_gti,gti_rts_s,gti_len_s,gti_t0_s,gti_t1_s = \
        make_sn(data_diced,mask=pi_mask,rate=pred_rate,usez=args.usez,
                min_gti=args.mingti)

    exposure = np.cumsum(gti_len_s)

    if args.minexp is not None:
        # depress S/N for values that do not satisfy eposure cuts
        mask = exposure/exposure[-1] < args.minexp
        sn[mask] = 0
        sn0[mask] = 0
        hs[mask] = 0

    if args.usez:
        hsig = sig2sigma(chi2.sf(hs,4))
    else:
        hsig = local_h2sig(hs)
    
    # scale exposure to the expected S/N
    amax = np.argmax(sn)
    exposure_scale = sn[amax]/exposure[amax]**0.5
    
    if args.nopulsetest:
        Hmax = amax
    else:
        Hmax = np.argmax(hsig)

    plt.figure(5); plt.clf()
    if args.usez:
        plt.plot(gti_rts_s,hsig,label='Z-test significance')
    else:
        plt.plot(gti_rts_s,hsig,label='H-test significance')
    plt.axvline(gti_rts_s[amax],color='k',ls='--',label='No  H-test (sig={:0.3f})'.format(hsig[amax]))
    if not args.nopulsetest:
        plt.axvline(gti_rts_s[Hmax],color='r',ls='--',label='Max H-test (sig={:0.3f})'.format(hsig[Hmax]))
    plt.xlabel('Background Rate (ct/s)')
    plt.ylabel('Significance (sigma)')
    plt.title('{} - [{:0.2f},{:0.2f}]'.format(args.name,eminbest,emaxbest))
    plt.legend(loc='lower right')
    plt.savefig('{}_sig_bestrange.png'.format(args.outfile))

    plt.clf()
    plt.scatter(eminlist,emaxlist, c=hgrid, s=10, edgecolor='')
    cbar = plt.colorbar()
    if args.usez:
        cbar.set_label('Z-test')
    else:
        cbar.set_label('H-test')
    plt.xlabel('Low Energy Cut (keV)')
    plt.ylabel('High Energy Cut (keV)')
    plt.savefig('{}_grid.png'.format(args.outfile))
    plt.clf()
    
    nbins=args.nbins
    select_ph = np.concatenate(ph_gti[:Hmax+1]).ravel()
    profbins = np.linspace(0.0,1.0,nbins+1,endpoint=True)
    profile, edges = np.histogram(select_ph,bins=profbins)
    bbins = np.concatenate((profbins, profbins[1:]+1.0, profbins[1:]+2.0))
    fullprof = np.concatenate((profile,profile,profile,np.array([profile[0]])))
    plt.errorbar(bbins-(0.5/nbins),fullprof,
                 yerr=fullprof**0.5,
                 marker ='',
                 drawstyle='steps-mid',
                 linewidth=1.5,
                 color='k',
                 label= '{:0.2f}-{:0.2f} keV'.format(eminbest,emaxbest)
    )
    
    plt.xlim((0.0,2.0))
    plt.ylabel('Photons')
    plt.xlabel('Phase')
    plt.title(args.name)
    plt.savefig('{}_profile.png'.format(args.outfile))

    print("Maximum significance: {:0.3f} sigma".format(hbest))
    print("   obtained in {:0.2f} (out of {:0.2f} ksec)".format(
        exposure[Hmax]/1000,exposure[-1]/1000))
    print("   between {:0.2f} and {:0.2f} keV".format(eminbest,emaxbest))
    print("   for {} events".format(len(select_ph)))

    plt.savefig('{}_profile.png'.format(args.outfile))

    if args.savefile:
        out_ph = np.concatenate(ph_gti[:Hmax+1])
        out_pi = np.concatenate(pi_gti[:Hmax+1])
        output = np.asarray([out_pi,out_ph]).transpose()
        np.savetxt('%s_optimzed.evt.gz'%(args.outfile),
                output,fmt='%04d %1.6f')

    if args.writegti:
        write_gtis(gti_t0_s[:Hmax+1],gti_t1_s[:Hmax+1],args.outfile)
    
else:
    
    print("Maximum significance: {:0.3f} sigma".format(hsig[Hmax]))
    print("   obtained in {:0.2f} ksec (out of {:0.2f} ksec)".format(exposure[Hmax]/1000,exposure[-1]/1000))
    print("   for {} events".format(len(select_ph)))

