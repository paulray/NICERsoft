"""
Calculate TOAs with a variety of methods.

Authors: Paul S. Ray <paul.ray@nrl.navy.mil>
         Matthew Kerr <matthew.kerr@gmail.com>
"""

import sys
import numpy as np
import pylab as pl
import scipy.stats
from scipy.optimize import fmin,golden,brentq
from stats import hm,hmw,sf_hm
from collections import deque

SECSPERDAY = 86400.

class TOAGenerator(object):
    """Manage a data set and set of options to produce the required TOAs from photon events."""

    def init(self):
        pass

    def __init__(self,data,polyco,template,**kwargs):

        self.data = data; self.polyco = polyco; self.template = template
        self.init()
        self.__dict__.update(**kwargs)
        self.mean_err = 1
        # keep track of errors in phase shift for later analysis
        self.phases = deque()
        self.phase_errs = deque()

    def get_toas(self,binner,use_midpoint=True):
        """ Calculate the TOAs specified by the binner."""

        self.phases.clear(); self.phase_errs.clear()

        # Compute observation duration for each TOA
        toas = np.empty(binner.ntoa)
        err_toas = np.empty(binner.ntoa)
        tim_strings = ['FORMAT 1']
        self.counter = 0
        tmp_tim_strings = ['FORMAT 1']

        for ii,(mjdstart,mjdstop) in enumerate(binner):

            # Compute freq and period at middle of observation
            tmid = (mjdstop + mjdstart)/2.
            pe   = self.polyco.getentry(tmid)
            freq = pe.evalfreq(tmid)
            period = 1.0/freq
            if period < 0:
                raise ValueError(
                    'Something went horribly wrong with the folding period.')

            # Compute phase at start of observation or at midpoint
            phase_time = tmid if use_midpoint else mjdstart
            pe = self.polyco.getentry(phase_time)
            polyco_phase0 = pe.evalphase(phase_time)

            # Select phases
            phases,weights = self.data.toa_data(mjdstart,mjdstop)
            if len(phases) == 0: continue

            tau,tau_err,prob,logl = \
                self.get_phase_shift(phases,weights,polyco_phase0)
            self.mean_err = (self.mean_err*ii + tau_err)/(ii+1)

            # compute RMS between linear assumption and actual polyco
            # over the integration period
            #dom = np.linspace(mjdstart,mjdstop,100)
            #ph = self.polyco.vec_evalabsphase(dom)
            #ph -= ph.min()
            #x = dom - tmid
            #p = np.polyfit(x,ph,1)
            #rms = (ph-np.polyval(p,x)).std()*(1e6/freq) # mus
            rms = 0.

            # Prepare a string to write to a .tim file or to send to STDOUT
            #print period
            toa = phase_time + (tau*period)/SECSPERDAY
            if tau_err < 100:
                toa_err = tau_err*period*1.0e6
            else:
                toa_err = 1e7 # hard code to 10s errors for nondetections
            frac_err = tau_err
            frame_label = 'BAT' if self.data.bary else 'GEO'
            weight_string = '' if (weights is None) else '-nwp %.2f'%(weights.sum())
            duration_string = '-tstart %s -tstop %s'%(mjdstart,mjdstop)
            rms_string = '-drms %.2f'%(rms)
            logl_string = '-logl %.2f'%(logl)
            s = " %s 0.0 %.12f %.2f %s -i LAT %s -np %d %s -chanceprob %.2e -pherr %.3f %s %s" % (frame_label,toa,toa_err,pe.obs,duration_string,len(phases),weight_string,prob,frac_err,rms_string,logl_string)
            toas[ii] = toa
            err_toas[ii] = toa_err
            tim_strings.append(s)
            self.counter += 1

        # Note TOAS in MJD, err_toas in microseconds, tim_strings a line for a FORMAT 1 .tim file
        return toas,err_toas,tim_strings

class UnbinnedTOAGenerator(TOAGenerator):

    def init(self):
        self.seeds = np.arange(0.00,0.991,0.01)
        self.good_ephemeris = True
        self.phi0 = self.template.get_location()
        self.prev_peak = self.phi0
        self.phase_drift = False
        self.weights = None
        self.plot_stem = None
        self.display = True
        self.likelihood_threshold = 5

    def __toa_error__(self,val,*args):
        f      = self.__toa_loglikelihood__
        delta  = 0.01
        d2 = (f( [val + delta], *args) - 2*f([val], *args) + f( [val - delta], *args))/(delta)**2
        if d2 < 0:
            print 'WARNING! Could not estimate error properly.  Setting to a large value...'
            return 0.2
        return d2**-0.5

    def __toa_loglikelihood__(self,p,*args):
        self.template.set_overall_phase(p[0])
        if args[1] is None:
            return -np.log(self.template(args[0])).sum()
        #return -np.log(1+args[1]*(self.template(args[0],suppress_bg=True)-1)).sum()
        return -np.log(1+args[1]*(self.template(args[0])-1)).sum()

    def get_phase_shift(self,phases,weights,polyco_phase0):

        f   = self.__toa_loglikelihood__
        jump = False
        if self.plot_stem is not None:
            name = ('%3d'%(self.counter+1)).replace(' ','0')
            plot_output = '%s%s.eps'%(self.plot_stem,name)
        else:
            plot_output = None
        # if the ephemeris is deemed good, we use it to track the
        # solution such that we don't adopt a fluctuation for the TOA
        if (self.good_ephemeris):
            seed = self.prev_peak
        else:
            seed = None
        x0,x0_err,best_ll = profile_analysis(
            f,(phases,weights),pred_phase=seed,plot_output=plot_output,
            thresh=self.likelihood_threshold)
        if x0_err < 1e2:
            self.prev_peak = x0
        else:
            x0 = self.prev_peak # track solution with "fake" TOA

         # the phase of profile relative to fiducal phase of the template
        peak_shift = x0-self.phi0

        # check for phase wraps
        if peak_shift < 0:
            # check for TOA wrapped higher
            if abs(peak_shift+1) < abs(peak_shift):
                peak_shift += 1
        elif peak_shift > 0:
            # check for predicted phase wrapped higher
            if abs(peak_shift-1) < abs(peak_shift):
                peak_shift -= 1

        # the phase relative to the reference phase, computed from the
        # reference time within the integration
        tau = peak_shift - polyco_phase0
        tau_err = x0_err
        if self.display:
            print 'Peak Shift: %.5f +/- %.5f'%(peak_shift,tau_err)
        self.phases.append(peak_shift)
        self.phase_errs.append(tau_err)
        h = hm(phases) if (weights is None) else hmw(phases,weights)
        return tau,tau_err,sf_hm(h),best_ll

def profile_analysis(logl,logl_args,pred_phase=None,nsamp=100,thresh=5,
    plot_output=None,max_jump=0.25):

    # (0) establish profile
    f = lambda x: logl([x],*logl_args)
    dom = np.linspace(0,1,nsamp+1)[:-1]
    cod = np.asarray(map(f,dom))

    # (1) find all local minima
    mask = (cod < np.roll(cod,1)) & (cod < np.roll(cod,-1))

    # (2) require that all local minima surpass a likelihood threshold
    m2 = mask & (cod < -abs(thresh))

    # (3) if no peaks satisfy conditions, allow global minimum
    if m2.sum() > 0:
        mask = m2

    # (4) if given a predicted phase, choose the peak closest to it as
    # TOA; otherwise, the global minimum
    if (mask.sum() > 1) and (pred_phase is not None):
        diffs = np.abs(dom-pred_phase)
        diffs = np.minimum(diffs,1-diffs)
        #d1 = np.abs(dom-pred_phase)
        #d2 = np.abs(pred_phase+(1-dom))
        #d3 = np.abs(dom-(pred_phase-1))
        #diffs = np.minimum(np.minimum(d1,d2),d3)
        idx = np.argmin(diffs[mask])
    else:
        idx = np.argmin(cod[mask])
    idx = np.arange(nsamp)[mask][idx] # index into main array

    # TODO -- something more sophisticated for 0.5 aliasing -- perhaps
    # allow less significant peaks... or possibly "search" at half the
    # frequency...

    # (5) find the minimum
    # define a shifted likelihood function to avoid phase wraps
    phi0 = dom[idx]
    g = lambda x: f(x+phi0)
    xmin = golden(g,brack=[-1./nsamp,0,1./nsamp])
    fmin = g(xmin)
    phi0 += xmin

    # (6) find the error bounds, slowly but surely
    ldiffs = cod - cod[idx] - 2 # cod[idx] vs. fmin is conservative
    h = lambda x: f(x+phi0) - fmin - 2

    rt_mask = np.roll(ldiffs,-idx) > 0
    if not np.any(rt_mask):
        rt = float(nsamp-1)/2
    else:
        rt_diff = np.arange(nsamp)[rt_mask][0]
        rt_brack = float(rt_diff)/nsamp-xmin
        if h(rt_brack) < 0:
            rt_brack *= 1.1 # hopefully deal with numerical slop
        rt = brentq(h,0,rt_brack)

    lt_mask = np.roll(ldiffs,nsamp-idx-1)[::-1] > 0
    if not np.any(lt_mask):
        lt = float(nsamp-1)/2
    else:
        lt_diff = np.arange(nsamp)[lt_mask][0]
        lt_brack = -float(lt_diff)/nsamp-xmin
        if h(lt_brack) < 0:
            lt_brack *= 1.1 # hopefully deal with numerical slop
        lt = brentq(h,lt_brack,0)

    # (7) construct TOA as the mean of the error positions
    tau = (rt+lt)/2 + phi0
    tau_err = (rt-lt)/4 # by 2 for average, by 2 again for 2->1 sigma

    if plot_output is not None:
        pl.rcParams['xtick.labelsize'] = 'large'
        pl.rcParams['ytick.labelsize'] = 'large'
        dom1 = dom
        cod1 = cod
        inset_delta = max(abs(rt),abs(lt))+0.01
        dom2 = np.linspace(tau-inset_delta,tau+inset_delta,30)
        cod2 = np.asarray(map(f,dom2))
        pl.figure(10); pl.clf();
        ax1 = pl.axes([0.10,0.10,0.88,0.87])
        #ax1.axhline(0,color='k')
        ax1.axhline(0,0,0.2,color='k') # hard hack for publn
        # calculate coordinates for inset; should use transAxes?
        #ax2 = pl.axes([0.18+0.4*(tau<0.5),0.58,0.28,0.28])
        ax2 = pl.axes([0.05+0.4*(tau<0.5),0.18,0.40,0.40])
        for i,(dom,cod,ax) in enumerate(zip([dom1,dom2],[cod1,cod2],[ax1,ax2])):
            if i != 1:
                ax.plot(dom,cod,color='gray',lw=2)
            # omit these lines for publication
            #ax.axvline(phi0,color='red',ls='--')
            #ax.axvline(tau,color='red')
            #ax.axvline(tau-tau_err,color='green',ls='-')
            #ax.axvline(tau+tau_err,color='green',ls='-')
            #if pred_phase is not None:
            #    ax.axvline(pred_phase,color='k',ls='-')
            if i==1:
                y = cod-cod.min()
                ax.plot(dom,y,color='gray',lw=2)
                ax.axvline(phi0,color='salmon',ls='--',lw=2)
                ax.axvline(tau,color='k',lw=2)
                ax.xaxis.set_major_locator(pl.matplotlib.ticker.MaxNLocator(4))
                #ax.axis([dom[0],dom[-1],cod.min(),cod.min()+5])
                ax.axis([dom[0],dom[-1],0,5])
                #ax.axhline(fmin+0.5,color='blue',ls='--')
                #ax.axhline(fmin+2.0,color='blue',ls='-.')
                #ax.axhline(0.5,color='k',ls='-')
                ax.axhline(2.0,color='k',ls='-')
                ax.axvline(tau-tau_err,color='k',ls='-')
                ax.axvline(tau+tau_err,color='k',ls='-')
                #ax.axvline(xmin-2*tau_err,color='green',ls='-.')
                #ax.axvline(xmin+2*tau_err,color='green',ls='-.')
                #pl.axvline(np.mod(rt+phi0,1),color='purple',ls='-')
                #pl.axvline(np.mod(lt+phi0,1),color='purple',ls='-')
        ax1.set_xlabel('Rotational Phase',size='x-large')
        ax1.set_ylabel('Negative Log Likelihood',size='x-large')
        ax2.set_ylabel('Rel. Log Likelihood',size='large')
        pl.savefig(plot_output)

    if (m2.sum() == 0):
        tau_err = 100
    # this is to catch aliases and prevent following TOAs from having
    # incorrect seed phase
    if pred_phase is not None:
        diff = abs(tau-pred_phase)
        diff = min(diff,1.-diff)
        if diff > max_jump:
            tau_err = 100
            tau = pred_phase
    return tau,tau_err,fmin
