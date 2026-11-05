import numpy as np
from skymaps import SkyDir
from uw.utilities.fitstools import rad_extract
from uw.like.pypsf import PsfOverlap
from scipy.integrate import simps
from cPickle import dump

class EdepAperture(object):
    def init(self):
        self.flaperture = np.radians(1)
    def __init__(self): pass
    def flat_aperture(en,ct): return self.flaperture

def calc_ratios(roi,edict,which=0):
    """ edict is a dictionary containing the event information.
        Calculate for each photon the ratio of the source to the background.
    """

    print """WARNING: This code doesn't seem to play well with pointlike in
        it's current state.  In particular, a proper way of interacting
        with diffuse sources needs to be established.  Recommend using the
        subset method."""
    en = edict['ENERGY']
    ct = edict['CONVERSION_TYPE']
    psm = roi.psm; bgm = roi.bgm
    dirs = map(SkyDir,edict['RA'],edict['DEC'])
    signals = np.empty(len(en))
    backs   = np.empty(len(en))
    psf = roi.sa.psf
    skydir = psm.point_sources[which].skydir
    if not hasattr(roi.sa,'ltcube'):
        psf.set_weights(roi.sa.dataspec.ltcube,skydir)
    else:
        psf.set_weights(roi.sa.ltcube,skydir)

    for i,(e,c,d) in enumerate(zip(en,ct,dirs)):
        psfluxes  = np.asarray([m(e) for m in psm.models])
        psdiffs   = np.asarray([d.difference(ps.skydir) for ps in psm.point_sources])
        psrates   = psfluxes * psf(e,c,psdiffs,density=True)
        if np.any(np.isnan(psrates)):
            print np.arange(len(psrates))[np.isnan(psrates)]
            raise ValueError('Sources got problems.')
        # NB preconvolved diffuse models -- something of a kluge
        bgrates   = np.asarray([m.get_dmodel(c).value(d,e)*m.smodel(e) for m in bgm.bgmodels])
        if np.any(np.isnan(bgrates)):
            print np.arange(len(bgrates))[np.isnan(bgrates)]
            raise ValueError('Sources got problems.')
        signals[i] = psrates[which]
        backs[i]   = (psrates.sum() - psrates[which] + bgrates.sum())

    return signals,backs

def calc_ratios2(roi,edict,which=0):
    """ edict is a dictionary containing the event information.
        Calculate for each photon the ratio of the source to the background.
    """

    #ens = [[b.emax for b in roi.bands if b.ct==0],[b.emax for b in roi.bands if b.ct==1]]
    ens = [[b.emax for b in roi.bands if b.ct==ct] for ct in [0,1]]
    bands = [[b for b in roi.bands if b.ct==ct] for ct in [0,1]]
    #pix_indices = [[np.asarray([b.b.index(x) for x in b.wsdl]) for b in roi.bands if b.ct==0],
                   #[np.asarray([b.b.index(x) for x in b.wsdl]) for b in roi.bands if b.ct==1]]
    pix_indices = [[np.asarray([b.b.index(x) for x in b.wsdl]) for b in roi.bands if b.ct==ct] for ct in [0,1]]
    #pix_helpers = [[np.arange(len(x)) for x in pix_indices[0]],[np.arange(len(x)) for x in pix_indices[1]]]
    pix_helpers = [[np.arange(len(x)) for x in pix_indices[ct]] for ct in [0,1]]

    def find(e,c,d):
        idx = np.searchsorted(ens[c],e)
        if idx == len(ens[c]):
            raise IndexError('I did not find the right energy')
        band = bands[c][idx]
        mask = band.b.index(d)==pix_indices[c][idx]
        if not np.any(mask):
            raise IndexError('Requested photon position not in ROI.')
        return band,(pix_helpers[c][idx])[mask]

    en = edict['ENERGY']
    ct = edict['CONVERSION_TYPE']
    dirs = map(SkyDir,edict['RA'],edict['DEC'])
    signals = np.empty(len(en))
    backs   = np.empty(len(en))

    for i,(e,c,d) in enumerate(zip(en,ct,dirs)):
        try:
            band,pix_idx = find(e,c,d)
        except IndexError: # this will happen if using an FT1 file with a superset of photons for ROI
            signals[i] = 0
            backs[i] = -1
            continue
        bgrates = band.bg_all_pix_counts[pix_idx] + band.ps_all_pix_counts[pix_idx]
        signals[i] = band.ps_counts[which]*band.ps_pix_counts[pix_idx,which]
        backs[i] = bgrates - signals[i]

    return signals,backs

def write_weights(roi,ft1files,colnames=['WEIGHT','TOTAL'],emin=100,which=0,subset_method=True):
    """ Write the signal and the total rates using the current state."""
    import pyfits
    if not hasattr(ft1files,'__iter__'):  ft1files = [ft1files]
    d = dict()
    for ft1 in ft1files:
        f = pyfits.open(ft1)
        keys = ['ENERGY','CONVERSION_TYPE','RA','DEC']
        for key in keys:
            d[key] = np.asarray(f[1].data.field(key),dtype=float if key!='CONVERSION_TYPE' else int)
        mask = d['ENERGY'] > emin
        for key in keys: d[key] = d[key][mask]
        if not subset_method:
            pre_signals,pre_backs = calc_ratios(roi,d,which=which)
        else:
            pre_signals,pre_backs = calc_ratios2(roi,d,which=which)
        signals = np.empty(len(mask))
        backs = np.empty(len(mask))
        signals[~mask] = 0; backs[~mask] = 1
        signals[mask] = pre_signals; backs[mask] = pre_backs
        d = f[1].data.names
        cols = []
        if (colnames[0] in d):
            f[1].data.field(colnames[0])[:] = signals/(signals+backs)
        else:
            cols += [pyfits.Column(name=colnames[0],format='E',array=signals/(signals+backs))]
        if (colnames[1] in d):
            f[1].data.field(colnames[1])[:] = signals+backs
        else:
            cols += [pyfits.Column(name=colnames[1],format='E',array=signals+backs)]
        if len(cols) > 0:
            tbhdu = pyfits.new_table(f[1].columns.data+cols,header=f[1].header)
            f[1] = tbhdu
        f.writeto(ft1,clobber=True)
        f.close()

class LCLikelihood(object):
    """ Implement a likelihood test using all available information about a set of photons.
        Uses the spectral fit from an ROIAnalysis to calculate the source and background rates
        for each photon.  Then supports a search in phase space for 1- and 2-peak light curves
        to find a maximum likelihood.
    """

    def init(self):
        self.phase_col_name  = 'PULSE_PHASE'
        self.scan_width      = 0.03 # width of gaussian to use in phase scan
        self.scan_resolution = 0.01  # steps in phase to take in phase scan
        self.period          = None
        self.which   = 0
        self.flaperture      = 1
        self.edepaperture    = None
        self.psrcataperture  = None
        self.prms            = 0
        self.eventclass      = 3
        self.cuts            = None # cuts on FT1 columns in radial extraction

        self.emin = 1e2
        self.emax = 1e4

    def __init__(self,roi,ft1files=None,**kwargs):
        self.init()
        self.__dict__.update(kwargs)
        self.roi = roi

        # if ft1file is None, try to get it from roi
        ft1files = ft1files if ft1files is not None else roi.sa.pixeldata.ft1files
        self.ft1files = ft1files

        skydir = roi.psm.point_sources[self.which].skydir

        # one can specify a pulsar catalog-like aperture which takes precedent
        if self.psrcataperture is not None:
            def r(e,event_class=0):
                return np.maximum(np.minimum(self.psrcataperture,0.8*(1000./e)**0.75),0.35)
        # one can pass a function taking as arguments energy and event conversion type
        elif self.edepaperture is not None: r = self.edepaperture
        # finally, the default is a cookie cutter aperture
        else: r = self.flaperture

        return_cols = ['TIME','CONVERSION_TYPE','EVENT_CLASS'] + ([self.phase_col_name] if self.period is None else [])
        edict = rad_extract(ft1files,skydir,r,return_cols=return_cols,cuts=self.cuts)
        mask = np.logical_and(edict['ENERGY'] > self.emin,edict['ENERGY'] < self.emax)
        mask = np.logical_and(edict['EVENT_CLASS'] >= self.eventclass, mask)
        for k,v in edict.iteritems(): edict[k] = v[mask]
        edict['CONVERSION_TYPE'] = edict['CONVERSION_TYPE'].astype(int)

        self.edict  = edict # save values for now for debugging; don't need in likelihood algorithm
        self.signals,self.backs = calc_ratios(roi,edict,which=self.which)
        self.weights = self.signals / (self.signals + self.backs)
        self.phases = edict[self.phase_col_name] if self.period is None else np.mod((edict['TIME']-edict['TIME'].min())/self.period,1)
        self.times  = edict['TIME']

    def pickle(self,filename):
        self.roi = None
        dump(self,file(filename,'wb'))

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

           
"""
def logchi2sf(test_statistic,dof):

   from scipy.stats import chi2
   from scipy.special import gamma
   dof = np.asarray(dof).astype(float)
   ts  = np.asarray(test_statistic).astype(float)

   p = chi2.sf(ts,dof)
   try:
      mask = p>0
      p[mask]= np.log(p)
      nmask = np.logical_not(p)
      p[nmask] = -ts[nmask]/2 + (dof[nmask]/2-1)*np.log(dof[nmask]/2) - np.log(gamma(dof[nmask]/2))
      return p
   except:
      if p >0: return np.log(p)
      else: return -ts/2 + (dof/2-1)*np.log(dof/2) - np.log(gamma(dof/2))

   # underflow conditions - use asymptotic expansion
   return p

def ts2sig(ts,dof):
    from scipy.stats import chi2
    from scipy.special import gamma

   if ts < 500:
      from uw.utilities.phasetools import sig2sigma
      from scipy.stats import chi2
      return sig2sigma(chi2.sf(ts,dof))

   from scipy.optimize import fsolve

   ts = float(ts)
   dof= float(dof)
   x0sq = ts - (dof - 2)*np.log(dof/2) + 2*np.log(gamma(dof/2))
   eta  = np.log(np.pi/2) - x0sq

   f    = lambda xsq: xsq + np.log(xsq) + eta

   return fsolve(f,x0sq)**0.5
"""
