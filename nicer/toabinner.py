"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/toabinner.py,v 1.4 2013/03/30 21:14:41 kerrm Exp $

handle binning of data for extraction of TOAs

Authors: Paul S. Ray <paul.ray@nrl.navy.mil>
         Matthew Kerr <matthew.kerr@gmail.com>
"""
import numpy as np
import sys

class TOABinner(object):
    """ Implement an interface to calculate bin edges for a span of LAT
        observations.  Note the default implementation is uniform in time."""

    def __init__(self,ntoa,phase_data,template=None):
        self.counter = 0
        self.ntoa = ntoa
        self.data = phase_data
        self.template = template # not needed for standard calculation
        self.make_edges()

    def __iter__(self):
        self.counter = 0
        return self

    def next(self):
        if self.counter >= self.ntoa: raise StopIteration
        idx,self.counter = self.counter,self.counter+1
        return self.starts[idx],self.stops[idx]

    def make_edges(self):
        bin_edges = np.linspace(self.data.mjd_start,self.data.mjd_stop,self.ntoa+1)
        self.starts = bin_edges[:-1]
        self.stops  = bin_edges[1:]

class UniformLogBinner(TOABinner):
    """ Uses the information in the profile to bin the data in uniform "log drops".
        This is based on a likelihood ratio test for the pulsation significance
        and should do a good job of dividing the data up into regions that give
        approximately equal power in estimating a TOA.
        
        Note on the formulation below --
        The null hypothesis is uniform emission, so the log likelihood for
        this hypothesis is just n_photons*log(constant) = constant.  So the choice
        of constant just sets the 0 point, and it's convenient to have
        everything positive.  For an unweighted profile, this means taking the
        log of the minimum of the alternate hypothesis, or 1 - norm.
        
        For the weighted likelihood, there isn't a closed form for the min,
        so we just subtract off the empirical min."""
    
    def make_edges(self):
        phases,weights,mjds = self.data.ph,self.data.weights,self.data.mjds
        if (weights is None): 
            logs = np.log(self.template(phases))
            if self.template.norm() < 1:
                logs -= np.log(1-self.template.norm())
            else: # for KD etc. where the DC component isn't analytic
                logs -= logs.min()
        else:
            logs = np.log(1+weights*(self.template(phases)-1))
            logs -= logs.min()
        delta = logs.sum()/self.ntoa
        my_delta = 0
        index = 1
        indices = np.asarray([-1]*(self.ntoa+1))
        for i in xrange(len(phases)):
            if my_delta > delta:
                indices[index] = i
                index += 1
                my_delta = my_delta - delta
            my_delta += logs[i]
        indices[-1] = len(phases)-1
        indices[0] = 0
        assert(not np.any(indices==-1))
        self.starts = np.empty(self.ntoa); self.stops = np.empty_like(self.starts)
        for i in xrange(self.ntoa):
            if i == 0:
                self.starts[i] = mjds[0] - 1
            else:
                self.starts[i] = (mjds[indices[i]] + mjds[indices[i]-1])/2
            if i == (self.ntoa-1):
                self.stops[i] = mjds[-1] + 1
            else:
                self.stops[i] = (mjds[indices[i+1]] + mjds[indices[i+1]-1])/2
            
class PrebinnedBinner(TOABinner):
    def __init__(self,starts,stops):
        self.counter = 0
        self.ntoa = len(starts)
        self.starts = starts
        self.stops = stops
