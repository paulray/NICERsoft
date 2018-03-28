#!/usr/bin/env python
from __future__ import division, print_function
import argparse
import numpy as np
from pylab import *
import astropy.io.fits as pyfits
from nicer.values import PI_TO_KEV

parser = argparse.ArgumentParser(description="Read event file with PULSE_PHASE and PI columns and make 2-D plot")
parser.add_argument("evfile", help="Name of event files to process")
parser.add_argument("--outfile", help="Name of output file for plot")
parser.add_argument("--nbins", help="Number of phase bins", type=int, default=32)
parser.add_argument("--emin", help="Min energy to plot", type=float, default=0.3)
parser.add_argument("--emax", help="Max energy to plot", type=float, default=2.0)
parser.add_argument("--nebins", help="Number of phase bins", type=int, default=35)
args = parser.parse_args()

pi = pyfits.getdata(args.evfile,'EVENTS').field('PI')
phase = pyfits.getdata(args.evfile,'EVENTS').field('PULSE_PHASE')

pi2 = np.append(pi,pi)
pimin = int(args.emin/PI_TO_KEV)
pimax = int(args.emax/PI_TO_KEV)
downsamp = int(1+(pimax-pimin)/args.nebins)
pirange=pimin+arange(args.nebins)*downsamp

phase2 = np.append(phase,phase+1.0)
phrange = arange(args.nbins*2,dtype=np.float)/float(args.nbins)

H,piedges,phedges = np.histogram2d(pi2,phase2,bins=[pirange,phrange])
eedges = piedges*PI_TO_KEV

pcolormesh(phedges,eedges,H,cmap="jet")
xlabel("Pulse Phase")
ylabel("Energy (keV)")
cbar = colorbar()
cbar.set_label("Photon count")

if args.outfile is not None:
    savefig(args.outfile)
else:
    show()
