#!/usr/bin/env python
import argparse
import numpy as np
from os import path
from pylab import *
import astropy.io.fits as pyfits
from nicer.values import PI_TO_KEV
from nicer.values import datadir

parser = argparse.ArgumentParser(
    description="Read event file with PULSE_PHASE and PI columns and make 2-D plot"
)
parser.add_argument("evfile", help="Name of event files to process")
parser.add_argument("--outfile", help="Name of output file for plot")
parser.add_argument("--nbins", help="Number of phase bins", type=int, default=32)
parser.add_argument("--emin", help="Min energy to plot", type=float, default=0.3)
parser.add_argument("--emax", help="Max energy to plot", type=float, default=2.0)
parser.add_argument("--nebins", help="Number of energy bins", type=int, default=35)
parser.add_argument("--log", help="Log scale histogram", action="store_true")
parser.add_argument("--sqrt", help="Log scale histogram", action="store_true")
parser.add_argument("--norm_effarea", help="Normalized by effective area", action="store_true")
parser.add_argument("--arfile", help="ARF file for effective area normalization", action="store_true")
args = parser.parse_args()

pi = pyfits.getdata(args.evfile, "EVENTS").field("PI")
phase = pyfits.getdata(args.evfile, "EVENTS").field("PULSE_PHASE")

pi2 = np.append(pi, pi)
pimin = int(args.emin / PI_TO_KEV)
pimax = int(args.emax / PI_TO_KEV)
downsamp = int(1 + (pimax - pimin) / args.nebins)
pirange = pimin + arange(args.nebins) * downsamp
print("PI bins ", pirange)

phase2 = np.append(phase, phase + 1.0)
phrange = arange(args.nbins * 2 + 1, dtype=float) / float(args.nbins)

H, piedges, phedges = np.histogram2d(pi2, phase2, bins=[pirange, phrange])
eedges = piedges * PI_TO_KEV

if args.norm_effarea:
    ## TO DO -- read from user provided ARF, other wise take ARF from /data/
    if args.arfile:
        arfile = args.arfile
    else:
        arfile = path.join(datadir, "nicer.arf")

    effarea = pyfits.getdata(arfile, "SPECRESP").field("SPECRESP")
    try:
        effarea_energy = pyfits.getdata(arfile, "SPECRESP").field("ENERGY")
    except:
        effarea_elow = pyfits.getdata(arfile, "SPECRESP").field("ENERG_LO")
        effarea_ehigh = pyfits.getdata(arfile, "SPECRESP").field("ENERG_HI")
        effarea_energy = (effarea_ehigh + effarea_elow) / 2.

    for e, ebin in enumerate(eedges[:-1]):
        # Calculate energy bin central energy
        energy = (ebin+eedges[e+1])/2
        # linear interp with the effective area curve
        eff = np.interp(energy, effarea_energy, effarea)
        # Normalize each "row" of the 2D histogram by the effective area at the central energy of the bin
        H[e,:] = H[e,:]/eff

if args.log:
    pcolormesh(phedges, eedges, np.log10(H), cmap="jet")
    cbar = colorbar()
    cbar.set_label("Log10(Photon count)")
elif args.sqrt:
    pcolormesh(phedges, eedges, np.sqrt(H), cmap="jet")
    cbar = colorbar()
    cbar.set_label("Sqrt(Photon count)")
else:
    pcolormesh(phedges, eedges, H, cmap="jet")
    cbar = colorbar()
    cbar.set_label("Photon count")

xlabel("Pulse Phase")
ylabel("Energy (keV)")

if args.outfile is not None:
    savefig(args.outfile)
else:
    show()
