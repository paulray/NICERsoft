#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as pyfits
import argparse
from pint.eventstats import z2m, hm, sf_z2m, sf_hm, sig2sigma, h2sig
from nicer.values import *

parser = argparse.ArgumentParser(
    description="Plot a binned light curve from a PULSE_PHASE column in a FITS file"
)
parser.add_argument(
    "evtfile", help="Input event file to process. Must have PULSE_PHASE column"
)
parser.add_argument(
    "--nbins", help="Number of profile bins (default 32)", default=32, type=int
)
parser.add_argument(
     "--emin", type=float, default=0.25, help="Minimum energy to include."
)
parser.add_argument(
    "--emax", type=float, default=12.0, help="Maximum energy to include."
)
parser.add_argument(
    "--outfile",
    help="Output file for plot (type determined by extension)",
    default=None,
)
args = parser.parse_args()


hdulist = pyfits.open(args.evtfile)
dat = hdulist[1].data
hdr = hdulist[1].header

ph = dat["PULSE_PHASE"]

# FILTER ON ENERGY
en = dat["PI"] * PI_TO_KEV
idx = np.where(np.logical_and((en >= args.emin), (en <= args.emax)))
ph = ph[idx]
print(f"Energy cuts left {len(ph)} out of {len(en)} events in {args.emin} to {args.emax} keV.")

print(
    "Z^2_2 test = {}".format(z2m(ph)[-1]), end=" "
)
try:
    print("({} sig)".format(sig2sigma(sf_z2m(z2m(ph)[-1]))))

except:
    print("")


print("H test = {} ({} sig)".format(hm(ph), h2sig(hm(ph))))

fig, ax = plt.subplots(figsize=(8, 4.5))


h, edges = np.histogram(ph, bins=np.linspace(0.0, 1.0, args.nbins, endpoint=True))
try:
    if hdr["EXPOSURE"] > 0:
        h = np.array(h, dtype=float) / (float(hdr["EXPOSURE"]) / args.nbins)
        rates = True
except:
    rates = False
ax.step(
    np.concatenate((edges[:-1], 1.0 + edges)),
    np.concatenate((h, h, np.array(h[-1:]))),
    where="post",
)
ax.set_xlabel("Phase")
if rates:
    ax.set_ylabel("Rate (c/s)")
else:
    ax.set_ylabel("Counts")

try:
    tstart = (hdr["TSTART"] + hdr["TIMEZERO"]) / 86400 + hdr["MJDREFF"] + hdr["MJDREFI"]
    ax.set_title("{0} : {1:.3f}".format(hdr["OBS_ID"], tstart))
except:
    ax.set_title("{0}".format(hdr["DATE-OBS"]))
ax.set_xlim((0.0, 2.0))
ax.grid(True)
if args.outfile is not None:
    plt.savefig(args.outfile)

plt.show()
