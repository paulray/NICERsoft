#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as pyfits
import astropy.stats
import astropy.units as u
import argparse
from pint.eventstats import z2m, hm, sf_z2m, sf_hm, sig2sigma, h2sig
from nicer.values import *
from nicer.fourier import *
import yaml
import os.path

outdict = dict()


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="Do RPP catalog pulse profile analysis on a merged event file.",
)
parser.add_argument(
    "evtfile", help="Input event file to process. Must have PULSE_PHASE column"
)
parser.add_argument("--nbins", help="Number of profile bins", default=32, type=int)
parser.add_argument(
    "--numharm",
    help="Number of Fourier harmonics to use (-1=auto)",
    default=4,
    type=int,
)
parser.add_argument(
    "--optemin", type=float, default=0.4, help="Minimum energy to include."
)
parser.add_argument(
    "--optemax", type=float, default=2.0, help="Maximum energy to include."
)
parser.add_argument(
    "--outfile",
    help="Output file for plot (type determined by extension)",
    default=None,
)
parser.add_argument(
    "--fermi",
    help="Fermi LAT event file (FT1) name, with PULSE_PHASE and MODEL_WEIGHT columns",
    default=None,
)
args = parser.parse_args()


hdulist = pyfits.open(args.evtfile)
dat = hdulist[1].data
hdr = hdulist[1].header
basename = os.path.splitext(args.evtfile)[0]

exp = float(hdr["EXPOSURE"])
objname = hdr["OBJECT"]

print(f"Object {objname}")
print(f"Exposure = {exp/1000.0} ks")

# Get phases and sort them in order
phases = dat["PULSE_PHASE"]
sortidx = np.argsort(phases)
ph = phases[sortidx]

# Get photon energies, with same sort
energies = dat["PI"] * PI_TO_KEV
en = energies[sortidx]

# Make 3 sets of events: optimal band, soft band, hard band
opt_idx = np.where(np.logical_and((en >= args.optemin), (en <= args.optemax)))
ph_opt = ph[opt_idx]
soft_idx = np.where(np.logical_and((en >= SOFT_EMIN), (en <= SOFT_EMAX)))
ph_soft = ph[soft_idx]
hard_idx = np.where(np.logical_and((en >= HARD_EMIN), (en <= HARD_EMAX)))
ph_hard = ph[hard_idx]


def band_analysis(ph_band, bandemin, bandemax, ax=None):
    "Perform analysis for a specific energy band and return dict of results"
    resdict = dict()
    print(
        f"Band  {bandemin} - {bandemax}: {len(ph_band)} photons, {len(ph_band)/exp:.3f} c/s"
    )
    z = z2m(ph_band)
    h = hm(ph_band)
    resdict['Htest'] = float(h)
    resdict['Ztest'] = float(z[-1])
    print(f"    Z^2_2 test = {z[-1]:.2f}", end=" ")
    try:
        print(f"({sig2sigma(sf_z2m(z[-1])):.2f} sig)")

    except:
        print("")
    print(f"    H test = {h:.2f} ({h2sig(h):.2f} sig)")
    n, c, s = compute_fourier(ph_band, nh=args.numharm)
    model = evaluate_fourier(n, c, s, args.nbins)
    model_bins = np.arange(args.nbins) / args.nbins
    pcounts = (model - model.min()).sum()
    pcounts_err = np.sqrt(model.sum() + model.min() * len(model))
    resdict['pulsed_rate'] = float(pcounts/exp)
    resdict['pulsed_rate_err'] = float(pcounts_err/exp)
    
    print(
        "    Pulsed counts = {0:.3f}, pulsed count rate = {1:.3f}+/-{2:.4f} c/s".format(
            pcounts, pcounts / exp, pcounts_err / exp
        )
    )

    prof, edges = np.histogram(
        ph_band, bins=np.linspace(0.0, 1.0, args.nbins, endpoint=True)
    )
    prof = np.array(prof, dtype=float)
    rates = prof / (exp / args.nbins)

    # Compute Fvar, which is the fractional RMS variability amplitude (excess
    # above Poisson)
    # Equation 10 of Vaughan et al (2003, MNRAS, 345, 1271)
    if prof.var() - prof.mean() >= 0.0:
        fracrms = np.sqrt(prof.var() - prof.mean()) / prof.mean()
    else:
        fracrms = -1
    print("    Fractional RMS is {0:.4f}".format(fracrms))
    resdict['Fvar'] = float(fracrms)

    # Compute the Bayesian Block histogram
    ##### WARNING, this probably doesn't handle wrapping through 1.0 correctly. Need to fix that!!!
    bb_hist, bb_edges = astropy.stats.histogram(
        ph_band, bins="blocks", range=[0.0, 1.0]
    )
    bb_widths = bb_edges[1:] - bb_edges[:-1]
    # print(f"{len(bb_hist)} {len(bb_edges)} {len(bb_widths)}")
    bb_rates = bb_hist / (exp * bb_widths)
    minidx = np.argmin(bb_rates)
    print(f"    BB Edges : {bb_edges}")
    print(f"    BB Rates : {bb_rates}")
    resdict['BB_OFFPULSE'] = [float(bb_edges[minidx]), float(bb_edges[minidx+1])]

    if ax:
        ax.step(
            np.concatenate((edges[:-1], 1.0 + edges)),
            np.concatenate((rates, rates, np.array(rates[-1:]))),
            where="post",
        )
        ax.step(
            bb_edges,
            np.concatenate((bb_rates, np.array(bb_rates[-1:]))),
            where="post",
            color="r",
        )
        ax.plot(
            model_bins + model_bins[1] / 2.0,
            model / (exp / args.nbins),
            color="g",
            lw=1,
        )
        ax.set_ylabel("Rate (c/s)")
        ax.set_xlim((0.0, 2.0))
    return resdict


fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(9, 16), sharex=True)

# Analyze all 3 bands
optres = band_analysis(ph_opt, args.optemin, args.optemax, ax=axs[0])
softres = band_analysis(ph_soft, SOFT_EMIN, SOFT_EMAX, ax=axs[1])
hardres = band_analysis(ph_hard, HARD_EMIN, HARD_EMAX, ax=axs[2])

with open(f"{basename}_profinfo.yml", 'w') as outfile:
    outdict["name"] = objname
    outdict["exposure"] = exp
    outdict["optres"] = optres
    outdict["softres"] = softres
    outdict["hardres"] = hardres
    yaml.dump(outdict, outfile, default_flow_style=False)


axs[0].set_title(f"{objname} Exp={exp/1000.0:.3f} ks")
axs[2].set_xlabel("Phase")

if args.outfile is not None:
    plt.savefig(args.outfile)

if args.fermi:
    from pint.fits_utils import read_fits_event_mjds

    latfig, latax = plt.subplots(1, 1)
    f = pyfits.open(args.fermi)
    latph = f["EVENTS"].data.field("PULSE_PHASE")
    laten = f["EVENTS"].data.field("ENERGY")
    latwt = f["EVENTS"].data.field("MODEL_WEIGHT")
    latmjds = read_fits_event_mjds(f["EVENTS"]) * u.d
    nfermi, fedges = np.histogram(
        latph, bins=np.linspace(0.0, 1.0, 64, endpoint=True), weights=latwt
    )
    x = np.concatenate((fedges[:-1], fedges[:-1] + 1.0))
    latax.step(x, np.concatenate((nfermi, nfermi)), where="post")
    latax.set_xlim((0.0, 2.0))
    latax.set_title("LAT Profile")


plt.show()
