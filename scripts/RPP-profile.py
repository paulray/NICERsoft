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

outdict = {}


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
    "--p0", type=float, default=0.005, help="p0 probability for Bayesian block analysis"
)
parser.add_argument(
    "--outbase",
    help="Output file for plot (type determined by extension)",
    default=None,
)
parser.add_argument(
    "--srcname",
    help="Source name (used in output .yml file).",
    default=None,
)
parser.add_argument(
    "--fermi",
    help="Fermi LAT event file (FT1) name, with PULSE_PHASE and MODEL_WEIGHT columns",
    default=None,
)
parser.add_argument(
    "--lat3pc",
    help="Fermi LAT 3PC data file for this pulsar, e.g. J0218+4232_3PC_data.fits",
    default=None,
)
args = parser.parse_args()


hdulist = pyfits.open(args.evtfile)
dat = hdulist[1].data
hdr = hdulist[1].header
basename = os.path.splitext(args.evtfile)[0]

exp = float(hdr["EXPOSURE"])
objname = hdr["OBJECT"]

outdict[objname] = {}

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


def compute_blocks(ph,p0=0.005):
    "Compute Bayesian Block edges for a set of photon phases."

    # Replicate photons before and after so they go from [-1.0, 2.0)
    bbphases = np.concatenate([ph - 1.0, ph, ph + 1.0])

    if len(ph) < 50000:
        # Compute unbinned Bayesian Block edges
        print(f"Unbinned BB ({p0})")
        edges = astropy.stats.bayesian_blocks(bbphases, fitness="events", p0=p0)
    else:
        # Compute binned to speed up process
        h, hedges = np.histogram(
            bbphases,
            range=(-1.0, 2.0),
            bins=np.linspace(-1.0, 2.0, 301, endpoint=True),
        )
        print(f"Binned BB ({p0})")
        edges = astropy.stats.bayesian_blocks(
            t=hedges[:-1],
            x=h,
            sigma=np.sqrt(h + 1),
            fitness="events",
            p0=p0,
        )
    # Select out only the bin edges between 0 and 1
    idx = np.logical_and(edges >= 0.0, edges < 1.0)
    mybins = np.concatenate((np.zeros(1), edges[idx], np.ones(1)))
    bb_hist, bb_edges = np.histogram(ph, range=(0.0, 1.0), bins=mybins)
    bb_widths = bb_edges[1:] - bb_edges[:-1]
    # Note that the first and last bin should be merged, i.e. it is (mybins[0],mybins[1]) U (mybins[-2],mybins[-1])
    c1 = bb_hist[0]
    c2 = bb_hist[-1]
    w1 = bb_widths[0]
    w2 = bb_widths[-1]
    cnts = c1 + c2
    width = w1 + w2
    rate = cnts / width
    c1new = rate * w1
    c2new = rate * w2
    bb_hist[0] = c1new
    bb_hist[-1] = c2new
    return bb_hist, bb_edges


def band_analysis(ph_band, bandemin, bandemax, ax=None, plotoffpulse=False):
    "Perform analysis for a specific energy band and return dict of results"
    print(
        f"Band  {bandemin} - {bandemax}: {len(ph_band)} photons, {len(ph_band)/exp:.3f} c/s"
    )
    z = z2m(ph_band)
    h = hm(ph_band)
    resdict = {
        "Htest": float(h),
        "Ztest": float(z[-1]),
        "Emin": float(bandemin),
        "Emax": float(bandemax),
    }
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
    resdict["pulsed_rate"] = float(pcounts / exp)
    resdict["pulsed_rate_err"] = float(pcounts_err / exp)

    print(
        "    Pulsed counts = {0:.3f}, pulsed count rate = {1:.3f}+/-{2:.4f} c/s".format(
            pcounts, pcounts / exp, pcounts_err / exp
        )
    )

    prof, edges = np.histogram(
        ph_band,
        range=(0.0, 1.0),
        bins=np.linspace(0.0, 1.0, args.nbins + 1, endpoint=True),
    )
    # print(f"Edges {edges}")
    prof = np.array(prof, dtype=float)
    rates = prof / (exp / args.nbins)
    # print(f"sum of prof {prof.sum()}, mean rate {rates.mean()} {len(rates)/args.nbins}")

    # Compute Fvar, which is the fractional RMS variability amplitude (excess
    # above Poisson)
    # Equation 10 of Vaughan et al (2003, MNRAS, 345, 1271)
    if prof.var() - prof.mean() >= 0.0:
        fracrms = np.sqrt(prof.var() - prof.mean()) / prof.mean()
    else:
        fracrms = -1
    print("    Fractional RMS (Fvar) is {0:.4f}".format(fracrms))
    resdict["Fvar"] = float(fracrms)

    # Compute the Bayesian Block histogram, handling wrapping at 1.0
    bb_hist, bb_edges = compute_blocks(ph_band,p0=args.p0)

    # Convert histogram to rates
    bb_widths = bb_edges[1:] - bb_edges[:-1]
    bb_rates = bb_hist / (exp * bb_widths)
    # Fix rate calc for first and last bins, which are merged.
    wraprate = (bb_hist[0] + bb_hist[-1])/(exp*(bb_widths[0] + bb_widths[-1]))
    bb_rates[0] = wraprate
    bb_rates[-1] = wraprate

    minidx = np.argmin(bb_rates)
    # print(f" Sum of bb_hist {bb_hist.sum()}, sum of bb_widths {bb_widths.sum()}")
    print(f"    BB Edges : {bb_edges}")
    print(f"    BB Rates : {bb_rates}")
    # If minidx is 0 or the last block, then we need to merge them!
    if (minidx == 0):
        resdict["BB_OFFPULSE"] = [float(bb_edges[-2]), float(bb_edges[1])]
    elif minidx == len(bb_rates)-1:
        resdict["BB_OFFPULSE"] = [float(bb_edges[-2]), float(bb_edges[1])]
    else:
        resdict["BB_OFFPULSE"] = [float(bb_edges[minidx]), float(bb_edges[minidx + 1])]
    print(f"BB_OFFPULSE: {resdict['BB_OFFPULSE']}")
        
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
        if plotoffpulse:
            ax.axvline(float(resdict["BB_OFFPULSE"][0]), linestyle="--", color="k")
            ax.axvline(float(resdict["BB_OFFPULSE"][1]), linestyle="--", color="k")
        ax.set_ylabel(f"{bandemin}-{bandemax} keV Rate (c/s)")
        ax.set_xlim((0.0, 2.0))
        ax.grid(True)
        ax.set_title(f"H-test: {float(h):.2f} ({h2sig(float(h)):.2f} sigma)")
    return resdict


# fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(9, 16), sharex=True)
fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(6, 8), sharex=True)

# Analyze all 3 bands
optres = band_analysis(ph_opt, args.optemin, args.optemax, ax=axs[0], plotoffpulse=True)
softres = band_analysis(ph_soft, SOFT_EMIN, SOFT_EMAX, ax=axs[1])
hardres = band_analysis(ph_hard, HARD_EMIN, HARD_EMAX, ax=axs[2])

if args.srcname is None:
    outfile = open(f"{basename}_profinfo.yml", "w")
else:
    outfile = open(args.srcname + "_profinfo.yml", "w")
outdict[objname]["name"] = objname
outdict[objname]["exposure"] = exp
outdict[objname]["optres"] = optres
outdict[objname]["softres"] = softres
outdict[objname]["hardres"] = hardres
yaml.dump(outdict, outfile, default_flow_style=False)


if args.srcname is not None:
    fig.suptitle(f"{args.srcname} Exp={exp/1000.0:.3f} ks", y=0.95)
else:
    fig.suptitle(f"{objname} Exp={exp/1000.0:.3f} ks", y=0.95)
axs[2].set_xlabel("Phase")

if args.outbase is not None:
    fig.savefig(f"{args.outbase}.png")

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
    if args.outbase is not None:
        latfig.savefig(f"{args.outbase}_LAT.png")

if args.lat3pc:
    print(f"3PC {args.lat3pc}")
    multifig, multiax = plt.subplots(1, 1)
    prof, edges = np.histogram(
        ph_opt,
        range=(0.0, 1.0),
        bins=np.linspace(0.0, 1.0, args.nbins + 1, endpoint=True),
    )
    prof = np.array(prof, dtype=float)
    rates = prof / (exp / args.nbins)

    hasFit = False
    hasRadio = False
    with pyfits.open(args.lat3pc) as f:
        dataphmin = f["GAMMA_LC"].data["Ph_Min"]
        dataphmax = f["GAMMA_LC"].data["Ph_Max"]
        dataCounts = f["GAMMA_LC"].data["GT100_WtCnt"]
        dataCountsUnc = f["GAMMA_LC"].data["Unc_GT100_WtCnt"]

        if "BEST_FIT_LC" in f:
            fitphmin = f["BEST_FIT_LC"].data["Ph_Min"]
            fitphmax = f["BEST_FIT_LC"].data["Ph_Max"]
            fit = f["BEST_FIT_LC"].data["Intensity"]
            hasFit = True

        if "RADIO_PROFILE" in f:
            radiophmin = f["RADIO_PROFILE"].data["Ph_Min"]
            radiophmax = f["RADIO_PROFILE"].data["Ph_Max"]
            radio = f["RADIO_PROFILE"].data["Norm_Intensity"]
            hasRadio = True

            # PLOT FERMI LAT
            middle = (dataphmin + dataphmax) / 2.0
            dmin = dataCounts.min()
            dmax = dataCounts.max()
            scdata = (dataCounts-dmin)/(dmax-dmin)
            multiax.step(middle, scdata + 2.0, "k-", where="post", label="Fermi/LAT")

            # PLOT RADIO
            if hasRadio:
                middle = (radiophmin + radiophmax) / 2.0
                rmin = radio.min()
                rmax = radio.max()
                multiax.plot(middle, (radio-rmin)/(rmax-rmin), "r-", label="Radio",alpha=0.7)

            # if hasFit:
            #     middle = (fitphmin + fitphmax) / 2.0
            #     multiax.plot(middle, fit, "b-", label="LATile Fit")

            # Now plot NICER
            rmin = rates.min()
            rmax = rates.max()
            scrates = (rates-rmin)/(rmax-rmin)
            multiax.step(
                np.concatenate((edges[:-1], 1.0 + edges)),
                np.concatenate((scrates, scrates, np.array(scrates[-1:])))+1.0,
                where="post",
                label="NICER"
            )

            multiax.set_xlabel("Phase")
            multiax.set_ylabel("Relative Intensity")
            multiax.set_title(f"Multiband Profile for {args.srcname}")
            multiax.grid(True)
            multiax.legend()
            multiax.set_ylim(0.0,3.05)
            multiax.set_xlim((0.0, 2.0))
            if args.outbase is not None:
                plt.savefig(f"{args.outbase}_multi.png")

if args.outbase is None:
    plt.show()
