#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import tqdm

# from astropy import log
from loguru import logger as log
from nicer.values import *
from astropy.io import fits
from astropy.table import Table, vstack
import astropy.units as u
from pint.eventstats import hm, z2m, h2sig


def cached_hm(mask):
    nph = mask.sum()
    if nph == 0:
        return 0
    s = (cached_hm._cache[..., mask].sum(axis=2) ** 2).sum(axis=1)
    return (
        (2.0 / nph) * np.cumsum(s) - 4 * np.arange(0, cached_hm._cache.shape[0])
    ).max()


def cached_zm(mask):
    nph = mask.sum()
    if nph == 0:
        return 0
    s = (cached_zm._cache[..., mask].sum(axis=2) ** 2).sum()
    return (2.0 / nph) * s


parser = argparse.ArgumentParser(
    description="Read event file with phase and optimize cuts"
)
parser.add_argument("evfiles", help="Name of event files to process", nargs="+")
parser.add_argument(
    "--minlowE", help="minimum E for the low cut.", default=0.25, type=float
)
parser.add_argument(
    "--maxlowE", help="maximum E for the low cut.", default=4.00, type=float
)
parser.add_argument(
    "--minhighE", help="minimum E for the high cut.", default=None, type=float
)
parser.add_argument(
    "--maxhighE", help="maximum E for the high cut.", default=12.0, type=float
)
parser.add_argument(
    "--steplowE", help="step in energy for the low cut.", default=0.05, type=float
)
parser.add_argument(
    "--stephighE", help="step in energy for the high cut.", default=0.1, type=float
)
parser.add_argument(
    "--noplot", help="Suppress plotting.", action="store_true", default=False
)
parser.add_argument(
    "--grid", help="Shows grid of H-test", action="store_true", default=False
)
parser.add_argument(
    "-n", "--nbins", help="Number of bins of outout profiles", default=16, type=int
)
parser.add_argument(
    "-m",
    "--maxharm",
    help="Search up to harmonic m; [default=20, H-test]",
    default=20,
    type=int,
)
parser.add_argument(
    "-z",
    "--ztest",
    help="Use Z-test rather than H-test; this searches with exactly the number of harmonics specified.  m=1 is Rayleigh, m=2 is Z^2_2, etc.",
    action="store_true",
    default=False,
)
parser.add_argument(
    "--min_good_time",
    help="Ignore files with GTI sum below specified value [default None; unit sec].",
    default=None,
    type=float,
)
args = parser.parse_args()

object = None
tlist = []

for fn in args.evfiles:
    log.info("Reading file {0}".format(fn))
    if object is None:
        try:
            object = fits.open(fn)[0].header["OBJECT"].replace("_", " ")
        except:
            object = "unknown object"
        log.info("Opening events for {}".format(object))
    t = Table.read(fn, hdu=1)
    if len(t) == 0:
        continue
    if args.min_good_time is not None:
        tgti = Table.read(fn, hdu="gti")
        good_time = np.sum(tgti["STOP"] - tgti["START"])
        if good_time < args.min_good_time:
            log.info("Ignoring file {0} with good time {1:0.2f}.".format(fn, good_time))
            continue
    tlist.append(t)

log.info("Concatenating files")
if len(tlist) == 1:
    etable = tlist[0]
else:
    etable = vstack(tlist, metadata_conflicts="silent")
del tlist

m = args.maxharm
ts_func = cached_zm if args.ztest else cached_hm
phasesinitial = etable["PULSE_PHASE"].astype(float)
h_ini = z2m(phasesinitial, m=m)[-1] if args.ztest else hm(phasesinitial, m=m)
hbest = h_ini
eminbest = 0.0
emaxbest = 100.0
ts_name = "Ztest" if args.ztest else "Htest"
if args.ztest:
    log.info(
        "Initial {0} = {1} in range {2:.2f}-{3:.2f} keV ".format(
            ts_name,
            np.round(hbest, 3),
            min(etable["PI"] * PI_TO_KEV),
            max(etable["PI"] * PI_TO_KEV),
        )
    )
else:
    log.info(
        "Initial {0} = {1} ({2} sigma) in range {3:.2f}-{4:.2f} keV".format(
            ts_name,
            np.round(hbest, 3),
            np.round(h2sig(hbest), 3),
            min(etable["PI"] * PI_TO_KEV),
            max(etable["PI"] * PI_TO_KEV),
        )
    )

# assemble cache
cache = np.empty([m, 2, len(phasesinitial)], dtype=float)
for i in range(m):
    cache[i, 0] = np.cos(phasesinitial * (2 * np.pi * (i + 1)))
    cache[i, 1] = np.sin(phasesinitial * (2 * np.pi * (i + 1)))
cached_hm._cache = cached_zm._cache = cache

eminlist = []
emaxlist = []
hgrid = []

emins = np.arange(args.minlowE, args.maxlowE + args.steplowE, args.steplowE)

for emin in tqdm.tqdm(emins):

    if args.minhighE is not None:
        minhighE = max(emin + args.steplowE, args.minhighE)
    else:
        minhighE = emin + args.steplowE
    emaxs = np.arange(minhighE, args.maxhighE + args.stephighE, args.stephighE)

    for emax in emaxs:
        mask = np.logical_and(
            etable["PI"] * PI_TO_KEV > emin, etable["PI"] * PI_TO_KEV < emax
        )
        h = ts_func(mask)
        eminlist.append(emin)
        emaxlist.append(emax)
        hgrid.append(h)
        if h >= hbest:
            hbest = h
            eminbest = emin
            emaxbest = emax


if args.ztest:
    log.info("Final {0} = {1}".format(ts_name, np.round(hbest, 3)))
else:
    log.info(
        "Final {0} = {1} ({2} sigma)".format(
            ts_name, np.round(hbest, 3), np.round(h2sig(hbest), 3)
        )
    )
log.info(
    "Best range:  emin {0} keV -- emax {1} keV".format(
        np.round(eminbest, 3), np.round(emaxbest, 3)
    )
)

if not args.noplot:
    profbins = np.linspace(0.0, 1.0, args.nbins + 1, endpoint=True)
    bbins = np.concatenate((profbins, profbins[1:] + 1.0, profbins[1:] + 2.0))

    # Initial pulse profile
    prof_ini, edges_ini = np.histogram(phasesinitial, bins=profbins)
    fullprof_ini = np.concatenate(
        (prof_ini, prof_ini, prof_ini, np.array([prof_ini[0]]))
    )

    # Optimal pulse profile
    idx = np.where(
        np.logical_and(
            etable["PI"] * PI_TO_KEV > eminbest, etable["PI"] * PI_TO_KEV < emaxbest
        )
    )[0]
    phases = etable["PULSE_PHASE"][idx]
    prof_opt, edges_ini = np.histogram(phases, bins=profbins)
    fullprof_opt = np.concatenate(
        (prof_opt, prof_opt, prof_opt, np.array([prof_opt[0]]))
    )

    fig, ax = plt.subplots()

    ax.errorbar(
        bbins - (0.5 / args.nbins),
        fullprof_ini,
        yerr=fullprof_ini**0.5,
        marker="",
        linewidth=1.5,
        color="r",
        label="Init. E range ({0} = {1:.2f})".format(ts_name, np.round(h_ini, 3)),
        drawstyle="steps-mid",
    )

    ax.errorbar(
        bbins - (0.5 / args.nbins),
        fullprof_opt,
        yerr=fullprof_opt**0.5,
        marker="",
        linewidth=1.5,
        color="b",
        label="{:0.2f}-{:0.2f} keV ({} = {:.2f})".format(
            eminbest, emaxbest, ts_name, np.round(hbest, 3)
        ),
        drawstyle="steps-mid",
    )

    # ax.text(
    #     0.1,
    #     0.1,
    #     "{0} = {1:.2f}".format(ts_name, np.round(hbest, 3)),
    #     transform=ax.transAxes,
    # )

    ax.set_ylabel("Counts")
    ax.set_xlabel("Pulse Phase")
    ax.set_title("Pulse Profile for {}".format(object))
    ax.set_xlim(0.0, 2.0)
    ax.set_ylim(ax.get_ylim()[0], 1.1 * ax.get_ylim()[1])
    ax.legend(loc="best")
    plt.savefig(
        "OptimalProfile_{:.2f}keV_{:.2f}keV.png".format(eminbest, emaxbest),
        bbox_inches="tight",
    )
    plt.clf()

if args.grid:
    plt.scatter(eminlist, emaxlist, c=hgrid, s=5)
    cbar = plt.colorbar()
    cbar.set_label("H-test")
    plt.title("{} vs energy range for {}".format(ts_name, object))
    plt.xlabel("Low Energy Cut (keV)")
    plt.ylabel("High Energy Cut (keV)")
    plt.savefig("GridSearch.png")
