#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse

# from astropy import log
from loguru import logger as log
from nicer.values import *
from astropy.io import fits
from astropy.table import Table, vstack
import astropy.units as u
from pint.eventstats import hm, z2m, h2sig


def channel_harmonic_sums(phases, pi, m):
    """Cumulative (over PI channel) cos and sin sums of harmonics 1..m.

    The sums for the photons in channels lo <= PI < hi are then
    cum[..., hi] - cum[..., lo], so each energy band costs O(m) instead of
    a pass over every photon.

    Returns (cumN, cumC, cumS) with shapes (nchan+1,), (m, nchan+1), (m, nchan+1).
    """
    nchan = int(pi.max()) + 1
    cumN = np.zeros(nchan + 1)
    cumC = np.zeros((m, nchan + 1))
    cumS = np.zeros((m, nchan + 1))
    cumN[1:] = np.cumsum(np.bincount(pi, minlength=nchan))
    for i in range(m):
        arg = phases * (2 * np.pi * (i + 1))
        cumC[i, 1:] = np.cumsum(np.bincount(pi, weights=np.cos(arg), minlength=nchan))
        cumS[i, 1:] = np.cumsum(np.bincount(pi, weights=np.sin(arg), minlength=nchan))
    return cumN, cumC, cumS


def band_ts(cum, lo, hi, ztest=False):
    """H-test (or Z^2_m) for the photons in channels lo <= PI < hi, for arrays of bands."""
    cumN, cumC, cumS = cum
    lo, hi = np.broadcast_arrays(lo, hi)
    nph = cumN[hi] - cumN[lo]
    s = (cumC[:, hi] - cumC[:, lo]) ** 2 + (cumS[:, hi] - cumS[:, lo]) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        if ztest:
            ts = (2.0 / nph) * s.sum(axis=0)
        else:
            m = s.shape[0]
            ts = (
                (2.0 / nph) * np.cumsum(s, axis=0) - 4 * np.arange(0, m)[:, None]
            ).max(axis=0)
    return np.where(nph > 0, ts, 0.0)


def get_parser():
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
    return parser


def main(argv=None):
    args = get_parser().parse_args(argv)

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
        t.keep_columns(["PI", "PULSE_PHASE"])
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
    phasesinitial = etable["PULSE_PHASE"].astype(float)
    pi = np.asarray(etable["PI"]).astype(np.int64)
    # Negative PI values never pass an energy cut (emin >= 0); bincount needs PI >= 0
    goodpi = pi >= 0
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

    # Harmonic sums per PI channel, accumulated over channel
    cum = channel_harmonic_sums(phasesinitial[goodpi], pi[goodpi], m)
    # Channel energies, computed exactly as PI * PI_TO_KEV is for each event,
    # so that searchsorted reproduces the cut  emin < PI * PI_TO_KEV < emax
    chan_keV = np.arange(len(cum[0]) - 1) * PI_TO_KEV

    eminlist = []
    emaxlist = []
    hgrid = []

    emins = np.arange(args.minlowE, args.maxlowE + args.steplowE, args.steplowE)

    for emin in emins:

        if args.minhighE is not None:
            minhighE = max(emin + args.steplowE, args.minhighE)
        else:
            minhighE = emin + args.steplowE
        emaxs = np.arange(minhighE, args.maxhighE + args.stephighE, args.stephighE)

        lo = np.searchsorted(chan_keV, emin, side="right")  # first channel > emin
        hi = np.searchsorted(chan_keV, emaxs, side="left")  # first channel >= emax
        hi = np.maximum(hi, lo)
        eminlist.extend([emin] * len(emaxs))
        emaxlist.extend(emaxs)
        hgrid.extend(band_ts(cum, lo, hi, ztest=args.ztest))

    # Same choice as a sequential scan updating on h >= hbest: the last band
    # with the maximum value, if that is at least the full-band value
    hgrid = np.asarray(hgrid)
    if len(hgrid) > 0 and hgrid.max() >= hbest:
        ibest = len(hgrid) - 1 - np.argmax(hgrid[::-1])
        hbest = hgrid[ibest]
        eminbest = eminlist[ibest]
        emaxbest = emaxlist[ibest]

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

    # Everything in log.info() is written to stderr. Write the emin,emax to stdout so that users can both see the progress and time remaining while the code is running, and capture the emin,emax values programmatically from the output:
    print("Best range:  emin {0} keV -- emax {1} keV".format(np.round(eminbest, 3), np.round(emaxbest, 3)))

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


if __name__ == "__main__":
    main()
