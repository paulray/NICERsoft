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

# --- TOA accuracy (Cramer-Rao lower bound) settings ---
TOA_EXP = 5000.0  # fiducial integration time (s) for the quoted TOA accuracy
TOA_MAXHARM = 50  # highest harmonic considered
TOA_ZJ_THRESH = 15.0  # harmonic j is significant if its Z^2 power z_j > this (p ~ 5e-4)
TOA_MIN_COUNTS = 100  # bands with fewer photons are skipped in the band search
TOA_NGRID = 4096  # phase grid for template integrals and the Monte Carlo check
# Energy grid for the TOA band search (same as the nioptcuts.py defaults)
TOA_EMIN_RANGE = (0.25, 4.0, 0.05)  # (start, stop, step) for the low edge, keV
TOA_EMAX_MAX, TOA_EMAX_STEP = 12.0, 0.1  # high edge runs from Emin+0.05 to 12 keV


def get_parser():
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
    parser.add_argument(
        "--toa-mc",
        action="store_true",
        help="Check the TOA Cramer-Rao bound with a Monte Carlo of simulated 5 ks observations",
    )
    parser.add_argument(
        "--toa-mc-trials",
        type=int,
        default=500,
        help="Number of simulated observations for --toa-mc",
    )
    return parser


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


def band_analysis(ph_band, bandemin, bandemax, ax=None, plotoffpulse=False, *, exp, args):
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


############################################################################
# TOA accuracy
#
# For a profile f(phi) (total rate vs phase, including background and
# unpulsed emission, normalized to mean 1), the Cramer-Rao lower bound on the
# phase of a TOA from N photons is sigma_phi^2 = 1 / (N * Int f'^2/f dphi).
# In the weak-signal limit (f ~ 1) this is the Hartman formula,
# 0.1125 P / (X sqrt(N)), and it can be written with the per-harmonic Z^2
# powers z_j = (2/N)|sum_k exp(2 pi i j phi_k)|^2, which have E[z_j] = 2 + signal:
#     weak-signal phase information over the exposure = (2 pi)^2 sum_j j^2 (z_j - 2)
# The "- 2" removes the noise bias. The weak-signal value is optimistic for
# high-contrast profiles, so the quoted accuracy (sigma_phase_5ks) is the exact
# CRLB, computed from a bias-corrected Fourier template of the significant
# harmonics, and the energy band is chosen to minimize it. All quantities here
# are in phase units; multiply by the period to get time.
############################################################################


def band_channels(emin, emax):
    "PI channel range [lo, hi] (inclusive) matching emin <= PI*PI_TO_KEV <= emax."
    lo = np.ceil(np.asarray(emin) / PI_TO_KEV - 1e-6).astype(int)
    hi = np.floor(np.asarray(emax) / PI_TO_KEV + 1e-6).astype(int)
    return lo, hi


def channel_harmonic_sums(ph, pi, maxharm=TOA_MAXHARM):
    """Cumulative (over PI channel) sums of cos and sin(2 pi j phi), j=1..maxharm.

    With these, the harmonic sums for any energy band are a difference of two
    entries, so a full grid of bands costs almost nothing.

    Returns (cumN, cumC, cumS) with shapes (nchan+1,), (maxharm, nchan+1), (maxharm, nchan+1).
    """
    nchan = int(pi.max()) + 1
    cumN = np.zeros(nchan + 1)
    cumC = np.zeros((maxharm, nchan + 1))
    cumS = np.zeros((maxharm, nchan + 1))
    cumN[1:] = np.cumsum(np.bincount(pi, minlength=nchan))
    z1 = np.exp(2j * np.pi * ph)
    zj = np.ones_like(z1)
    for j in range(maxharm):
        zj = zj * z1  # exp(2 pi i (j+1) phi)
        cumC[j, 1:] = np.cumsum(np.bincount(pi, weights=zj.real, minlength=nchan))
        cumS[j, 1:] = np.cumsum(np.bincount(pi, weights=zj.imag, minlength=nchan))
    return cumN, cumC, cumS


def band_zpowers(cum, lo, hi):
    """Photon count, cos/sin sums and Z^2 powers z_j for band(s) lo <= PI <= hi.

    lo and hi can be scalars or arrays of bands; z has shape (maxharm,) or (maxharm, nbands).
    """
    cumN, cumC, cumS = cum
    lo = np.clip(lo, 0, len(cumN) - 1)
    hi = np.clip(hi, -1, len(cumN) - 2)
    n = cumN[hi + 1] - cumN[lo]
    c = cumC[:, hi + 1] - cumC[:, lo]
    s = cumS[:, hi + 1] - cumS[:, lo]
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(n > 0, 2.0 * (c**2 + s**2) / n, 0.0)
    return n, c, s, z


def select_nharm(z):
    "Number of harmonics to use: the highest harmonic with z_j > TOA_ZJ_THRESH (at least 1)."
    sig = np.nonzero(z > TOA_ZJ_THRESH)[0]
    return int(sig[-1] + 1) if len(sig) else 1


def weak_info(z, nharm):
    "Bias-corrected weak-signal phase information over the exposure, (2pi)^2 sum j^2 (z_j - 2)."
    j2 = np.arange(1, nharm + 1) ** 2
    if z.ndim == 2:
        j2 = j2[:, None]
    return (2 * np.pi) ** 2 * (j2 * (z[:nharm] - 2.0)).sum(axis=0)


def toa_template(n, c, s, z, nharm, ngrid=TOA_NGRID):
    """Bias-corrected Fourier template f(phi) (mean 1) and f'(phi) on a phase grid.

    Each harmonic's amplitude is scaled by sqrt(1 - 2/z_j), which removes the
    noise bias in its power (|c_j|^2 -> |c_j|^2 - 1/N), so that the
    weak-signal information of this template equals weak_info().

    n can be a scalar (c, s, z of shape (maxharm,)) or an array of bands
    (c, s, z of shape (maxharm, nbands)); f and fp then have shape (nbands, ngrid).
    """
    j = np.arange(1, nharm + 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        shrink = np.sqrt(np.clip(1.0 - 2.0 / z[:nharm], 0.0, None))
    a = 2.0 * c[:nharm] * shrink / n
    b = 2.0 * s[:nharm] * shrink / n
    jj = j if a.ndim == 1 else j[:, None]
    phi = (np.arange(ngrid) + 0.5) / ngrid
    theta = 2 * np.pi * np.outer(j, phi)
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    f = 1.0 + a.T @ cos_t + b.T @ sin_t
    fp = 2 * np.pi * ((jj * b).T @ cos_t - (jj * a).T @ sin_t)
    return f, fp


def exact_info(n, f, fp):
    """Exact phase information over the exposure, N * Int f'^2/f dphi (NaN where f <= 0).

    Note the weak-signal value N * Int f'^2 dphi is *larger* than this for
    high-contrast profiles (f > 1 at the peaks, where f' is large), so the
    weak-signal (Hartman) approximation is optimistic for sharp, bright pulses.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        info = n * np.mean(fp**2 / f, axis=-1)
    return np.where(f.min(axis=-1) > 0, info, np.nan)


def toa_search_band(cum, nharm, ngrid=1024, chunk=1000):
    "Grid search for the energy band that maximizes the exact (CRLB) phase information."
    emins = np.arange(TOA_EMIN_RANGE[0], TOA_EMIN_RANGE[1] + 1e-9, TOA_EMIN_RANGE[2])
    eminlist, emaxlist = [], []
    for emin in emins:
        emaxs = np.arange(emin + TOA_EMIN_RANGE[2], TOA_EMAX_MAX + 1e-9, TOA_EMAX_STEP)
        eminlist.append(np.full(len(emaxs), emin))
        emaxlist.append(emaxs)
    eminlist = np.concatenate(eminlist)
    emaxlist = np.concatenate(emaxlist)
    lo, hi = band_channels(eminlist, emaxlist)
    n, c, s, z = band_zpowers(cum, lo, hi)
    winfo = weak_info(z, nharm)
    ok = np.nonzero((n >= TOA_MIN_COUNTS) & (winfo > 0))[0]
    if len(ok) == 0:
        # No signal anywhere: return the band with the most weak-signal information
        winfo[n < TOA_MIN_COUNTS] = -np.inf
        best = np.argmax(winfo)
        return int(lo[best]), int(hi[best])
    info = np.full(len(lo), -np.inf)
    for idx in np.array_split(ok, max(1, len(ok) // chunk)):
        f, fp = toa_template(n[idx], c[:, idx], s[:, idx], z[:, idx], nharm, ngrid=ngrid)
        info[idx] = np.nan_to_num(exact_info(n[idx], f, fp), nan=-np.inf)
    if not np.isfinite(info).any():
        # Every template rings below zero: fall back to the weak-signal information
        info = np.where(n >= TOA_MIN_COUNTS, winfo, -np.inf)
    best = np.argmax(info)
    return int(lo[best]), int(hi[best])


def toa_mc(f, nphot_mean, ntrials, sigma_crlb, seed=1):
    """Monte Carlo check of the CRLB at the fiducial exposure.

    Draws ntrials simulated observations (Poisson number of photons with mean
    nphot_mean) from the template f, measures each one's phase by maximum
    likelihood against the same template, and returns the RMS phase error and
    the fraction of trials with errors > 5 sigma_crlb.
    """
    rng = np.random.default_rng(seed)
    nbins = len(f)
    cdf = np.cumsum(f) / f.sum()
    conj_G = np.conj(np.fft.rfft(np.log(f)))
    errs = np.empty(ntrials)
    for i in range(ntrials):
        nph = rng.poisson(nphot_mean)
        tau = rng.random()
        ph = (np.searchsorted(cdf, rng.random(nph)) + rng.random(nph)) / nbins + tau
        h = np.bincount((np.mod(ph, 1.0) * nbins).astype(int) % nbins, minlength=nbins)
        # Log-likelihood vs trial shift: L[k] = sum_b h[b] log f[b - k]
        L = np.fft.irfft(np.fft.rfft(h) * conj_G, n=nbins)
        k = np.argmax(L)
        lm, l0, lp = L[(k - 1) % nbins], L[k], L[(k + 1) % nbins]
        denom = lm - 2 * l0 + lp
        dk = 0.5 * (lm - lp) / denom if denom < 0 else 0.0
        # Bin k holds phases [k, k+1)/nbins and the template is at bin centers
        est = (k + dk) / nbins
        errs[i] = np.mod(est - tau + 0.5, 1.0) - 0.5
    rms = float(np.sqrt(np.mean(errs**2)))
    outliers = float(np.mean(np.abs(errs) > 5 * sigma_crlb))
    return rms, outliers


def toa_band_result(cum, lo, hi, nharm, exp):
    "TOA quantities (phase units, at the fiducial exposure) for one band and number of harmonics."
    n, c, s, z = band_zpowers(cum, lo, hi)
    scale = TOA_EXP / exp
    winfo = weak_info(z, nharm)
    res = {
        "nharm": int(nharm),
        "counts": float(n),
        "rate": float(n / exp),
        "zj": [float(x) for x in z[:nharm]],
        "Z2_signal_5ks": float((z[:nharm] - 2.0).sum() * scale),
    }
    res["sigma_phase_5ks_weak"] = (
        float(1.0 / np.sqrt(winfo * scale)) if winfo > 0 else float("nan")
    )
    f, fp = toa_template(n, c, s, z, nharm)
    res["template_fmin"] = float(f.min())
    einfo = exact_info(n, f, fp)
    if winfo > 0 and np.isfinite(einfo) and einfo > 0:
        res["sigma_phase_5ks"] = float(1.0 / np.sqrt(einfo * scale))
        res["exact_ok"] = True
    else:
        # Template rings below zero (or no signal): fall back to the weak-signal
        # value, which is optimistic for high-contrast profiles. Flagged by exact_ok.
        res["sigma_phase_5ks"] = res["sigma_phase_5ks_weak"]
        res["exact_ok"] = False
    return res, f


def toa_info(ph, pi, exp, optemin, optemax, args):
    """TOA accuracy for the H-test optimal band and for the band that minimizes the TOA error.

    Returns (optdict, toares): optdict holds toa_* keys to add to optres;
    toares is the result for the TOA-optimal band, with its Emin/Emax.
    """
    good = pi >= 0
    cum = channel_harmonic_sums(ph[good], pi[good].astype(int))

    # H-test optimal band (from nioptcuts.py), for comparison
    lo0, hi0 = band_channels(optemin, optemax)
    z0 = band_zpowers(cum, lo0, hi0)[3]
    m0 = select_nharm(z0)
    optband, _ = toa_band_result(cum, lo0, hi0, m0, exp)
    optdict = {
        "toa_nharm": optband["nharm"],
        "toa_Z2_signal_5ks": optband["Z2_signal_5ks"],
        "toa_sigma_phase_5ks": optband["sigma_phase_5ks"],
        "toa_sigma_phase_5ks_weak": optband["sigma_phase_5ks_weak"],
    }

    # Band search. The number of harmonics is set in the current band and
    # held fixed during the search, then re-chosen in the winning band.
    m = m0
    for _ in range(5):
        lo, hi = toa_search_band(cum, m)
        m_new = select_nharm(band_zpowers(cum, lo, hi)[3])
        if m_new == m:
            break
        m = m_new
    toares, f = toa_band_result(cum, lo, hi, m_new, exp)
    toares = {
        "Emin": float(round(lo * PI_TO_KEV, 3)),
        "Emax": float(round(hi * PI_TO_KEV, 3)),
        "fiducial_exposure": TOA_EXP,
        **toares,
    }
    band = np.logical_and(pi >= lo, pi <= hi)
    toares["Htest"] = float(hm(ph[band]))

    if args.toa_mc:
        if toares["exact_ok"]:
            rms, outl = toa_mc(
                f, toares["counts"] * TOA_EXP / exp, args.toa_mc_trials, toares["sigma_phase_5ks"]
            )
            toares["mc_rms_phase_5ks"] = rms
            toares["mc_ratio"] = float(rms / toares["sigma_phase_5ks"])
            toares["mc_outlier_frac"] = outl
            toares["mc_trials"] = int(args.toa_mc_trials)
        else:
            print("    TOA MC skipped: template not positive")

    print(
        f"TOA (H-test band {optemin}-{optemax} keV): nharm={m0} "
        f"Z2_5ks={optband['Z2_signal_5ks']:.1f} sigma_phase(5ks)={optband['sigma_phase_5ks']:.3g}"
    )
    print(
        f"TOA band {toares['Emin']}-{toares['Emax']} keV: nharm={toares['nharm']} "
        f"Z2_5ks={toares['Z2_signal_5ks']:.1f} sigma_phase(5ks)={toares['sigma_phase_5ks']:.3g} "
        f"(weak {toares['sigma_phase_5ks_weak']:.3g}, exact_ok={toares['exact_ok']})"
    )
    if "mc_ratio" in toares:
        print(
            f"    MC: rms={toares['mc_rms_phase_5ks']:.3g} ratio to CRLB={toares['mc_ratio']:.2f} "
            f"outliers(>5 sigma)={toares['mc_outlier_frac']:.3f}"
        )
    return optdict, toares


def main(argv=None):
    args = get_parser().parse_args(argv)


    hdulist = pyfits.open(args.evtfile)
    dat = hdulist[1].data
    hdr = hdulist[1].header
    basename = os.path.splitext(args.evtfile)[0]
    mjdstart = (dat["TIME"].min() + float(hdr["TIMEZERO"]))/86400.0 + float(hdr["MJDREFI"]) + float(hdr["MJDREFF"])
    mjdstop = (dat["TIME"].max() + float(hdr["TIMEZERO"]))/86400.0 + float(hdr["MJDREFI"]) + float(hdr["MJDREFF"])

    print(f"mjdstart = {mjdstart}")

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


    # fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(9, 16), sharex=True)
    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(6, 8), sharex=True)

    # Analyze all 3 bands
    optres = band_analysis(
        ph_opt, args.optemin, args.optemax, ax=axs[0], plotoffpulse=True, exp=exp, args=args
    )
    softres = band_analysis(ph_soft, SOFT_EMIN, SOFT_EMAX, ax=axs[1], exp=exp, args=args)
    hardres = band_analysis(ph_hard, HARD_EMIN, HARD_EMAX, ax=axs[2], exp=exp, args=args)

    # TOA accuracy: in the H-test band (added to optres) and in the band that
    # minimizes the TOA uncertainty (toares). Phase units, for TOA_EXP seconds.
    pi = np.asarray(dat["PI"])[sortidx]
    toa_opt, toares = toa_info(ph, pi, exp, args.optemin, args.optemax, args)
    optres.update(toa_opt)

    if args.srcname is None:
        outfile = open(f"{basename}_profinfo.yml", "w")
    else:
        outfile = open(args.srcname + "_profinfo.yml", "w")
    outdict[objname]["name"] = objname
    outdict[objname]["exposure"] = exp
    outdict[objname]["optres"] = optres
    outdict[objname]["softres"] = softres
    outdict[objname]["hardres"] = hardres
    outdict[objname]["toares"] = toares
    outdict[objname]["mjdstart"] =  float(mjdstart)
    outdict[objname]["mjdstop"] =  float(mjdstop)
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

        # --- NICER profile (always available): compute, but don't plot yet
        prof, edges = np.histogram(
            ph_opt,
            range=(0.0, 1.0),
            bins=np.linspace(0.0, 1.0, args.nbins + 1, endpoint=True),
        )
        prof = np.array(prof, dtype=float)
        rates = prof / (exp / args.nbins)
        rmin = rates.min()
        rmax = rates.max()
        scrates = (rates - rmin) / (rmax - rmin + 1e-12)

        # Flags and data holders
        has_lat = False
        has_radio = False
        lat_x = lat_y = None
        radio_x = radio_y = None

        # --- Try to read LAT and RADIO from the 3PC file
        with pyfits.open(args.lat3pc) as f:
            if "GAMMA_LC" in f:
                has_lat = True
                dataphmin = f["GAMMA_LC"].data["Ph_Min"]
                dataphmax = f["GAMMA_LC"].data["Ph_Max"]
                dataCounts = f["GAMMA_LC"].data["GT100_WtCnt"]
                mid = (dataphmin + dataphmax) / 2.0
                dmin = np.nanmin(dataCounts)
                dmax = np.nanmax(dataCounts)
                scdata = (dataCounts - dmin) / (dmax - dmin + 1e-12)
                lat_x, lat_y = mid, scdata

            if "RADIO_PROFILE" in f:
                has_radio = True
                radiophmin = f["RADIO_PROFILE"].data["Ph_Min"]
                radiophmax = f["RADIO_PROFILE"].data["Ph_Max"]
                radio = f["RADIO_PROFILE"].data["Norm_Intensity"]
                midr = (radiophmin + radiophmax) / 2.0
                rmin_r = np.nanmin(radio)
                rmax_r = np.nanmax(radio)
                sradio = (radio - rmin_r) / (rmax_r - rmin_r + 1e-12)
                radio_x, radio_y = midr, sradio

        # --- Dynamic stacking: assign vertical levels only to available series
        level = 0.0
        if has_radio:
            multiax.plot(radio_x, radio_y, "r-", label="Radio", alpha=0.7)
            level += 1.0

        # NICER goes next
        multiax.step(
            np.concatenate((edges[:-1], 1.0 + edges)),
            np.concatenate((scrates, scrates, np.array(scrates[-1:]))) + (level if level else 0.0),
            where="post",
            label="NICER"
        )
        level += 1.0

        # LAT on top if present
        if has_lat:
            multiax.step(lat_x, lat_y + level, "k-", where="post", label="Fermi/LAT")
            level += 1.0


        # Off-pulse region
        multiax.plot([float(optres["BB_OFFPULSE"][0]),float(optres["BB_OFFPULSE"][0])],[1,2], linestyle="--", color="k",label='Off-pulse')
        if float(optres["BB_OFFPULSE"][1]) > float(optres["BB_OFFPULSE"][0]):
            multiax.plot([float(optres["BB_OFFPULSE"][1]),float(optres["BB_OFFPULSE"][1])],[1,2], linestyle="--", color="k")
        else:
            # It wraps over phase=1; move it so that it doesn't look weird/discontinuous when 2 rotations are shown:
            multiax.plot([float(optres["BB_OFFPULSE"][1])+1,float(optres["BB_OFFPULSE"][1])+1],[1,2], linestyle="--", color="k")

        # Layout
        multiax.set_xlabel("Phase")
        multiax.set_ylabel("Relative Intensity")
        multiax.set_title(f"Multiband Profile for {args.srcname}")
        multiax.grid(True)
        multiax.legend(loc='upper center',ncol=4)
        multiax.set_xlim((0.0, 2.0))
        multiax.set_ylim(0.0, level + 0.35) # enough room on top for legend

        if args.outbase is not None:
            plt.savefig(f"{args.outbase}_multi.png")

    if args.outbase is None:
        plt.show()


if __name__ == "__main__":
    main()
