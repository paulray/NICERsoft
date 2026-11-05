#!/usr/bin/env python
"""Fourier (harmonic) analysis of pulse phases."""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from pint.templates import lctemplate, lcprimitives, lcfitters
from pint.eventstats import z2m, sf_z2m, hm, sf_hm, sig2sigma
import sys
import scipy.stats


def compute_fourier(phases, nh=10, pow_phase=False, with_phase_errors=False):
    """Compute the Fourier components of a set of pulse phases.

    The DC term is not computed or returned.

    Parameters
    ----------
    phases : array_like
        Pulse phases, in the range [0, 1).
    nh : int, optional
        Number of harmonics (1 = fundamental only).
    pow_phase : bool, optional
        If True, return powers and phases instead of cosine and sine
        components.
    with_phase_errors : bool, optional
        If True (and ``pow_phase`` is True), also return the phase
        uncertainties.

    Returns
    -------
    n : int
        Number of phases.
    c, s : numpy.ndarray
        Cosine and sine components of each harmonic, normalized so that a
        pure sinusoid with 100% modulation has amplitude 1. Returned if
        ``pow_phase`` is False.
    fourier_pow : numpy.ndarray
        Leahy-normalized power of each harmonic. Returned if ``pow_phase``
        is True.
    fourier_phases : numpy.ndarray
        Phase of each harmonic, in radians. Returned if ``pow_phase`` is True.
    phase_errs : numpy.ndarray
        1-sigma uncertainty of each phase, in radians. Returned if
        ``with_phase_errors`` is also True. The Gaussian approximation used
        breaks down below a signal-to-noise ratio of about 3.
    """
    phis = 2.0 * np.pi * phases  # Convert phases to radians
    n = len(phis)
    c = np.asarray([(np.cos(k * phis)).sum() for k in range(1, nh + 1)]) / n
    s = np.asarray([(np.sin(k * phis)).sum() for k in range(1, nh + 1)]) / n
    c *= 2.0
    s *= 2.0

    if not pow_phase:
        return n, c, s
    # These should be Leahy normalized powers
    fourier_pow = (n / 2) * (c**2 + s**2)
    # These phases are in radians
    fourier_phases = np.arctan2(s, c)

    if not with_phase_errors:
        return n, fourier_pow, fourier_phases

    # Here we compute the Fourier amplitudes (not powers, and not Leahy normalized)
    # This allows us to compute the uncertainties, but note:
    # * When SNR<3, the phase is poorly constrained and Gaussian approx breaks down
    # * Also at low SNR, the amplitude is biased high
    # Here the amplitude is normalized such that for a pure sinusoid, A_1 = 1 means 100% modulation,
    # i.e. represents the peak-to-peak variation
    fourier_amp = np.sqrt(c**2 + s**2)
    sigma_component = np.sqrt(2.0 / n)
    amplitude_errs = np.full(nh, sigma_component)
    # Note that SNR = fourier_amp / amplitude_errs
    phase_errs = sigma_component / fourier_amp  # in radians
    return n, fourier_pow, fourier_phases, phase_errs


def evaluate_fourier(n, c, s, nbins, k=None):
    """Evaluate a Fourier model of a pulse profile at the bin centers.

    Parameters
    ----------
    n : int
        Total number of counts in the profile.
    c, s : array_like
        Cosine and sine components, as returned by :func:`compute_fourier`.
    nbins : int
        Number of profile bins.
    k : int, optional
        If given, include only harmonic ``k + 1`` (plus the constant level).

    Returns
    -------
    numpy.ndarray
        Model counts in each bin.
    """
    # This should be updated to do a little integral over each bin.
    # Currently evaluates the model at the center of each bin
    model = np.zeros(nbins) + n / nbins
    theta = 2.0 * np.pi * np.arange(nbins, dtype=float) / nbins
    theta += theta[1] / 2.0
    if k is not None:
        model += (n / nbins) * (
            c[k] * np.cos((k + 1) * theta) + s[k] * np.sin((k + 1) * theta)
        )
    else:
        for k in range(len(c)):
            model += (n / nbins) * (
                c[k] * np.cos((k + 1) * theta) + s[k] * np.sin((k + 1) * theta)
            )

    return model


def evaluate_chi2(hist, model):
    """Compute the chi-squared of a binned profile against a model.

    The model is used for the variance in each bin.

    Parameters
    ----------
    hist : array_like
        Counts in each bin.
    model : array_like
        Model counts in each bin.

    Returns
    -------
    float
        The chi-squared value.
    """
    # Question here is whether error should be sqrt(data) or sqrt(model)
    return ((hist - model) ** 2 / model).sum()
