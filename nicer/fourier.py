#!/usr/bin/env python

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from pint.templates import lctemplate, lcprimitives, lcfitters
from pint.eventstats import z2m, sf_z2m, hm, sf_hm, sig2sigma
import sys
import scipy.stats


def compute_fourier(phases, nh=10, pow_phase=False):
    """Compute Fourier amplitudes from an array of pulse phases
    phases should be [0,1.0)
    nh is the number of harmonics (1 = fundamental only)
    Returns: cos and sin component arrays, unless pow_phase is True
    then returns Fourier power (Leahy normalized) and phase arrays
    DC bin is not computed or returned
    """
    phis = 2.0 * np.pi * phases  # Convert phases to radians
    n = len(phis)
    c = np.asarray([(np.cos(k * phis)).sum() for k in range(1, nh + 1)]) / n
    s = np.asarray([(np.sin(k * phis)).sum() for k in range(1, nh + 1)]) / n
    c *= 2.0
    s *= 2.0

    if not pow_phase:
        return n, c, s
    # CHECK!  There could be errors here!
    # These should be Leahy normalized powers
    fourier_pow = (n / 2) * (c**2 + s**2)
    fourier_phases = np.arctan2(s, c)
    return n, fourier_pow, fourier_phases


def evaluate_fourier(n, c, s, nbins, k=None):
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
    # Question here is whether error should be sqrt(data) or sqrt(model)
    return ((hist - model) ** 2 / model).sum()
