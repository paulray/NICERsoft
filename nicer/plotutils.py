"""Plotting and event-selection helpers used by the quick-look scripts.

Most functions take an event table (``etable``, an
:class:`astropy.table.Table` read from a NICER event file, usually with the
``TIME`` column renamed to ``MET``), a filter-file table (``mktable``, read
from the ``.mkf`` file), and/or a GTI table (``gtitable``, with ``START``,
``STOP`` and, for the elapsed-time helpers, ``CUMTIME`` columns).

The plotting functions draw into the current :mod:`matplotlib` axes, and
the page-layout functions in :mod:`nicer.sci_plots`,
:mod:`nicer.eng_plots` and :mod:`nicer.bkg_plots` combine them into
figures.
"""

import numpy as np
import matplotlib.pyplot as plot
import matplotlib as mpl
import copy
import scipy
from scipy import ndimage
from loguru import logger as log
from glob import glob
from nicer.values import *
from os import path
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
from astropy.stats import mad_std, sigma_clipped_stats
import astropy.io.fits as pyfits
import astropy.units as u


# ------------------------THIS MAKES THE TOTAL COUNT HISTOGRAM---------------------------
def event_counter(etable):
    """Count the events from each detector.

    Parameters
    ----------
    etable : astropy.table.Table
        Event table with a ``DET_ID`` column.

    Returns
    -------
    numpy.ndarray
        Number of events for each detector, in the order of
        :data:`nicer.values.IDS`.
    """
    IDevents = np.zeros_like(IDS)

    for i, det_id in enumerate(IDS):
        IDevents[i] = np.count_nonzero(etable["DET_ID"] == det_id)
    return IDevents


def find_hot_detectors(etable):
    """Identify detectors with an anomalously high number of events.

    The per-detector counts from :func:`event_counter` are sigma-clipped
    (3 sigma, 3 iterations), and any detector with more than the clipped mean
    plus 3 times the clipped standard deviation is flagged as hot.

    Parameters
    ----------
    etable : astropy.table.Table
        Event table with a ``DET_ID`` column.

    Returns
    -------
    numpy.ndarray or None
        ``DET_ID`` values of the hot detectors, or None if none were found.
    """
    # Compute which detectors to mask on the fly
    det_events = event_counter(etable)
    # Remove detector with no events
    # temp_events = np.delete(det_events, np.where(det_events == 0))
    stats = sigma_clipped_stats(det_events, maxiters=3, sigma_lower=3, sigma_upper=3)
    bad_dets = IDS[det_events > stats[0] + 3.0 * stats[2]]
    log.info("Detector Count Mean {0}, std {1}".format(stats[0], stats[2]))
    if len(bad_dets) > 0:
        log.warning("!!! Found hot detectors {0}".format(bad_dets))
        return bad_dets
    return None


def hist_use(etable):
    """Count events per detector and pick a bar color for each detector.

    Hot detectors (see :func:`find_hot_detectors`) are colored red, all
    others black.

    Parameters
    ----------
    etable : astropy.table.Table
        Event table with a ``DET_ID`` column.

    Returns
    -------
    IDevents : numpy.ndarray
        Number of events for each detector in :data:`nicer.values.IDS`.
    colors : numpy.ndarray
        Matplotlib color for each detector.
    """
    # Make array of event counts by DET_ID
    IDevents = event_counter(etable)
    colors = np.array(["k"] * len(IDevents))

    # Color using auto identified hot detectors
    bad_dets = find_hot_detectors(etable)
    if bad_dets is not None:
        for i, det_id in enumerate(IDS):
            if det_id in bad_dets:
                colors[i] = "r"

    # Remove any that have 0 counts and compute std dev
    # temp = np.delete(IDevents, np.where(IDevents == 0))
    # stdev = np.std(temp)

    # Set points that are off by more than 2 sigma to red
    # diff = np.array(IDevents,dtype=float)-np.mean(temp)
    # idx = np.where(diff>2.0*stdev)[0]
    # colors[idx] = 'r'

    return IDevents, colors


def plot_total_count_hist(etable, ax_rate, ax_counts):
    """Plot a bar chart of the number of events from each detector.

    Hot detectors are drawn in red. The twin axis ``ax_rate`` is scaled to show
    the equivalent count rate, using ``etable.meta["EXPOSURE"]``.

    Parameters
    ----------
    etable : astropy.table.Table
        Event table with a ``DET_ID`` column and ``EXPOSURE`` in its metadata.
    ax_rate : matplotlib.axes.Axes
        Axes for the count-rate scale (typically ``ax_counts.twinx()``).
    ax_counts : matplotlib.axes.Axes
        Axes to draw the bar chart in.

    Returns
    -------
    numpy.ndarray
        Number of events for each detector in :data:`nicer.values.IDS`.
    """

    num_events, colors = hist_use(etable)

    tc = ax_counts.bar(IDS, num_events, color=colors)

    ax_rate.set_ylabel("c/s")
    cntmin, cntmax = ax_counts.get_ylim()
    ax_rate.set_ylim(
        (cntmin / etable.meta["EXPOSURE"], cntmax / etable.meta["EXPOSURE"])
    )
    # countrate.set_ylim([np.min(rate)-20,np.max(rate)+20])

    ax_counts.set_xlabel("DET_ID")
    ax_counts.set_ylabel("# of Events")
    plot.locator_params(nbins=20)
    plot.title("Total (Filtered) Event Count by Detector")
    # total_counts.set_ylim([np.min(num_events)-20, np.max(num_events)+20])

    return num_events


# ----------------------THIS MAKES THE GRAYSCALE ID/EVENT COUNT CHART---------------------
def structure(etable, num_events):
    """Arrange per-detector values on the focal-plane grid.

    The ``RAWX``, ``RAWY`` position of each detector is taken from its first
    event in ``etable``. Detectors with no events are left at zero.

    Parameters
    ----------
    etable : astropy.table.Table
        Event table with ``DET_ID``, ``RAWX`` and ``RAWY`` columns.
    num_events : array_like
        Value for each detector, in the order of :data:`nicer.values.IDS`.

    Returns
    -------
    numpy.ndarray
        Array of shape (7, 8), indexed as ``[RAWY, RAWX]``.
    """
    rawx = np.zeros_like(IDS, dtype=int)
    rawy = np.zeros_like(IDS, dtype=int)

    # getting RAWX and RAWY vals for each ID
    for count, detid in enumerate(IDS):
        idx = np.where(etable["DET_ID"] == detid)[0]
        if len(idx) > 0:
            rawx[count] = etable["RAWX"][idx][0]
            rawy[count] = etable["RAWY"][idx][0]
        else:
            log.info("No counts for det {0}".format(detid))
            rawx[count] = -1
            rawy[count] = -1

    # In STRUCTURE, each element corresponds to the geometric position
    # of each detector, while the value is the # of counts
    structure = np.zeros(shape=(7, 8))

    for i in range(len(rawx)):
        if rawx[i] >= 0 and rawy[i] >= 0:
            structure[rawy[i]][rawx[i]] = num_events[i]

    return structure


def plot_detector_chart(etable, num_events, ax_map):
    """Plot per-detector values as an image of the focal plane.

    The image is drawn in the current axes, using :func:`structure` to place
    the values.

    Parameters
    ----------
    etable : astropy.table.Table
        Event table with ``DET_ID``, ``RAWX`` and ``RAWY`` columns.
    num_events : array_like
        Value for each detector, in the order of :data:`nicer.values.IDS`.
    ax_map : matplotlib.axes.Axes
        Not used; the image is drawn in the current axes.
    """
    # WANT TO GET THE ORIGIN IN THE TOP RIGHT HAND CORNER
    struct = structure(etable, num_events)
    # plot.style.use('grayscale')
    ax_img = plot.imshow(struct, origin="lower")
    plot.gca().invert_yaxis()
    plot.gca().invert_xaxis()
    plot.rcParams.update({"font.size": 8})
    plot.title("Filtered Event Count by Detector Location")
    plot.xlabel("Raw X")
    plot.ylabel("Raw Y")
    plot.colorbar(ax_img, orientation="horizontal")

    return


# ----------------------THIS MAKES THE LIGHT CURVE---------------------------
def light_curve(etable, startmet, stopmet, binsize):
    """Bin event times into a light curve.

    Parameters
    ----------
    etable : astropy.table.Table
        Event table with a ``MET`` column, sorted by time.
    startmet, stopmet : float or None
        Time range (MET seconds) to include. If both are None, the full span of
        ``etable`` is used.
    binsize : float
        Bin size in seconds.

    Returns
    -------
    bins : numpy.ndarray
        Left edge of each bin, in seconds since ``startmet``.
    sums : numpy.ndarray
        Number of events in each bin.
    """
    if startmet is None and stopmet is None:
        startmet = etable["MET"][0]
        t = etable["MET"] - startmet
        stopmet = etable["MET"][-1]
    else:
        t = (
            etable["MET"][
                np.where(
                    np.logical_and(etable["MET"] < stopmet, etable["MET"] > startmet)
                )
            ]
            - startmet
        )

    duration = stopmet - startmet
    # Add 1 bin to make sure last bin covers last events
    bins = np.arange(0.0, duration + binsize, binsize)
    sums, edges = np.histogram(t, bins=bins, range=(0.0, duration))

    # Chop off last bin edge, which is only for computing histogram, not plotting
    return bins[:-1], sums


def gti_colormap():
    """Colors used to tell GTI segments apart in time-series plots.

    Returns
    -------
    colornames : list of str
        The color names, in order.
    cmap : matplotlib.colors.Colormap
        Colormap mapping GTI index to color.
    norm : matplotlib.colors.Normalize
        Normalization to use with ``cmap``.
    """
    colornames = [
        "black",
        "green",
        "red",
        "blue",
        "magenta",
        "orange",
        "cyan",
        "yellow",
        "gray",
    ]
    colorlevels = np.arange(len(colornames))
    cmap, norm = mpl.colors.from_levels_and_colors(
        levels=colorlevels, colors=colornames, extend="max"
    )
    return colornames, cmap, norm


def plot_light_curve(etable, lclog, gtitable, binsize=1.0, noplot=False, plot_pos=None):
    """Compute and plot a light curve over a set of GTIs.

    Each GTI is binned separately with :func:`light_curve`, and the segments
    are concatenated in elapsed good time, so gaps between GTIs are removed.
    Points are colored by GTI segment and the mean rate is drawn as a dashed
    line.

    Parameters
    ----------
    etable : astropy.table.Table
        Event table with a ``MET`` column.
    lclog : bool
        If True, use a logarithmic y axis.
    gtitable : astropy.table.Table
        GTI table with ``START`` and ``STOP`` columns.
    binsize : float, optional
        Bin size in seconds.
    noplot : bool, optional
        If True, only compute the light curve.
    plot_pos : {None, "corner", "left", "center"}, optional
        Position of this panel in a grid of plots (see
        :func:`nicer.eng_plots.plot_all_lc`), which controls the axis labels.

    Returns
    -------
    mean_rate : float
        Mean count rate (counts/s).
    sums : numpy.ndarray
        Number of events in each bin.
    """
    #'Compute binned light curve of events and return mean rate,plots light curve'
    # EDGE CASE FOR FIRST INSTANCE
    bins, sums = light_curve(
        etable, gtitable["START"][0], gtitable["STOP"][0], binsize=binsize
    )
    cc = np.zeros_like(bins, dtype=float)
    cumtime = bins[-1] + binsize
    # THE REST OF THE GOOD INTERVALS
    for i in range(1, len(gtitable["START"])):
        mybins, mysums = light_curve(
            etable, gtitable["START"][i], gtitable["STOP"][i], binsize=binsize
        )
        bins = np.append(bins, mybins + cumtime)
        cumtime += mybins[-1] + binsize
        sums = np.append(sums, mysums)
        mycolors = np.zeros_like(mybins, dtype=float) + float(i)
        cc = np.append(cc, mycolors)

    # Compute mean rate
    rate = sums / binsize
    mean_rate = rate.mean()
    if not noplot:
        colornames, cmap, norm = gti_colormap()
        plot.scatter(
            bins,
            rate,
            c=np.fmod(cc, len(colornames)),
            cmap=cmap,
            norm=norm,
            marker="+",
            label="Light Curve",
        )
        label = "Mean Rate: {0:.3f} c/s".format(mean_rate)
        # Plot line at mean counts per bin
        plot.axhline(
            y=mean_rate, xmin=bins[0], xmax=bins[-1], linestyle="dashed", label=label
        )

        bininfo = "bin size = {} sec".format(binsize)
        plot.annotate(bininfo, xy=(0.02, 0.85), xycoords="axes fraction")

        # plot.legend(loc = 4)
        if lclog:
            plot.yscale("log")
            plot.ylim(ymin=0.1)

        # Default option
        if plot_pos is None:
            plot.ylabel("c/s")

        ## Options for plot_all_spec
        if plot_pos == "corner":
            plot.ylabel("c/s")
            plot.tick_params(axis="x", which="both", bottom="on", labelbottom="on")
        if plot_pos == "left":
            plot.ylabel("c/s")
            plot.tick_params(axis="x", which="both", labelbottom="off")
        if plot_pos == "center":
            plot.tick_params(axis="x", which="both", labelbottom="off")

    return mean_rate, sums


# -------------------------------THIS PLOTS THE FAST TO SLOW___------------------
def plot_slowfast(etable, args):
    """Plot the ratio of slow to fast PI against energy.

    Only events that triggered both the slow and fast chains are used. Events
    above the "trumpet" cut (``PI/PI_FAST > 1.1 + 120/PI``, the same cut as
    :func:`filt_ratio_trumpet`) are drawn in red. If there are more than 50000
    events, a subset is plotted. The numbers of fast-only, slow-only, both and
    rejected events are written on the plot.

    Parameters
    ----------
    etable : astropy.table.Table
        Event table with ``PI``, ``PI_FAST`` and ``EVENT_FLAGS`` columns.
    args : argparse.Namespace
        Not used.
    """
    log.info("Counting slow and fast")
    # First do some counts
    nfastonly = np.count_nonzero(
        np.logical_and(
            etable["EVENT_FLAGS"][:, FLAG_FAST],
            np.logical_not(etable["EVENT_FLAGS"][:, FLAG_SLOW]),
        )
    )
    nslowonly = np.count_nonzero(
        np.logical_and(
            etable["EVENT_FLAGS"][:, FLAG_SLOW],
            np.logical_not(etable["EVENT_FLAGS"][:, FLAG_FAST]),
        )
    )
    nboth = np.count_nonzero(
        np.logical_and(
            etable["EVENT_FLAGS"][:, FLAG_SLOW], etable["EVENT_FLAGS"][:, FLAG_FAST]
        )
    )
    log.info("Using only SLOW+FAST events for ratio plot")
    # Only compute ratio for events with both triggers
    etable = etable[
        np.logical_and(
            etable["EVENT_FLAGS"][:, FLAG_SLOW], etable["EVENT_FLAGS"][:, FLAG_FAST]
        )
    ]
    downsampfac = None
    if len(etable) > 50000:
        log.warning("Too many events for ratio plot. Plotting subset of points")
        downsampfac = len(etable) // 50000
        etable = etable[::downsampfac]
    log.info("Computing ratio")
    # Ratio is SLOW to FAST. Edge events should have ratio bigger than cut
    ratio = np.array(etable["PI"], dtype=float) / np.array(
        etable["PI_FAST"], dtype=float
    )

    #    fastconst = 1.0
    #    fastsig = 250.0
    #    fastquart = 4.0e-11
    fastconst = 1.1
    fastsig = 1200.0
    fastquart = 0.0

    x = np.array(etable["PI"], dtype=float)
    ratio_cut = fastconst + (fastsig / 10.0) / x + fastquart * x**3

    colors = np.array(["k"] * len(ratio))
    idx = np.where(ratio > ratio_cut)[0]
    colors[idx] = "r"
    log.debug("Plotting the points")
    plot.scatter(etable["PI"] * PI_TO_KEV, ratio, s=0.4, c=colors)

    if downsampfac is None:
        plot.title("PI Slow to Fast Ratio vs Energy")
    else:
        plot.title("PI Slow to Fast Ratio vs Energy (SUBSET)")
    plot.xlabel("Energy")
    plot.ylabel("PI Ratio")
    plot.ylim([0.5, 2.5])
    fast_str = "# of fast only  : " + str(nfastonly)
    slow_str = "# of slow only  : " + str(nslowonly)
    total = "# of both       : " + str(nboth)
    bad = "# of bad points : " + str(len(idx))
    plot.annotate(fast_str, xy=(0.03, 0.85), xycoords="axes fraction")
    plot.annotate(slow_str, xy=(0.03, 0.8), xycoords="axes fraction")
    plot.annotate(total, xy=(0.03, 0.75), xycoords="axes fraction")
    plot.annotate(bad, xy=(0.03, 0.7), xycoords="axes fraction")
    # plot.annotate("Ratio cut = {0:.2f}".format(ratio_cut),xy=(0.65,0.85),xycoords='axes fraction')
    x = np.linspace(min(etable["PI"]), max(etable["PI"]), 100, dtype=float)
    cutline = fastconst + (fastsig / 10.0) / x + fastquart * x**3
    plot.plot(x * PI_TO_KEV, cutline, "g--", linewidth=1.5)
    return


# -------------------------------THIS PLOTS THE ENERGY SPECTRUM------------------
def calc_pi(etable, calfile):
    """Compute PI from the slow-chain PHA, using a linear calibration per detector.

    Parameters
    ----------
    etable : astropy.table.Table
        Event table with ``PHA`` and ``DET_ID`` columns.
    calfile : str
        Text file with three columns: ``DET_ID``, offset ``E0`` (keV) and
        gain (keV per PHA channel).

    Returns
    -------
    numpy.ndarray
        PI value for each event (integer, units of 10 eV).
    """

    # Load calibration file
    det_ids, e0s, gains = np.loadtxt(calfile, unpack=True)
    det_ids = np.array(det_ids, dtype=int)

    e_keV = np.zeros_like(etable["PHA"], dtype=float)

    for d, e0, g in zip(det_ids, e0s, gains):
        idx = np.where(etable["DET_ID"] == d)[0]
        e_keV[idx] = e0 + g * etable["PHA"][idx]

    pi = np.array(e_keV / PI_TO_KEV, dtype=int)
    return pi


def plot_energy_spec(etable, binscale=1.0, plot_pos=None):
    """Plot the PI energy spectrum in counts/keV on log-log axes.

    Bins are 0.02 keV wide below 2 keV and 0.1 keV wide from 2 to 15 keV, both
    multiplied by ``binscale``.

    Parameters
    ----------
    etable : astropy.table.Table
        Event table with a ``PI`` column.
    binscale : float, optional
        Factor to multiply the bin widths by.
    plot_pos : {None, "corner", "left", "bottom", "center"}, optional
        Position of this panel in a grid of plots (see
        :func:`nicer.eng_plots.plot_all_spectra`), which controls the axis
        labels.
    """
    bb = np.concatenate(
        (np.arange(0.0, 2.0, 0.02 * binscale), np.arange(2.0, 15.0, 0.1 * binscale))
    )
    hh, hh_bins = np.histogram(etable["PI"] * PI_TO_KEV, bins=bb)
    widths = bb[1:] - bb[:-1]
    plot.step(bb[:-1], hh / widths, where="post")
    #    plot.hist(etable['PI']*PI_TO_KEV, bins=bb,
    #        histtype='step',log=False,density=False)
    plot.yscale("log")
    plot.xscale("log")
    plot.xlim((0.1, 20.0))

    # Default option
    if plot_pos is None:
        plot.title("PI Spectrum")
        plot.xlabel("Energy (keV)")
        plot.ylabel("Counts/keV")

    ## Options for plot_all_spec
    if plot_pos == "corner":
        plot.xlabel("Energy (keV)")
        plot.ylabel("Counts/keV")
        plot.tick_params(axis="x", which="both", bottom="on", labelbottom="on")
    if plot_pos == "left":
        plot.ylabel("Counts/keV")
        plot.tick_params(axis="x", which="both", labelbottom="off")
    if plot_pos == "bottom":
        plot.xlabel("Energy (keV)")
    if plot_pos == "center":
        plot.tick_params(axis="x", which="both", labelbottom="off")

    return


# -------------------------------THIS PLOTS THE POWER SPECTRUM (FFT)--------------
def choose_N(orig_N):
    """Choose an efficient FFT length.

    Returns a value at least as large as ``orig_N`` that is either a power of
    two or a highly factorable number, whichever is smaller. The returned value
    is divisible by at least 16 (twice the maximum downsample factor of 8).

    Parameters
    ----------
    orig_N : int or float
        Minimum number of points.

    Returns
    -------
    int
        The chosen length (at least 1024).
    """
    # A list of 4-digit numbers that are highly factorable by small primes
    goodfactors = [
        1008,
        1024,
        1056,
        1120,
        1152,
        1200,
        1232,
        1280,
        1296,
        1344,
        1408,
        1440,
        1536,
        1568,
        1584,
        1600,
        1680,
        1728,
        1760,
        1792,
        1920,
        1936,
        2000,
        2016,
        2048,
        2112,
        2160,
        2240,
        2304,
        2352,
        2400,
        2464,
        2560,
        2592,
        2640,
        2688,
        2800,
        2816,
        2880,
        3024,
        3072,
        3136,
        3168,
        3200,
        3360,
        3456,
        3520,
        3584,
        3600,
        3696,
        3840,
        3872,
        3888,
        3920,
        4000,
        4032,
        4096,
        4224,
        4320,
        4400,
        4480,
        4608,
        4704,
        4752,
        4800,
        4928,
        5040,
        5120,
        5184,
        5280,
        5376,
        5488,
        5600,
        5632,
        5760,
        5808,
        6000,
        6048,
        6144,
        6160,
        6272,
        6336,
        6400,
        6480,
        6720,
        6912,
        7040,
        7056,
        7168,
        7200,
        7392,
        7680,
        7744,
        7776,
        7840,
        7920,
        8000,
        8064,
        8192,
        8400,
        8448,
        8624,
        8640,
        8800,
        8960,
        9072,
        9216,
        9408,
        9504,
        9600,
        9680,
        9856,
    ]
    if orig_N < 1024:
        return 1024
    # Get the number represented by the first 4 digits of orig_N
    first4 = int(str(orig_N)[:4])
    # Now get the number that is just bigger than orig_N
    # that has its first 4 digits equal to "factor"
    for factor in goodfactors:
        if factor > first4:
            break
    new_N = factor
    while new_N < orig_N:
        new_N *= 10
    # Finally, compare new_N to the closest power_of_two
    # greater than orig_N.  Take the closest.
    two_N = 2
    while two_N < orig_N:
        two_N *= 2
    if two_N < new_N:
        return two_N
    else:
        return new_N


def plot_fft_of_power(etable, nyquist, pslog, writeps):
    """Plot the power spectrum of the event times.

    Events are binned at ``0.5/nyquist`` seconds, with the number of bins
    chosen by :func:`choose_N`. The power is normalized by the number of
    events, and the lowest 50 frequency bins are set to zero. The frequency and
    power of the highest peak are printed.

    Parameters
    ----------
    etable : astropy.table.Table
        Event table with a ``MET`` column.
    nyquist : float
        Nyquist frequency (Hz).
    pslog : bool
        If True, use a logarithmic y axis.
    writeps : bool
        If True, also write the power spectrum to ``powspec.txt``.
    """

    dt = 0.5 / nyquist
    METmin = etable["MET"].min()
    T = etable["MET"].max() - etable["MET"].min()

    # Choose good number of bins for efficient FFT
    n = choose_N(T / float(dt))
    bins = np.arange(n) * dt
    log.info("{0} {1}".format(T / dt, n))
    log.info(
        "Computing FFT with {0} bins of {1} s, covering {2} total time (Nyquist = {3})".format(
            n, dt, T, nyquist
        )
    )
    ts, edges = np.histogram(etable["MET"] - METmin, bins)

    ft = np.fft.rfft(ts)
    power = (ft * ft.conj()).real
    power /= len(etable["MET"])
    power[0:50] = 0.0
    x = np.fft.rfftfreq(len(ts), dt)
    # idx = np.where(power>20)
    idx = np.argmax(power)
    print(x[idx], power[idx])
    if pslog:
        plot.semilogy(x, power)
    else:
        plot.plot(x, power)
    if writeps:
        data = np.array([x, power])
        data = data.T
        np.savetxt(file("powspec.txt", "w"), data, fmt=["%f", "%f"])

    plot.title("Power Spectrum")
    plot.xlabel("Frequency")
    plot.ylabel("Power")
    return


# -------------------------------THIS PLOTS THE DEADTIME HISTOGRAM------------------
def plot_deadtime(etable):
    """Plot a histogram of the per-event deadtime, in 1 microsecond bins.

    Parameters
    ----------
    etable : astropy.table.Table
        Event table with a ``DEADTIME`` column (seconds).
    """
    us = etable["DEADTIME"] * 1.0e6
    # Bin at 1 us resolution
    max = np.floor(us.max()) + 1
    plot.hist(us, range=(0.0, max), bins=int(max), log=True)
    # plot.ticklabel_format(style = 'sci', axis='x', scilimits = (0,0))
    plot.title("Histogram of deadtime")
    plot.xlabel("Deadtime [microseconds]")
    plot.ylabel("Frequency [occurrences]")

    return


# -------------------------PULSE PROFILE----------------------------------


def pulse_profile_fixed(etable, F0):
    """Fold the events at a fixed frequency and plot a 32-bin pulse profile.

    Phases are measured from the first event in ``etable``; no barycentering
    is done.

    Parameters
    ----------
    etable : astropy.table.Table
        Event table with a ``MET`` column.
    F0 : float
        Folding frequency (Hz).
    """
    phase = np.fmod((etable["MET"] - etable["MET"][0]) * F0, 1.0)
    plot.hist(phase, bins=32)
    plot.ylabel("Counts")
    plot.xlabel("Pulse Phase")
    plot.title("Pulse Profile (F0={0:.6f})".format(F0))


def pulse_profile(ax, etable, args):
    """Compute pulse phases with PINT and plot a 32-bin pulse profile.

    The H-test value is written on the plot. Nothing is done unless both
    ``args.orb`` and ``args.par`` are given.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw the profile in.
    etable : astropy.table.Table
        Event table with a ``T`` column of event times (as accepted by
        :class:`pint.toa.TOA`).
    args : argparse.Namespace
        Must have ``orb`` (NICER orbit file) and ``par`` (timing model) attributes.
    """
    if (args.orb is None) or (args.par is None):
        log.warning("You did not specify orbfile or parfile")
        log.info("Please input files for orb and par with --orb and --par")
        return
    import pint
    import astropy.io.fits as pyfits
    from astropy.time import TimeDelta
    import pint.toa, pint.models
    from pint.plot_utils import phaseogram_binned
    from pint.observatory.satellite_obs import get_satellite_observatory
    from pint.eventstats import hm

    ### Make arguments for parfile and orbfile and only do this if both are present
    log.info(
        "Event file TELESCOPE = {0}, INSTRUMENT = {1}".format(
            etable.meta["TELESCOP"], etable.meta["INSTRUME"]
        )
    )
    # Instantiate NICERObs once so it gets added to the observatory registry
    log.info("Setting up NICER observatory")
    get_satellite_observatory("NICER", args.orb)

    log.info("Reading model from PARFILE")
    # Load PINT model objects
    modelin = pint.models.get_model(args.par)
    log.info(str(modelin))

    # Read event file and return list of TOA objects
    log.info("doing the load_toas thing")
    # tl  = load_NICER_TOAs(pulsefilename[0])

    # Create TOA list
    tl = []
    for t in etable["T"]:
        tl.append(pint.toa.TOA(t, obs="NICER"))

    planets = False
    if "PLANET_SHAPIRO" in modelin.params:
        if modelin.PLANET_SHAPIRO.value:
            planets = True
    ts = pint.toa.get_TOAs_list(
        tl,
        planets=planets,
        include_bipm=False,
    )
    # No longer needed, since Feb 28 reprocessing
    #    log.warning('Applying -1.0s time correction to event time TOAs for pulse phase plot')
    #    ts.adjust_TOAs(TimeDelta(np.ones(len(ts.table))*-1.0*u.s,scale='tt'))
    # Note: adjust_TOAs recomputes TDBs and posvels so no need to do again.
    #    ts.compute_TDBs()
    #    ts.compute_posvels(ephem='DE421',planets=True)

    # Compute phases
    phss = modelin.phase(ts, abs_phase=True).frac
    # Strip the units, because PINT may return u.cycle
    phss = np.array(phss)
    # ensure all postive
    phases = np.where(phss < 0.0, phss + 1.0, phss)
    mjds = ts.get_mjds()

    h = hm(phases)
    if not np.isfinite(h):
        log.error("H not finite, using {0} phases!".format(len(phases)))
        print("Phases from {0} to {1}\n".format(h.min(), h.max()))
    else:
        log.info("H = {0} from {1} phases".format(h, len(phases)))
    ax.hist(phases, bins=32)
    ax.text(0.1, 0.1, "H = {0:.2f}".format(h), transform=ax.transAxes)

    # np.savetxt('{0}.phases'.format(args.basename),np.transpose([etable['MET'], etable['PI'],phases]))

    plot.ylabel("Counts")
    plot.xlabel("Pulse Phase")
    plot.title("Pulse Profile")
    return


# -------------------------OVERSHOOT RATE FOR RATIO------------------------------


def apply_gti(etable, gtitable):
    """Select the events that fall inside a set of GTIs.

    Parameters
    ----------
    etable : astropy.table.Table
        Event table with a ``MET`` column.
    gtitable : astropy.table.Table
        GTI table with ``START`` and ``STOP`` columns (MET seconds).

    Returns
    -------
    astropy.table.Table
        The events with ``START < MET < STOP`` for some GTI.
    """
    mets = etable["MET"]
    idx = np.where(
        np.logical_and(mets > gtitable["START"][0], mets < gtitable["STOP"][0])
    )
    goodlist = [etable[idx]]
    for ii in range(1, len(gtitable["START"])):
        idx = np.where(
            np.logical_and(mets > gtitable["START"][ii], mets < gtitable["STOP"][ii])
        )
        goodlist.append(etable[idx])

    return vstack(goodlist, metadata_conflicts="silent")


def convert_from_elapsed_goodtime(elapsedstarts, elapsedstops, gtitable):
    """Convert intervals in elapsed good time back to MET.

    This is the inverse of :func:`convert_to_elapsed_goodtime`, used for
    selecting time ranges interactively on a light curve plot.

    Parameters
    ----------
    elapsedstarts, elapsedstops : array_like
        Interval start and stop times in elapsed good time (seconds).
    gtitable : astropy.table.Table
        GTI table with ``START`` and ``CUMTIME`` columns.

    Returns
    -------
    startmets, stopmets : numpy.ndarray
        The intervals in MET seconds.
    """
    startmets = np.array([])
    stopmets = np.array([])
    for timex in range(0, len(elapsedstarts)):
        starttime = elapsedstarts[timex]
        stoptime = elapsedstarts[timex]
        for idx in range(0, len(gtitable)):
            if idx < len(gtitable) - 1:
                if np.logical_and(
                    starttime >= gtitable["CUMTIME"][idx],
                    stoptime < gtitable["CUMTIME"][idx + 1],
                ):
                    startmets = np.append(startmets, starttime + gtitable["START"][idx])
                    stopmets = np.append(startmets, stoptime + gtitable["START"][idx])
                else:
                    continue
            else:
                if starttime > gtitable["CUMTIME"][idx]:
                    startmets = np.append(startmets, starttime + gtitable["START"][idx])
                    stopmets = np.append(startmets, stoptime + gtitable["START"][idx])
    return startmets, stopmets


def convert_to_elapsed_goodtime(mets, vals, gtitable):
    """Select values that fall inside GTIs and convert their times to elapsed good time.

    Elapsed good time is the time since the start of the GTI plus the total
    duration of all earlier GTIs (the ``CUMTIME`` column), so plots made
    against it have no gaps between GTIs.

    Parameters
    ----------
    mets : array_like
        Times of the values (MET seconds).
    vals : array_like
        Values at those times.
    gtitable : astropy.table.Table
        GTI table with ``START``, ``STOP`` and ``CUMTIME`` columns.

    Returns
    -------
    etimes : numpy.ndarray
        Elapsed good time of each selected value (seconds).
    goodvals : numpy.ndarray
        The selected values.
    cc : numpy.ndarray
        Index of the GTI each value falls in, for coloring plots by GTI.
    """
    mets = np.asarray(mets)
    vals = np.asarray(vals)
    idx = np.where(
        np.logical_and(mets > gtitable["START"][0], mets < gtitable["STOP"][0])
    )
    goodvals = vals[idx]
    etimes = mets[idx] - gtitable["START"][0]
    cc = np.zeros_like(goodvals, dtype=float)
    for ii in range(1, len(gtitable["START"])):
        idx = np.where(
            np.logical_and(mets > gtitable["START"][ii], mets < gtitable["STOP"][ii])
        )
        goodvals = np.append(goodvals, vals[idx])
        etimes = np.append(
            etimes, mets[idx] - gtitable["START"][ii] + gtitable["CUMTIME"][ii]
        )
        cc = np.append(cc, np.zeros_like(vals[idx], dtype=float) + float(ii))

    # Returns the arrays of elapsed times, values, and an array of what segment it is in, used for setting plot colors by GTI segment
    return etimes, goodvals, cc


# def plot_overshoot(etable, overshootrate, gtitable, args, hkmet, bothrate, mktable):
def plot_overshoot(mktable, ovbintable, gtitable, args):
    """Plot the overshoot rate against elapsed good time.

    Plots ``FPM_OVERONLY_COUNT`` from the filter file, colored by GTI, along
    with the expected overshoot rate ``1.52 * COR_SAX**(-0.633)``. If
    ``ovbintable`` is given, the binned overshoot rate is also plotted;
    otherwise ``FPM_DOUBLE_COUNT`` (events with both overshoot and undershoot
    flags) is plotted.

    Parameters
    ----------
    mktable : astropy.table.Table
        Filter (``.mkf``) table.
    ovbintable : astropy.table.Table or None
        Binned overshoot table from
        :meth:`nicer.NicerFileSet.NicerFileSet.getbinnedovershoots`.
    gtitable : astropy.table.Table
        GTI table with ``START``, ``STOP`` and ``CUMTIME`` columns.
    args : argparse.Namespace
        Not used.
    """

    # etime, overshoot, cc = convert_to_elapsed_goodtime(hkmet, overshootrate, gtitable)
    # etime, overshoot, cc = convert_to_elapsed_goodtime(
    #     mktable["TIME"], mktable["NUM_FPM_ON"] * mktable["FPM_OVERONLY_COUNT"], gtitable
    # )
    etime, overshoot, cc = convert_to_elapsed_goodtime(
        mktable["TIME"], mktable["FPM_OVERONLY_COUNT"], gtitable
    )

    cortime, cor, corcc = convert_to_elapsed_goodtime(
        mktable["TIME"], mktable["COR_SAX"], gtitable
    )

    # Delete the nan values from undershoot and cc so they don't plot
    etimeNan = np.where(np.isnan(overshoot))
    ovtime = np.delete(etime, etimeNan)
    overshoot = np.delete(overshoot, etimeNan)
    cor = np.delete(cor, etimeNan)
    cc = np.delete(cc, etimeNan)

    colornames, cmap, norm = gti_colormap()

    plot.scatter(
        ovtime,
        overshoot,
        c=np.fmod(cc, len(colornames)),
        cmap=cmap,
        norm=norm,
        marker="+",
        label="Overshoots/FPM",
    )

    # Add this to plot overonly_expr="1.52*OVERONLY_NORM*COR_SAX**(-0.633)"
    OVERONLY_NORM = 1.0
    plot.plot(
        ovtime,
        1.52 * OVERONLY_NORM * cor ** (-0.633),
        color="r",
        linewidth=3.0,
        label="OVERONLY_EXPR",
    )

    if ovbintable is not None:
        etime, binnedOV, cc = convert_to_elapsed_goodtime(
            ovbintable["TIME"],
            ovbintable["FPM_OVERONLY_COUNT"],
            gtitable,
        )
        plot.plot(etime, binnedOV, linewidth=2.0)
    else:
        # if bothrate is not None:
        #      etime, both, cc = convert_to_elapsed_goodtime(hkmet, bothrate, gtitable)
        etime, both, cc = convert_to_elapsed_goodtime(
            mktable["TIME"],
            mktable["FPM_DOUBLE_COUNT"],
            gtitable,
        )
        plot.scatter(
            etime, both, color="c", marker=".", label="Both Under and Over Flags"
        )
        # plot.yscale("symlog", linthreshy=10.0)

    plot.legend(loc=2)
    plot.ylabel("Overshoots/FPM")
    plot.grid(True)

    return


# def plot_SAA(mktable, gtitable, overshootrate):
def plot_SAA(mktable, gtitable):
    """Mark the times the spacecraft was in the SAA.

    SAA times (``NICER_SAA`` flag in the filter file) are drawn at the height
    of the maximum overshoot rate, so they can be overlaid on
    :func:`plot_overshoot`.

    Parameters
    ----------
    mktable : astropy.table.Table
        Filter (``.mkf``) table.
    gtitable : astropy.table.Table
        GTI table with ``START``, ``STOP`` and ``CUMTIME`` columns.
    """
    time, insaa, colors = convert_to_elapsed_goodtime(
        mktable["TIME"], mktable["NICER_SAA"], gtitable
    )
    time = np.delete(time, np.where(insaa == 0))
    insaa = np.delete(insaa, np.where(insaa == 0))

    # Get Overshoot Max, just to put the SAA on the same plot
    ovtime, overshoot, cc = convert_to_elapsed_goodtime(
        mktable["TIME"], mktable["NUM_FPM_ON"] * mktable["FPM_OVERONLY_COUNT"], gtitable
    )
    # insaa[np.where(insaa == 1)] = max(overshootrate)
    insaa[np.where(insaa == 1)] = np.nanmax(overshoot)

    plot.scatter(time, insaa, color="y", label="In the SAA", marker="_")
    plot.legend(loc=2)
    return


# -------------------------UNDERSHOOT RATE FOR RATIO------------------------------
# def plot_undershoot(etable, undershootrate, gtitable, args, hkmet, mktable):
def plot_undershoot(mktable, gtitable, args):
    """Plot the undershoot rate against elapsed good time, marking times in sunshine.

    Parameters
    ----------
    mktable : astropy.table.Table
        Filter (``.mkf``) table with ``FPM_UNDERONLY_COUNT`` and ``SUNSHINE``
        columns.
    gtitable : astropy.table.Table
        GTI table with ``START``, ``STOP`` and ``CUMTIME`` columns.
    args : argparse.Namespace
        If ``args.lclog`` is True, a logarithmic y axis is used.
    """

    # etime, undershoot, cc = convert_to_elapsed_goodtime(hkmet, undershootrate, gtitable)
    # etime, undershoot, cc = convert_to_elapsed_goodtime(mktable['TIME'], mktable['NUM_FPM_ON']*mktable['FPM_UNDERONLY_COUNT'], gtitable)
    etime, undershoot, cc = convert_to_elapsed_goodtime(
        mktable["TIME"], mktable["FPM_UNDERONLY_COUNT"], gtitable
    )

    # Delete the nan values from undershoot and cc so they don't plot
    etimeNan = np.where(np.isnan(undershoot))
    udtime = np.delete(etime, etimeNan)
    undershoot = np.delete(undershoot, etimeNan)
    cc = np.delete(cc, etimeNan)

    colornames, cmap, norm = gti_colormap()
    plot.scatter(
        udtime,
        undershoot,
        c=np.fmod(cc, len(colornames)),
        cmap=cmap,
        norm=norm,
        marker="+",
    )

    # Add sunshine
    sunmet = mktable["TIME"]
    sunshine = mktable["SUNSHINE"]
    sunt, suny, suncc = convert_to_elapsed_goodtime(sunmet, sunshine, gtitable)
    sidx = np.where(suny == 0)
    # Delete the 0 values so they don't plot
    sunt = np.delete(sunt, sidx)
    suny = np.delete(suny, sidx)

    plot.scatter(
        sunt, 1.1 * suny * undershoot.max(), color="y", label="Sunshine", marker="_"
    )
    plot.legend(loc=2)
    plot.grid(True)
    plot.ylabel("Undershoots/FPM")
    if args.lclog:
        plot.yscale("log")

    return


# -------------------------SUN / EARTH / MOON ANGLES-----------------------------
def plot_angles(mktable, gtitable):
    """Plot the Sun, bright Earth and Moon angles and the elevation against elapsed good time.

    Parameters
    ----------
    mktable : astropy.table.Table
        Filter (``.mkf``) table with ``SUN_ANGLE``, ``BR_EARTH``,
        ``MOON_ANGLE`` and ``ELV`` columns.
    gtitable : astropy.table.Table
        GTI table with ``START``, ``STOP`` and ``CUMTIME`` columns.
    """
    sun = mktable["SUN_ANGLE"]
    earth = mktable["BR_EARTH"]
    moon = mktable["MOON_ANGLE"]
    elv = mktable["ELV"]
    met = mktable["TIME"]

    goodtime, sunangle, cc = convert_to_elapsed_goodtime(met, sun, gtitable)
    goodtime, earthangle, cc = convert_to_elapsed_goodtime(met, earth, gtitable)
    goodtime, moonangle, cc = convert_to_elapsed_goodtime(met, moon, gtitable)
    goodtime, elvangle, cc = convert_to_elapsed_goodtime(met, elv, gtitable)

    plot.scatter(goodtime, sunangle, marker=".", s=3, color="y", alpha=0.5, label="Sun")
    plot.scatter(
        goodtime,
        earthangle,
        marker=".",
        s=3,
        color="b",
        alpha=0.5,
        label="Bright Earth",
    )
    plot.scatter(
        goodtime, moonangle, marker=".", s=3, color="grey", alpha=0.5, label="Moon"
    )
    plot.scatter(goodtime, elvangle, marker=".", s=3, color="m", alpha=0.5, label="ELV")
    plot.legend(loc=2)
    plot.ylim((0.0, 180.0))
    plot.grid(True)
    # plot.yticks([0.0, 45.0, 90.0, 135.0, 180.0])
    # plot.yticks(np.arange(0, 180, 15.0), minor=True)
    frame1 = plot.gca()
    frame1.axes.set_yticks([0.0, 45.0, 90.0, 135.0, 180.0], minor=False)
    frame1.axes.set_yticks(np.arange(0, 180, 15.0), minor=True)
    # plot.minorticks_on()
    plot.grid(which="major", axis="y", linewidth=0.5, alpha=1.0)
    plot.grid(which="minor", axis="y", linewidth=0.5, alpha=0.4)
    plot.ylabel("Angle (deg)")

    return


# --------------------------POINTING---------------------------------------------
def plot_pointing(mktable, gtitable):
    """Plot the pointing offset (``ANG_DIST``) against elapsed good time.

    A red line marks 1 arcminute.

    Parameters
    ----------
    mktable : astropy.table.Table
        Filter (``.mkf``) table.
    gtitable : astropy.table.Table
        GTI table with ``START``, ``STOP`` and ``CUMTIME`` columns.
    """
    time, pointing, cc = convert_to_elapsed_goodtime(
        mktable["TIME"], mktable["ANG_DIST"], gtitable
    )

    colornames, cmap, norm = gti_colormap()

    plot.scatter(
        time,
        pointing,
        c=np.fmod(cc, len(colornames)),
        cmap=cmap,
        norm=norm,
        marker="+",
        label="Pointing Offset",
    )
    plot.legend(loc=2)
    plot.ylabel("Angle (deg)")
    plot.yscale("log")
    plot.axhline(1.0 / 60.0, c="r")
    plot.ylim((0.0001, 100.0))

    return


# --------------------------LAT / LON---------------------------------------------
def plot_latlon(mktable, gtitable):
    """Plot the spacecraft latitude against elapsed good time.

    Parameters
    ----------
    mktable : astropy.table.Table
        Filter (``.mkf``) table with a ``SAT_LAT`` column.
    gtitable : astropy.table.Table
        GTI table with ``START``, ``STOP`` and ``CUMTIME`` columns.
    """
    time, lat, cc = convert_to_elapsed_goodtime(
        mktable["TIME"], mktable["SAT_LAT"], gtitable
    )
    # time, lon, cc2 = convert_to_elapsed_goodtime(mktable['TIME'], mktable['SAT_LON'], gtitable)

    colornames, cmap, norm = gti_colormap()

    plot.scatter(
        time,
        lat,
        c=np.fmod(cc, len(colornames)),
        norm=norm,
        cmap=cmap,
        marker="^",
        label="Latitude",
    )
    # plot.scatter(time, lon, c = colors, cmap = cmap, marker = '_', label = 'Longitude')
    plot.legend(loc=2)
    plot.ylim((-60.0, 60.0))
    plot.xlabel("Elapsed Time (s)", labelpad=1)
    plot.grid(True)
    plot.ylabel("Degrees")
    return


def plot_cor(mktable, gtitable):
    """Plot the cutoff rigidity (``COR_SAX``, GeV) against elapsed good time.

    A dashed line marks 1.5 GeV.

    Parameters
    ----------
    mktable : astropy.table.Table
        Filter (``.mkf``) table.
    gtitable : astropy.table.Table
        GTI table with ``START``, ``STOP`` and ``CUMTIME`` columns.
    """
    time, cor, cc = convert_to_elapsed_goodtime(
        mktable["TIME"], mktable["COR_SAX"], gtitable
    )
    # time, lon, cc2 = convert_to_elapsed_goodtime(mktable['TIME'], mktable['SAT_LON'], gtitable)

    colornames, cmap, norm = gti_colormap()

    plot.scatter(
        time,
        cor,
        c=np.fmod(cc, len(colornames)),
        norm=norm,
        cmap=cmap,
        marker=".",
        s=3,
        label="COR_SAX",
    )
    # plot.scatter(time, lon, c = colors, cmap = cmap, marker = '_', label = 'Longitude')
    plot.legend(loc=2)
    plot.ylim((0.0, 20.0))
    plot.axhline(1.5, linestyle="--", color="r")
    plot.xlabel("Elapsed Time (s)", labelpad=1)
    plot.grid(True)
    plot.minorticks_on()
    plot.grid(which="minor", axis="y", linewidth=0.5, alpha=0.3)
    plot.ylabel("GeV")
    return


# -------------------------THIS PLOTS USEFUL TEXT AT THE TOP OF THE SUPLOT-------
# def calc_nresets(etable, IDS):
def calc_nresets(mktable, IDS):
    """Count the undershoot-only resets of each detector.

    Parameters
    ----------
    mktable : astropy.table.Table
        Filter (``.mkf``) table with an ``MPU_UNDERONLY_COUNT`` column.
    IDS : numpy.ndarray
        Detector IDs (normally :data:`nicer.values.IDS`).

    Returns
    -------
    numpy.ndarray
        Total number of resets for each detector, summed over the table.
    """

    nresets = np.zeros_like(IDS)

    #     # For each DET_ID count the number of events with undershoot flag set
    #     for i in range(len(IDS)):
    #         idx = np.where(np.logical_and(etable['DET_ID'] == IDS[i],
    #                         etable['EVENT_FLAGS'][:,FLAG_UNDERSHOOT]))[0]
    #         nresets[i] = len(idx)

    MPU_ud_only = mktable["MPU_UNDERONLY_COUNT"].reshape(
        (len(mktable["MPU_UNDERONLY_COUNT"]), 56)
    )
    # Reshaping is necessary in the case of MPU_UNDERONLY_COUNT
    # given as a 7x8 matrix instead of 56 array in the mktable

    for i in range(len(IDS)):
        nresets[i] = MPU_ud_only[:, i].sum()

    return nresets


def plot_resetrate(IDS, reset_rates):
    """Plot a bar chart of the reset rate of each detector.

    Parameters
    ----------
    IDS : array_like
        Detector IDs.
    reset_rates : array_like
        Reset rate (Hz) for each detector.
    """
    reset = plot.bar(IDS, reset_rates, width=0.85)
    plot.title("Undershoot-only reset rate by Detector")
    plot.ylabel("Reset Rate [Hz]")
    plot.xlabel("DET_ID")
    # plot.ylim([0, np.max(reset_rates)+2])


# ----------------------Filter Ratio Cut-----------------------------------------
def filt_ratio(etable, ratiocut):
    """Remove events with a slow/fast PHA ratio at or above a cut.

    Events whose ``PHA`` or ``PHA_FAST`` is not positive are kept.

    Parameters
    ----------
    etable : astropy.table.Table
        Event table with ``PHA`` and ``PHA_FAST`` columns.
    ratiocut : float
        Maximum allowed ``PHA/PHA_FAST``.

    Returns
    -------
    astropy.table.Table
        The events that pass the cut.
    """
    # Filters out the points > filtratio
    ratio = np.zeros_like(etable["PI"], dtype=float)
    idx = np.where(np.logical_and(etable["PHA"] > 0, etable["PHA_FAST"] > 0))[0]
    ratio[idx] = np.asarray(etable["PHA"][idx], dtype=float) / np.asarray(
        etable["PHA_FAST"][idx], dtype=float
    )
    etable = etable[np.where(ratio < ratiocut)[0]]
    return etable


def filt_ratio_trumpet(etable):
    """Remove events above the energy-dependent "trumpet" cut on the slow/fast PI ratio.

    Events are kept if ``PI/PI_FAST < 1.1 + 120/PI``. Events whose ``PI`` or
    ``PI_FAST`` is not positive are also kept.

    Parameters
    ----------
    etable : astropy.table.Table
        Event table with ``PI`` and ``PI_FAST`` columns.

    Returns
    -------
    astropy.table.Table
        The events that pass the cut.
    """
    # Filters out the points > filtratio
    ratio = np.zeros_like(etable["PI"], dtype=float)
    idx = np.where(np.logical_and(etable["PI"] > 0, etable["PI_FAST"] > 0))[0]
    ratio[idx] = np.asarray(etable["PI"][idx], dtype=float) / np.asarray(
        etable["PI_FAST"][idx], dtype=float
    )
    # Old version
    # fastconst = 1.0
    # fastsig = 250.0
    # fastquart = 4.0e-11
    # New version 2018 Feb
    fastconst = 1.1
    fastsig = 1200.0
    fastquart = 0.0

    x = np.array(etable["PI"], dtype=float)
    ratio_cut = fastconst + (fastsig / 10.0) / x + fastquart * x**3

    etable = etable[np.where(ratio < ratio_cut)[0]]
    return etable


def format_gti_longstring(gticolumn, nCut=3):
    """Format one GTI column as fixed-width text for annotating a figure.

    If there are more than ``2*nCut`` rows, only the first and last ``nCut``
    are shown, and the omitted rows are replaced by ``...`` (or, for the
    ``DURATION`` column, by their total duration in parentheses).

    Parameters
    ----------
    gticolumn : astropy.table.Column
        A GTI table column, with a unit.
    nCut : int, optional
        Number of rows to show at each end.

    Returns
    -------
    str
        The formatted column, with a header.
    int
        Number of omitted rows.
    """

    gtidata = gticolumn.data
    units = "({})".format(gticolumn.unit.name)

    if gticolumn.name == "DURATION":
        omitted_gtis_duration = np.sum(gtidata[nCut:-nCut])
        sumgti = "({:.1f})".format(omitted_gtis_duration)
    else:
        sumgti = "..."

    # Find column width
    width = max(len(str(np.max(gtidata))) + 1, len(gticolumn.name), len(sumgti))

    # String header
    string_gtidata = "{}\n{}\n{}\n".format(
        gticolumn.name.center(width, " "), units.center(width, " "), "_" * width
    )

    # String body (without omitted GTIs is too many GTIs)
    if len(gtidata) > nCut * 2:
        short_gtidata = np.delete(gtidata, np.arange(nCut, len(gtidata) - nCut, 1))
        string_gtidata += "\n".join(
            ["{val: {width}}".format(val=i, width=width) for i in short_gtidata[:nCut]]
        )
        string_gtidata += "\n{}\n".format(sumgti.rjust(width, " "))
        string_gtidata += "\n".join(
            ["{val: {width}}".format(val=i, width=width) for i in short_gtidata[-nCut:]]
        )
    else:
        short_gtidata = gtidata
        string_gtidata += "\n".join(
            ["{val: {width}}".format(val=i, width=width) for i in short_gtidata]
        )

    return string_gtidata, len(gtidata[nCut:-nCut])
