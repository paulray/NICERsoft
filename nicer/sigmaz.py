#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gammaincinv
import pint.models, pint.toa
import astropy.units as u
from pint.residuals import Residuals
import traceback

# np.seterr(all='raise')


def sigmaz(t, y, err, nseg, diagplot=False):
    """Compute sigma_z from lists of measurement times and values.

    Input:
    ------
    t: array of floats
      The measurement times (days).
    y: array of floats
      The measurement values (seconds).
    err: array of floats (1D)
      Error bars of the measurements (seconds).
    nseg : array of ints
      In each iteration, the total time span of the measurements will be split into Nseg segments. This array contains all the values of Nseg we want to use.
    diagplot: bool
      Make a diagnostic plot of the polynomial fit to the full set of measurements.

    Output:
    -------
    sz_corr : array of floats
      Values of bias-corrected sigma-z for different segment length tau.
    szerr_lower : array of floats
      Lower error bars for the sigma-z values.
    szerr_upper : array of floats
      Upper error bars for the sigma-z values.
    tz : array of floats
      The values of the segment lengths, tau (days), for which sigma-z was calculated.
    nsegments : array of ints
      How many segments of each recorded length passed all criteria for being used when calculating tau and sigma-z statistics.
    """

    # The length of the output arrays depends on how many segments meet our criteria of more than 6 points, and longer than T/sqrt(2)
    sz = []
    tz = []
    ngood = []  # How many good segments went into each tz,sz point

    toas = t
    toaerr = err
    toares = y

    # Total span of the TOAs
    durationday = toas[-1] - toas[0]  # days
    durationsec = durationday * 86400.0  # seconds

    # The amount of wiggle room for the TOAs to fall on the other side of the segment range and still be included. It's really only supposed to account for roundoff error.  We have no reason to expect TOAs to fall on the border except for the first and last TOAs of the whole batch, so I don't believe we're in danger of double counting any TOAs.
    wiggle = 1e-5

    # Polynomial order to fit (a cubic may fail to produce a good fit for a long TOA span for pulsars with a lot of red noise; it fails for the NANOGrav data set on B1937+21).
    polyorder = 3

    for iseg in nseg:
        # For each duration of length durationday/iseg compute sz.
        dur_oneseg = durationday / iseg  # how long is one segment
        ngoodsegs = 0  # Reset the counter for good segments
        C3sqr = 0  # This will accumulate values of C3sqr
        C3un_sum = 0  # This will accumulate the sum of 1/C3_sigma^2 to normalize the C3^2 weights at the end

        n_sing_matrix = (
            0  # how many segments make polyfit fail with a singular matrix error
        )
        n_few_points = 0  # how many segments have too few points
        n_short_dataspan = 0  # in how many segments the points are clustered within too small a portion of the selected time span
        n_C3_neg_var = 0  # for how many segments the C3 coefficient has a negative variance in the covariance matrix
        for jseg in range(0, iseg):  # Now loop through each segment of this length
            # for iseq > 1 there are multiple segments we need to analyze
            segrange = (toas[0] + dur_oneseg * jseg, toas[0] + dur_oneseg * (jseg + 1))
            centertime = (
                segrange[1] + segrange[0]
            ) / 2.0  # Midpoint of observation interval
            # Fit the polynomial using only the toas in the interval
            desind = np.where(
                (toas > (segrange[0] - wiggle)) & (toas < (segrange[1] + wiggle))
            )

            if (
                np.size(desind)
            ) > polyorder + 3:  # if cov. matrix needed for error estimates on fitted params
                # if (np.size(desind))>polyorder: # if cov. matrix not needed
                dataspan = np.max(toas[desind]) - np.min(toas[desind])
            else:
                n_few_points = n_few_points + 1
                continue

            # Matsakis recommends segment be longer than dur_oneseg/sqrt(2)
            if dataspan <= (dur_oneseg / np.sqrt(2)):  # xAL added this criterion
                n_short_dataspan = n_short_dataspan + 1
                continue
            else:
                res = toares[desind]
                toaerrs = toaerr[desind]

                try:
                    # NOTE: polyfit needs 1/sigma, not 1/sigma^2 weights. Times and residuals need to be in the same units, here are in seconds
                    p, pcov = np.polyfit(
                        (toas[desind] - centertime) * 86400.0,
                        res.astype(float),
                        polyorder,
                        cov=True,
                        full=False,
                        w=np.abs(1.0 / toaerrs),
                    )
                    # p = np.polyfit((toas[desind]-centertime)*86400.0,
                    #    res.astype(float),polyorder, cov=False, full=False, w = np.abs(1./toaerrs) )
                except:
                    # print('Polyfit failed!')
                    # traceback.print_exc()
                    n_sing_matrix = n_sing_matrix + 1
                    continue

                # Get C3 coefficient uncertainty from the covariance matrix
                C3variance = np.diag(pcov)[-4]
                if C3variance < 0:
                    n_C3_neg_var = n_C3_neg_var + 1
                    # print('C3variance = %e' % C3variance)
                    continue

                C3un = np.sqrt(C3variance)
                C3un_sum = (
                    C3un_sum + 1.0 / C3un**2
                )  # for normalizing weights at the end
                C3sqr = C3sqr + p[-4] ** 2 / C3un**2
                # C3sqr=C3sqr+p[0]**2    # Accumulate to eventually find avg C3^2
                ngoodsegs += (
                    1  # the number of good segments (with at least 6 TOAs in them)
                )

            # Plot data and fit for case where the full set of resids is treated as one segment
            if iseg == 1 and diagplot:
                fig = plt.figure()
                ax = fig.add_subplot(1, 1, 1)
                toas_secs = (toas[desind] - centertime) * 86400.0
                ax.plot(toas[desind], res.astype(float) * 1.0e6, "ko")
                ax.errorbar(
                    toas[desind],
                    res.astype(float) * 1.0e6,
                    yerr=toaerr[desind] * 1.0e6,
                    fmt="none",
                    color="k",
                    capsize=2.0,
                )
                ax.plot(toas[desind], np.polyval(p, toas_secs) * 1.0e6, "r")
                ax.set_xlabel("MJD")
                ax.set_ylabel("Res (us)")
                plt.title("Order-%d polynomial fit to full TOA set" % polyorder)
                # plt.savefig("sigmaz-diagnostic.png", dpi=300, format='png', bbox_inches='tight')

        print(
            "Divided data into %d segments of length %.1f days. Number of good segments: %d"
            % (iseg, dur_oneseg, ngoodsegs)
        )
        if n_few_points > 0:
            print("--->Segments with too few TOAs: %d" % n_few_points)
        if n_short_dataspan > 0:
            print("--->Segments with too short TOA span: %d" % n_short_dataspan)
        if n_sing_matrix > 0:
            print(
                "--->Segments causing singular matrix error in polyfit: %d"
                % n_sing_matrix
            )
        if n_C3_neg_var > 0:
            print("--->Segments with C3 variance <0: %d" % n_C3_neg_var)

        if ngoodsegs != 0:
            # C3sqr=C3sqr/ngoodsegs # unweighted average
            C3sqr = (
                C3sqr / C3un_sum
            )  # average weighted by the uncertainties in fitted C3 values
            sz.append(
                (dur_oneseg * 86400) ** 2 * np.sqrt(C3sqr) / (2.0 * np.sqrt(5.0))
            )  # sigma_z formula
            tz.append(dur_oneseg)  # days
            ngood.append(ngoodsegs)

    # Sigma-z bias correction and error bars
    nsegments = np.array(ngood)
    x16 = np.divide(gammaincinv(0.5 * nsegments, 0.16), 0.5 * nsegments)
    x50 = np.divide(gammaincinv(0.5 * nsegments, 0.50), 0.5 * nsegments)
    x84 = np.divide(gammaincinv(0.5 * nsegments, 0.84), 0.5 * nsegments)
    sz_corr = np.divide(sz, np.sqrt(x50))
    szerr_upper = np.multiply(sz_corr, np.sqrt(np.divide(x50, x16)) - 1.0)
    szerr_lower = np.multiply(sz_corr, 1.0 - np.sqrt(np.divide(x50, x84)))

    return sz_corr, szerr_lower, szerr_upper, tz, nsegments


def psr_sigmaz(parfile, timfile, nseg, diagplot=False):
    """Compute sigma-z from a pulsar par and tim file.

    Input:
    ------
    timfile : string
      The file containing TOAs, in Tempo/Tempo2 format.
    parfile : string
      The file containing the timing model.
    nseg : array of ints
      In each iteration, the total time span of the TOAs will be split into Nseg segments. This array contains all the values of Nseg we want to use.
    diagplot: bool
      Make a diagnostic plot of the polynomial fit to the full set of TOAs.

    Output:
    -------
    sz_corr : array of floats
      Values of bias-corrected sigma-z for different segment length tau.
    szerr_lower : array of floats
      Lower error bars for the sigma-z values.
    szerr_upper : array of floats
      Upper error bars for the sigma-z values.
    tz : array of floats
      The values of the segment lengths, tau (days), for which sigma-z was calculated.
    nsegments : array of ints
      How many segments of each recorded length passed all criteria for being used when calculating tau and sigma-z statistics.
    """

    # Read in the TOAs and timing model and compute the residuals
    m = pint.models.get_model(parfile)
    t = pint.toa.get_TOAs(timfile, ephem="DE430")
    t.compute_pulse_numbers(m)
    toas = t.get_mjds().value  # MJD
    toaerr = t.get_errors().to(u.s).value  # s
    toares = Residuals(t, m).time_resids.to(u.s).value  # s

    return sigmaz(toas, toares, toaerr, nseg, diagplot)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: sigmaz.py [parfile] [timfile]")
        sys.exit()

    par = sys.argv[1]
    tim = sys.argv[2]

    seggrid = np.logspace(0.1, 3, num=50)
    seggrid = seggrid.astype(int)
    seggrid = np.unique(seggrid)
    # seggrid = 2**np.arange(0,10)

    sz, errlo, errhi, tz, ns = psr_sigmaz(par, tim, seggrid, diagplot=True)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(tz, sz, "ko", label="Sigma-z")
    ax.errorbar(tz, sz, yerr=[errlo, errhi], fmt="none", color="k", capsize=2.0)
    ax.grid(which="both")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Sigma-z")
    plt.title("%s\n%s" % (par, tim))
    plt.show()
