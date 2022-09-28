import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from loguru import logger as log
from astropy.time import Time, TimeDelta

from nicer.plotutils import *


def sci_plots(etable, gtitable, args):
    # GRID SET UP
    figure2 = plt.figure(figsize=(11, 8.5), facecolor="white")
    sci_grid = gridspec.GridSpec(5, 7)

    # Build PHA Fast/Slow ratio plot before filtering by ratio
    # Only do this if powerspectrum not requested

    if not args.powspec:
        log.info("Building fast/slow subplot")
        plt.subplot(sci_grid[1:3, 2:5])
        plot_slowfast(etable, args)

    # Now, filter out the points above the ratio cut, if requested
    if args.filtratio:
        log.info("Applying ratio filter using trumpet")
        etable = filt_ratio_trumpet(etable)

    # Light Curve
    log.info("Building light curve")
    plt.subplot(sci_grid[3:5, :7])
    meanrate, a = plot_light_curve(etable, args.lclog, gtitable, binsize=args.lcbinsize)
    plot.title("Light Curve")
    plot.xlabel("Time Elapsed (s)")

    # Energy Spectrum
    log.info("Building energy spectrum")
    plt.subplot(sci_grid[1:3, :2])
    plot_energy_spec(etable)

    # Power Spectrum
    if args.powspec:
        log.info("Looking at power spectrum")
        plt.subplot(sci_grid[1:3, 2:5])
        # plot_fft_of_power(etable, args.nyquist, args.pslog, args.writeps)

    # PULSE PROFILE
    log.info("Building pulse profile")
    axprofile = plt.subplot(sci_grid[1:3, 5:7])
    if (args.orb is not None) and (args.par is not None):
        log.info("Calling pulse profile using PINT")
        pulse_profile(axprofile, etable, args)
    elif args.foldfreq > 0.0:
        log.info("Calling pulse profile with fixed frequency")
        pulse_profile_fixed(etable, args.foldfreq)
    else:
        pass

    # Making the plot all nice and stuff
    plt.subplots_adjust(
        left=0.07, right=0.99, bottom=0.05, top=0.9, wspace=0.8, hspace=0.8
    )

    figure2.suptitle(
        "ObsID {0}: {1} on {2}".format(
            etable.meta["OBS_ID"],
            etable.meta["OBJECT"],
            etable.meta["DATE-OBS"].replace("T", " at "),
        ),
        fontsize=14,
    )

    # tstart, tstop, exposure
    exposure = float(etable.meta["EXPOSURE"])
    # tstart = etable.meta['DATE-OBS'].replace('T',' at ')
    # tend = etable.meta['DATE-END'].replace('T', ' at ')
    tstart = TimeDelta(etable.meta["TSTART"], format="sec", scale="tt") + Time(
        etable.meta["MJDREFI"] + etable.meta["MJDREFF"], format="mjd", scale="tt"
    )
    tend = TimeDelta(etable.meta["TSTOP"], format="sec", scale="tt") + Time(
        etable.meta["MJDREFI"] + etable.meta["MJDREFF"], format="mjd", scale="tt"
    )
    fraction = exposure / (float(etable.meta["TSTOP"]) - float(etable.meta["TSTART"]))

    # Add text info here:
    plt.figtext(
        0.07,
        0.93,
        "Start time: {0} - End time: {1}".format(tstart.iso, tend.iso),
        fontsize=10,
    )
    plt.figtext(
        0.07,
        0.90,
        "Total clock time between start and stop: {0:.1f} s".format(
            float(etable.meta["TSTOP"]) - float(etable.meta["TSTART"])
        ),
        fontsize=10,
    )
    plt.figtext(
        0.07,
        0.87,
        "Exposure: {0:.1f} s   -->   Coverage fraction is {1:.3f}".format(
            exposure, fraction
        ),
        fontsize=10,
    )
    plt.figtext(0.07, 0.84, "Mean count rate:  {0:.3f} c/s".format(meanrate), fontsize=10)
    # plt.figtext(.07, .84, etable.meta['FILT_STR'], fontsize=10)
    if args.mask:
        plt.figtext(0.07, 0.81, "IDS {0} are masked".format(args.mask), fontsize=10)

    stringtable_start = str(gtitable["START"][:]).split('\n')
    stringtable_stop = str(gtitable["STOP"][:]).split('\n')
    stringtable_duration = str(gtitable["DURATION"][:]).split('\n')

    # In case the duration_table contains '...'
    bad_durations = np.char.endswith(stringtable_duration, '...')
    if bad_durations.any():
        idx_to_remove = np.argwhere(bad_durations).flatten()
        for idx in np.flip(idx_to_remove):
            del stringtable_start[idx]
            del stringtable_stop[idx]
            del stringtable_duration[idx]

    # Bit of formatting of the duration table
    stringtable_duration[0:2] = ["  {}".format(i) for i in stringtable_duration[0:2]]
    stringtable_duration[2] = "---{}".format(stringtable_duration[2])
    stringtable_duration[3:] = ["  {}".format(i) for i in stringtable_duration[3:]]

    # omits GTIs in the display table if there are more than 7 (the first 3 lines are for the table header)
    if len(stringtable_start)> 10:
        stringtable_start = stringtable_start[0:6] + ['...'] + stringtable_start[-3:]
        stringtable_stop = stringtable_stop[0:6] + ['...'] + stringtable_stop[-3:]
        stringtable_duration = stringtable_duration[0:6] + ['({:9.1f})'.format(np.sum(list(map(float, stringtable_duration[7:-3]))))] + stringtable_duration[-3:]
    stringtable_start ='\n'.join(stringtable_start)
    stringtable_stop ='\n'.join(stringtable_stop)
    stringtable_duration ='\n'.join(stringtable_duration)

    plt.figtext(0.5, 0.77, stringtable_start, fontsize=10, fontname='Courier')
    plt.figtext(0.58, 0.77, stringtable_stop, fontsize=10, fontname='Courier')
    plt.figtext(0.66, 0.77, stringtable_duration, fontsize=10, fontname='Courier')

    return figure2
