import numpy as np
import matplotlib.pyplot as plot
import copy
import scipy
from scipy import ndimage
from astropy import log

from nicer.values import *

'''
CTRL + F WORDS FOR QUICK TROUBLESHOOTING:

*event_counter(data1)
    counts events / ID

*plot_total_count_hist(data1, countrate, total_counts)
    plots the histogram of event count / ID. highlights those > 1 std dev from the mean.

*structure(data1, IDS, num_events)
    creates an array of rawx rawy coordinates with the total count per rawx,rawy pair.

*plot_detector_chart(data1, IDS, num_events, sci_grid, detector_map)
    plots the structure (created above) as a grayscale map of event count intensities.

*slow_fast(data1)
    calculates number of slow only, fast only, and both from pha fast and slow. This is used to create text on the slow vs fast plot.

*plot_fft_of_power(time)
    plots the power spectrum. bins data (filtered by no_more_resets) at .5 ms, takes the rfft, and plots.

*plot_deadtime(data1, sci_grid, dead)
    plots a histogram of the deadtime.

*plot_pulseprofile()
    plots pulse profile

*reset_rate(data1, event_flags)
    calculates the reset rate. Is used by plot_resetrate and the text maker on both plots.

*plot_resetrate(IDS, reset, rate)
    plots a histogram of the reset rate per ID

*plot_nicetextinfo(num_events, sci_grid, figure1, reset_rate, info)
    creates text strings that are included in both figures.

'''

#------------------------THIS MAKES THE TOTAL COUNT HISTOGRAM---------------------------
def event_counter(data1):
    IDS = np.array([0, 1, 2, 3, 4, 5,6, 7, 10, 11, 12, 13, 14, 15, 16, 17, 20, 21, 22, 23, 24, 25, 26, 27, 30, 31, 32, 33, 34, 35, 36, 37, 40, 41, 42, 43, 44, 45, 46, 47, 50, 51, 52, 53, 54, 55, 56, 57, 60, 61, 62, 63, 64, 65, 66, 67])
    IDevents = np.zeros(shape = (2, 56))
    IDevents[0] = IDS[:]

        #finding all the instances of each ID
    for id in IDS:
        IDevents[1][np.where(IDS == id)] = np.count_nonzero(data1[5] == (id))

    return IDevents

def hist_use(data1):
    IDevents = event_counter(data1)
    IDS = copy.deepcopy(IDevents[0])
    temp = np.delete(IDevents[1], np.where(IDevents[1] == 0))
    stdev = np.std(temp)
    colors = range(IDevents[0].size)

        #creating colors for the histogram
    for i in xrange(0, IDevents[1].size):
        if abs(IDevents[1][i] - np.mean(temp))> stdev:
            colors[i] = 'r'
        else:
            colors[i] = 'k'

    return IDS, IDevents[1], stdev, colors


def plot_total_count_hist(data1, countrate, total_counts):
    IDS, num_events, stdev, colors = hist_use(data1)

    tc = total_counts.bar((IDS), num_events, color = colors)
    rate = num_events / 60

    ##cr = countrate.bar(IDS, rate, color = 'y')
    countrate.set_ylabel('Counts per second')
    countrate.set_ylim([np.min(rate)-20,np.max(rate)+20])

    total_counts.set_xlabel('DET_ID')
    total_counts.set_ylabel('Total # of Events')
    plot.locator_params(nticks = 20)
    plot.title('Total Event Count by Detector')
    total_counts.set_ylim([np.min(num_events)-20, np.max(num_events)+20])

    return tc, IDS, num_events, rate

#----------------------THIS MAKES THE GRAYSCALE ID/EVENT COUNT CHART---------------------
def structure(data1, IDS, num_events):
    rawx = np.zeros(IDS.size)
    rawy = np.zeros(IDS.size)
    count = 0

    #getting x and y vals for each ID
    for id in IDS:
        x = np.where(data1[5]==(id))
        if len(x[0]) > 0:
            rawx[count] = data1[1,x[0][0]]
            rawy[count] = data1[2,x[0][0]]

        else:
            rawx[count] = 0
            rawy[count] = 0

        x = []
        count += 1

        #In STRUCTURE, each element corresponds to the geometric position of each detector, while the value is the # of counts
    structure = np.zeros(shape = (7,8))

    for i in xrange(len(rawx)):
        structure[int(rawy[i])][int(rawx[i])] = num_events[i]

    return structure

def plot_detector_chart(data1, IDS, num_events, sci_grid, detector_map):
    #WANT TO GET THE ORIGIN IN THE TOP RIGHT HAND CORNER
    struct = structure(data1, IDS, num_events)
    plot.style.use('grayscale')
    detector_map = plot.imshow(struct, origin = 'lower')
    plot.gca().invert_yaxis()
    plot.gca().invert_xaxis()
    plot.rcParams.update({'font.size' : 8})
    plot.title('Total Event Count by Detector Location')
    plot.xlabel('Raw X')
    plot.ylabel('Raw Y')
    plot.colorbar(detector_map, orientation = 'horizontal')

    return detector_map
#----------------------THIS MAKES THE LIGHT CURVE---------------------------
def light_curve(etable,binsize):
    met0 = etable['MET'].min()
    t = etable['MET']-met0

    # Add 1 bin to make sure last bin covers last events
    bins = np.arange(0.0,t.max()+binsize,binsize)
    sums, edges = np.histogram(t, bins=bins)

    # Chop off last bin edge, which is only for computing histogram, not plotting
    return bins[:-1], sums

def plot_light_curve(etable,binsize=1.0):
    'Compute binned light curve of events and return mean rate'
    bins, sums = light_curve(etable, binsize=binsize)
    light_curve_plot = plot.plot(bins, sums, linewidth = .6)


    # Compute mean rate
    rate = sums/binsize
    mean_rate = rate.mean()

    label = 'Mean Rate: {0:.3f} c/s'.format(rate.mean())
    # Plot line at mean counts per bin
    mean_counts = sums.mean()
    plot.plot([bins[0],bins[-1]], [mean_counts,mean_counts], 'r--', label = label)

    plot.legend(loc = 4)
    plot.title('Light Curve')
    plot.xlabel('Time Elapsed (s)')
    plot.ylabel('Counts')

    # Compute the mean rate
    return mean_rate

#-------------------------------THIS PLOTS THE FAST TO SLOW AND SLOW TO FAST------------------

def plot_slowfast(etable):
    'Scatter plot of slow and fast PHA, highlighting points above ratio cut'

    # Ratio is SLOW to FAST. Edge events should have ratio bigger than cut
    ratio = etable['PHA']/etable['PHA_FAST']

    ratio_cut = 1.8
    colors = np.array(['k']*len(ratio))
    idx = np.where(ratio>ratio_cut)
    colors[idx] = 'r'

    fastslow_ratio = plot.scatter(etable['PHA'], etable['PHA_FAST'], s=.4, c = colors)

    phax = np.arange(0,etable['PHA'].max())
    plot.plot(phax, phax/ratio_cut, 'g--', linewidth = 0.3)

    plot.title('PHA Fast vs. PHA Slow')
    plot.xlabel('PHA')
    plot.ylabel('PHA_FAST')

    # fast_str = "# of fast only : " + str(slow)
    # slow_str = "# of slow only : " + str(fast)
    # total =    "# of both      : " + str(total)
    # plot.annotate(fast_str, xy=(0.03, 0.85), xycoords='axes fraction')
    # plot.annotate(slow_str, xy=(0.03, 0.8), xycoords='axes fraction')
    # plot.annotate(total, xy=(0.03, 0.75), xycoords='axes fraction')

    return fastslow_ratio

#-------------------------------THIS PLOTS THE ENERGY SPECTRUM------------------
def calc_pi(etable, calfile):
    'Compute PI from PHA (slow) using approximate linear calibration'

    # Load calibration file
    det_ids, e0s, gains = np.loadtxt(calfile,unpack=True)
    det_ids = np.array(det_ids,dtype=int)

    e_keV = np.zeros_like(etable['PHA'],dtype=np.float)

    for d, e0, g  in zip(det_ids,e0s,gains):
        idx = np.where(etable['DET_ID'] == d)[0]
        e_keV[idx] = e0 + g*etable['PHA'][idx]

    pi = np.array(e_keV/PI_TO_KEV,dtype=np.int)
    return pi

def plot_energy_spec(etable):
    plot.hist(etable['PI']*PI_TO_KEV, bins=200, range=(0.0,15.0),
        histtype='step')
    plot.yscale('log')
    #plot.xscale('log')
    plot.title('PI Spectrum')
    plot.xlabel('Energy (keV)')
    plot.ylabel('Counts')

    return
#-------------------------------THIS PLOTS THE POWER SPECTRUM (FFT)--------------

def plot_fft_of_power(etable):
    #taking out the event flags
    bins = np.arange(time[0], time[-1], .0005)
    stuff, edges = np.histogram(time,bins)

    ft = np.fft.rfft(stuff)
    power = ft * ft.conj().real
    power[0:10] = 0
    power = power / len(data1[0])
    x = np.fft.rfftfreq(len(stuff), .0005)
    noflags = plot.plot(x,power)
    plot.title('Power Spectrum')
    plot.xlabel('Frequency')
    plot.ylabel('Power')
    return noflags
#-------------------------------THIS PLOTS THE DEADTIME HISTOGRAM------------------
def plot_deadtime(data1, sci_grid, dead):
    ms = data1[6]*1000000
    dead = plot.hist(ms, 50)
    plot.ticklabel_format(style = 'sci', axis='x', scilimits = (0,0))
    plot.title('Histogram of deadtime')
    plot.xlabel('Deadtime [microseconds]')
    plot.ylabel('Frequency [occurrences]')
    plot.ticklabel_format(style = 'plain')

    return dead

#-------------------------PULSE PROFILE----------------------------------
def pulse_profile(data1):

##    from astropy import log
##    import astropy.io.fits as pyfits
##    import pint
##
##    ### Make arguments for parfile and orbfile and only do this if both are present
##
##    hdr = pyfits.getheader(args.eventname,ext=1)
##
##    log.info('Event file TELESCOPE = {0}, INSTRUMENT = {1}'.format(hdr['TELESCOP'],
##        hdr['INSTRUME']))
##    if hdr['TELESCOP'] == 'NICER':
##        # Instantiate NICERObs once so it gets added to the observatory registry
##        if args.orbfile is not None:
##            log.info('Setting up NICER observatory')
##            NICERObs(name='NICER',FPorbname=args.orbfile,tt2tdb_mode='none')
##        # Read event file and return list of TOA objects
##        tl  = load_NICER_TOAs(args.eventname)
##
##    elif hdr['TELESCOP'] == 'XTE':
##        # Instantiate RXTEObs once so it gets added to the observatory registry
##        if args.orbfile is not None:
##            # Determine what observatory type is.
##            log.info('Setting up RXTE observatory')
##            RXTEObs(name='RXTE',FPorbname=args.orbfile,tt2tdb_mode='none')
##        # Read event file and return list of TOA objects
##        tl  = load_RXTE_TOAs(args.eventname)
##
##    ts = pint.toa.TOAs(toalist=tl)
##    ts.filename = args.eventname
##    ts.compute_TDBs()
##    ts.compute_posvels(ephem=args.ephem,planets=args.planets)
##
##
##    # Load PINT model objects
##    modelin = pint.models.get_model(args.parname)
##    log.info(str(modelin))
##
##    phss = modelin.phase(ts.table)[1]
##    # ensure all postive
##    phases = np.where(phss < 0.0, phss + 1.0, phss)
##
##    mjds = ts.get_mjds()

    pulse = plot.plot(data1[0],data1[3])
    return pulse
#-------------------------THIS PLOTS USEFUL TEXT AT THE TOP OF THE SUPLOT-----------
def reset_rate(etable, IDS):
    'Count resets (detector undershoots) for each detector'

    nresets = np.zeros_like(IDS)

    # For each DET_ID count the number of events with undershoot flag set
    for i in range(len(IDS)):
        idx = np.where(np.logical_and(etable['DET_ID'] == IDS[i],
                        etable['EVENT_FLAGS'][:,FLAG_UNDERSHOOT]))[0]
        nresets[i] = len(idx)

    return nresets

def plot_resetrate(IDS, reset_rates):
    reset = plot.bar(IDS, reset_rates.value, width = .85)
    plot.title('Reset Rate by Detector')
    plot.ylabel('Reset Rate [Hz]')
    plot.xlabel('DET_ID')
    #plot.ylim([0, np.max(reset_rates)+2])
