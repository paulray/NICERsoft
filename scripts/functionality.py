import numpy as np
import matplotlib.pyplot as plot
import copy
import scipy
from scipy import ndimage

from nicer.values import *

#CURRENTLY RESET RATE IS NOT FUNCTIONAL
'''
CTRL + F WORDS FOR QUICK TROUBLESHOOTING:

*get_FITS(filename)
    gets data from 1 FITS file and puts it into an array called data1, along with event_flags and the identifying information in info

*smush()
    takes data from multiple FITS files for the same time scale and organizes them into one data1 array, an event_flags array, and 1 info array.

*event_counter(data1)
    counts events / ID

*plot_total_count_hist(data1, countrate, total_counts)
    plots the histogram of event count / ID. highlights those > 1 std dev from the mean.

*structure(data1, IDS, num_events)
    creates an array of rawx rawy coordinates with the total count per rawx,rawy pair.

*plot_detector_chart(data1, IDS, num_events, sci_grid, detector_map)
    plots the structure (created above) as a grayscale map of event count intensities.

*light_curve(data1,event_flags)
    creates a histogram of total event count / second. Calls no_more_resets which filters out all RESETS

*plot_light_curve(data1, sci_grid, light_curve_plot, event_flags)
    plots the light curve based on the histogram created by light_curve

*slow_fast(data1)
    calculates number of slow only, fast only, and both from pha fast and slow. This is used to create text on the slow vs fast plot.

*plot_slowx_fasty(data1, sci_grid, fastslow_ratio, fast, slow, total)
    plots slow on x, fast on y. Includes an expected maximum ratio of 1.8. Highlights fast data > ratio = 1.8. Includes text created in slow_fast and a legend.

*PI(data1)
    maps pha_slow values to a PI value based on E0 and gains which are defined here.

*plot_power_spec(pha_slow, sci_grid, power_spec, event_flags)
    Plots the energy spectrum [sorry the name is wrong]

*no_more_resets(data, event_flags)
    This pulls out all data with the same indices as a RESET as flagged by event_flags. THIS IS CALLED BY MULTIPLE OTHER FUNCTIONS AND IS A GOOD PLACE TO TROUBLESHOOT.

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
def light_curve(data1, event_flags):
    time = no_more_resets(data1[0], event_flags)
    bins = np.arange(time[0],time[-1],1)
    sums, edges = np.histogram(time, bins)
    average_counts = np.mean(sums)

    return sums, average_counts

def plot_light_curve(data1, sci_grid, light_curve_plot, event_flags):
    sums, average_counts = light_curve(data1, event_flags)
    time_elapsed = range(1,len(sums)+1)
    light_curve_plot = plot.plot(time_elapsed, sums, linewidth = .6)

    #Plotting the Mean line
    mean_shown = np.array([average_counts for i in xrange(len(time_elapsed))])
    label = 'Mean Rate: ' + str(int(average_counts)) + ' counts / sec'
    plot.plot(time_elapsed, mean_shown, 'r--', label = label)
    plot.legend(loc = 4)
    plot.title('Light Curve')
    plot.xlabel('Time Elapsed (s)')
    plot.ylabel('Counts')

    return light_curve_plot

#-------------------------------THIS PLOTS THE FAST TO SLOW AND SLOW TO FAST------------------
def slow_fast(data1):
    fast = len(np.where(data1[4] == 0)[0])
    slow = len(np.where(data1[3] ==0)[0])
    total = len(data1[3]) - fast - slow

    return fast, slow, total

def plot_slowx_fasty(data1, sci_grid, fastslow_ratio, fast, slow, total):
    acceptable_range = data1[3] * 1.8
    colors = range(len(acceptable_range))

    for i in xrange(0,len(data1[4])-1):
        if data1[4][i] > acceptable_range[i]:
            colors[i] = 'r'
        else:
            colors[i] = 'k'

    fastslow_ratio = plot.scatter(data1[3], data1[4], s=.4, c = colors)
    label = 'Expected ratio of 1.8'
    plot.plot(data1[3], acceptable_range, 'g', linewidth = '.3', label = label)

    plot.legend(loc = 0)
    plot.title('PHA Slow vs. PHA Fast')
    plot.xlabel('PHA Slow')
    plot.ylabel('PHA Fast')

    fast_str = "# of fast only: " + str(slow)
    slow_str = "# of slow only: " + str(fast)
    total = "# of both: " + str(total)
    plot.annotate(fast_str, xy=(0.03, 0.85), xycoords='axes fraction')
    plot.annotate(slow_str, xy=(0.03, 0.8), xycoords='axes fraction')
    plot.annotate(total, xy=(0.03, 0.75), xycoords='axes fraction')

    return fastslow_ratio

#-------------------------------THIS PLOTS THE ENERGY SPECTRUM------------------
def PI(data1):

    '''
    gains[0] = DET_IDS
    gains[1] = E0
    gains[2] = gains
    '''

    a = np.array([0, 1, 2, 3, 4, 5,6, 7, 10, 11, 12, 13, 14, 15, 16, 17, 20, 21, 22, 23, 24, 25, 26, 27, 30, 31, 32, 33, 34, 35, 36, 37, 40, 41, 42, 43, 44, 45, 46, 47, 50, 51, 52, 53, 54, 55, 56, 57, 60, 61, 62, 63, 64, 65, 66, 67])
    b = -1 * np.array([1.192, 1.099, 1.042, 1.147, 1.279, 1.277, 1.215, 1.311, .865, 1.190, 1.044,1.221,1.240,1.212,1.557,1.356,1.294,1.281,1.315,1.023,1.260, 1.076,1.093,1.195,1.058,1.293,1.344,1.386,1.046,1.204,1.143,1.316,1.247,1.148,1.423,1.191,1.244,1.297,1.298,1.227, 1.280,1.351,1.191,1.379,1.288,1.238,1.128,1.110,1.035,1.060,1.228,1.231,1.261,1.314,1.110,1.247])
    c = .001* np.array([3.6645, 3.6818,3.6945,3.886,3.7975,3.7161,3.6516,3.7166,3.7661,3.7229,3.5702,3.6318,3.7039,3.7458,3.7408,3.7955,3.7282,3.6350,3.7843,3.6560,3.6763,3.6945,3.7050,3.6872,3.6768,3.6293,3.7926,3.6682,3.7413,3.7667,3.7721,3.5859,3.6732,3.7951,3.8908,3.6876,3.7517,3.6848,3.7658,3.6901,3.7326,3.6779,3.7166,3.7133,3.6560,3.7387,3.7423,3.6756,3.7016,3.7115,3.7597,3.6493,3.7406,3.6480,3.7113,3.7657])
    converter = np.vstack((a,b,c))
    PI_list = np.zeros(len(data1[0])-1)

    count = 0
    for i in xrange(0,len(data1[0])-1):
        E0 = converter[1][np.where(converter[0] == data1[5][i])[0][0]]
        gain = converter[2][np.where(converter[0] == data1[5][i])[0][0]]
        PI_list[count] = E0 + (gain*data1[3][i])
        count += 1

    return PI_list

def plot_power_spec(data1, sci_grid, power_spec, event_flags, PI_flag):
    if PI_flag:
        PIshit = PI(data1)
        pha_slow = no_more_resets(PIshit, event_flags)
    else:
        pha_slow = no_more_resets(data1[3],event_flags)
    power_spec = plot.hist((pha_slow / (1000)), 1000, histtype='step')
    plot.yscale('log')
    plot.title('Energy Spectrum')
    plot.xlabel('Energy (keV)')
    plot.ylabel('Counts')

    return power_spec
#-------------------------------THIS PLOTS THE POWER SPECTRUM (FFT)--------------
def no_more_resets(data, event_flags):
    temp = np.array(copy.deepcopy(data))

    for i in xrange(0,len(event_flags)-1):
        if event_flags[i][5]:
            temp[i] = 0
        elif event_flags[i][6]:
            temp[i] = 0
        elif event_flags[i][7]:
            temp[i] = 0

    filtered = temp[np.nonzero(temp)]

    return filtered

def plot_fft_of_power(time, data1):
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
