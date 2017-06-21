import numpy as np
import matplotlib.pyplot as plot
import copy
import scipy
from scipy import ndimage
from astropy import log

from values import *

#------------------------THIS MAKES THE TOTAL COUNT HISTOGRAM---------------------------
def event_counter(etable):
    'Count events by DET_ID'
    IDevents = np.zeros_like(IDS)

    for i, id in enumerate(IDS):
        IDevents[i] = np.count_nonzero(etable['DET_ID'] == id)

    return IDevents

def hist_use(etable):
    'Creates array of event count per ID and colors to those > 1 sigma from mean red'
    # Make array of event counts by DET_ID
    IDevents = event_counter(etable)

    # Remove any that have 0 counts and compute std dev
    temp = np.delete(IDevents, np.where(IDevents == 0))
    stdev = np.std(temp)
    colors = np.array(['k']*len(IDevents))

    # Set points that are off by more that 1 sigma to red
    diff = np.array(IDevents,dtype=np.float)-np.mean(temp)
    idx = np.where(diff>stdev)[0]
    colors[idx] = 'r'

    return IDevents, colors


def plot_total_count_hist(etable, ax_rate, ax_counts):
    'Plots event count per ID as a histogram with event count and countrate on y axes'

    num_events, colors = hist_use(etable)

    tc = ax_counts.bar(IDS, num_events, color = colors)
    rate = np.array(num_events,dtype=np.float) / etable.meta['EXPOSURE']

    ax_rate.set_ylabel('c/s')
    #countrate.set_ylim([np.min(rate)-20,np.max(rate)+20])

    ax_counts.set_xlabel('DET_ID')
    ax_counts.set_ylabel('# of Events')
    plot.locator_params(nticks = 20)
    plot.title('Total Event Count by Detector')
    #total_counts.set_ylim([np.min(num_events)-20, np.max(num_events)+20])

    return num_events

#----------------------THIS MAKES THE GRAYSCALE ID/EVENT COUNT CHART---------------------
def structure(etable, num_events):
    'Creates a grid where the xy pair corresponds to RAWX,RAWY and value at each entry is event count'
    rawx = np.zeros_like(IDS)
    rawy = np.zeros_like(IDS)


    #getting RAWX and RAWY vals for each ID
    for count, id in enumerate(IDS):
        idx = np.where(etable['DET_ID']==(id))[0]
        if len(idx) > 0:
            rawx[count] = etable['RAWX'][idx][0]
            rawy[count] = etable['RAWY'][idx][0]

        else:
            rawx[count] = 0
            rawy[count] = 0

    #In STRUCTURE, each element corresponds to the geometric position
    # of each detector, while the value is the # of counts
    structure = np.zeros(shape = (7,8))

    for i in xrange(len(rawx)):
        structure[int(rawy[i])][int(rawx[i])] = num_events[i]

    return structure

def plot_detector_chart(etable, num_events,  ax_map):
    'Plots the structure created in structure() above as a grayscale grid'
    #WANT TO GET THE ORIGIN IN THE TOP RIGHT HAND CORNER
    struct = structure(etable, num_events)
    #plot.style.use('grayscale')
    ax_img = plot.imshow(struct, origin = 'lower')
    plot.gca().invert_yaxis()
    plot.gca().invert_xaxis()
    plot.rcParams.update({'font.size' : 8})
    plot.title('Total Event Count by Detector Location')
    plot.xlabel('Raw X')
    plot.ylabel('Raw Y')
    plot.colorbar(ax_img, orientation = 'horizontal')

    return

#----------------------THIS MAKES THE LIGHT CURVE---------------------------
def light_curve(etable,binsize):
    'Bins events as a histogram to be plotted as the light curve. returns bins and the histogram'
    met0 = etable.meta['TSTART']
    t = etable['MET']-met0

    # Add 1 bin to make sure last bin covers last events
    bins = np.arange(0.0,t.max()+binsize,binsize)
    sums, edges = np.histogram(t, bins=bins)

    # Chop off last bin edge, which is only for computing histogram, not plotting
    return bins[:-1], sums

def plot_light_curve(etable, lclog, binsize=1.0):
    'Compute binned light curve of events and return mean rate,plots light curve'
    bins, sums = light_curve(etable, binsize=binsize)


    # Compute mean rate
    rate = sums/binsize
    mean_rate = rate.mean()

    plot.plot(bins, rate, linewidth = .6)

    label = 'Mean Rate: {0:.3f} c/s'.format(mean_rate)
    # Plot line at mean counts per bin
    plot.plot([bins[0],bins[-1]], [mean_rate,mean_rate], 'r--', label = label)
    #plot.legend(loc = 4)
    plot.title('Light Curve')
    plot.xlabel('Time Elapsed (s)')
    plot.ylabel('c/s')
    if lclog:
    	plot.yscale('log')
    #Plot the counts / second on the other y axis

    return mean_rate

#-------------------------------THIS PLOTS THE FAST TO SLOW AND SLOW TO FAST------------------
def plot_slowfast(etable):
    'Scatter plot of PI and fast PHA, highlighting points above ratio cut'
    
    # First do some counts
    nfastonly = np.count_nonzero(np.logical_and(etable['EVENT_FLAGS'][:,FLAG_FAST],
                                            np.logical_not(etable['EVENT_FLAGS'][:,FLAG_SLOW])))
    nslowonly = np.count_nonzero(np.logical_and(etable['EVENT_FLAGS'][:,FLAG_SLOW],
                                            np.logical_not(etable['EVENT_FLAGS'][:,FLAG_FAST])))
    nboth = np.count_nonzero(np.logical_and(etable['EVENT_FLAGS'][:,FLAG_SLOW],
                                            etable['EVENT_FLAGS'][:,FLAG_FAST]))

    # Only compute ratio for events with both triggers
    etable = etable[np.logical_and(etable['EVENT_FLAGS'][:,FLAG_SLOW],etable['EVENT_FLAGS'][:,FLAG_FAST])]

    # Ratio is SLOW to FAST. Edge events should have ratio bigger than cut
    ratio = np.array(etable['PHA'],dtype=np.float)/np.array(etable['PHA_FAST'],dtype=np.float)

    ratio_cut = 1.4
    colors = np.array(['k']*len(ratio))
    idx = np.where(ratio>ratio_cut)[0]
    colors[idx] = 'r'

    plot.scatter(etable['PI'],ratio, s=.4, c = colors)

    x = np.arange(min(etable['PI']), max(etable['PI']))
    phax = np.ones_like(x)*ratio_cut
    
    plot.plot(x, phax , 'g--', linewidth = 0.5)

    plot.title('PHA Slow / Fast vs PI')
    plot.xlabel('PI')
    plot.ylabel('PHA Ratio')

    fast_str = "# of fast only  : " + str(nfastonly)
    slow_str = "# of PHA only : " + str(nslowonly)
    total =    "# of both        : " + str(nboth)
    bad = "# of bad points: " + str(len(idx))
    plot.annotate(fast_str, xy=(0.03, 0.85), xycoords='axes fraction')
    plot.annotate(slow_str, xy=(0.03, 0.8), xycoords='axes fraction')
    plot.annotate(total, xy=(0.03, 0.75), xycoords='axes fraction')
    plot.annotate(bad, xy = (.03, .7), xycoords='axes fraction')
    plot.annotate("Ratio cut = {0:.2f}".format(ratio_cut),xy=(0.65,0.1),xycoords='axes fraction')
    return

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
    'plots the energy spectrum of PI'
    plot.hist(etable['PI']*PI_TO_KEV, bins=200, range=(0.0,15.0),
        histtype='step',log=False)
    plot.yscale('log')
    plot.xscale('log')
    plot.xlim((0.1,20.0))
    plot.title('PI Spectrum')
    plot.xlabel('Energy (keV)')
    plot.ylabel('Counts')

    return
#-------------------------------THIS PLOTS THE POWER SPECTRUM (FFT)--------------
def choose_N(orig_N):
    """
    choose_N(orig_N):
        Choose a time series length that is larger than
            the input value but that is highly factorable.
            Note that the returned value must be divisible
            by at least the maximum downsample factor * 2.
            Currently, this is 8 * 2 = 16.
    """
    # A list of 4-digit numbers that are highly factorable by small primes
    goodfactors = [1008, 1024, 1056, 1120, 1152, 1200, 1232, 1280, 1296,
                   1344, 1408, 1440, 1536, 1568, 1584, 1600, 1680, 1728,
                   1760, 1792, 1920, 1936, 2000, 2016, 2048, 2112, 2160,
                   2240, 2304, 2352, 2400, 2464, 2560, 2592, 2640, 2688,
                   2800, 2816, 2880, 3024, 3072, 3136, 3168, 3200, 3360,
                   3456, 3520, 3584, 3600, 3696, 3840, 3872, 3888, 3920,
                   4000, 4032, 4096, 4224, 4320, 4400, 4480, 4608, 4704,
                   4752, 4800, 4928, 5040, 5120, 5184, 5280, 5376, 5488,
                   5600, 5632, 5760, 5808, 6000, 6048, 6144, 6160, 6272,
                   6336, 6400, 6480, 6720, 6912, 7040, 7056, 7168, 7200,
                   7392, 7680, 7744, 7776, 7840, 7920, 8000, 8064, 8192,
                   8400, 8448, 8624, 8640, 8800, 8960, 9072, 9216, 9408,
                   9504, 9600, 9680, 9856]
    if orig_N < 1024:
        return 1024
    # Get the number represented by the first 4 digits of orig_N
    first4 = int(str(orig_N)[:4])
    # Now get the number that is just bigger than orig_N
    # that has its first 4 digits equal to "factor"
    for factor in goodfactors:
        if factor > first4: break
    new_N = factor
    while new_N < orig_N:
        new_N *= 10
    # Finally, compare new_N to the closest power_of_two
    # greater than orig_N.  Take the closest.
    two_N = 2
    while two_N < orig_N:
        two_N *= 2
    if two_N < new_N: return two_N
    else: return new_N

def plot_fft_of_power(etable,nyquist):
    'plots the power spectrum'

    dt = 2.0/nyquist
    METmin = etable['MET'].min()
    T = etable['MET'].max() - etable['MET'].min()

    # Choose good number of bins for efficient FFT
    n = choose_N(T/float(dt))
    bins = np.arange(n)*dt
    log.info('{0} {1}'.format(T/dt,n))
    log.info('Computing FFT with {0} bins of {1} s, covering {2} total time'.format(n,dt,T))
    ts, edges = np.histogram(etable['MET']-METmin,bins)

    ft = np.fft.rfft(ts)
    power = (ft * ft.conj()).real
    power /= len(etable['MET'])
    power[0:100] = 0.0
    x = np.fft.rfftfreq(len(ts), dt)
    #idx = np.where(power>20)
    idx = np.argmax(power)
    print(x[idx], power[idx])
    plot.plot(x,power)
    plot.title('Power Spectrum')
    plot.xlabel('Frequency')
    plot.ylabel('Power')
    return
#-------------------------------THIS PLOTS THE DEADTIME HISTOGRAM------------------
def plot_deadtime(etable):
    'Plot histogram of detector deadtime in microseconds.'
    us = etable['DEADTIME']*1.0e6
    # Bin at 1 us resolution
    max = np.floor(us.max())+1
    plot.hist(us,range=(0.0,max), bins=int(max),log=True)
    #plot.ticklabel_format(style = 'sci', axis='x', scilimits = (0,0))
    plot.title('Histogram of deadtime')
    plot.xlabel('Deadtime [microseconds]')
    plot.ylabel('Frequency [occurrences]')

    return

#-------------------------PULSE PROFILE----------------------------------
def pulse_profile_fixed(etable, F0):
    phase = np.fmod((etable['MET']-etable['MET'][0])*F0,1.0)
    plot.hist(phase,bins=32)
    plot.ylabel('Counts')
    plot.xlabel('Pulse Phase')
    plot.title('Pulse Profile (F0={0:.6f})'.format(F0))

'''
def pulse_profile(etable, orbfile, parfile):

    import astropy.io.fits as pyfits
    import pint.toa, pint.models
    from pint.event_toas import load_NICER_TOAs
    from pint.event_toas import load_RXTE_TOAs
    from pint.plot_utils import phaseogram_binned
    from pint.observatory.nicer_obs import NICERObs
    from pint.observatory.rxte_obs import RXTEObs

   ### Make arguments for parfile and orbfile and only do this if both are present

    log.info('Event file TELESCOPE = {0}, INSTRUMENT = {1}'.format(etable.meta['TELESCOP'],
       etable.meta['INSTRUME']))
    if etable.meta['TELESCOP'] == 'NICER':
        # Instantiate NICERObs once so it gets added to the observatory registry
       if orbfile is not None:
           log.info('Setting up NICER observatory')
           NICERObs(name='NICER',FPorbname=orbfile,tt2tdb_mode='none')
       # Read event file and return list of TOA objects
       tl  = load_NICER_TOAs(eventname)

   elif hdr['TELESCOP'] == 'XTE':
       # Instantiate RXTEObs once so it gets added to the observatory registry
       if orbfile is not None:
           # Determine what observatory type is.
           log.info('Setting up RXTE observatory')
           RXTEObs(name='RXTE',FPorbname=orbfile,tt2tdb_mode='none')
       # Read event file and return list of TOA objects
       tl  = load_RXTE_TOAs(args.eventname)

   ts = pint.toa.TOAs(toalist=tl)
   ts.compute_TDBs()
   ts.compute_posvels(ephem='DE421',planets=False)


   # Load PINT model objects
   modelin = pint.models.get_model(parfile)
   log.info(str(modelin))

   phss = modelin.phase(ts.table)[1]
   # ensure all postive
   phases = np.where(phss < 0.0, phss + 1.0, phss)

    mjds = ts.get_mjds()

    # Histogram phases to make pulse profile

    return
'''
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
    'Plots reset rates'
    reset = plot.bar(IDS, reset_rates, width = .85)
    plot.title('Reset Rate by Detector')
    plot.ylabel('Reset Rate [Hz]')
    plot.xlabel('DET_ID')
    #plot.ylim([0, np.max(reset_rates)+2])
