import numpy as np
import matplotlib.pyplot as plot
import matplotlib as mpl
import copy
import scipy
from scipy import ndimage
from astropy import log
from glob import glob
from values import *
from os import path
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
import astropy.io.fits as pyfits
import astropy.units as u
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

    # Set points that are off by more that 2 sigma to red
    diff = np.array(IDevents,dtype=np.float)-np.mean(temp)
    idx = np.where(diff>2.0*stdev)[0]
    colors[idx] = 'r'

    return IDevents, colors


def plot_total_count_hist(etable, ax_rate, ax_counts):
    'Plots event count per ID as a histogram with event count and countrate on y axes'

    num_events, colors = hist_use(etable)

    tc = ax_counts.bar(IDS, num_events, color = colors)

    ax_rate.set_ylabel('c/s')
    cntmin, cntmax = ax_counts.get_ylim()
    ax_rate.set_ylim((cntmin/etable.meta['EXPOSURE'],cntmax/etable.meta['EXPOSURE']))
    #countrate.set_ylim([np.min(rate)-20,np.max(rate)+20])

    ax_counts.set_xlabel('DET_ID')
    ax_counts.set_ylabel('# of Events')
    plot.locator_params(nticks = 20)
    plot.title('Total (Filtered) Event Count by Detector')
    #total_counts.set_ylim([np.min(num_events)-20, np.max(num_events)+20])

    return num_events

#----------------------THIS MAKES THE GRAYSCALE ID/EVENT COUNT CHART---------------------
def structure(etable, num_events):
    'Creates a grid where the xy pair corresponds to RAWX,RAWY and value at each entry is event count'
    rawx = np.zeros_like(IDS,dtype=np.int)
    rawy = np.zeros_like(IDS,dtype=np.int)

    #getting RAWX and RAWY vals for each ID
    for count, detid in enumerate(IDS):
        idx = np.where(etable['DET_ID']==detid)[0]
        if len(idx) > 0:
            rawx[count] = etable['RAWX'][idx][0]
            rawy[count] = etable['RAWY'][idx][0]
        else:
            print("No counts for det ",detid)
            rawx[count] = -1
            rawy[count] = -1

    #In STRUCTURE, each element corresponds to the geometric position
    # of each detector, while the value is the # of counts
    structure = np.zeros(shape = (7,8))

    for i in xrange(len(rawx)):
        if rawx[i] >= 0 and rawy[i] >= 0:
            structure[rawy[i]][rawx[i]] = num_events[i]

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
    plot.title('Filtered Event Count by Detector Location')
    plot.xlabel('Raw X')
    plot.ylabel('Raw Y')
    plot.colorbar(ax_img, orientation = 'horizontal')

    return

#----------------------THIS MAKES THE LIGHT CURVE---------------------------
def light_curve(etable, startmet, stopmet, binsize):
    'Bins events as a histogram to be plotted as the light curve. returns bins and the histogram'
    if startmet is None and stopmet is None:
        startmet = etable['MET'][0]
        t = etable['MET'] - startmet
        stopmet = etable['MET'][-1]
    else:
        t = etable['MET'][np.where(np.logical_and(etable['MET'] < stopmet, etable['MET'] > startmet))] - startmet

    duration = stopmet-startmet
    # Add 1 bin to make sure last bin covers last events
    bins = np.arange(0.0,duration+binsize,binsize)
    sums, edges = np.histogram(t, bins=bins, range=(0.0,duration))

    # Chop off last bin edge, which is only for computing histogram, not plotting
    return bins[:-1], sums

def plot_light_curve(etable, lclog, gtitable, binsize=1.0):
   #'Compute binned light curve of events and return mean rate,plots light curve'
    #EDGE CASE FOR FIRST INSTANCE
    bins, sums = light_curve(etable, gtitable['START'][0], gtitable['STOP'][0], binsize=binsize)
    cc = np.zeros_like(bins,dtype=np.float)
    cumtime = bins[-1]+binsize
	#THE REST OF THE GOOD INTERVALS
    for i in xrange(1,len(gtitable['START'])):
        mybins, mysums = light_curve(etable, gtitable['START'][i], gtitable['STOP'][i], binsize=binsize)
        bins = np.append(bins, mybins+cumtime)
        cumtime += mybins[-1]+binsize
        sums = np.append(sums, mysums)
        mycolors = np.zeros_like(mybins,dtype=np.float)+np.float(i)
        cc = np.append(cc,mycolors)

    #Compute mean rate
    rate = sums/binsize
    mean_rate = rate.mean()
    colornames = ['black','green','red','blue','magenta']
    colorlevels = np.arange(len(colornames))
    cmap, norm = mpl.colors.from_levels_and_colors(levels=colorlevels, colors=colornames, extend='max')
    lc = plot.scatter(bins, rate, c=np.fmod(cc,len(colornames)), cmap=cmap,norm=norm,marker='+', label='Light Curve')
    label = 'Mean Rate: {0:.3f} c/s'.format(mean_rate)
    # Plot line at mean counts per bin
    plot.axhline(y=mean_rate, xmin=bins[0], xmax=bins[-1], linestyle='dashed', label = label)
    #plot.legend(loc = 4)
    plot.ylabel('c/s')
    if lclog:
    	plot.yscale('log')
        plot.ylim(ymin=0.1)

    return mean_rate, lc

#-------------------------------THIS PLOTS THE FAST TO SLOW___------------------
def plot_slowfast(etable,args):
    'Scatter plot of PI and fast PHA, highlighting points above ratio cut'
    log.info('Counting slow and fast')
   # First do some counts
    nfastonly = np.count_nonzero(np.logical_and(etable['EVENT_FLAGS'][:,FLAG_FAST],
                                            np.logical_not(etable['EVENT_FLAGS'][:,FLAG_SLOW])))
    nslowonly = np.count_nonzero(np.logical_and(etable['EVENT_FLAGS'][:,FLAG_SLOW],
                                            np.logical_not(etable['EVENT_FLAGS'][:,FLAG_FAST])))
    nboth = np.count_nonzero(np.logical_and(etable['EVENT_FLAGS'][:,FLAG_SLOW],
                                            etable['EVENT_FLAGS'][:,FLAG_FAST]))
    log.info('Using only SLOW+FAST events for ratio plot')
    # Only compute ratio for events with both triggers
    etable = etable[np.logical_and(etable['EVENT_FLAGS'][:,FLAG_SLOW],etable['EVENT_FLAGS'][:,FLAG_FAST])]
    downsampfac = None
    if len(etable) > 50000:
        log.warning('Too many events for ratio plot. Plotting subset of points')
        downsampfac = len(etable)//50000
        etable = etable[::downsampfac]
    log.info('Computing ratio')
    # Ratio is SLOW to FAST. Edge events should have ratio bigger than cut
    ratio = np.array(etable['PHA'],dtype=np.float)/np.array(etable['PHA_FAST'],dtype=np.float)

    ratio_cut = args.filtratio
    colors = np.array(['k']*len(ratio))
    idx = np.where(ratio>ratio_cut)[0]
    colors[idx] = 'r'
    log.info('Plotting the points')
    plot.scatter(etable['PI']*PI_TO_KEV,ratio, s=.4, c = colors)
    x = np.arange(min(etable['PI']),max(etable['PI']))
    phax = np.ones_like(x)*ratio_cut

    if downsampfac is None:
        plot.title('PHA Slow to Fast Ratio vs Energy')
    else:
        plot.title('PHA Slow to Fast Ratio vs Energy (SUBSET)')
    plot.xlabel('Energy')
    plot.ylabel('PHA Ratio')
    plot.ylim([min(ratio)-.5,ratio_cut + 1.5])
    fast_str = "# of fast only  : " + str(nfastonly)
    slow_str = "# of slow only  : " + str(nslowonly)
    total =    "# of both       : " + str(nboth)
    bad =      "# of bad points : " + str(len(idx))
    plot.annotate(fast_str, xy=(0.03, 0.85), xycoords='axes fraction')
    plot.annotate(slow_str, xy=(0.03, 0.8), xycoords='axes fraction')
    plot.annotate(total, xy=(0.03, 0.75), xycoords='axes fraction')
    plot.annotate(bad, xy = (.03, .7), xycoords='axes fraction')
    plot.annotate("Ratio cut = {0:.2f}".format(ratio_cut),xy=(0.65,0.85),xycoords='axes fraction')
    plot.plot(x*PI_TO_KEV, phax, 'g--', linewidth = 0.5)
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

def plot_fft_of_power(etable,nyquist, pslog, writeps):
    'plots the power spectrum'

    dt = 0.5/nyquist
    METmin = etable['MET'].min()
    T = etable['MET'].max() - etable['MET'].min()

    # Choose good number of bins for efficient FFT
    n = choose_N(T/float(dt))
    bins = np.arange(n)*dt
    log.info('{0} {1}'.format(T/dt,n))
    log.info('Computing FFT with {0} bins of {1} s, covering {2} total time (Nyquist = {3})'.format(n,dt,T, nyquist))
    ts, edges = np.histogram(etable['MET']-METmin,bins)

    ft = np.fft.rfft(ts)
    power = (ft * ft.conj()).real
    power /= len(etable['MET'])
    power[0:50] = 0.0
    x = np.fft.rfftfreq(len(ts), dt)
    #idx = np.where(power>20)
    idx = np.argmax(power)
    print(x[idx], power[idx])
    if pslog:
        plot.semilogy(x,power)
    else:
        plot.plot(x,power)
    if writeps:
        data = np.array([x,power])
        data = data.T
        np.savetxt(file('powspec.txt','w'), data, fmt=['%f','%f'])

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

def pulse_profile(ax, etable, args):
    if (args.orb is None) or (args.par is None):
        log.warning('You did not specify orbfile or parfile')
        log.info('Please input files for orb and par with --orb and --par')
        return
    import pint
    import astropy.io.fits as pyfits
    import pint.toa, pint.models
    from pint.plot_utils import phaseogram_binned
    from pint.observatory.nicer_obs import NICERObs
    from pint.eventstats import hm

    ### Make arguments for parfile and orbfile and only do this if both are present
    log.info('Event file TELESCOPE = {0}, INSTRUMENT = {1}'.format(etable.meta['TELESCOP'],
                                                                  etable.meta['INSTRUME']))
    # Instantiate NICERObs once so it gets added to the observatory registry
    log.info('Setting up NICER observatory')
    NICERObs(name='NICER',FPorbname=args.orb,tt2tdb_mode='none')

    # Read event file and return list of TOA objects
    log.info('doing the load_toas thing')
    #tl  = load_NICER_TOAs(pulsefilename[0])

    # Create TOA list
    tl = []
    for t in etable['T']:
        tl.append(pint.toa.TOA(t, obs='NICER'))

    ts = pint.toa.TOAs(toalist=tl)
    ts.compute_TDBs()
    ts.compute_posvels(ephem='DE421',planets=True)

    log.info('Did all the stuff, now to PARFILE')
    # Load PINT model objects
    modelin = pint.models.get_model(args.par)
    log.info(str(modelin))

    # Compute phases
    phss = modelin.phase(ts.table)[1]
    # ensure all postive
    phases = np.where(phss < 0.0, phss + 1.0, phss)
    mjds = ts.get_mjds()

    h = hm(phases)
    ax.hist(phases, bins = 32)
    ax.text(0.1, 0.1, 'H = {0:.2f}'.format(h), transform=ax.transAxes)

    np.savetxt('{0}.phases'.format(args.basename),np.transpose([etable['MET'], etable['PI'],phases]))

    plot.ylabel('Counts')
    plot.xlabel('Pulse Phase')
    plot.title('Pulse Profile')
    return
#-------------------------OVERSHOOT RATE FOR RATIO------------------------------

def apply_gti(etable, gtitable):
    mets = etable['MET']
    idx = np.where(np.logical_and(mets>gtitable['START'][0],mets<gtitable['STOP'][0]))
    goodlist = [ etable[idx] ]
    for ii in range(1,len(gtitable['START'])):
        idx = np.where(np.logical_and(mets>gtitable['START'][ii],mets<gtitable['STOP'][ii]))
        goodlist.append(etable[idx])

    return vstack(goodlist,metadata_conflicts='silent')



def convert_to_elapsed_goodtime(mets, vals, gtitable):
    'Given a set of values at METs, extract the values during the GTIs and return times that are in elapsed good time for plotting'
    mets = np.asarray(mets)
    vals = np.asarray(vals)
    idx = np.where(np.logical_and(mets>gtitable['START'][0],mets<gtitable['STOP'][0]))
    goodvals = vals[idx]
    etimes = mets[idx] - gtitable['START'][0]
    cc = np.zeros_like(goodvals, dtype=np.float)
    for ii in range(1,len(gtitable['START'])):
        idx = np.where(np.logical_and(mets>gtitable['START'][ii],mets<gtitable['STOP'][ii]))
        goodvals = np.append(goodvals, vals[idx])
        etimes = np.append(etimes, mets[idx]-gtitable['START'][ii] + gtitable['CUMTIME'][ii])
        cc = np.append(cc, np.zeros_like(vals[idx],dtype=np.float)+np.float(ii))

    # Returns the arrays of elapsed times, values, and an array of what segment it is in, used for setting plot colors by GTI segment
    return etimes, goodvals, cc

def plot_overshoot(etable, overshootrate, gtitable, args, hkmet, bothrate):

    etime, overshoot, cc = convert_to_elapsed_goodtime(hkmet, overshootrate, gtitable)
    colornames = ['black','green','red','blue','magenta']
    colorlevels = np.arange(len(colornames))
    cmap, norm = mpl.colors.from_levels_and_colors(levels=colorlevels, colors=colornames, extend='max')

    plot.scatter(etime, overshoot, c=np.fmod(cc,len(colornames)), cmap=cmap, norm=norm, marker='+')
    if bothrate is not None:
        etime, both, cc = convert_to_elapsed_goodtime(hkmet, bothrate, gtitable)
        plot.scatter(etime, both, color = 'c', marker='.', label='Both Under and Over Flags')
        plot.legend(loc = 2)
    plot.ylabel('Overshoot rate')
    plot.grid(True)

    if args.lclog:
        plot.yscale('log')
        plot.ylim(ymin=1.0)
    return

def plot_SAA(mktable, gtitable, overshootrate):
    time, insaa, colors = convert_to_elapsed_goodtime(mktable['TIME'], mktable['SAA'], gtitable)
    time = np.delete(time, np.where(insaa == 0))
    insaa = np.delete(insaa, np.where(insaa == 0))
    colornames = ['black','green','red','blue','magenta']
    colorlevels = np.arange(len(colornames))
    cmap, norm = mpl.colors.from_levels_and_colors(levels=colorlevels, colors=colornames, extend='max')
    insaa[np.where(insaa == 1)] = max(overshootrate) + 100
    plot.scatter(time, insaa, color = 'y', label = 'In the SAA',marker = '_')
    plot.legend(loc = 2)
    return

#-------------------------UNDERSHOOT RATE FOR RATIO------------------------------
def plot_undershoot(etable, undershootrate, gtitable, args, hkmet, mktable):

    etime, undershoot, cc = convert_to_elapsed_goodtime(hkmet, undershootrate, gtitable)

    colornames = ['black','green','red','blue','magenta']
    colorlevels = np.arange(len(colornames))
    cmap, norm = mpl.colors.from_levels_and_colors(levels=colorlevels, colors=colornames, extend='max')

    plot.scatter(etime, undershoot, c=np.fmod(cc,len(colornames)), cmap=cmap, norm=norm, marker='+')

    # Add sunshine
    sunmet = mktable['TIME']
    sunshine = mktable['SUNSHINE']
    sunt, suny, suncc = convert_to_elapsed_goodtime(sunmet, sunshine, gtitable)
    sidx = np.where(suny==0)
    # Delete the 0 values so they don't plot
    sunt=np.delete(sunt,sidx)
    suny=np.delete(suny,sidx)
    plot.scatter(sunt,suny*undershoot.max(), color = 'y', label = 'Sunshine', marker = '_')
    plot.legend(loc = 2)

    plot.ylabel('Undershoot rate')
    if args.lclog:
        plot.yscale('log')

    return


#-------------------------SUN / EARTH / MOON ANGLES-----------------------------
def plot_angles(mktable, gtitable):
    sun = mktable['SUN_ANGLE']
    earth = mktable['BR_EARTH']
    moon = mktable['MOON_ANGLE']
    elv = mktable['ELV']
    met = mktable['TIME']


    goodtime, sunangle, cc = convert_to_elapsed_goodtime(met, sun, gtitable)
    goodtime, earthangle, cc = convert_to_elapsed_goodtime(met, earth, gtitable)
    goodtime, moonangle, cc = convert_to_elapsed_goodtime(met, moon, gtitable)
    goodtime, elvangle, cc = convert_to_elapsed_goodtime(met, elv, gtitable)

    plot.scatter(goodtime, sunangle, marker = '.', color = 'y', alpha=0.5, label = 'Sun')
    plot.scatter(goodtime, earthangle, marker ='.', color = 'b', alpha=0.5, label = 'Bright Earth')
    plot.scatter(goodtime, moonangle, marker = '.', color = 'grey', alpha = 0.5, label = 'Moon')
    plot.scatter(goodtime, elvangle, marker = '.', color = 'm', alpha = 0.5, label = 'ELV')
    plot.legend(loc = 2)
    plot.ylim((0.0,180.0))
    plot.grid(True)
    plot.yticks([0.0, 45.0, 90.0, 135.0, 180.0])
    plot.ylabel('Angle (deg)')

    return

#--------------------------POINTING---------------------------------------------
def plot_pointing(mktable, gtitable):
    time, pointing, colors = convert_to_elapsed_goodtime(mktable['TIME'], mktable['ANG_DIST'], gtitable)

    colornames = ['black','green','red','blue','magenta']
    colorlevels = np.arange(len(colornames))
    cmap, norm = mpl.colors.from_levels_and_colors(levels=colorlevels, colors=colornames, extend='max')

    plot.scatter(time, pointing, c = colors, cmap = cmap,marker = '+', label='Pointing Offset')
    plot.legend(loc = 2)
    plot.ylabel('Angle (deg)')
    plot.yscale('log')
    plot.axhline(1.0/60.0,c='r')
    plot.ylim((0.0001,100.0))

    return
#--------------------------LAT / LON---------------------------------------------
def plot_latlon(mktable, gtitable):
    time, lat, colors = convert_to_elapsed_goodtime(mktable['TIME'], mktable['SAT_LAT'], gtitable)
    time, lon, colors = convert_to_elapsed_goodtime(mktable['TIME'], mktable['SAT_LON'], gtitable)

    colornames = ['black','green','red','blue','magenta']
    colorlevels = np.arange(len(colornames))
    cmap, norm = mpl.colors.from_levels_and_colors(levels=colorlevels, colors=colornames, extend='max')

    plot.scatter(time, lat, c = colors, cmap = cmap,marker = '^', label='Latitude')
    #plot.scatter(time, lon, c = colors, cmap = cmap, marker = '_', label = 'Longitude')
    plot.legend(loc = 2)
    plot.ylim((-60.0,60.0))
    plot.xlabel('Elapsed Time (s)', labelpad = 1)
    plot.grid(True)
    plot.ylabel('Degrees')
    return
#-------------------------THIS PLOTS USEFUL TEXT AT THE TOP OF THE SUPLOT-------
def calc_nresets(etable, IDS):
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

#----------------------Filter Ratio Cut-----------------------------------------
def filt_ratio(etable, ratiocut):
    #Filters out the points > filtratio
    ratio = np.zeros_like(etable['PI'],dtype=np.float)
    idx = np.where(np.logical_and(etable['PHA']>0, etable['PHA_FAST']>0))[0]
    ratio[idx] = np.asarray(etable['PHA'][idx],dtype=np.float)/np.asarray(etable['PHA_FAST'][idx],dtype=np.float)
    etable = etable[np.where(ratio < ratiocut)[0]]
    return etable
