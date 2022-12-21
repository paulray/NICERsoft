#!/usr/bin/env python
import numpy as np
import astropy.units as u
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.ticker import FixedLocator
from astropy.coordinates import SkyCoord, get_sun, get_moon, ICRS, Angle
from astropy.coordinates.name_resolve import get_icrs_coordinates
from pyorbital import tlefile
import argparse

#issorb.py
def PulsarVis(sourcename,times):
    
    if sourcename in {'PSRJ0030+0451','PSRJ0437-4715'}:
    	SunAvoidance = 55.0*u.deg
    else:
        SunAvoidance = 45.0*u.deg
    print("Sun Avoidance: {0:.3f}".format(SunAvoidance))	
    MoonAvoidance = 15.0*u.deg
    # Each pulsar's feasible observation range based on ISS hardware data
    if sourcename in 'PSRB1821-24': #Analysis ticket 890 
        HardwareUpperAvoidance = 120*u.deg #360*u.deg 
	HardwareLowerAvoidance = 88*u.deg #0*u.deg
    elif sourcename in 'PSRB1937+21':
        HardwareUpperAvoidance = 142.0*u.deg #360*u.deg 
	HardwareLowerAvoidance = 91.0*u.deg #0*u.deg
    elif sourcename in 'PSRJ0218+4232':
        HardwareUpperAvoidance = 156.0*u.deg #360*u.deg 
	HardwareLowerAvoidance = 88.0*u.deg #0*u.deg
    elif sourcename in 'PSRJ0030+0451':
        HardwareUpperAvoidance = 148.0*u.deg #360*u.deg 
	HardwareLowerAvoidance = 95.0*u.deg #0*u.deg
    elif sourcename in 'PSRJ0437-4715':
	HardwareUpperAvoidance = 95.0*u.deg #360*u.deg 
	HardwareLowerAvoidance = 35.0*u.deg #0*u.deg 
    else:	     
    	HardwareUpperAvoidance = 360.0*u.deg #360*u.deg 
    	HardwareLowerAvoidance = 0.0*u.deg #0*u.deg  
    
    print("ISS Upper:",HardwareUpperAvoidance,"Lower:",HardwareLowerAvoidance,"\n".format(HardwareLowerAvoidance,HardwareUpperAvoidance))

    tle160lines = ['1 25544U 98067A   17160.91338884 +.00001442 +00000-0 +29152-4 0  9993',
                   '2 25544 051.6425 074.5823 0004493 253.3640 193.9362 15.54003243060621']
        
    tle167lines = ['1 25544U 98067A   17167.53403196 +.00002711 +00000-0 +48329-4 0  9994',
                  '2 25544 051.6431 041.5846 0004445 283.0899 147.8207 15.54043876061656']

    platform = 'ISS (ZARYA)'

    # Set up two TLEs so we can compute the precession rate from the
    # change of RA of ascending node with time
    tle1 = tlefile.read(platform,line1=tle160lines[0],line2=tle160lines[1])
    tle2 = tlefile.read(platform,line1=tle167lines[0],line2=tle167lines[1])

#    print(platform)
    ISSInclination = tle2.inclination*u.deg
#    print("Inclination = {0:.3f}".format(ISSInclination))

    StarboardPoleDec0 = -1*(90.0*u.deg-ISSInclination)
    PortPoleDec0 = -1*StarboardPoleDec0

#print("Starboard Orbit Pole Declination {0:.2f}".format(StarboardPoleDec0))
#    print("Port Orbit Pole Declination {0:.2f}".format(PortPoleDec0))

    # Compute ISS precession rate in degrees per day
    ISSPrecessionRate = (tle2.right_ascension-tle1.right_ascension)/(tle2.epoch_day-tle1.epoch_day)*u.deg/u.d
#    print("ISS Precession Rate {0:.3f} ({1:.3f} period)".format(ISSPrecessionRate,
#                                                               np.abs(360.0*u.deg/ISSPrecessionRate)))

    ttle2 = Time("{0:4d}-01-01T00:00:00".format(int(tle2.epoch_year)+2000),
                format='isot',scale='utc') + tle2.epoch_day*u.d
#    print("ttle2 = ",ttle2.isot)
    StarboardPoleRA0 = np.fmod(tle2.right_ascension + 90.0,360.0)*u.deg
    PortPoleRA0 = np.fmod(StarboardPoleRA0 + 180.0*u.deg, 360.0*u.deg)

#   print("Starboard Pole RA @ ttle2 = {0:.3f}".format(StarboardPoleRA0))
#   print("Port Pole RA @ ttle2 = {0:.3f}".format(PortPoleRA0))

    def StarboardPoleDec(t):
       return np.ones_like(t)*StarboardPoleDec0

    def PortPoleDec(t):
        return np.ones_like(t)*StarboardPoleDec0
        
    def StarboardPoleRA(t):
            return np.fmod(StarboardPoleRA0 + (t-ttle2).to(u.d)*ISSPrecessionRate, 360.0*u.deg)
        
    def PortPoleRA(t):
            return np.fmod(StarboardPoleRA(t)+180.0*u.deg, 360.0*u.deg)
        
    def StarboardPoleCoord(t):
            return SkyCoord(StarboardPoleRA(t).value,StarboardPoleDec(t).value,unit=u.deg,frame="icrs")
        
    def PortPoleCoord(t):
            return SkyCoord(PortPoleRA(t).value,PortPoleDec(t).value,unit=u.deg,frame="icrs")
        
    now = Time.now()
    doy_now = float(now.yday.split(':')[1])
    #   print("Current DOY = {0}".format(np.int(doy_now)))
    #print("StarboardPoleRA (now) = {0:.3f}".format(StarboardPoleRA(now)))
    #print("PortPoleRA (now) = {0:.3f}".format(PortPoleRA(now)))

    if sourcename[0:2] != 'PSR':
    	SourcePos = get_icrs_coordinates(sourcename)
    else:
    	splitstr  = sourcename.split(',')
    	SourcePos = ICRS(ra=Angle(double(splitstr[0])),dec=Angle(double(splitstr[1])))
    #print("\nSource: {0} at {1}, {2}".format(sourcename,SourcePos.ra, SourcePos.dec))
    #print("Separation from Starboard Pole = {0:.3f}".format(SourcePos.separation(StarboardPoleCoord(now))))
    #print("Separation from Port Pole = {0:.3f}".format(SourcePos.separation(PortPoleCoord(now))))

    #Calculate terms and references to do check
    #doy2XXX = np.arange(365.0)
    #times = doy2XXX*u.d + startepoch

    issseps = SourcePos.separation(StarboardPoleCoord(times)).to(u.deg)

    # The Sun and Moon positions are returned in the GCRS frame
    # Convert them to ICRS so .separation doesn't go insane when
    # comparing different frames with different obstimes.
    SunPos = get_sun(times)
    SunPos = SkyCoord(SunPos.ra, SunPos.dec, frame='icrs')
    MoonPos = get_moon(times)
    MoonPos = SkyCoord(MoonPos.ra, MoonPos.dec, frame='icrs')
    # Hold indicies when Sun and Moon cause constraint violation
    sunseps = SourcePos.separation(SunPos).to(u.deg)
    idxsun = sunseps<SunAvoidance
    moonseps = SourcePos.separation(MoonPos).to(u.deg)
    idxmoon = moonseps<MoonAvoidance
    # Hold indicies when sep is outside of 91 to 142 deg due to hardware violation (KWood analysis 10/27/2017)
    idxang = ~((issseps>HardwareLowerAvoidance) & (issseps<HardwareUpperAvoidance))
    # Apply all vis constraints
    indxall = idxsun | idxmoon | idxang #true when constraint is violated
    # Generate and populate signal output
    vis = np.ones(len(doy2XXX))
    vis[indxall] = float('NaN')

    #Return result
    return vis

if __name__ == '__main__':
#Do the same thing as issorb.py, but for a list of targets
    parser = argparse.ArgumentParser(description="Forecast pulsar visibility for a list of targets")
    parser.add_argument("-source", dest='sourcefile', help="File with Source names. Will default to SEXTANT targets")
    parser.add_argument("-doy_in", dest='doy_in', help="Start doy epoch",default=1)
    parser.add_argument("-year_in",dest='year_in', help="Start year epoch",default=2018)
    args = parser.parse_args()

    if args.sourcefile is not None:
        sourcelist = [line.rstrip('\n') for line in open(args.sourcefile)]
    else:
        sourcelist = ['PSRB1821-24','PSRB1937+21','PSRJ0030+0451','PSRJ0218+4232','PSRJ0437-4715','PSRJ0534+2200']

#Plot each vis as it is generated for each pulsar
    fig, ax = plt.subplots()
    doy2XXX = np.arange(365.0)
    #print("{:4}:{:03d}".format(args.year_in,args.doy_in))
    startepoch = Time('{:4}:{:03d}'.format(int(args.year_in),int(args.doy_in)),format='yday',scale='utc')
    times = doy2XXX*u.d + startepoch
    row = 0
    for x in sourcelist:
        print("Processing {:20s}".format(sourcelist[row]))
        vis = PulsarVis(x,times) #assumes Nx1 array from PulsarVis(x)
        ax.plot(doy2XXX,vis*row,linewidth=5)
        row += 1

    #Formatting time labels
    monthnames = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','']
    monthstartdoy = [0,31,59,90,120,151,181,212,243,273,304,334,364]

    ax.set_title('Ideal Visibility Forecast for the SEXTANT Experiment')
    ax.set_xlabel('Elasped Time from {}\n'.format(startepoch))
    ax.set_xticks(monthstartdoy)
    ax.set_xticklabels(monthnames)
    ax.xaxis.set_minor_locator(FixedLocator(np.arange(0,365,1).tolist()))
    ax.set_yticks(range(len(sourcelist)))
    ax.set_yticklabels(sourcelist)
    ax.grid(True)
    plt.show()

    print("...done")
