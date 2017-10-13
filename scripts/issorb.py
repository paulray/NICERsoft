#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
import astropy.units as u
from astropy.time import Time
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, get_sun, get_moon, ICRS, Angle
from astropy.coordinates.name_resolve import get_icrs_coordinates
from pyorbital import tlefile
import argparse

SunAvoidance = 45.0*u.deg
MoonAvoidance = 15.0*u.deg

tle160lines = ['1 25544U 98067A   17160.91338884 +.00001442 +00000-0 +29152-4 0  9993',
        '2 25544 051.6425 074.5823 0004493 253.3640 193.9362 15.54003243060621']

tle167lines = ['1 25544U 98067A   17167.53403196 +.00002711 +00000-0 +48329-4 0  9994',
        '2 25544 051.6431 041.5846 0004445 283.0899 147.8207 15.54043876061656']

platform = 'ISS (ZARYA)'

parser = argparse.ArgumentParser(description="Compute ISS orbit geometry vs time")
parser.add_argument("sourcename",help="Source name to look up coordinates")
parser.add_argument("-ra ",  dest='ra', help="Source RA  [default: lookup]", default=None)
parser.add_argument("-dec", dest='dec',help="Source DEC [default: lookup]", default=None)
parser.add_argument("--vis", help="Visibility file from Wayne", default=None)
args = parser.parse_args()

# Set up two TLEs so we can compute the precession rate from the
# change of RA of ascending node with time
tle1 = tlefile.read(platform,line1=tle160lines[0],line2=tle160lines[1])
tle2 = tlefile.read(platform,line1=tle167lines[0],line2=tle167lines[1])

print(platform)
ISSInclination = tle2.inclination*u.deg
print("Inclination = {0:.3f}".format(ISSInclination))

StarboardPoleDec0 = -1*(90.0*u.deg-ISSInclination)
PortPoleDec0 = -1*StarboardPoleDec0

print("Starboard Orbit Pole Declination {0:.2f}".format(StarboardPoleDec0))
print("Port Orbit Pole Declination {0:.2f}".format(PortPoleDec0))

# Compute ISS precession rate in degrees per day
ISSPrecessionRate = (tle2.right_ascension-tle1.right_ascension)/(tle2.epoch_day-tle1.epoch_day)*u.deg/u.d
print("ISS Precession Rate {0:.3f} ({1:.3f} period)".format(ISSPrecessionRate,
    np.abs(360.0*u.deg/ISSPrecessionRate)))

ttle2 = Time("{0:4d}-01-01T00:00:00".format(int(tle2.epoch_year)+2000),
    format='isot',scale='utc') + tle2.epoch_day*u.d
print("ttle2 = ",ttle2.isot)
StarboardPoleRA0 = np.fmod(tle2.right_ascension + 90.0,360.0)*u.deg
PortPoleRA0 = np.fmod(StarboardPoleRA0 + 180.0*u.deg, 360.0*u.deg)

print("Starboard Pole RA @ ttle2 = {0:.3f}".format(StarboardPoleRA0))
print("Port Pole RA @ ttle2 = {0:.3f}".format(PortPoleRA0))

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
doy_now = np.float(now.yday.split(':')[1])
print("Current DOY = {0}".format(np.int(doy_now)))
print("StarboardPoleRA (now) = {0:.3f}".format(StarboardPoleRA(now)))
print("PortPoleRA (now) = {0:.3f}".format(PortPoleRA(now)))

if args.ra is None or args.dec is None:
    SourcePos = get_icrs_coordinates(args.sourcename)
else:
    SourcePos = ICRS(ra=Angle(args.ra),dec=Angle(args.dec))
print("\nSource: {0} at {1}, {2}".format(args.sourcename,SourcePos.ra, SourcePos.dec))
print("Separation from Starboard Pole = {0:.3f}".format(SourcePos.separation(StarboardPoleCoord(now))))
print("Separation from Port Pole = {0:.3f}".format(SourcePos.separation(PortPoleCoord(now))))

# Plot separation from orbit pole for all of 2017
doy2017 = np.arange(365.0)
times = doy2017*u.d + Time("2017-01-01T00:00:00",format='isot',scale='utc')

fig, ax = plt.subplots()
seps = SourcePos.separation(StarboardPoleCoord(times))
ax.plot(doy2017, seps.to(u.deg))

# The Sun and Moon positions are returned in the GCRS frame
# Convert them to ICRS so .separation doesn't go insane when
# comparing different frames with different obstimes.
SunPos = get_sun(times)
SunPos = SkyCoord(SunPos.ra, SunPos.dec, frame='icrs')
MoonPos = get_moon(times)
MoonPos = SkyCoord(MoonPos.ra, MoonPos.dec, frame='icrs')
# Plot dots for Sun and Moon separation when they violate the constraint
sunseps = SourcePos.separation(SunPos).to(u.deg)
idx = np.where(sunseps<SunAvoidance)[0]
ax.plot(doy2017[idx], sunseps[idx],'o',c='y',label='Sun')
moonseps = SourcePos.separation(MoonPos).to(u.deg)
idx = np.where(moonseps<MoonAvoidance)[0]
ax.plot(doy2017[idx], moonseps[idx],'o',c='k',label='Moon')

# Plot orange dot for now
ax.plot(doy_now,SourcePos.separation(StarboardPoleCoord(now)),'o',label='Today')

if args.vis is not None:
    doy, fr = np.loadtxt(args.vis,unpack=True)
ax.plot(doy,fr*100)
ax.set_title(args.sourcename)
ax.legend()
ax.set_xlabel('DOY 2017')
ax.set_ylim((0.0,180.0))
ax.grid(True)
ax.set_ylabel('Angle from Starboard Pole (deg)')
plt.show()
