#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
import astropy.units as u
from astropy.time import Time
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.coordinates.name_resolve import get_icrs_coordinates
from pyorbital import tlefile

tle160lines = ['1 25544U 98067A   17160.91338884 +.00001442 +00000-0 +29152-4 0  9993',
        '2 25544 051.6425 074.5823 0004493 253.3640 193.9362 15.54003243060621']

tle167lines = ['1 25544U 98067A   17167.53403196 +.00002711 +00000-0 +48329-4 0  9994',
        '2 25544 051.6431 041.5846 0004445 283.0899 147.8207 15.54043876061656']

platform = 'ISS (ZARYA)'

tle160 = tlefile.read(platform,line1=tle160lines[0],line2=tle160lines[1])
tle167 = tlefile.read(platform,line1=tle167lines[0],line2=tle167lines[1])

# Most recent TLE
tle = tle167

print(platform)
ISSInclination = tle.inclination*u.deg
print("Inclination = {0:.3f}".format(ISSInclination))

StarboardPoleDec0 = -1*(90.0*u.deg-ISSInclination)
PortPoleDec0 = -1*StarboardPoleDec0

print("Starboard Orbit Pole Declination {0:.2f}".format(StarboardPoleDec0))
print("Port Orbit Pole Declination {0:.2f}".format(PortPoleDec0))

# Compute ISS precession rate in degrees per day
ISSPrecessionRate = (tle167.right_ascension-tle160.right_ascension)/(tle167.epoch_day-tle160.epoch_day)*u.deg/u.d
print("ISS Precession Rate {0:.3f} ({1:.3f} period)".format(ISSPrecessionRate,
    np.abs(360.0*u.deg/ISSPrecessionRate)))

t0 = Time("{0:4d}-01-01T00:00:00".format(int(tle167.epoch_year)+2000),
    format='isot',scale='utc') + tle167.epoch_day*u.d
print("t0 = ",t0.isot)
StarboardPoleRA0 = np.fmod(tle167.right_ascension + 90.0,360.0)*u.deg
PortPoleRA0 = np.fmod(StarboardPoleRA0 + 180.0*u.deg, 360.0*u.deg)

print("Starboard Pole RA @ t0 = {0:.3f}".format(StarboardPoleRA0))
print("Port Pole RA @ t0 = {0:.3f}".format(PortPoleRA0))

def StarboardPoleDec(t):
    return np.ones_like(t)*StarboardPoleDec0

def PortPoleDec(t):
    return np.ones_like(t)*StarboardPoleDec0

def StarboardPoleRA(t):
    return np.fmod(StarboardPoleRA0 + (t-t0).to(u.d)*ISSPrecessionRate, 360.0*u.deg)

def PortPoleRA(t):
    return np.fmod(StarboardPoleRA(t)+180.0*u.deg, 360.0*u.deg)

def StarboardPoleCoord(t):
    return SkyCoord(StarboardPoleRA(t).value,StarboardPoleDec(t).value,unit=u.deg,frame="icrs")

def PortPoleCoord(t):
    return SkyCoord(PortPoleRA(t).value,PortPoleDec(t).value,unit=u.deg,frame="icrs")

now = Time.now()

print("StarboardPoleRA (now) = {0:.3f}".format(StarboardPoleRA(now)))
print("PortPoleRA (now) = {0:.3f}".format(PortPoleRA(now)))

PosJ0437 = get_icrs_coordinates('PSRJ0437-4715')
print("Separation from Starboard Pole = {0:.3f}".format(PosJ0437.separation(StarboardPoleCoord(now))))
print("Separation from Port Pole = {0:.3f}".format(PosJ0437.separation(PortPoleCoord(now))))

# Plot separation from orbit pole for all of 2017
times = np.arange(360)*u.d + t0
fig, ax = plt.subplots()
seps = PosJ0437.separation(StarboardPoleCoord(times))
ax.plot((times-t0).to(u.d), seps.to(u.deg))
ax.set_title('PSR J0437-4715')
ax.set_xlabel('DOY 2017')
ax.set_ylabel('Angle from Starboard Pole (deg)')
plt.show()
