#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
import os.path as path
import argparse
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy import log
from astropy.coordinates import SkyCoord
import sys
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib.colors import Normalize, LogNorm
from mpl_toolkits.basemap import Basemap

from nicer.mcc import MCC
from nicer.values import *
import nicer.yday_custom

header = """                                                                         18 Oct 2017 12:33:47
Satellite-NICER-To-AreaTarget-SAA:  Access Summary Report


NICER-To-SAA
------------
                  Access     Start Time (ISO-YD)       Stop Time (ISO-YD)     Duration (sec)
                  ------    ---------------------    ---------------------    --------------
"""

parser = argparse.ArgumentParser(description = "Compute SAA and polar horn entry/exit times.")
parser.add_argument("mcc", help="Name of MCC ephemeris file to use")
args = parser.parse_args()

log.info('Reading MCC ephemeris')
eph = MCC(args.mcc)
dur_s = (eph.t[-1]-eph.t[0]).to(u.s).value

# Make time spacing be 10 seconds over range of ephemeris
myt = eph.t[0] + np.arange(0.0,dur_s,10.0)*u.s
mymet = (myt-MET0).to(u.s).value

#llpoints =  zip(eph.lon.value,eph.lat.value)
mylat,mylon  = eph.latlon(mymet)
llpoints = zip(mylon.value,mylat.value)

log.info('Getting SAA data')
saa_poly = np.loadtxt(path.join(datadir,'saa_lonlat.txt'))
saa_path = mplPath.Path(saa_poly,closed=True)
in_saa = saa_path.contains_points(llpoints)

nph_poly = np.loadtxt(path.join(datadir,'nph_lonlat.txt'))
nph_path = mplPath.Path(nph_poly)
in_nph = nph_path.contains_points(llpoints)

neph_poly = np.loadtxt(path.join(datadir,'neph_lonlat.txt'))
neph_path = mplPath.Path(neph_poly)
in_neph = neph_path.contains_points(llpoints)

sph_poly = np.loadtxt(path.join(datadir,'sph_lonlat.txt'))
sph_path = mplPath.Path(sph_poly)
in_sph = sph_path.contains_points(llpoints)

nicer_poly = np.loadtxt(path.join(datadir,'nicer_saa.txt'))
nicer_path = mplPath.Path(nicer_poly)
in_nicer = nicer_path.contains_points(llpoints)

in_any = in_saa | in_nph | in_neph | in_sph | in_nicer

tstart = None
naccess=1
print(header)
prev = in_any[0]
for i in range(len(in_any)):
    if in_any[i] != prev:
        if in_any[i]:
            tstart = myt[i]
        else:
            tstop = myt[i]
            if tstart is not None:
                print("                  {0:6d}     {1}     {2}  {3:14.3f}".format(
                naccess,tstart.yday_custom,tstop.yday_custom,(tstop-tstart).to(u.s).value))
                naccess += 1
    prev = in_any[i]

map = Basemap(projection='cyl', resolution = 'l',  llcrnrlon=-180, llcrnrlat=-61,
    urcrnrlon=180, urcrnrlat=61, lat_0 = 0, lon_0 = 0)
map.drawcoastlines()
map.scatter(mylon[in_any],mylat[in_any],s=1.0)
plt.show()
