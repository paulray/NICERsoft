#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os.path as path
import argparse
import astropy.units as u
from astropy.time import Time
from astropy.table import Table
import astropy.io.fits as pyfits
from astropy import log
import sys
from scipy.interpolate import InterpolatedUnivariateSpline

from mpl_toolkits.basemap import Basemap

from nicer.mcc import MCC
from nicer.values import *


parser = argparse.ArgumentParser(description = "Plot HKP data on map")
parser.add_argument("mpuhkfiles", help="Input file names (should end in mpu?.hk)", nargs='+')
parser.add_argument("--mcc", help="Specify MCC file to use for positions", default=None)
parser.add_argument("--sps", help="Specify SPS file (_apid0260.hk) to use for positions", default=None)
args = parser.parse_args()

log.info('Reading '+args.mpuhkfiles[0])
hdulist = pyfits.open(args.mpuhkfiles[0])
td = hdulist[1].data
met = td['TIME']
log.info("MET Range {0} to {1}".format(met.min(),met.max()))
t = MET0+met*u.s
overshootrate = td['MPU_OVER_COUNT'].sum(axis=1)

for fn in args.mpuhkfiles[1:]:
    log.info('Reading '+fn)
    hdulist = pyfits.open(fn)
    mytd = hdulist[1].data
    mymet = td['TIME']
    myt = MET0+met*u.s
    myovershootrate = td['MPU_OVER_COUNT'].sum(axis=1)
    if not np.all(mymet == met):
        log.error('TIME axes are not compatible')
        sys.exit(1)
    overshootrate += myovershootrate

if args.mcc is not None:
    eph = MCC(args.mcc)
    lat, lon = eph.latlon(met)
elif args.sps is not None:
    spshdu = pyfits.open(args.sps)
    spsmet = spshdu[1].data['TIME']
    spslat = np.rad2deg(spshdu[1].data['GPS_SPS_LAT'])
    spslon = np.rad2deg(spshdu[1].data['GPS_SPS_LON'])
    log.info('SPS Range {0} to {1}'.format(spsmet.min(),spsmet.max()))
    latinterp = InterpolatedUnivariateSpline(spsmet,spslat,ext='extrapolate')
    # WARNING: Interpolating longitude is a bad idea since it has discontinuities!
    loninterp = InterpolatedUnivariateSpline(spsmet,spslon,ext='extrapolate')
    # Check that we aren't trying to extrapolate too far
    if (met.min() < spsmet.min()-10.0) or (met.max() > spsmet.max()+10):
        log.error('Trying to extrapolate too far!')
    lat = latinterp(met)
    lon = loninterp(met)
else:
    log.error('Must specify --sps or --mcc')
    sys.exit(2)

saa_lon, saa_lat = np.loadtxt(path.join(datadir,'saa_lonlat.txt'),unpack=True)

fig, ax = plt.subplots(figsize=(16,9))
map = Basemap(projection='cyl', resolution = 'l',
              lat_0=0, lon_0=0)
map.drawcoastlines()
map.scatter(lon, lat,c=overshootrate,cmap='jet')
map.plot(saa_lon,saa_lat,'g',lw=2)
plt.show()
