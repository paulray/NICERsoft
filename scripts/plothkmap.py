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
from matplotlib.colors import Normalize, LogNorm
from mpl_toolkits.basemap import Basemap

from nicer.mcc import MCC
from nicer.sps import SPS
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
    eph = SPS(args.sps)
    lat, lon = eph.latlon(met)
else:
    log.error('Must specify --sps or --mcc')
    sys.exit(2)

saa_lon, saa_lat = np.loadtxt(path.join(datadir,'saa_lonlat.txt'),unpack=True)
nph_lon, nph_lat = np.loadtxt(path.join(datadir,'nph_lonlat.txt'),unpack=True)
neph_lon, neph_lat = np.loadtxt(path.join(datadir,'neph_lonlat.txt'),unpack=True)
sph_lon, sph_lat = np.loadtxt(path.join(datadir,'sph_lonlat.txt'),unpack=True)

fig, ax = plt.subplots(figsize=(16,9))
map = Basemap(projection='cyl', resolution = 'l',
              lat_0=0, lon_0=0)
map.drawcoastlines()
sc = map.scatter(lon, lat,c=overshootrate,norm=LogNorm(vmin=10.0,vmax=1000.0),cmap='jet')
map.plot(saa_lon,saa_lat,'r',lw=2)
map.plot(nph_lon,nph_lat,color='orange',marker='o',markersize=10.0,linestyle='-')
map.plot(neph_lon,neph_lat,color='orange',marker='o',markersize=10.0,linestyle='-')
map.plot(sph_lon,sph_lat,'orange',marker='o',markersize=10.0,linestyle='-')
cbar = map.colorbar(sc, location='bottom',pad='5%')
cbar.set_label('Overshoot Rate')
plt.show()
