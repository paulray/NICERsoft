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
parser.add_argument("obsdirs", help="Input directories", nargs='+')
args = parser.parse_args()

saa_lon, saa_lat = np.loadtxt(path.join(datadir,'saa_lonlat.txt'),unpack=True)
nph_lon, nph_lat = np.loadtxt(path.join(datadir,'nph_lonlat.txt'),unpack=True)
sph_lon, sph_lat = np.loadtxt(path.join(datadir,'sph_lonlat.txt'),unpack=True)
fig, ax = plt.subplots(figsize=(16,9))
map = Basemap(projection='cyl', resolution = 'l',
              lat_0=0, lon_0=0)
map.drawcoastlines()
map.plot(saa_lon,saa_lat,'r',lw=2)
map.plot(nph_lon,nph_lat,color='orange',marker='o',markersize=10.0,linestyle='-')
map.plot(sph_lon,sph_lat,'orange',marker='o',markersize=10.0,linestyle='-')

for obsdir in args.obsdirs:
    log.info('Processing '+obsdir)
    mpuhkfiles = glob(path.join(obsdir,'xti/hk/*mpu*.hk'))
    hdulist = pyfits.open(mpuhkfiles[0])
    td = hdulist[1].data
    met = td['TIME']
    log.info("MET Range {0} to {1}".format(met.min(),met.max()))
    t = MET0+met*u.s
    overshootrate = td['MPU_OVER_COUNT'].sum(axis=1)

    for fn in mpuhkfiles[1:]:
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

    spsname = glob(path.join(obsdir,'auxil/*apid0260.hk'))[0]
    log.info('Reading SPS file {0}'.format(spsname))
    eph = SPS(spsname)
    lat, lon = eph.latlon(met)

    sc = map.scatter(lon, lat, c=overshootrate,
        norm=LogNorm(vmin=10.0,vmax=1000.0),cmap='jet')

cbar = map.colorbar(sc, location='bottom',pad='5%')
cbar.set_label('Overshoot Rate')
plt.show()
