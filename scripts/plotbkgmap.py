#!/usr/bin/env python
from __future__ import print_function, division
import argparse
import numpy as np
import matplotlib.pyplot as plot
from astropy import log
import astropy.units as u
from astropy.time import Time
from astropy.table import Table
from os import path
from matplotlib.colors import Normalize, LogNorm
from cartopy.crs import PlateCarree

from nicer.values import *
from nicer.mcc import MCC
#from nicer.sps import SPS
from nicer.latloninterp import LatLonInterp

parser = argparse.ArgumentParser(description="Plot background info on a map")
parser.add_argument("bkffiles", help="Name of bkf files to process", nargs='+')
parser.add_argument("--column", help="Which bkf column to plot", default="EV_OVERSHOOT")
args = parser.parse_args()

log.info('Getting SAA data')
saa_lon, saa_lat = np.loadtxt(path.join(datadir,'saa_lonlat.txt'),unpack=True)
nph_lon, nph_lat = np.loadtxt(path.join(datadir,'nph_lonlat.txt'),unpack=True)
neph_lon, neph_lat = np.loadtxt(path.join(datadir,'neph_lonlat.txt'),unpack=True)
sph_lon, sph_lat = np.loadtxt(path.join(datadir,'sph_lonlat.txt'),unpack=True)


#Creating the plots and figure
log.info('plotting map')
fig = plot.figure(figsize = (11,8), facecolor = 'white')

#fig, ax = plot.subplots(figsize=(16,9))
ax = plot.subplot(1,1,1,projection=PlateCarree())
ax.coastlines()
ax.set_extent([-180, 180, -61, 61], crs=PlateCarree())
ax.set_xticks([])
ax.set_yticks([])

if args.column == 'BAD_RATIO':
    vmin = 0.1
    vmax = 100.0
else:
    vmin = 5.0
    vmax = 100.0
for bk in args.bkffiles:
    bkftable = Table.read(bk,hdu=1)
    overshootrate = bkftable[args.column]
    sc = ax.scatter(bkftable['LON'], bkftable['LAT'],c=overshootrate,
                     norm=LogNorm(vmin=vmin,vmax=vmax),cmap='jet',s=2.0)

ax.plot(saa_lon,saa_lat,color='k',linestyle='dashed')
ax.plot(nph_lon,nph_lat,color='k',linestyle='dotted')
ax.plot(neph_lon,neph_lat,color='k',linestyle='dotted')
ax.plot(sph_lon,sph_lat,color='k',linestyle='dotted')
cbar = plot.colorbar(sc, location='bottom',pad=0.05,aspect=60)
plot.title(args.column)

plot.show()
