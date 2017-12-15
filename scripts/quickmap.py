#!/usr/bin/env python
from __future__ import print_function, division
import argparse
import numpy as np
import matplotlib.pyplot as plot
import matplotlib as mpl
from astropy import log
import astropy.units as u
from astropy.time import Time
from astropy.table import Table
from os import path
from matplotlib.colors import Normalize, LogNorm
from mpl_toolkits.basemap import Basemap
from nicer.plotutils import *

from nicer.values import *
from nicer.mcc import MCC
#from nicer.sps import SPS
from nicer.latloninterp import LatLonInterp

parser = argparse.ArgumentParser(description="Plot points on a map, filtering by tot.gti")
parser.add_argument("bkffiles", help="Name of bkf (or mkf) files to process", nargs='+')
args = parser.parse_args()

log.info('Getting SAA data')
saa_lon, saa_lat = np.loadtxt(path.join(datadir,'saa_lonlat.txt'),unpack=True)
nph_lon, nph_lat = np.loadtxt(path.join(datadir,'nph_lonlat.txt'),unpack=True)
neph_lon, neph_lat = np.loadtxt(path.join(datadir,'neph_lonlat.txt'),unpack=True)
sph_lon, sph_lat = np.loadtxt(path.join(datadir,'sph_lonlat.txt'),unpack=True)
nicer_lon, nicer_lat = np.loadtxt(path.join(datadir,'nicer_saa.txt'),unpack=True)

bkftable = Table.read(args.bkffiles[0],hdu=1)
gtitable = Table.read("tot.gti",hdu=2)

if gtitable is not None:
    startmet = gtitable['START'][0]
    stopmet = gtitable['STOP'][0]
    duration = stopmet-startmet
    # Add 1 bin to make sure last bin covers last events
    lcbinsize=1.0
    lc_elapsed_bins = np.arange(0.0,duration+lcbinsize,lcbinsize)
    lc_met_bins = startmet+lc_elapsed_bins
    cumtimes = [ 0.0 ]
    cumtime = lc_elapsed_bins[-1]+lcbinsize
    for i in range(1,len(gtitable['START'])):
        startmet = gtitable['START'][i]
        stopmet = gtitable['STOP'][i]
        duration = stopmet-startmet
        myelapsedbins = np.arange(0.0,duration+lcbinsize,lcbinsize)
        lc_elapsed_bins = np.append(lc_elapsed_bins,cumtime+myelapsedbins)
        lc_met_bins = np.append(lc_met_bins,np.arange(startmet,stopmet+lcbinsize,lcbinsize))
        mylcduration = myelapsedbins[-1]+lcbinsize
        cumtimes.append(cumtime)
        cumtime += mylcduration
    gtitable['CUMTIME'] = np.array(cumtimes)

#Creating the plots and figure
log.info('plotting map')
fig = plot.figure(figsize = (11,8), facecolor = 'white')

plot.subplot(1,1,1)

map = Basemap(projection='cyl', resolution = 'l',  llcrnrlon=-180, llcrnrlat=-61,
urcrnrlon=180, urcrnrlat=61, lat_0 = 0, lon_0 = 0)
map.drawcoastlines()

lon = bkftable['SAT_LON']
# Convert longitudes to -180, 180
lon[lon>180] -= 360.0
lat = bkftable['SAT_LAT']

colornames, cmap, norm = gti_colormap()

etime, goodlon, cc = convert_to_elapsed_goodtime(bkftable['TIME'], lon, gtitable)
etime, goodlat, cc = convert_to_elapsed_goodtime(bkftable['TIME'], lat, gtitable)

print(len(lon), len(lat), len(cc))
sc = map.scatter(goodlon, goodlat,c=np.fmod(cc,len(colornames)),s=2.0,cmap=cmap,
    norm=norm)

map.plot(saa_lon,saa_lat,'r',lw=2)
map.plot(nph_lon,nph_lat,color='orange',linestyle='-')
map.plot(neph_lon,neph_lat,color='orange',linestyle='-')
map.plot(sph_lon,sph_lat,'orange',linestyle='-')
map.plot(nicer_lon,nicer_lat,'black',linestyle='-')
cbar = map.colorbar(sc, location='bottom',pad='5%')
plot.title("Map")

plot.show()
