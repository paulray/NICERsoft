#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plot
from astropy import log
import astropy.units as u
from astropy.time import Time
from os import path
from matplotlib.colors import Normalize, LogNorm
import matplotlib.gridspec as gridspec

from nicer.values import *
from nicer.mcc import MCC
#from nicer.sps import SPS
from nicer.latloninterp import LatLonInterp

def cartography(hkmet, overshootrate, args, undershootrate, etable, mktable):
    from mpl_toolkits.basemap import Basemap
#Getting the data
    log.info('Getting SAA data')
    saa_lon, saa_lat = np.loadtxt(path.join(datadir,'saa_lonlat.txt'),unpack=True)
    nph_lon, nph_lat = np.loadtxt(path.join(datadir,'nph_lonlat.txt'),unpack=True)
    sph_lon, sph_lat = np.loadtxt(path.join(datadir,'sph_lonlat.txt'),unpack=True)


    #eph = MCC(path.join(datadir,'MCC1_On_Console_20171631440_V01.txt'))

    # Use GPS SPS data for ISS location

    llinterp = LatLonInterp(mktable['TIME'], mktable['SAT_LAT'], mktable['SAT_LON'])

    #eph = SPS(args.sps)
    #log.info('got SPS ephemeris')
    #lat, lon = eph.latlon(hkmet)

    # Interpolate lon, lat from MKF housekeeping data
    lat, lon = llinterp.latlon(hkmet)

    #Creating the plots and figure
    log.info('plotting map')
    fig = plot.figure(figsize = (9,12), facecolor = 'white')
    mapper = gridspec.GridSpec(19,4)

    #fig, ax = plot.subplots(figsize=(16,9))
    overshoot = plot.subplot(mapper[1:10,0:4])

    map = Basemap(projection='cyl', resolution = 'l',  llcrnrlon=-180, llcrnrlat=-61,
    urcrnrlon=180, urcrnrlat=61, lat_0 = 0, lon_0 = 0)
    map.drawcoastlines()
    sc = map.scatter(lon, lat,c=overshootrate,norm=LogNorm(vmin=10.0,vmax=1000.0),cmap='jet')
    #sctest1 = map.scatter(lon, lat,c=np.ones(len(lon))*1000.0,norm=LogNorm(vmin=10.0,vmax=1000.0),cmap='jet')
    #sctest2 = map.scatter(mktable['SAT_LON'], mktable['SAT_LAT'],c=np.ones(len(mktable['SAT_LON']))*100.0,norm=LogNorm(vmin=10.0,vmax=1000.0),cmap='jet',alpha=0.5)
    map.plot(saa_lon,saa_lat,'r',lw=2)
    map.plot(nph_lon,nph_lat,color='orange',linestyle='-')
    map.plot(sph_lon,sph_lat,'orange',linestyle='-')
    cbar = map.colorbar(sc, location='bottom',pad='5%')
    plot.title('Overshoot Rate')
    #cbar.set_label('Overshoot Rate')

    undershoot = plot.subplot(mapper[10:19,0:4])

    map = Basemap(projection='cyl', resolution = 'l', llcrnrlon=-180, llcrnrlat=-61,
    urcrnrlon=180, urcrnrlat=61,lat_0 = 0, lon_0 = 0)
    map.drawcoastlines()
    sc = map.scatter(lon, lat,c=undershootrate,norm=LogNorm(vmin=10.0,vmax=1000.0),cmap='jet')
    map.plot(saa_lon,saa_lat,'r',lw=2)
    map.plot(nph_lon,nph_lat,color='orange',marker='o',markersize=5,linestyle='-')
    map.plot(sph_lon,sph_lat,'orange',marker='o',markersize=5,linestyle='-')
    cbar = map.colorbar(sc, location='bottom',pad='5%')
    plot.title('Undershoot Rate')
    cbar.set_label('Undershoot Rate')
    fig.suptitle('ObsID {0}: {1} on {2}'.format(etable.meta['OBS_ID'],
            etable.meta['OBJECT'],etable.meta['DATE-OBS'].replace('T',' at ')),
            fontsize=18)
    #plot.tight_layout()
    return fig
