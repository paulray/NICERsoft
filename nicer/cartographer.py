#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plot
from astropy import log
import astropy.units as u
from astropy.time import Time
from os import path
from matplotlib.colors import Normalize, LogNorm

from nicer.values import *
from nicer.plotutils import *
# These next three are different ways of computing ISS lat, lon vs time
#from nicer.mcc import MCC
#from nicer.sps import SPS
from nicer.latloninterp import LatLonInterp

def cartography(hkmet, overshootrate, args, undershootrate, etable, mktable, gtitable):
    from mpl_toolkits.basemap import Basemap
#Getting the data
    log.info('Getting SAA data')
    saa_lon, saa_lat = np.loadtxt(path.join(datadir,'saa_lonlat.txt'),unpack=True)
    nph_lon, nph_lat = np.loadtxt(path.join(datadir,'nph_lonlat.txt'),unpack=True)
    neph_lon, neph_lat = np.loadtxt(path.join(datadir,'neph_lonlat.txt'),unpack=True)
    sph_lon, sph_lat = np.loadtxt(path.join(datadir,'sph_lonlat.txt'),unpack=True)
    nicer_lon, nicer_lat = np.loadtxt(path.join(datadir,'nicer_saa.txt'),unpack=True)


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
    fig = plot.figure(figsize = (8,11), facecolor = 'white')

    #fig, ax = plot.subplots(figsize=(16,9))
    overshoot = plot.subplot(3,1,1)

    map = Basemap(projection='cyl', resolution = 'l',  llcrnrlon=-180, llcrnrlat=-61,
    urcrnrlon=180, urcrnrlat=61, lat_0 = 0, lon_0 = 0)
    map.drawcoastlines()
    sc = map.scatter(lon, lat,c=overshootrate, s=2.0,
        norm=LogNorm(vmin=10.0,vmax=1000.0), cmap='jet')
    #sctest1 = map.scatter(lon, lat,c=np.ones(len(lon))*1000.0,norm=LogNorm(vmin=10.0,vmax=1000.0),cmap='jet')
    #sctest2 = map.scatter(mktable['SAT_LON'], mktable['SAT_LAT'],c=np.ones(len(mktable['SAT_LON']))*100.0,norm=LogNorm(vmin=10.0,vmax=1000.0),cmap='jet',alpha=0.5)
    map.plot(saa_lon,saa_lat,'orange',linestyle='-')
    map.plot(nph_lon,nph_lat,color='orange',linestyle='-')
    map.plot(neph_lon,neph_lat,color='orange',linestyle='-')
    map.plot(sph_lon,sph_lat,'orange',linestyle='-')
    map.plot(nicer_lon,nicer_lat,'red',linestyle='-')
    cbar = map.colorbar(sc, location='bottom',pad='5%')
    plot.ylabel('Overshoot Rate')
    #cbar.set_label('Overshoot Rate')

    undershoot = plot.subplot(3,1,2)

    map = Basemap(projection='cyl', resolution = 'l', llcrnrlon=-180, llcrnrlat=-61,
    urcrnrlon=180, urcrnrlat=61,lat_0 = 0, lon_0 = 0)
    map.drawcoastlines()
    sc = map.scatter(lon, lat,c=undershootrate, s=2.0,
        norm=LogNorm(vmin=10.0,vmax=1000.0), cmap='jet')
    map.plot(saa_lon,saa_lat,color='orange',marker='o',markersize=2,linestyle='-')
    map.plot(nph_lon,nph_lat,color='orange',marker='o',markersize=2,linestyle='-')
    map.plot(neph_lon,neph_lat,color='orange',marker='o',markersize=2,linestyle='-')
    map.plot(sph_lon,sph_lat,'orange',marker='o',markersize=2,linestyle='-')
    map.plot(nicer_lon,nicer_lat,'red',linestyle='-')
    cbar = map.colorbar(sc, location='bottom',pad='5%')
    plot.ylabel('Undershoot Rate')
    #cbar.set_label('Undershoot Rate')

    gtiplot = plot.subplot(3,1,3)

    # Plot map with just colors for GTIs
    map = Basemap(projection='cyl', resolution = 'l', llcrnrlon=-180, llcrnrlat=-61,
    urcrnrlon=180, urcrnrlat=61,lat_0 = 0, lon_0 = 0)
    map.drawcoastlines()
    etime, goodlon, cc = convert_to_elapsed_goodtime(hkmet, lon, gtitable)
    etime, goodlat, cc = convert_to_elapsed_goodtime(hkmet, lat, gtitable)
    colornames, cmap, norm = gti_colormap()
    sc = map.scatter(goodlon, goodlat,c=np.fmod(cc,len(colornames)),s=2.0,cmap=cmap,
        norm=norm)
    map.plot(saa_lon,saa_lat,color='orange',marker='o',markersize=2,linestyle='-')
    map.plot(nph_lon,nph_lat,color='orange',marker='o',markersize=2,linestyle='-')
    map.plot(neph_lon,neph_lat,color='orange',marker='o',markersize=2,linestyle='-')
    map.plot(sph_lon,sph_lat,'orange',marker='o',markersize=2,linestyle='-')
    map.plot(nicer_lon,nicer_lat,'red',linestyle='-')
    plot.ylabel('GTI Colors')

    fig.suptitle('ObsID {0}: {1} on {2}'.format(etable.meta['OBS_ID'],
            etable.meta['OBJECT'],etable.meta['DATE-OBS'].replace('T',' at ')),
            fontsize=14)

    #plot.tight_layout()


    return fig
