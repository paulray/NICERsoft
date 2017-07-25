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
from nicer.sps import SPS

def cartography(hkmet, overshootrate, args, undershootrate, etable, mktable):
    from mpl_toolkits.basemap import Basemap
#Getting the data
    log.info('Getting SAA data')
    saa_lon, saa_lat = np.loadtxt(path.join(datadir,'saa_lonlat.txt'),unpack=True)
    nph_lon, nph_lat = np.loadtxt(path.join(datadir,'nph_lonlat.txt'),unpack=True)
    sph_lon, sph_lat = np.loadtxt(path.join(datadir,'sph_lonlat.txt'),unpack=True)

    #eph = MCC(path.join(datadir,'MCC1_On_Console_20171631440_V01.txt'))
   
    # Use GPS SPS data for ISS location
    
    lat = mktable['SAT_LAT']
    lon = mktable['SAT_LON']
    
    #Making the overshootrate the same length as lat/lon
    evenbins = np.histogram(overshootrate, bins = len(lat))[0]
    sum = []
    mysum = 0

    for num in evenbins:
        count = 0
        for i in xrange(count,count+num):
            mysum = mysum + overshootrate[i]
        sum = np.append(sum, (mysum/num))
        mysum = 0
        count = count + num
    
    overshootrate = sum

    #Making the undershootrate the same length as the lat/lon
    evenbins = np.histogram(undershootrate, bins = len(lat))[0]
    sum = []
    mysum = 0
        
    for num in evenbins:
        count = 0
        for i in xrange(count,count+num):
            mysum = mysum + undershootrate[i]
        sum = np.append(sum, (mysum/num))
        mysum = 0
        count = count + num

    undershootrate = sum

    #eph = SPS(args.sps)
    #log.info('got SPS ephemeris')
    #lat, lon = eph.latlon(hkmet)

#Creating the plots and figure
    log.info('plotting map')
    fig = plot.figure(figsize = (9,12), facecolor = 'white')
    mapper = gridspec.GridSpec(19,4)
    
    #fig, ax = plot.subplots(figsize=(16,9))
    overshoot = plot.subplot(mapper[1:10,0:4])

    map = Basemap(projection='cyl', resolution = 'l', lat_0=0, lon_0=0)
    map.drawcoastlines()
    sc = map.scatter(lon, lat,c=overshootrate,norm=LogNorm(vmin=10.0,vmax=1000.0),cmap='jet')
    #sc = map.scatter(lon,lat)
    map.plot(saa_lon,saa_lat,'r',lw=2)
    map.plot(nph_lon,nph_lat,color='orange',marker='o',markersize=10.0,linestyle='-')
    map.plot(sph_lon,sph_lat,'orange',marker='o',markersize=10.0,linestyle='-')
    cbar = map.colorbar(sc, location='bottom',pad='5%')
    plot.title('Overshoot Rate')
    #cbar.set_label('Overshoot Rate')

    undershoot = plot.subplot(mapper[10:19,0:4])

    map = Basemap(projection='cyl', resolution = 'l', lat_0=0, lon_0=0)
    map.drawcoastlines()
    sc = map.scatter(lon, lat,c=undershootrate,norm=LogNorm(vmin=10.0,vmax=1000.0),cmap='jet')
    map.plot(saa_lon,saa_lat,'r',lw=2)
    map.plot(nph_lon,nph_lat,color='orange',marker='o',markersize=10.0,linestyle='-')
    map.plot(sph_lon,sph_lat,'orange',marker='o',markersize=10.0,linestyle='-')
    cbar = map.colorbar(sc, location='bottom',pad='5%')
    plot.title('Undershoot Rate')
    cbar.set_label('Undershoot Rate')
    fig.suptitle('ObsID {0}: {1} on {2}'.format(etable.meta['OBS_ID'],
            etable.meta['OBJECT'],etable.meta['DATE-OBS'].replace('T',' at ')),
            fontsize=18)
    plot.tight_layout()
    return fig
