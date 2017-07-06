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
from nicer.mcc import MCC
from nicer.sps import SPS

def cartography(hkmet, overshootrate, args):
    from mpl_toolkits.basemap import Basemap

    log.info('Getting SAA data')
    saa_lon, saa_lat = np.loadtxt(path.join(datadir,'saa_lonlat.txt'),unpack=True)
    nph_lon, nph_lat = np.loadtxt(path.join(datadir,'nph_lonlat.txt'),unpack=True)
    sph_lon, sph_lat = np.loadtxt(path.join(datadir,'sph_lonlat.txt'),unpack=True)

    #eph = MCC(path.join(datadir,'MCC1_On_Console_20171631440_V01.txt'))
    # Use GPS SPS data for ISS location
    eph = SPS(args.sps)
    log.info('got SPS ephemeris')
    lat, lon = eph.latlon(hkmet)


    log.info('plotting map')
    fig, ax = plot.subplots(figsize=(16,9))
    map = Basemap(projection='cyl', resolution = 'l', lat_0=0, lon_0=0)
    map.drawcoastlines()
    sc = map.scatter(lon, lat,c=overshootrate,norm=LogNorm(vmin=10.0,vmax=1000.0),cmap='jet')
    map.plot(saa_lon,saa_lat,'r',lw=2)
    map.plot(nph_lon,nph_lat,color='orange',marker='o',markersize=10.0,linestyle='-')
    map.plot(sph_lon,sph_lat,'orange',marker='o',markersize=10.0,linestyle='-')
    cbar = map.colorbar(sc, location='bottom',pad='5%')
    cbar.set_label('Overshoot Rate')

    return fig
