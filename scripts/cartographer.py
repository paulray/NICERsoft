#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plot
from astropy import log
import astropy.units as u
from astropy.time import Time
from os import path

from nicer.mcc import MCC
from nicer.values import *

def cartography(hkmet, overshootrate):
    from mpl_toolkits.basemap import Basemap

    log.info('Getting SAA data')
    saa_lon, saa_lat = np.loadtxt(path.join(datadir,'saa_lonlat.txt'),unpack=True)

    eph = MCC(path.join(datadir,'MCC1_On_Console_20171631440_V01.txt'))
    log.info('got eph')

    lat, lon = eph.latlon(hkmet)

    log.info('plotting map')
    fig, ax = plot.subplots(figsize=(16,9))
    map = Basemap(projection='cyl', resolution = 'l', lat_0=0, lon_0=0)
    map.drawcoastlines()
    map.scatter(lon, lat,c=overshootrate,cmap='jet')
    map.plot(saa_lon,saa_lat,'g',lw=2)

    return fig
