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

#def cartography(hkmet, overshootrate, args, undershootrate, etable, mktable, gtitable):
def cartography(etable, mktable, gtitable, args):

    from cartopy.crs import PlateCarree

    #Getting the data
    log.info('Getting SAA data')
    saa_lon, saa_lat = np.loadtxt(path.join(datadir,'saa_lonlat.txt'),unpack=True)
    nph_lon, nph_lat = np.loadtxt(path.join(datadir,'nph_lonlat.txt'),unpack=True)
    neph_lon, neph_lat = np.loadtxt(path.join(datadir,'neph_lonlat.txt'),unpack=True)
    sph_lon, sph_lat = np.loadtxt(path.join(datadir,'sph_lonlat.txt'),unpack=True)
    nicer_lon, nicer_lat = np.loadtxt(path.join(datadir,'nicer_saa.txt'),unpack=True)


    #eph = MCC(path.join(datadir,'MCC1_On_Console_20171631440_V01.txt'))

    # Use GPS SPS data for ISS location

    #llinterp = LatLonInterp(mktable['TIME'], mktable['SAT_LAT'], mktable['SAT_LON'])

    #eph = SPS(args.sps)
    #log.info('got SPS ephemeris')
    #lat, lon = eph.latlon(hkmet)

    # Interpolate lon, lat from MKF housekeeping data
    #lat, lon = llinterp.latlon(hkmet)

    #Creating the plots and figure
    log.info('plotting map')
    fig = plot.figure(figsize = (8,11), facecolor = 'white')

    # Get lat, lon, overshoots from the .mkf table
    ovtime, overshootrate, cc = convert_to_elapsed_goodtime(mktable['TIME'], mktable['NUM_FPM_ON']*mktable['FPM_OVERONLY_COUNT'], gtitable)
    ovtime, lat, cc = convert_to_elapsed_goodtime(mktable['TIME'], mktable['SAT_LAT'], gtitable)
    ovtime, lon, cc = convert_to_elapsed_goodtime(mktable['TIME'], mktable['SAT_LON'], gtitable)

    # remove "nan" from arrays
    ovtimeNan = np.where(np.isnan(overshootrate))
    ovtimeNew=np.delete(ovtime,ovtimeNan)
    overshootrate=np.delete(overshootrate,ovtimeNan)
    lat = np.delete(lat,ovtimeNan)
    lon = np.delete(lon,ovtimeNan)
    # Adjust lon to be -180 to 180 instead of 0 to 360
    lon[lon>180.0] -= 360.0

    # Plot Overshoot map
    ax_ov  = plot.subplot(3,1,1,projection=PlateCarree())
    ax_ov.coastlines()
    ax_ov.set_extent([-180, 180, -61, 61], crs=PlateCarree())
    ax_ov.set_xticks([])
    ax_ov.set_yticks([])

    sc = ax_ov.scatter(lon, lat,c=overshootrate, s=2.0, norm=LogNorm(vmin=10.0,vmax=1000.0), cmap='jet')
    ax_ov.plot(saa_lon,saa_lat,'orange',linestyle='-')
    ax_ov.plot(nph_lon,nph_lat,color='orange',linestyle='-')
    ax_ov.plot(neph_lon,neph_lat,color='orange',linestyle='-')
    ax_ov.plot(sph_lon,sph_lat,'orange',linestyle='-')
    ax_ov.plot(nicer_lon,nicer_lat,'red',linestyle='-')
    cbar = plot.colorbar(sc, location='bottom',pad=0.05,aspect=60)
    ax_ov.set_aspect('auto', adjustable=None)
    ax_ov.set_ylabel('Overshoot Rate')

    # Get lat, lon, undershoots from the .mkf table
    udtime, undershootrate, cc = convert_to_elapsed_goodtime(mktable['TIME'], 52*mktable['FPM_UNDERONLY_COUNT'], gtitable)
    udtime, lat, cc = convert_to_elapsed_goodtime(mktable['TIME'], mktable['SAT_LAT'], gtitable)
    udtime, lon, cc = convert_to_elapsed_goodtime(mktable['TIME'], mktable['SAT_LON'], gtitable)

    # remove "nan" from arrays
    udtimeNan = np.where(np.isnan(undershootrate))
    udtimeNew=np.delete(udtime,udtimeNan)
    undershootrate=np.delete(undershootrate,udtimeNan)
    lat = np.delete(lat,udtimeNan)
    lon = np.delete(lon,udtimeNan)
    # Adjust lon to be -180 to 180 instead of 0 to 360
    lon[lon>180.0] -= 360.0

    # Plot Undershoot map
    ax_ud = plot.subplot(3,1,2,projection=PlateCarree())
    ax_ud.coastlines()
    ax_ud.set_extent([-180, 180, -61, 61], crs=PlateCarree())
    ax_ud.set_xticks([])
    ax_ud.set_yticks([])

    sc = ax_ud.scatter(lon, lat,c=undershootrate, s=2.0, norm=LogNorm(vmin=10.0,vmax=1000.0), cmap='jet')
    ax_ud.plot(saa_lon,saa_lat,color='orange',marker='o',markersize=2,linestyle='-')
    ax_ud.plot(nph_lon,nph_lat,color='orange',marker='o',markersize=2,linestyle='-')
    ax_ud.plot(neph_lon,neph_lat,color='orange',marker='o',markersize=2,linestyle='-')
    ax_ud.plot(sph_lon,sph_lat,'orange',marker='o',markersize=2,linestyle='-')
    ax_ud.plot(nicer_lon,nicer_lat,'red',linestyle='-')
    cbar = plot.colorbar(sc, location='bottom',pad=0.05,aspect=60)
    ax_ud.set_aspect('auto', adjustable=None)
    ax_ud.set_ylabel('Undershoot Rate')

    ax_gti = plot.subplot(3,1,3,projection=PlateCarree())
    ax_gti.coastlines()
    ax_gti.set_extent([-180, 180, -61, 61], crs=PlateCarree())
    ax_gti.set_xticks([])
    ax_gti.set_yticks([])


    # Plot map with just colors for GTIs
    etime, goodlon, cc = convert_to_elapsed_goodtime(mktable['TIME'], mktable['SAT_LON'], gtitable)
    etime, goodlat, cc = convert_to_elapsed_goodtime(mktable['TIME'], mktable['SAT_LAT'], gtitable)
    # Adjust lon to be -180 to 180 instead of 0 to 360
    goodlon[goodlon>180.0] -= 360.0

    colornames, cmap, norm = gti_colormap()
    sc = ax_gti.scatter(goodlon, goodlat,c=np.fmod(cc,len(colornames)),s=2.0,cmap=cmap, norm=norm)
    ax_gti.plot(saa_lon,saa_lat,color='orange',marker='o',markersize=2,linestyle='-')
    ax_gti.plot(nph_lon,nph_lat,color='orange',marker='o',markersize=2,linestyle='-')
    ax_gti.plot(neph_lon,neph_lat,color='orange',marker='o',markersize=2,linestyle='-')
    ax_gti.plot(sph_lon,sph_lat,'orange',marker='o',markersize=2,linestyle='-')
    ax_gti.plot(nicer_lon,nicer_lat,'red',linestyle='-')
    ax_gti.set_ylabel('GTI Colors')

    fig.suptitle('ObsID {0}: {1} on {2}'.format(etable.meta['OBS_ID'],
            etable.meta['OBJECT'],etable.meta['DATE-OBS'].replace('T',' at ')),
            fontsize=14)

    #plot.tight_layout()


    return fig
