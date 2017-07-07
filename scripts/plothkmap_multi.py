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
from astropy.coordinates import SkyCoord
import sys
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib.colors import Normalize, LogNorm
from mpl_toolkits.basemap import Basemap

from nicer.mcc import MCC
from nicer.sps import SPS
from nicer.values import *

parser = argparse.ArgumentParser(description = "Plot HKP data on map")
parser.add_argument("obsdirs", help="Input directories", nargs='+')
parser.add_argument("--under", help="Plot undershoot instead",action='store_true')
args = parser.parse_args()

if args.under:
    colname = 'MPU_UNDER_COUNT'
    desc = 'Undershoot (reset) Rate'
    vmin = 100.0
    vmax = 10000.0
else:
    colname = 'MPU_OVER_COUNT'
    desc = 'Overshoot Rate'
    vmin = 10.0
    vmax = 3000.0

saa_lon, saa_lat = np.loadtxt(path.join(datadir,'saa_lonlat.txt'),unpack=True)
nph_lon, nph_lat = np.loadtxt(path.join(datadir,'nph_lonlat.txt'),unpack=True)
sph_lon, sph_lat = np.loadtxt(path.join(datadir,'sph_lonlat.txt'),unpack=True)
fig, (ax1, ax2) = plt.subplots(2,1,figsize=(8.5,11))
map1 = Basemap(projection='cyl', resolution = 'l',
              lat_0=0, lon_0=0, ax=ax1)
map1.drawcoastlines()
map1.plot(saa_lon,saa_lat,'k',lw=2,linestyle='--')
map1.plot(nph_lon,nph_lat,color='orange',marker='o',markersize=10.0,linestyle='-')
map1.plot(sph_lon,sph_lat,'orange',marker='o',markersize=10.0,linestyle='-')

map2 = Basemap(projection='cyl', resolution = 'l',
              lat_0=0, lon_0=0, ax=ax2)
map2.drawcoastlines()
map2.plot(saa_lon,saa_lat,'k',lw=2,linestyle='--')
map2.plot(nph_lon,nph_lat,color='orange',marker='o',markersize=10.0,linestyle='-')
map2.plot(sph_lon,sph_lat,'orange',marker='o',markersize=10.0,linestyle='-')

#fig2, ax3 = plt.subplots()
#ax3.set_xlim((0.0,180.0))
#ax3.set_xlabel('Sun Angle (deg)')
#ax3.set_ylabel(desc)
#ax3.grid(True)

for obsdir in args.obsdirs:
    log.info('Processing '+obsdir)
    mkfname = glob(path.join(obsdir,'auxil/*.mkf'))[0]
    mkf = Table.read(mkfname,hdu=1)
    mkf.sort('TIME')
    sunmet = mkf['TIME']
    earthpos = SkyCoord(mkf['EARTH_RA'],mkf['EARTH_DEC'],frame='icrs')
    sunpos = SkyCoord(mkf['SUN_RA'],mkf['SUN_DEC'],frame='icrs')
    sunshine = mkf['SUNSHINE']

    spsname = glob(path.join(obsdir,'auxil/*apid0260.hk'))[0]
    log.info('Reading SPS file {0}'.format(spsname))
    eph = SPS(spsname)

    mpuhkfiles = glob(path.join(obsdir,'xti/hk/*mpu*.hk'))
    if len(mpuhkfiles) == 0:
        log.info('No files found')
        continue
    hdulist = pyfits.open(mpuhkfiles[0])
    td = hdulist[1].data
    met = td['TIME']
    log.info("MET Range {0} to {1}".format(met.min(),met.max()))
    t = MET0+met*u.s
    rate = td[colname].sum(axis=1)
    sep = earthpos.separation(sunpos)
    ESangle = np.interp(met,sunmet,sep)*sep.unit
    sunidx  = np.where(ESangle>70.0*u.deg)[0]
    darkidx = np.where(ESangle<70.0*u.deg)[0]

    for fn in mpuhkfiles[1:]:
        log.info('Reading '+fn)
        hdulist = pyfits.open(fn)
        mytd = hdulist[1].data
        mymet = td['TIME']
        myt = MET0+met*u.s
        myrate = td[colname].sum(axis=1)
        if not np.all(mymet == met):
            log.error('TIME axes are not compatible')
            sys.exit(1)
        rate += myrate

    if len(sunidx)>2:
        trimmet = met[sunidx]
        trimrate = rate[sunidx]
        lat, lon = eph.latlon(trimmet)
        sc1 = map1.scatter(lon, lat, s=10, c=trimrate, alpha=0.5,
            norm=LogNorm(vmin=vmin,vmax=vmax),cmap='jet')
    if len(darkidx)>2:
        trimmet = met[darkidx]
        trimrate = rate[darkidx]
        lat, lon = eph.latlon(trimmet)
        sc2 = map2.scatter(lon, lat, s=10, c=trimrate, alpha=0.5,
            norm=LogNorm(vmin=vmin,vmax=vmax),cmap='jet')

    #idx = np.isfinite(mkf['SUN_ANGLE'])
    #sa = mkf['SUN_ANGLE'][idx]
    #samet = sunmet[idx]
    #ax3.plot(sa,np.interp(samet,met,rate),'.')

ax1.set_title('In Sunshine')
cbar1 = map1.colorbar(sc1, location='bottom',pad='5%')
cbar1.set_label(desc)
ax2.set_title('In Darkness')
cbar2 = map2.colorbar(sc2, location='bottom',pad='5%')
cbar2.set_label(desc)

plt.show()
