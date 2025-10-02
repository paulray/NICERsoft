#!/usr/bin/env python
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
from cartopy.crs import PlateCarree

from nicer.mcc import MCC
from nicer.sps import SPS
from nicer.values import *

parser = argparse.ArgumentParser(description="Plot HKP data on map")
parser.add_argument("obsdirs", help="Input directories", nargs="+")
parser.add_argument("--under", help="Plot undershoot instead", action="store_true")
parser.add_argument("--mc", help="Plot MCILWAIN_L instead", action="store_true")
args = parser.parse_args()

if args.under:
    colname = "MPU_UNDER_COUNT"
    desc = "Undershoot (reset) Rate"
    vmin = 100.0
    vmax = 10000.0
else:
    colname = "COR_SAX"
    desc = "Cutoff Rigidity"
    vmin = 1.0
    vmax = 13.0

if args.mc:
    colname = "MCILWAIN_L"
    desc = "MCILWAIN_L"
    vmin = -90.0
    vmax = 90


saa_lon, saa_lat = np.loadtxt(path.join(datadir, "saa_lonlat.txt"), unpack=True)
nph_lon, nph_lat = np.loadtxt(path.join(datadir, "nph_lonlat.txt"), unpack=True)
neph_lon, neph_lat = np.loadtxt(path.join(datadir, "neph_lonlat.txt"), unpack=True)
sph_lon, sph_lat = np.loadtxt(path.join(datadir, "sph_lonlat.txt"), unpack=True)

fig, ax1 = plt.subplots(1, 1, figsize=(11, 8.5), projection=PlateCarree())
ax1.coastlines()
ax1.set_extent([-180, 180, -61, 61], crs=PlateCarree())
ax1.set_xticks([])
ax1.set_yticks([])

ax1.plot(saa_lon, saa_lat, color="orange", marker="o", markersize=2.0, linestyle="-")
ax1.plot(nph_lon, nph_lat, color="orange", marker="o", markersize=2.0, linestyle="-")
ax1.plot(neph_lon, neph_lat, color="orange", marker="o", markersize=2.0, linestyle="-")
ax1.plot(sph_lon, sph_lat, "orange", marker="o", markersize=2.0, linestyle="-")

for obsdir in args.obsdirs:
    log.info("Processing " + obsdir)
    mkfname = glob(path.join(obsdir, "auxil/*.mkf"))[0]
    mkf = Table.read(mkfname, hdu=1)
    mkf.sort("TIME")
    mkfmet = mkf["TIME"]
    earthpos = SkyCoord(mkf["EARTH_RA"], mkf["EARTH_DEC"], frame="icrs")
    sunpos = SkyCoord(mkf["SUN_RA"], mkf["SUN_DEC"], frame="icrs")
    sunshine = mkf["SUNSHINE"]
    lon = mkf["SAT_LON"]
    lat = mkf["SAT_LAT"]

    # Adjust lon to be -180 to 180 instead of 0 to 360
    lon[lon > 180.0] -= 360.0

    sc1 = ax1.scatter(
        lon,
        lat,
        s=5,
        c=mkf[colname],
        norm=Normalize(vmin=vmin, vmax=vmax),
        cmap="Paired",
    )


ax1.set_title("Map")
cbar1 = plot.colorbar(sc1, location="bottom", pad=0.05, aspect=60)
ax1.set_aspect("auto", adjustable=None)
cbar1.set_label(desc)

plt.show()
