#!/usr/bin/env python
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
from cartopy.crs import PlateCarree
from nicer.plotutils import *

from nicer.values import *
from nicer.mcc import MCC

# from nicer.sps import SPS
from nicer.latloninterp import LatLonInterp

parser = argparse.ArgumentParser(
    description="Plot points on a map, filtering by tot.gti"
)
parser.add_argument("bkffiles", help="Name of bkf (or mkf) files to process", nargs="+")
args = parser.parse_args()

log.info("Getting SAA data")
saa_lon, saa_lat = np.loadtxt(path.join(datadir, "saa_lonlat.txt"), unpack=True)
nph_lon, nph_lat = np.loadtxt(path.join(datadir, "nph_lonlat.txt"), unpack=True)
neph_lon, neph_lat = np.loadtxt(path.join(datadir, "neph_lonlat.txt"), unpack=True)
sph_lon, sph_lat = np.loadtxt(path.join(datadir, "sph_lonlat.txt"), unpack=True)
nicer_lon, nicer_lat = np.loadtxt(path.join(datadir, "nicer_saa.txt"), unpack=True)

bkftable = Table.read(args.bkffiles[0], hdu=1)
gtitable = Table.read("tot.gti", hdu=2)

if gtitable is not None:
    startmet = gtitable["START"][0]
    stopmet = gtitable["STOP"][0]
    duration = stopmet - startmet
    # Add 1 bin to make sure last bin covers last events
    lcbinsize = 1.0
    lc_elapsed_bins = np.arange(0.0, duration + lcbinsize, lcbinsize)
    lc_met_bins = startmet + lc_elapsed_bins
    cumtimes = [0.0]
    cumtime = lc_elapsed_bins[-1] + lcbinsize
    for i in range(1, len(gtitable["START"])):
        startmet = gtitable["START"][i]
        stopmet = gtitable["STOP"][i]
        duration = stopmet - startmet
        myelapsedbins = np.arange(0.0, duration + lcbinsize, lcbinsize)
        lc_elapsed_bins = np.append(lc_elapsed_bins, cumtime + myelapsedbins)
        lc_met_bins = np.append(
            lc_met_bins, np.arange(startmet, stopmet + lcbinsize, lcbinsize)
        )
        mylcduration = myelapsedbins[-1] + lcbinsize
        cumtimes.append(cumtime)
        cumtime += mylcduration
    gtitable["CUMTIME"] = np.array(cumtimes)

# Creating the plots and figure
log.info("plotting map")
fig = plot.figure(figsize=(11, 8), facecolor="white")

ax = plot.subplot(1, 1, 1, projection=PlateCarree())
ax.coastlines()
ax.set_extent([-180, 180, -61, 61], crs=PlateCarree())
ax.set_xticks([])
ax.set_yticks([])

lon = bkftable["SAT_LON"]
# Convert longitudes to -180, 180
lon[lon > 180] -= 360.0
lat = bkftable["SAT_LAT"]

colornames, cmap, norm = gti_colormap()

etime, goodlon, cc = convert_to_elapsed_goodtime(bkftable["TIME"], lon, gtitable)
etime, goodlat, cc = convert_to_elapsed_goodtime(bkftable["TIME"], lat, gtitable)

print(len(lon), len(lat), len(cc))
sc = ax.scatter(
    goodlon, goodlat, c=np.fmod(cc, len(colornames)), s=2.0, cmap=cmap, norm=norm
)

ax.plot(saa_lon, saa_lat, "r", lw=2)
ax.plot(nph_lon, nph_lat, color="orange", linestyle="-")
ax.plot(neph_lon, neph_lat, color="orange", linestyle="-")
ax.plot(sph_lon, sph_lat, "orange", linestyle="-")
ax.plot(nicer_lon, nicer_lat, "black", linestyle="-")
cbar = plot.colorbar(sc, location="bottom", pad=0.05, aspect=60)
ax.set_aspect("auto", adjustable=None)
plot.title("Map")

plot.show()
