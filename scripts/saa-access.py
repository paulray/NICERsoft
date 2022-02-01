#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
import os.path as path
import argparse
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy import log
from astropy.coordinates import SkyCoord
import sys
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib.colors import Normalize, LogNorm
from cartopy.crs import PlateCarree


from nicer.mcc import MCC
from nicer.values import *
import nicer.yday_custom

tnow = Time.now()
dt = tnow.datetime
mnames = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
# Date Format: 18 Oct 2017 12:33:47
header = """                                                                         {0:2d} {1} {2:4d} {3:02d}:{4:02d}:{5:02d}
Satellite-NICER-To-AreaTarget-SAA:  Access Summary Report


NICER-To-SAA
------------
                  Access     Start Time (ISO-YD)       Stop Time (ISO-YD)     Duration (sec)
                  ------    ---------------------    ---------------------    --------------""".format(dt.day, mnames[dt.month], dt.year, dt.hour, dt.minute, dt.second)

parser = argparse.ArgumentParser(description = "Compute SAA and polar horn entry/exit times.")
parser.add_argument("mcc", help="Name of MCC ephemeris file to use")
args = parser.parse_args()

log.info('Reading MCC ephemeris')
eph = MCC(args.mcc)
dur_s = (eph.t[-1]-eph.t[0]).to(u.s).value

# Make a unique filename version
doystr = eph.t[0].yday.split(':')[1]
ver=0
while True:
    ver += 1
    outname = 'STK1_SAA_2017{0:3s}_V{1:02d}.txt'.format(doystr,ver)
    if not path.exists(outname):
        break
outfile = open(outname,'w')

# Make time spacing be 10 seconds over range of ephemeris
myt = eph.t[0] + np.arange(0.0,dur_s,10.0)*u.s
mymet = (myt-MET0).to(u.s).value

#llpoints =  zip(eph.lon.value,eph.lat.value)
mylat,mylon  = eph.latlon(mymet)
llpoints = zip(mylon.value,mylat.value)

log.info('Getting SAA data')
saa_poly = np.loadtxt(path.join(datadir,'saa_lonlat.txt'))
saa_path = mplPath.Path(saa_poly,closed=True)
in_saa = saa_path.contains_points(llpoints)

nph_poly = np.loadtxt(path.join(datadir,'nph_lonlat.txt'))
nph_path = mplPath.Path(nph_poly)
in_nph = nph_path.contains_points(llpoints)

neph_poly = np.loadtxt(path.join(datadir,'neph_lonlat.txt'))
neph_path = mplPath.Path(neph_poly)
in_neph = neph_path.contains_points(llpoints)

sph_poly = np.loadtxt(path.join(datadir,'sph_lonlat.txt'))
sph_path = mplPath.Path(sph_poly)
in_sph = sph_path.contains_points(llpoints)

nicer_poly = np.loadtxt(path.join(datadir,'nicer_saa.txt'))
nicer_path = mplPath.Path(nicer_poly)
in_nicer = nicer_path.contains_points(llpoints)

in_any = in_saa | in_nph | in_neph | in_sph | in_nicer

tstart = None
naccess=1
print(header,file=outfile)
prev = in_any[0]
tstartlist = []
tstoplist = []
durlist = []
accesslist = []
for i in range(len(in_any)):
    if in_any[i] != prev:
        if in_any[i]:
            tstart = myt[i]
        else:
            tstop = myt[i]
            if tstart is not None:
                print("                  {0:6d}    {1}    {2}    {3:14.3f}".format(
                    naccess,tstart.yday_custom,tstop.yday_custom,
                    (tstop-tstart).to(u.s).value),
                    file=outfile)
                tstartlist.append(tstart)
                tstoplist.append(tstop)
                durlist.append((tstop-tstart).to(u.s).value)
                accesslist.append(naccess)
                naccess += 1
    prev = in_any[i]

tstarts = np.array(tstartlist)
tstops = np.array(tstoplist)
durs = np.array(durlist)*u.s
accesses = np.array(accesslist,dtype=np.int)

# Global Statistics
# -----------------
# Min Duration          52    2017-303T21:07:28.946    2017-303T21:07:49.194            20.248
# Max Duration          34    2017-299T21:20:26.105    2017-299T21:27:55.512           449.407
# Mean Duration                                                                        308.829
# Total Duration                                                                     19147.423

print("""
Global Statistics
-----------------
Min Duration        {0:4d}    {1}    {2}     {3:9.3f}
Max Duration        {4:4d}    {5}    {6}     {7:9.3f}
Mean Duration                                                                     {8:10.3f}
Total Duration                                                                    {9:10.3f}
""".format(accesses[durs.argmin()], tstarts[durs.argmin()], tstops[durs.argmin()], durs.min().value,
           accesses[durs.argmax()], tstarts[durs.argmax()], tstops[durs.argmax()], durs.max().value,
           durs.mean().value, durs.sum().value),
           file=outfile)

outfile.close()

fig = plot.figure(figsize=(8, 11), facecolor='white')

ax = plot.subplot(1,1,1,projection=PlateCarree())
ax.coastlines()
ax.set_extent([-180, 180, -61, 61], crs=PlateCarree())
ax.set_xticks(np.arange(-60,60,15))
ax.set_yticks(np.arange(-180,180,15))

ax.scatter(mylon[in_any],mylat[in_any],s=1.0)
plt.show()
