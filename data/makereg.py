#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np

saacols = np.loadtxt('saa_lonlat.txt')
# For .reg file need to convert lon to [0,360)
lon = saacols[:,0]
lon[lon<0]+= 360.0
saacols[:,0] = lon

outfile = file('saa.reg','w')
print("-polygon(",file=outfile,end="")
coordstring = ", ".join(["{0}".format(x) for x in saacols.flatten()])
print(coordstring,file=outfile,end="")
print(")",file=outfile)
outfile.close()

outfile = file('polarhorns.reg','w')
nphcols = np.loadtxt('nph_lonlat.txt')
# For .reg file need to convert lon to [0,360)
lon = nphcols[:,0]
lon[lon<0]+= 360.0
nphcols[:,0] = lon

print("-polygon(",file=outfile,end="")
coordstring = ", ".join(["{0}".format(x) for x in nphcols.flatten()])
print(coordstring,file=outfile,end="")
print(")",file=outfile)

sphcols = np.loadtxt('sph_lonlat.txt')
# For .reg file need to convert lon to [0,360)
lon = sphcols[:,0]
lon[lon<0]+= 360.0
sphcols[:,0] = lon

print("-polygon(",file=outfile,end="")
coordstring = ", ".join(["{0}".format(x) for x in sphcols.flatten()])
print(coordstring,file=outfile,end="")
print(")",file=outfile)
outfile.close()
