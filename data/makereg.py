#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np

saacols = np.loadtxt('saa_lonlat.txt')

outfile = file('saa.reg','w')
print("-polygon(",file=outfile,end="")
coordstring = ", ".join(["{0}".format(x) for x in saacols.flatten()])
print(coordstring,file=outfile,end="")
print(")",file=outfile)
outfile.close()

outfile = file('polarhorns.reg','w')
nphcols = np.loadtxt('nph_lonlat.txt')
print("-polygon(",file=outfile,end="")
coordstring = ", ".join(["{0}".format(x) for x in nphcols.flatten()])
print(coordstring,file=outfile,end="")
print(")",file=outfile)

sphcols = np.loadtxt('sph_lonlat.txt')
print("-polygon(",file=outfile,end="")
coordstring = ", ".join(["{0}".format(x) for x in sphcols.flatten()])
print(coordstring,file=outfile,end="")
print(")",file=outfile)
outfile.close()

