#!/usr/bin/env python
from __future__ import print_function, division
import argparse
import numpy as np
import astropy.units as u

from astropy.time import Time
MET0 = Time("2014-01-01T00:00:00.0",scale='utc')
GPS0 = Time("1980-01-06T00:00:00",scale='utc')

parser = argparse.ArgumentParser(description="Print NICER times in a variety of formats.")
parser.add_argument("--iso",help="ISO Date (e.g. 2017-05-16)",default=None)
parser.add_argument("--met",type=float,help="MET (e.g. 108893169.0)",default=None)
parser.add_argument("--gps",type=float,help="GPS Time (e.g. 1180440832.0)",default=None)
args = parser.parse_args()

if args.met is not None:
    x = MET0+args.met*u.s
elif args.iso is not None:
    x = Time(args.iso,scale='utc')
elif args.gps is not None:
    x = GPS0+args.gps
else:
    x = Time.now()

print('UTC        : ',x)
x.format = 'yday'
x.out_subfmt = 'date'
print('YEAR:DOY   : ',x)

print("MET        : ",(x-MET0).to(u.s).value)
print("GPS Time   : ",(x-GPS0).to(u.s).value)
print("MJD (UTC)  : ",x.utc.mjd)
print("MJD (TT)   : ",x.tt.mjd)
