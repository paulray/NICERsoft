#!/usr/bin/env python
from astropy.table import Table, vstack
import numpy as np
import sys

# This should be a list of _apid0260.hk files
filenames = sys.argv[1:]

for fn in filenames:
    tlist = []
    tlist.append(Table.read(fn, hdu=1))

    etable = vstack(tlist, metadata_conflicts="silent")
    # print(np.where(etable['GPS_SPS_TIME_VALID'] == 0))

    TIMEzerocount = len(np.where(etable["GPS_SPS_TIME_VALID"] == 0)[0])
    TIMEonecount = len(np.where(etable["GPS_SPS_TIME_VALID"] == 1)[0])
    TIMEdiffcount = len(np.where(etable["GPS_SPS_TIME_VALID"] > 1)[0])
    zerocount = len(np.where(etable["GPS_SPS_VALID"] == 0)[0])
    onecount = len(np.where(etable["GPS_SPS_VALID"] == 1)[0])

    idx = np.where(etable["GPS_SPS_GDOP"] != 0.0)[0]
    gdop = etable["GPS_SPS_GDOP"][idx].mean()

    sol = "{0}: SPS_VALID (0 = {1:5d}, 1= {2:5d}), SPS_TIME_VALID (0 = {3:5d}, 1 = {4:5d}, >1 = {5:5d}) <GDOP> = {6:.2f}".format(
        fn[48:71], zerocount, onecount, TIMEzerocount, TIMEonecount, TIMEdiffcount, gdop
    )

    print(sol)
