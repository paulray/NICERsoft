#!/usr/bin/env python
from astropy.table import Table, vstack
import numpy as np

filenames = [
'/data/NICER/preliminary/fits1/1706140400/auxil/ni1706140400_apid0260.hk',
'/data/NICER/preliminary/fits1/1706140637/auxil/ni1706140637_apid0260.hk',
'/data/NICER/preliminary/fits1/1706140806/auxil/ni1706140806_apid0260.hk',
'/data/NICER/preliminary/fits1/1706141302/auxil/ni1706141302_apid0260.hk',
'/data/NICER/preliminary/fits1/1706152137/auxil/ni1706152137_apid0260.hk',
'/data/NICER/preliminary/fits1/1706152233/auxil/ni1706152233_apid0260.hk',
'/data/NICER/preliminary/fits1/1706152346/auxil/ni1706152346_apid0260.hk',
'/data/NICER/preliminary/fits1/1706160735/auxil/ni1706160735_apid0260.hk',
'/data/NICER/preliminary/fits1/1706161001/auxil/ni1706161001_apid0260.hk',
'/data/NICER/preliminary/fits1/1706161405/auxil/ni1706161405_apid0260.hk',
'/data/NICER/preliminary/fits1/1706161524/auxil/ni1706161524_apid0260.hk',
'/data/NICER/preliminary/fits1/1706161601/auxil/ni1706161601_apid0260.hk',
'/data/NICER/preliminary/fits1/1706162128/auxil/ni1706162128_apid0260.hk',
'/data/NICER/preliminary/fits1/1706162208/auxil/ni1706162208_apid0260.hk',
'/data/NICER/preliminary/fits1/1706170609/auxil/ni1706170609_apid0260.hk',
'/data/NICER/preliminary/fits1/1706171350/auxil/ni1706171350_apid0260.hk'
]

for fn in filenames:
    tlist = []
    tlist.append(Table.read(fn,hdu=1))
    
    etable = vstack(tlist,metadata_conflicts='silent')
    #print(np.where(etable['GPS_SPS_TIME_VALID'] == 0))

    TIMEzerocount = len(np.where(etable['GPS_SPS_TIME_VALID'] == 0)[0])
    TIMEonecount = len(np.where(etable['GPS_SPS_TIME_VALID'] == 1)[0])
    TIMEdiffcount = len(np.where(etable['GPS_SPS_TIME_VALID'] > 1)[0])
    zerocount = len(np.where(etable['GPS_SPS_VALID'] == 0)[0])
    onecount = len(np.where(etable['GPS_SPS_VALID'] == 1)[0])

    idx = np.where(etable['GPS_SPS_GDOP'] != 0.0)[0]
    gdop = etable['GPS_SPS_GDOP'][idx].mean()

    sol = "{0}: SPS_VALID (0 = {1:5d}, 1= {2:5d}), SPS_TIME_VALID (0 = {3:5d}, 1 = {4:5d}, >1 = {5:5d}) <GDOP> = {6:.2f}".format(fn[48:71],  
zerocount, onecount, TIMEzerocount, TIMEonecount, TIMEdiffcount, gdop)
    
    print(sol)
    
