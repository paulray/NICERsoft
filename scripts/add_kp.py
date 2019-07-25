#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
import argparse
from astropy import log
from os import path, environ
from glob import glob
import astropy.io.fits as pyfits
from astropy.time import Time
import astropy.units as u

from nicer.values import *

def read_kpfiles_potsdam():
    kpfiles = glob(path.join(datadir,'KP-potsdam','kp*.tab'))
    kpfiles.sort()
    hours = np.array([0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0])*u.hour

    kpdict = {            '0o': 0.0, '0+': 0.33, 
              '1-': 0.66, '1o': 1.0, '1+': 1.33,
              '2-': 1.66, '2o': 2.0, '2+': 2.33,
              '3-': 2.66, '3o': 3.0, '3+': 3.33,
              '4-': 3.66, '4o': 4.0, '4+': 4.33,
              '5-': 4.66, '5o': 5.0, '5+': 5.33,
              '6-': 5.66, '6o': 6.0, '6+': 6.33,
              '7-': 6.66, '7o': 7.0, '7+': 7.33,
              '8-': 7.66, '8o': 8.0, '8+': 8.33,
              '9-': 8.66, '9o': 9.0, '9+': 9.33
    }

    kpmets = []
    kpvals = []
    for kpname in kpfiles:
        with open(kpname,'r') as file:
            for line in file:
                if len(line.strip()) <= 0:
                    continue
                cols = line.split()
                if len(cols[0]) == 6:
                    t0 = Time("20{0}-{1}-{2}T00:00:00".format(cols[0][0:2],cols[0][2:4],cols[0][4:6]),format='isot',scale='utc')
                    for kpstring,h in zip(cols[1:9],hours):
                        met = ((t0+h)-MET0).to(u.s).value
                        kpval = kpdict[kpstring]
                        kpmets.append(met)
                        kpvals.append(kpval)
    kpmets = np.array(kpmets)
    kpvals = np.array(kpvals)
    log.info('KP table goes from {0} to {1}'.format(kpmets[0],kpmets[-1]))

    return (kpmets,kpvals)

def read_kp_GOF(fname):
    dat = pyfits.getdata(fname,extname='KP')
    return (dat['TIME'], dat['KP'])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Add Kp values as a column to an existing MKF file')
    parser.add_argument("mkf", help="Input MKF file process (will be modified)", nargs='+')
    parser.add_argument("--potsdam",help="Use Potsdam Kp values from NICERsoft data (default is NOA data from NICER GOF)",action="store_true",default=False)
    parser.add_argument("--nearest",help="Use nearest Kp value, instead of the default linear interpolation",action="store_true",default=False)
    parser.add_argument("--kp", help="Path to NICER GOF kp.fits file (default $CALDB/data/gen/pcf/kp.fits). Only used if not --potsdam.", default=None)
    args = parser.parse_args()

    # Read the Pottsdam Kp index data files
    if args.potsdam:
        log.info('Reading Potsdam Kp files...')
        kpmets, kpvals = read_kpfiles_potsdam()
    else:
        log.info('Reading NICER GOF KP files')
        if args.kp is None:
            args.kp = path.join(environ['CALDB'],"data/gen/pcf/kp.fits")
            if not path.isfile(args.kp):
                log.error("KP file {0} does not exist.  Please download https://heasarc.gsfc.nasa.gov/FTP/caldb/data/gen/pcf/kp.fits",args.kp)
                raise RuntimeError
        kpmets, kpvals = read_kp_GOF(args.kp)
                

    for mkf in args.mkf:
        # Read the METs from the MKF file
        with pyfits.open(mkf) as hdulist:
            log.info('Processing {0}'.format(mkf))
            #mkfdat = pyfits.getdata(mkf,extname='PREFILTER')
            mkfmets = hdulist[1].data['TIME']

            # Look up, or interpolate, all the Kp values
            if args.nearest:
                kp = np.array([kpvals[kpmets.searchsorted(m)] for m in mkfmets])
            else:   
                kp = np.interp(mkfmets, kpmets, kpvals)

            names     = hdulist[1].data.names

            if 'KP' in names:
                log.info('Found KP column, overwriting...')
                hdulist[1].data.field('KP')[:] = kp
            else:
                log.info('Adding KP column to MKF file...')
                col   = pyfits.Column(name='KP',format='E',array=kp)
                tbhdu = pyfits.BinTableHDU.from_columns(hdulist[1].data.columns+col,
                        header=hdulist[1].header)
                hdulist[1] = tbhdu

            hdulist.writeto(mkf,overwrite=True)
            hdulist.close()
