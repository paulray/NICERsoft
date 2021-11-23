#!/usr/bin/env python
# Program: nicerfits2bin.py
# Authors: Paul S. Ray <paul.ray@nrl.navy.mil>
# Description:
# Reads a NICER clean FITS file of photon event times
# and generates a .bin and .inf file suitable for PRESTO analysis
# The .bin is just a binary file of double precision MJD times
from __future__ import division, print_function
import os,sys
import argparse
import numpy as np
from astropy import log
import astropy.units as u
import astropy.io.fits as pyfits
from astropy.time import Time
from astropy.table import Table
from nicer.values import *
from nicer.plotutils import choose_N
from subprocess import check_call

def runcmd(cmd):
    # CMD should be a list of strings since it is not processed by a shell
    log.info('CMD: '+" ".join(cmd))
    #log.info(cmd)
    check_call(cmd,env=os.environ)

desc="""Convert NICER events to PRESTO bin+inf files. """

parser = argparse.ArgumentParser(description=desc)
parser.add_argument("evfile", help="Name of event file to convert (must be barycentered)")
parser.add_argument("--dt", help="Time series bin size", type=float, default=1.0/128.0)
parser.add_argument("--ra",help="RA in HH:MM:SS.s",default="00:00:00")
parser.add_argument("--dec",help="DEC in DD:MM:SS.s", default="00:00:00")

parser.add_argument("--observer",help="Name of person analyzing data",default=None)
parser.add_argument("--force",help="Force processing even if not barycentered",default=False,action="store_true")
parser.add_argument("--search",help="Do an FFT and search the data after creating it",default=False,action="store_true")
parser.add_argument("--numharm",help="Number of harmonics to sum in search",default="2")
parser.add_argument("--zmax",help="Max z for acceleration search",default="50")
parser.add_argument("--flo",help="Lowest freq (of highest harmonic) to search",default="0.02")
args = parser.parse_args()
if args.observer is None:
    import getpass
    args.observer = getpass.getuser()

# Grab base filename for output
base = os.path.splitext(os.path.basename(args.evfile))[0]

etable = Table.read(args.evfile,hdu=1)
if etable.meta['TIMESYS'] != 'TDB' and not args.force:
    log.error('Event file must be barycentered!')
    sys.exit(1)
gtitable = Table.read(args.evfile,hdu='GTI')

epoch_met = gtitable['START'][0]
# WARNING: This loses precision! Should be done with astropy times
# Should make utility routine to convert FITS TIME column to astropy times properly
epoch_mjd = (etable.meta['MJDREFI'] + etable.meta['MJDREFF']
    + (etable.meta['TIMEZERO'] + epoch_met)/86400.0)

# Write event times to bin file
eventtimes = np.array(etable['TIME'],dtype=float)-epoch_met
log.info('Event times: {0} to {1}'.format(eventtimes.min(),eventtimes.max()))
eventtimes.tofile('{0}.events'.format(base))
#
nbins = choose_N((gtitable['STOP'][-1]-epoch_met)/args.dt)
log.info('Using {0} bins of width {1} seconds for total of {2} s'.format(nbins,
    args.dt,nbins*args.dt))

# Now write out .inf file
inf = open('%s.inf' % base,'w')
infstr = """ Data file name without suffix          =  {0}
 Telescope used                         =  NICER
 Instrument used                        =  XTI
 Object being observed                  =  {1}
 J2000 Right Ascension (hh:mm:ss.ssss)  =  {2}
 J2000 Declination     (dd:mm:ss.ssss)  =  {3}
 Data observed by                       =  NICER
 Epoch of observation (MJD)             =  {4:.9f}
 Barycentered?           (1=yes, 0=no)  =  1
 Number of bins in the time series      =  {5}
 Width of each time series bin (sec)    =  {6}
 Any breaks in the data? (1=yes, 0=no)  =  0
 Type of observation (EM band)          =  X-ray
 Field-of-view diameter (arcsec)        =  180.00
 Central energy (kev)                   =  4.2
 Energy bandpass (kev)                  =  3.8
 Data analyzed by                       =  {7}
 Any additional notes:
    None
""".format(base,etable.meta['OBJECT'], args.ra, args.dec, epoch_mjd, nbins,
    args.dt,args.observer)

inf.write(infstr)
inf.close()

# Now bin data and write .dat file
bins = np.arange(nbins+1,dtype=float)*args.dt
sums, edges = np.histogram(eventtimes, bins=bins)
dat = np.array(sums,np.float32)
dat.tofile('{0}.dat'.format(base))


if args.search:
    cmd = ['realfft', '{0}.dat'.format(base)]
    runcmd(cmd)
    cmd = ['accelsearch', '-zmax', args.zmax, '-flo',args.flo, '-numharm', args.numharm,'-sigma','1.5','-photon','{0}.fft'.format(base)]
    runcmd(cmd)
    # If candidates are found, fold the most significant one
    if os.path.isfile('{0}_ACCEL_{1}.cand'.format(base,args.zmax)):
        cmd = ['prepfold', '-noxwin', '-events', '-double', '-accelcand', '1', 
            '-accelfile', '{0}_ACCEL_{1}.cand'.format(base,args.zmax), 
            '{0}.events'.format(base)]
        runcmd(cmd)
