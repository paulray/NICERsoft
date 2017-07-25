import matplotlib.pyplot as plt
import numpy as np
import argparse
from astropy import log
from astropy.table import Table, vstack
import astropy.io.fits as pyfits
import astropy.units as u
from astropy.time import Time
from os import path
from glob import glob
#from nicer.plotutils import *
from nicer.fitsutils import *
from nicer.values import *
import sys

class NicerFileSet:

    def __init__(self, args):
        log.info('Initializing the data object')
        self.args = args

        #Getting the names of the event files from obsdir
        self.evfiles = glob(path.join(self.args.obsdir,'xti/event_cl/ni*mpu?_cl.evt'))
        self.evfiles.sort()
        if len(self.evfiles) == 0:
            log.error("No event files found!")
            sys.exit(1)
        log.info('Found event files: {0}'.format("\n" + "    \n".join(self.evfiles)))

        # Get name of orbit file from obsdir
        try:
            self.args.orb = glob(path.join(self.args.obsdir,'auxil/ni*.orb'))[0]
        except:
            log.error("Orbit file not found!")
        log.info('Found the orbit file: {0}'.format(self.args.orb))

        # Get name of SPS HK file (apid0260)
        if self.args.sps is None:
            try:
                self.args.sps = glob(path.join(self.args.obsdir,'auxil/ni*_apid0260.hk'))[0]
            except:
                self.args.sps = None

        # Get name of MPU housekeeping files
        self.hkfiles = glob(path.join(self.args.obsdir,'xti/hk/ni*.hk'))
        self.hkfiles.sort()
        log.info('Found the MPU housekeeping files: {0}'.format("\n"+"\t\n".join(self.hkfiles)))

        # Get name of filter (.mkf) file
        self.mkfiles = glob(path.join(args.obsdir,'auxil/ni*.mkf'))[0]

        #Compiling Event Data
        self.getgti()
        if self.args.useftools:
            self.etable = filtallandmerge_ftools(self.evfiles,workdir=None)
        else:
            self.etable = self.createetable()
        self.sortmet()
        self.makebasename()

        #--------------------Editing / Filtering the event data-----------------
        #filtering out chosen IDS
        if self.args.mask is not None:
            log.info('Masking IDS')
            for id in self.args.mask:
                self.etable = self.etable[np.where(self.etable['DET_ID'] != id)]
        # If there are no PI columns, add them with approximate calibration
        if self.args.pi or not ('PI' in self.etable.colnames):
            log.info('Adding PI')
            calfile = path.join(datadir,'gaincal_linear.txt')
            pi = calc_pi(self.etable,calfile)
            self.etable['PI'] = pi

        if self.args.applygti is not None:
            g = Table.read(self.args.applygti)
            log.info('Applying external GTI from {0}'.format(self.args.applygti))
            g['DURATION'] = g['STOP']-g['START']
            # Only keep GTIs longer than 16 seconds
            g = g[np.where(g['DURATION']>16.0)]
            log.info('Applying external GTI')
            print(g)
            self.etable = apply_gti(self.etable,g)
            # Replacing this GTI does not work. It needs to be ANDed with the existing GTI
            self.etable.meta['EXPOSURE'] = g['DURATION'].sum()
            self.gtitable = g
        log.info('Exposure : {0:.2f}'.format(self.etable.meta['EXPOSURE']))
        #-----------------------------------------------------------------------

        #Compiling HK Data
        self.getmk()
        self.hkshootrate()
        self.eventundershootrate()
        #self.eventovershootrate()


    def createetable(self):
        log.info('Reading files')
        tlist = []
        for fn in self.evfiles:
            log.info('Reading file {0}'.format(fn))
            tlist.append(Table.read(fn,hdu=1))
            if len(tlist[0]) > 3000000:
                log.error('There is too much data to handle. Not processing...')
                sys.exit(3)
        log.info('Concatenating files')
        if len(tlist) == 1:
            self.etable = tlist[0]
        else:
            self.etable = vstack(tlist,metadata_conflicts='silent')
        del tlist

        return self.etable

    def sortmet(self):
        # Change TIME column name to MET to reflect what it really is
        self.etable.columns['TIME'].name = 'MET'
        # Update exposure to be sum of GTI durations
        self.etable.meta['EXPOSURE'] = self.gtitable['DURATION'].sum()

        # Sort table by MET
        self.etable.sort('MET')
        log.info("Event MET Range : {0} to {1}".format(self.etable['MET'].min(),
            self.etable['MET'].max(), self.etable['MET'].max()-self.etable['MET'].min()))
        log.info("TSTART {0}  TSTOP {1} (Span {2} seconds)".format(self.etable.meta['TSTART'],
            self.etable.meta['TSTOP'], self.etable.meta['TSTOP']-self.etable.meta['TSTART'] ))
        log.info("DATE Range {0} to {1}".format(self.etable.meta['DATE-OBS'],
            self.etable.meta['DATE-END']))

        if self.args.object is not None:
            self.etable.meta['OBJECT'] = self.args.object

    def getmk(self):
        #Making the MK Table
        log.info('Getting MKTable')
        if len(self.mkfiles) > 0:
            self.mktable = Table.read(self.mkfiles,hdu=1)
        else:
            self.mktable = None

    def hkshootrate(self):
        # getting the overshoot and undershoot rate from HK files.  Times are hkmet
        log.info('Getting HKP overshoot and undershoot rates')
        if len(self.hkfiles) > 0:
            log.info('Reading '+self.hkfiles[0])
            hdulist = pyfits.open(self.hkfiles[0])
            td = hdulist[1].data
            self.hkmet = td['TIME']
            log.info("HK MET Range {0} to {1} (Span = {2:.1f} seconds)".format(self.hkmet.min(),
                self.hkmet.max(),self.hkmet.max()-self.hkmet.min()))
            self.overshootrate = td['MPU_OVER_COUNT'].sum(axis=1)
            self.undershootrate = td['MPU_UNDER_COUNT'].sum(axis=1)
            nresets = td['MPU_UNDER_COUNT'].sum(axis=0)

            for fn in self.hkfiles[1:]:
                log.info('Reading '+fn)
                hdulist = pyfits.open(fn)
                mytd = hdulist[1].data
                mymet = mytd['TIME']
                myovershootrate = mytd['MPU_OVER_COUNT'].sum(axis=1)
                myundershootrate = mytd['MPU_UNDER_COUNT'].sum(axis=1)
                # If time axis is bad, skip this MPU.
                # Should fix this!
                if not np.all(mymet == self.hkmet):
                    log.error('TIME axes are not compatible')
                else:
                    self.overshootrate += myovershootrate
                    self.undershootrate += myundershootrate
                myreset = mytd['MPU_UNDER_COUNT'].sum(axis=0)
                nresets = np.append(nresets,myreset)
            del hdulist
            self.reset_rates = nresets / np.float(self.etable.meta['EXPOSURE'])

        else:
            hkmet = None
            overshootrate = None
            undershootrate = None
            nresets = calc_nresets(self.etable, IDS)
            reset_rates = nresets/self.etable.meta['EXPOSURE']

    def eventovershootrate(self):
        #Get a table of OVERSHOOT ONLY events
        log.info('Getting overshoot only events')
        self.ovstable = get_eventovershoots_ftools(self.evfiles,workdir=None)
        hkmetbins = np.append(self.hkmet,(self.hkmet[-1]+self.hkmet[1]-self.hkmet[0]))
        bins = np.arange(self.ovstable['TIME'][0], self.ovstable['TIME'][-1]+self.ovstable['TIME'][1]-self.ovstable['TIME'][0], self.args.lcbinsize)
        self.eventovershoot, edges = np.histogram(self.ovstable['TIME'],bins)

    def eventundershootrate(self):
        log.info('Getting event undershoot rates')
        # Define bins for hkmet histogram
        hkmetbins = np.append(self.hkmet,(self.hkmet[-1]+self.hkmet[1]-self.hkmet[0]))
        #Both under and overshoot
        idx = np.logical_and(self.etable['EVENT_FLAGS'][:,FLAG_UNDERSHOOT]==True,
                             self.etable['EVENT_FLAGS'][:,FLAG_OVERSHOOT]==True)
        self.bothrate, edges = np.histogram(self.etable['MET'][idx],hkmetbins)

        #Just undershoot
        idx = np.logical_and(self.etable['EVENT_FLAGS'][:,FLAG_UNDERSHOOT]==True,
                             self.etable['EVENT_FLAGS'][:,FLAG_OVERSHOOT]==False)
        self.eventundershoot, edges = np.histogram(self.etable['MET'][idx],hkmetbins)

        #Just overshoot
        idx = np.logical_and(self.etable['EVENT_FLAGS'][:,FLAG_UNDERSHOOT]==False,
                             self.etable['EVENT_FLAGS'][:,FLAG_OVERSHOOT]==True)
        self.eventovershoot, edges = np.histogram(self.etable['MET'][idx],hkmetbins)

        del idx, edges

    def getlatlon(self):
        self.lat = self.mktable['SAT_LAT']
        self.lon = self.mktable['SAT_LON']
        self.sun = self.mktable['SUNSHINE']

    def writeovsfile(self):
        # Write overshoot and undershoot rates to file for filtering
        log.info('Writing over/undershoot rates')
        tcol = pyfits.Column(name='TIME',unit='S',array=self.hkmet,format='D')
        ocol = pyfits.Column(name='HK_OVERSHOOT',array=self.overshootrate,format='D')
        ucol = pyfits.Column(name='HK_UNDERSHOOT',array=self.undershootrate,format='D')
        eocol = pyfits.Column(name='EV_OVERSHOOT',array=self.eventovershoot,format='D')
        eucol = pyfits.Column(name='EV_UNDERSHOOT',array=self.eventundershoot,format='D')
        bothcol = pyfits.Column(name='EV_BOTH',array=self.bothrate,format='D')
        lat = pyfits.Column(name='SAT_LAT',array=self.lat,format='D')
        lon = pyfits.Column(name='SAT_LON',array=self.lon,format='D')
        sun = pyfits.Column(name='SUNSHINE',array=self.sun,format='D')
        ovhdu = pyfits.BinTableHDU.from_columns([tcol,ocol,ucol, eocol, eucol, bothcol, lat, lon, sun], name='HKP')
        ovhdu.writeto("{0}.ovs".format(self.basename),overwrite=True,checksum=True)

    def getgti(self):
            # Read the GTIs from the first event FITS file
        self.gtitable = Table.read(self.evfiles[0],hdu=2)
        log.info('Got the good times from GTI')
        self.gtitable['DURATION'] = self.gtitable['STOP']- self.gtitable['START']
        # Only keep GTIs longer than 16 seconds
        idx = np.where(self.gtitable['DURATION']>16.0)[0]
        self.gtitable = self.gtitable[idx]
        print(self.gtitable)

    def makebasename(self):
        bn = path.basename(self.evfiles[0]).split('_')[0]
        log.info('OBS_ID {0}'.format(self.etable.meta['OBS_ID']))
        if self.etable.meta['OBS_ID'].startswith('000000'):
            log.info('Overwriting OBS_ID with {0}'.format(bn))
            self.etable.meta['OBS_ID'] = bn
            self.etable.meta['OBS_ID'] = bn
        if self.args.basename is None:
            self.basename = '{0}'.format(bn)
        else:
            self.basename = self.args.basename
