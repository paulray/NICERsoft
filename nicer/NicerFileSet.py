from __future__ import (print_function, division, unicode_literals, absolute_import)
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
from nicer.plotutils import *
from nicer.fitsutils import *
from nicer.values import *

class NicerFileSet:
    def __init__(self, args):
        log.info('Initializing the data object')
        self.args = args
        
        #Getting the names of the event files from obsdir
        self.evfiles = glob(path.join(self.args.obsdir,'xti/event_cl/ni*mpu?_cl.evt*'))
        self.evfiles.sort()
        if len(self.evfiles) == 0:
            log.error("No event files found!")
            sys.exit(1)
        log.info('Found event files: {0}'.format("\n" + "    \n".join(self.evfiles)))

        # Get name of orbit file from obsdir
        try:
            self.args.orb = glob(path.join(self.args.obsdir,'auxil/ni*.orb*'))[0]
        except:
            log.error("Orbit file not found!")
        log.info('Found the orbit file: {0}'.format(self.args.orb))

        # Get name of SPS HK file (apid0260)
        if self.args.sps is None:
            try:
                self.args.sps = glob(path.join(self.args.obsdir,'auxil/ni*_apid0260.hk*'))[0]
            except:
                self.args.sps = None

        # Get name of MPU housekeeping files
        self.hkfiles = glob(path.join(self.args.obsdir,'xti/hk/ni*.hk*'))
        self.hkfiles.sort()
        log.info('Found the MPU housekeeping files: {0}'.format("\n"+"\t\n".join(self.hkfiles)))

        # Get name of filter (.mkf) file
        self.mkfiles = glob(path.join(args.obsdir,'auxil/ni*.mkf*'))[0]

        #Compiling Event Data
        self.getgti()
        if self.args.useftools:
            self.etable = filtallandmerge_ftools(self.evfiles,workdir=None)
        else:
            self.etable = self.createetable()
        self.sortmet()
        self.makebasename()

        if args.applygti is not None:
            g = Table.read(args.applygti)
            log.info('Applying external GTI from {0}'.format(args.applygti))
            g['DURATION'] = g['STOP']-g['START']
            # Only keep GTIs longer than 16 seconds
            g = g[np.where(g['DURATION']>16.0)]
            log.info('Applying external GTI')
            print(g)
            etable = apply_gti(etable,g)
            # Replacing this GTI does not work. It needs to be ANDed with the existing GTI
            etable.meta['EXPOSURE'] = g['DURATION'].sum()
            gtitable = g
    
        #Compiling HK Data
        self.getmk()
        self.hkshootrate()
        self.geteventshoots()
        #self.eventovershootrate()

        #Compiling other data
        self.getlatlon()
        
    def createetable(self):
        log.info('Reading files')
        tlist = []
        for fn in self.evfiles:
            log.info('Reading file {0}'.format(fn))
            tlist.append(Table.read(fn,hdu=1))
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
            self.hkovershoots = td['MPU_OVER_COUNT'].sum(axis=1)
            self.hkundershoots = td['MPU_UNDER_COUNT'].sum(axis=1)
            nresets = td['MPU_UNDER_COUNT'].sum(axis=0)

            for fn in self.hkfiles[1:]:
                log.info('Reading '+fn)
                hdulist = pyfits.open(fn)
                mytd = hdulist[1].data
                mymet = mytd['TIME']
                myhkovershoots= mytd['MPU_OVER_COUNT'].sum(axis=1)
                myhkundershoots = mytd['MPU_UNDER_COUNT'].sum(axis=1)
                # If time axis is bad, skip this MPU.
                # Should fix this!
                if not np.all(mymet == self.hkmet):
                    log.error('TIME axes are not compatible')
                else:
                    self.hkovershoots += myhkovershoots
                    self.hkundershoots += myhkundershoots
                myreset = mytd['MPU_UNDER_COUNT'].sum(axis=0)
                nresets = np.append(nresets,myreset)
            del hdulist
            self.reset_rates = nresets / np.float(self.etable.meta['EXPOSURE'])

        else:
            hkmet = None
            hkovershoots = None
            hkundershoots = None
            nresets = calc_nresets(self.etable, IDS)
            reset_rates = nresets/self.etable.meta['EXPOSURE']

    def eventovershootrate(self):
        #Get a table of OVERSHOOT ONLY events
        log.info('Getting overshoot only events')
        self.ovstable = get_eventovershoots_ftools(self.evfiles,workdir=None)
        hkmetbins = np.append(self.hkmet,(self.hkmet[-1]+self.hkmet[1]-self.hkmet[0]))
        bins = np.arange(self.ovstable['TIME'][0], self.ovstable['TIME'][-1]+self.ovstable['TIME'][1]-self.ovstable['TIME'][0], self.args.lcbinsize)
        self.eventovershoot, edges = np.histogram(self.ovstable['TIME'],bins)

    def geteventshoots(self):
        log.info('Getting event undershoot rates')
        # Define bins for hkmet histogram
        hkmetbins = np.append(self.hkmet,(self.hkmet[-1]+self.hkmet[1]-self.hkmet[0]))
        
        #Both under and overshoot
        idx = np.logical_and(self.etable['EVENT_FLAGS'][:,FLAG_UNDERSHOOT]==True,
                             self.etable['EVENT_FLAGS'][:,FLAG_OVERSHOOT]==True)
        self.eventbothshoots, edges = np.histogram(self.etable['MET'][idx],hkmetbins)

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

    def writeovsfile(self, badlightcurve):
        # Write overshoot and undershoot rates to file for filtering
        log.info('Writing over/undershoot rates')
        tcol = pyfits.Column(name='TIME',unit='S',array=self.hkmet,format='D')
        ocol = pyfits.Column(name='HK_OVERSHOOT',array=self.hkovershoots,format='D')
        ucol = pyfits.Column(name='HK_UNDERSHOOT',array=self.hkundershoots,format='D')
        eocol = pyfits.Column(name='EV_OVERSHOOT',array=self.eventovershoot,format='D')
        eucol = pyfits.Column(name='EV_UNDERSHOOT',array=self.eventundershoot,format='D')
        bothcol = pyfits.Column(name='EV_BOTH',array=self.eventbothshoots,format='D')

        if badlightcurve is not None:
            badcol = pyfits.Column(name='BAD_LC', array=badlightcurve, format='D')
            ovhdu = pyfits.BinTableHDU.from_columns([tcol,ocol,ucol, eocol, eucol, bothcol, badcol], name='HKP')
            ovhdu.writeto("{0}.ovs".format(self.basename),overwrite=True,checksum=True)
        else:
            ovhdu = pyfits.BinTableHDU.from_columns([tcol,ocol,ucol, eocol, eucol, bothcol], name='HKP')
            ovhdu.writeto("{0}.ovs".format(self.basename),overwrite=True,checksum=True)
        return

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
