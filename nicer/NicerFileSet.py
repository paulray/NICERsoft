from __future__ import (print_function, division, unicode_literals, absolute_import)
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
from nicer.latloninterp import LatLonInterp
import sys

class NicerFileSet:
    def __init__(self, args):
        log.info('Initializing the data object')
        self.args = args

        #Getting the names of the event files from obsdir
        self.evfiles = glob(path.join(self.args.obsdir,'xti/event_cl/ni*mpu?_cl.evt*'))
        self.evfiles.sort()
        if len(self.evfiles) == 0:
            log.error("No event files found!")
            raise Exception('No event files found!')
        log.info('Found clean event files: {0}'.format("\n" + "    \n".join(self.evfiles)))

        self.ufafiles = glob(path.join(self.args.obsdir,'xti/event_cl/ni*mpu?_ufa.evt*'))
        self.ufafiles.sort()
        log.info('Found unfiltered event files: {0}'.format("\n" + "    \n".join(self.ufafiles)))

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
        self.mkfile = glob(path.join(args.obsdir,'auxil/ni*.mkf*'))[0]
        self.mktable = Table.read(self.mkfile,hdu=1)

        # Make lat, lon interpolater from mktable
        self.llinterp = LatLonInterp(self.mktable['TIME'], self.mktable['SAT_LAT'], self.mktable['SAT_LON'])

        #Compiling Event Data
        self.getgti()
        if self.args.useftools:
            self.etable = filtallandmerge_ftools(self.evfiles,workdir=None)
        else:
            self.etable = self.createetable()
        if len(self.etable) == 0:
            log.error("No events in etable! Aborting")
            raise Exception('No events in etable!')
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
            self.etable = apply_gti(self.etable,g)
            # Replacing this GTI does not work. It needs to be ANDed with the existing GTI
            self.etable.meta['EXPOSURE'] = g['DURATION'].sum()
            self.gtitable = g

        if args.gtirows is not None:
            log.info('Apply gti rows {}'.format(args.gtirows))
            g = self.gtitable[args.gtirows]
            print(g)
            self.etable = apply_gti(self.etable,g)
            self.gtitable = g

        #Compiling HK Data
        if args.extraphkshootrate:
            self.quickhkshootrate()
        else:
            self.hkshootrate()
        self.geteventshoots()

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

    def quickhkshootrate(self):
        'Compute HK shoot rates from a single MPU*7 instead of trying to sum all MPUs. This avoids errors when time axes are not compatible'
        log.info('Getting HKP quick overshoot and undershoot rates from one MPU')
        if len(self.hkfiles) > 0:
            log.info('Reading '+self.hkfiles[4])
            hdulist = pyfits.open(self.hkfiles[4])
            td = hdulist[1].data
            self.hkmet = td['TIME']
            log.info("HK MET Range {0} to {1} (Span = {2:.1f} seconds)".format(self.hkmet.min(),
                self.hkmet.max(),self.hkmet.max()-self.hkmet.min()))
            # Only read one MPU so multiply by 7 to extrapolate to full rates
            self.hkovershoots = td['MPU_OVER_COUNT'].sum(axis=1)*7
            self.hkundershoots = td['MPU_UNDER_COUNT'].sum(axis=1)*7
            # Here we don't get reset rates for all MPUS
            nresets = td['MPU_UNDER_COUNT'].sum(axis=0)
            del hdulist
            self.reset_rates = None

            timecol = pyfits.Column(name='HKTIME',array=self.hkmet, format = 'D')
            ovscol = pyfits.Column(name = 'HK_OVERSHOOT', array = self.hkovershoots, format = 'D')
            undcol = pyfits.Column(name = 'HK_UNDERSHOOT',array = self.hkundershoots, format = 'D')

            self.hkshoottable = Table([self.hkmet, self.hkovershoots, self.hkundershoots],names = ('HKMET', 'HK_OVERSHOOTS', 'HK_UNDERSHOOTS'))
        else:
            hkmet = None
            hkovershoots = None
            hkundershoots = None
            nresets = calc_nresets(self.etable, IDS)
            reset_rates = nresets/self.etable.meta['EXPOSURE']

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

                if not np.array_equal(mymet, self.hkmet):
                    log.error('Time axes are not compatible')
                else:
                    self.hkovershoots += myhkovershoots
                    self.hkundershoots += myhkundershoots

                myreset = mytd['MPU_UNDER_COUNT'].sum(axis=0)
                nresets = np.append(nresets,myreset)
            del hdulist
            self.reset_rates = nresets / np.float(self.etable.meta['EXPOSURE'])

            timecol = pyfits.Column(name='HKTIME',array=self.hkmet, format = 'D')
            ovscol = pyfits.Column(name = 'HK_OVERSHOOT', array = self.hkovershoots, format = 'D')
            undcol = pyfits.Column(name = 'HK_UNDERSHOOT',array = self.hkundershoots, format = 'D')

            self.hkshoottable = Table([self.hkmet, self.hkovershoots, self.hkundershoots],names = ('HKMET', 'HK_OVERSHOOTS', 'HK_UNDERSHOOTS'))

        else:
            hkmet = None
            hkovershoots = None
            hkundershoots = None
            nresets = calc_nresets(self.etable, IDS)
            reset_rates = nresets/self.etable.meta['EXPOSURE']

    def geteventshoots(self):
        log.info('Getting event shoot rates')
        # Define bins for hkmet histogram
        hkmetbins = np.append(self.hkmet,(self.hkmet[-1]+self.hkmet[1]-self.hkmet[0]))

        if self.args.eventshootrate or self.args.writebkf:
            #Both under and overshoot
            if len(self.ufafiles) == 0:
                filelist = self.evfiles
            else:
                filelist = self.ufafiles
            etable = get_eventbothshoots_ftools(filelist,workdir=None)
            self.eventbothshoots, edges = np.histogram(etable['TIME'],hkmetbins)
            #Just overshoot
            etable = get_eventovershoots_ftools(filelist,workdir=None)
            self.eventovershoots, edges = np.histogram(etable['TIME'],hkmetbins)

            self.eventshoottable = Table([self.hkmet, self.eventovershoots, self.eventbothshoots],names = ('HKMET', 'EVENT_OVERSHOOTS', 'EVENT_BOTHSHOOTS'))
            del etable
        else:
            self.eventbothshoots = None
            self.eventovershoots = None

        # Don't compute this unless specifically requested, because it can be slow
        if self.args.eventshootrate:
            etable = get_eventundershoots_ftools(filelist,workdir=None)
            self.eventundershoots, edges = np.histogram(etable['TIME'],hkmetbins)
            self.eventshoottable = Table([self.hkmet, self.eventovershoots, self.eventundershoots, self.eventbothshoots],names = ('HKMET', 'EVENT_OVERSHOOTS', 'EVENT_UNDERSHOOTS', 'EVENT_BOTHSHOOTS'))
            del etable

        else:
            self.eventundershoots = None

    def writebkffile(self):
        # Write useful rates  to file for filtering

        # Define bins for hkmet histogram
        hkmetbins = np.append(self.hkmet,(self.hkmet[-1]+self.hkmet[1]-self.hkmet[0]))

        # Extract bad ratio events and bin onto hkmet bins
        if len(self.ufafiles) == 0:
            badtable = get_badratioevents_ftools(self.evfiles,workdir=None)
        else:
            badtable = get_badratioevents_ftools(self.ufafiles,workdir=None)
        badlightcurve = np.histogram(badtable['TIME'], hkmetbins)[0]
        badlightcurve = np.array(badlightcurve,dtype=np.float)
        # Really should convolve in GTI segments!
        kernel = np.ones(32)/32.0
        badlightcurve = np.convolve(badlightcurve,kernel,mode='same')

        log.info('Writing over/undershoot rates')
        tcol = pyfits.Column(name='TIME',unit='S',array=self.hkmet,format='D')
        ocol = pyfits.Column(name='HK_OVERSHOOT',array=self.hkovershoots,format='D')
        ucol = pyfits.Column(name='HK_UNDERSHOOT',array=self.hkundershoots,format='D')
        eocol = pyfits.Column(name='EV_OVERSHOOT',array=self.eventovershoots,format='D')
        bothcol = pyfits.Column(name='EV_BOTH',array=self.eventbothshoots,format='D')
        badcol = pyfits.Column(name='BAD_RATIO', array=badlightcurve, format='D')

        lat, lon = self.llinterp.latlon(self.hkmet)
        latcol = pyfits.Column(name='LAT', array=lat, format='D')
        loncol = pyfits.Column(name='LON', array=lon, format='D')

        ovhdu = pyfits.BinTableHDU.from_columns([tcol,ocol,ucol, eocol, bothcol,
            badcol,latcol,loncol], name='HKP')
        ovhdu.header['TIMESYS'] = self.etable.meta['TIMESYS']
        ovhdu.header['TIMEREF'] = self.etable.meta['TIMEREF']
        ovhdu.header['MJDREFI'] = self.etable.meta['MJDREFI']
        ovhdu.header['MJDREFF'] = self.etable.meta['MJDREFF']
        ovhdu.header['TIMEZERO'] = self.etable.meta['TIMEZERO']
        ovhdu.header['TIMEUNIT'] = self.etable.meta['TIMEUNIT']
        ovhdu.header['TSTART'] = self.etable.meta['TSTART']
        ovhdu.header['TSTOP'] = self.etable.meta['TSTOP']
        ovhdu.writeto("{0}.bkf".format(self.basename),overwrite=True,checksum=True)
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
