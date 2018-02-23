from __future__ import division, print_function
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.table import Table
from astropy import log
from scipy.interpolate import InterpolatedUnivariateSpline
from astropy.coordinates import GCRS, ITRS, EarthLocation, CartesianRepresentation, get_body_barycentric_posvel

from nicer.values import *

class LatLonInterp:
    '''LatLonInterp
    Class to allow interpolating lat, lon

    Initialize with a set of met, lat, lon

    Provides interpolated lat, long for any MET time during the ephemeris valid
    interval.
    '''
    def __init__(self, met,lat,lon):
        self.met = met
        log.info('MET Range {0} to {1}'.format(self.met.min(),self.met.max()))

        # Read lat and lon. These come with units (rad) provided by the FITS file
        self.lat = lat
        self.lon = lon

        self.earthloc = EarthLocation.from_geodetic(lon=self.lon,lat=self.lat)

        # Convert to cartesian for smooth interpolation, since the wrapping
        # of latitude gives interpolators fits.
        self.ecef_x = InterpolatedUnivariateSpline(self.met,self.earthloc.x,ext='const')
        self.ecef_y = InterpolatedUnivariateSpline(self.met,self.earthloc.y,ext='const')
        self.ecef_z = InterpolatedUnivariateSpline(self.met,self.earthloc.z,ext='const')

    def latlon(self, met):
        x = self.ecef_x(met)
        y = self.ecef_y(met)
        z = self.ecef_z(met)
        el = EarthLocation.from_geocentric(x,y,z,unit=u.m)
        return (el.lat, el.lon)
