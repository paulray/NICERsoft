"""Interpolate the spacecraft latitude and longitude from tabulated values."""

import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.table import Table
from loguru import logger as log
from scipy.interpolate import InterpolatedUnivariateSpline
from astropy.coordinates import (
    GCRS,
    ITRS,
    EarthLocation,
    CartesianRepresentation,
    get_body_barycentric_posvel,
)

from nicer.values import *


class LatLonInterp:
    """Interpolate the spacecraft latitude and longitude at any time.

    The positions are converted to Earth-fixed Cartesian coordinates and
    interpolated with splines, which avoids problems with longitude wrapping.
    Outside the tabulated range, the first or last value is returned.

    Parameters
    ----------
    met : array_like
        Times of the tabulated positions, in MET seconds (for example
        ``mktable["TIME"]``).
    lat, lon : astropy.units.Quantity
        Latitudes and longitudes at those times (for example
        ``mktable["SAT_LAT"]`` and ``mktable["SAT_LON"]``).
    """

    def __init__(self, met, lat, lon):
        self.met = met
        log.info("MET Range {0} to {1}".format(self.met.min(), self.met.max()))

        # Read lat and lon. These come with units (rad) provided by the FITS file
        self.lat = lat
        self.lon = lon

        self.earthloc = EarthLocation.from_geodetic(lon=self.lon, lat=self.lat)

        # Convert to cartesian for smooth interpolation, since the wrapping
        # of latitude gives interpolators fits.
        self.ecef_x = InterpolatedUnivariateSpline(
            self.met, self.earthloc.x, ext="const"
        )
        self.ecef_y = InterpolatedUnivariateSpline(
            self.met, self.earthloc.y, ext="const"
        )
        self.ecef_z = InterpolatedUnivariateSpline(
            self.met, self.earthloc.z, ext="const"
        )

    def latlon(self, met):
        """Return the interpolated latitude and longitude.

        Parameters
        ----------
        met : float or array_like
            Time(s), in MET seconds.

        Returns
        -------
        lat : astropy.coordinates.Latitude
        lon : astropy.coordinates.Longitude
        """
        x = self.ecef_x(met)
        y = self.ecef_y(met)
        z = self.ecef_z(met)
        el = EarthLocation.from_geocentric(x, y, z, unit=u.m)
        return (el.lat, el.lon)
