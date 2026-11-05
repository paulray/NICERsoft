"""ISS position from the GPS SPS housekeeping data."""

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


class SPS:
    """Interpolate the ISS latitude and longitude from GPS SPS housekeeping data.

    Rows where both latitude and longitude are zero are ignored. Outside the
    time range of the file, the first or last position is returned.

    Parameters
    ----------
    spsname : str
        SPS housekeeping file (``*_apid0260.hk``).
    """

    def __init__(self, spsname):
        spstab = Table.read(spsname, hdu=1)
        # Apply TIMEZERO if needed
        if "TIMEZERO" in spstab.meta:
            log.info(
                "Applying TIMEZERO of {0} to spstab".format(spstab.meta["TIMEZERO"])
            )
            spstab["TIME"] += spstab.meta["TIMEZERO"]
            spstab.meta["TIMEZERO"] = 0.0

        idx = np.where(
            np.logical_not(
                np.logical_and(
                    spstab["GPS_SPS_LAT"] == 0.0, spstab["GPS_SPS_LON"] == 0.0
                )
            )
        )[0]
        spstab = spstab[idx]
        spstab.sort("TIME")
        self.met = spstab["TIME"]
        log.info("SPS Range {0} to {1}".format(self.met.min(), self.met.max()))

        # Read lat and lon. These come with units (rad) provided by the FITS file
        self.spslat = spstab["GPS_SPS_LAT"]
        self.spslon = spstab["GPS_SPS_LON"]

        self.earthloc = EarthLocation.from_geodetic(lon=self.spslon, lat=self.spslat)

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
