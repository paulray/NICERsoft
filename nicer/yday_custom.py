"""An astropy time format for year and day-of-year strings.

Importing this module registers the ``yday_custom`` format with
:class:`astropy.time.Time`.
"""

from astropy.time import TimeISO


class TimeYearDayTimeCustom(TimeISO):
    """Year, day-of-year and time as ``<YYYY>-<DOY>T<HH>:<MM>:<SS>``.

    The day of year (DOY) goes from 001 to 365 (366 in leap years). For
    example, ``2000-001T00:00:00`` is midnight on January 1, 2000.
    The allowed subformats are:

    - ``date_hms``: date + hours, minutes, seconds
    - ``date_hm``: date + hours, minutes
    - ``date``: date only
    """

    name = "yday_custom"
    subfmts = (
        (
            "date_hms",
            "%Y-%jT%H:%M:%S",
            "{year:d}-{yday:03d}T{hour:02d}:{min:02d}:{sec:02d}",
        ),
        ("date_hm", "%Y-%jT%H:%M", "{year:d}-{yday:03d}T{hour:02d}:{min:02d}"),
        ("date", "%Y-%j", "{year:d}-{yday:03d}"),
    )
