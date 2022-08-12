from astropy.time import TimeISO


class TimeYearDayTimeCustom(TimeISO):
    """
    Year, day-of-year and time as "<YYYY>-<DOY>T<HH>:<MM>:<SS.sss...>".
    The day-of-year (DOY) goes from 001 to 365 (366 in leap years).
    For example, 2000-001T00:00:00.000 is midnight on January 1, 2000.
    The allowed subformats are:
    - 'date_hms': date + hours, mins, secs (and optional fractional secs)
    - 'date_hm': date + hours, mins
    - 'date': date
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
