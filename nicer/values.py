"""Constants used throughout NICERsoft: time epochs, energy conversions,
standard energy bands, ``EVENT_FLAGS`` bits and detector IDs."""

import numpy as np
from astropy.time import Time
from os import path

# Base values for NICER MET and GPS time
# Both are counted as seconds since these epochs
#: Epoch of NICER Mission Elapsed Time (MET), which counts seconds since 2014-01-01 UTC.
MET0 = Time("2014-01-01T00:00:00.0", scale="utc")
#: Epoch of GPS time, which counts seconds since 1980-01-06 UTC.
GPS0 = Time("1980-01-06T00:00:00", scale="utc")

# Conversion from PI to keV (PI is in units of 10 eV)
#: Multiply a PI channel by this to get energy in keV (PI is in units of 10 eV).
PI_TO_KEV = 10.0 / 1000.0
#: Multiply an energy in keV by this to get the PI channel.
KEV_TO_PI = 1000.0 / 10.0

# Standard energy bands for RPP catalog analysis
#: Lower edge (keV) of the soft band used in the RPP catalog analysis.
SOFT_EMIN = 0.4
#: Upper edge (keV) of the soft band.
SOFT_EMAX = 2.0
#: Lower edge (keV) of the hard band used in the RPP catalog analysis.
HARD_EMIN = 2.0
#: Upper edge (keV) of the hard band.
HARD_EMAX = 10.0

# Bit assignments in the EVENT_FLAGS housekeeping
#: Index of the "first MPU" bit in the ``EVENT_FLAGS`` column.
FLAG_FIRST_MPU = 2
#: Index of the bit set when the slow chain triggered.
FLAG_SLOW = 3
#: Index of the bit set when the fast chain triggered.
FLAG_FAST = 4
#: Index of the software (forced) trigger bit.
FLAG_SWTRIG = 5
#: Index of the overshoot bit.
FLAG_OVERSHOOT = 6
#: Index of the undershoot bit.
FLAG_UNDERSHOOT = 7

#: IDs of the 4 detectors that were dead at launch.
DEAD_DETS = {11, 20, 22, 60}

# Array of DET_IDs that are used
#: The 56 ``DET_ID`` values (MPU number times 10 plus detector number, for MPUs 0-6).
IDS = np.array(
    [
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        30,
        31,
        32,
        33,
        34,
        35,
        36,
        37,
        40,
        41,
        42,
        43,
        44,
        45,
        46,
        47,
        50,
        51,
        52,
        53,
        54,
        55,
        56,
        57,
        60,
        61,
        62,
        63,
        64,
        65,
        66,
        67,
    ]
)

#: Path to the NICERsoft ``data`` directory.
datadir = path.join(path.split(__file__)[0], "../data")
