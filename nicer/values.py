import numpy as np
from astropy.time import Time
from os import path

# Base values for NICER MET and GPS time
# Both are counted as seconds since these epochs
MET0 = Time("2014-01-01T00:00:00.0", scale="utc")
GPS0 = Time("1980-01-06T00:00:00", scale="utc")

# Conversion from PI to keV (PI is in units of 10 eV)
PI_TO_KEV = 10.0 / 1000.0
KEV_TO_PI = 1000.0 / 10.0

# Bit assignments in the EVENT_FLAGS housekeeping
FLAG_FIRST_MPU = 2
FLAG_SLOW = 3
FLAG_FAST = 4
FLAG_SWTRIG = 5
FLAG_OVERSHOOT = 6
FLAG_UNDERSHOOT = 7

# IDs of the 4 detectors that were dead at launch
DEAD_DETS = set([11, 20, 22, 60])

# Array of DET_IDs that are used
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

datadir = path.join(path.split(__file__)[0], "../data")
