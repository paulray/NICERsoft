#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as pyfits
import astropy.stats
import astropy.units as u
import argparse
from pint.eventstats import z2m, hm, sf_z2m, sf_hm, sig2sigma, h2sig
from nicer.values import *
from nicer.fourier import *
import yaml
import os.path


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="Do RPP catalog pulsed spectral analysis on a merged event file.",
)
parser.add_argument(
    "evtfile", help="Input event file to process. Must have PULSE_PHASE column"
)
parser.add_argument(
    "profres", help="YAML file of profile resultsn"
)

args = parser.parse_args()


fd = open(args.profres, 'r')
y = yaml.safe_load(fd)
fd.close()

offrange = y['optres']['BB_OFFPULSE']

if offrange[1] > offrange[0]:
    # 0.5 - 0.6 means off pulse is in range 0.5 to 0.6 and rest is onpulse
    print(f"OFF: PULSE_PHASE>{offrange[0]} && PULSE_PHASE<{offrange[1]}")
    print(f"ON : PULSE_PHASE<{offrange[0]} || PULSE_PHASE>{offrange[1]}")
else:
    # 0.9 - 0.1 (max<min) means off pulse is in range 0.9 to 1.0 or 0.0 to 1.0  and rest is onpulse
    print(f"OFF: PULSE_PHASE>{offrange[0]} || PULSE_PHASE<{offrange[1]}")
    print(f"ON : PULSE_PHASE>{offrange[1]} && PULSE_PHASE<{offrange[0]}")


# Now call niextspec to extract onpulse and offpulse spectra

    


