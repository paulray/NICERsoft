#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as pyfits
import astropy.stats
import astropy.units as u
import argparse
from pint.eventstats import z2m, hm, sf_z2m, sf_hm, sig2sigma, h2sig
#from nicer.values import *
#from nicer.fourier import *
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
    "profres", help="YAML file of profile results"
)
parser.add_argument(
    "resp", help="Response file from non-phase-resolved spectrum"
)
parser.add_argument(
    "arf", help="ARF file from non-phase-resolved spectrum"
)

args = parser.parse_args()

hdul = pyfits.open(args.evtfile)

fd = open(args.profres, 'r')
y = yaml.safe_load(fd)
fd.close()

offrange = y[list(y.items())[0][0]]['optres']['BB_OFFPULSE']

on_phase_range=0
off_phase_range=0
if offrange[1] > offrange[0]:
    # 0.5 - 0.6 means off pulse is in range 0.5 to 0.6 and rest is onpulse
    print(f"OFF: PULSE_PHASE>{offrange[0]} && PULSE_PHASE<{offrange[1]}")
    print(f"ON : PULSE_PHASE<{offrange[0]} || PULSE_PHASE>{offrange[1]}")
    off = f'"PULSE_PHASE>{offrange[0]} && PULSE_PHASE<{offrange[1]}"'
    on = f'"PULSE_PHASE<{offrange[0]} || PULSE_PHASE>{offrange[1]}"'
    on_phase_range=(offrange[0]+(1-offrange[1]))
    off_phase_range=(offrange[1]-offrange[0])
else:
    # 0.9 - 0.1 (max<min) means off pulse is in range 0.9 to 1.0 or 0.0 to 1.0  and rest is onpulse
    print(f"OFF: PULSE_PHASE>{offrange[0]} || PULSE_PHASE<{offrange[1]}")
    print(f"ON : PULSE_PHASE>{offrange[1]} && PULSE_PHASE<{offrange[0]}")
    off = f'"PULSE_PHASE>{offrange[0]} || PULSE_PHASE<{offrange[1]}"'
    on = f'"PULSE_PHASE<{offrange[0]} && PULSE_PHASE>{offrange[1]}"'
    off_phase_range=(offrange[1]+(1-offrange[0]))
    on_phase_range=(offrange[1]-offrange[0])
    
os.system(f"ftselect {args.evtfile} on_pulse.evt {on} clobber=yes")
os.system(f"ftselect {args.evtfile} off_pulse.evt {off} clobber=yes")

os.system(f"niextspect on_pulse.evt on_pulse.pha clobber=yes")
os.system(f"niextspect off_pulse.evt off_pulse.pha clobber=yes")

exp = hdul[1].header['EXPOSURE']
on_exp = exp*on_phase_range
off_exp = exp*off_phase_range

print(f"On phase range:{on_phase_range}, Off phase range:{off_phase_range}")
#exposure
print(f"On Exposure:{on_exp}, Off Exposure:{off_exp}")

os.system(f"fthedit 'on_pulse.pha[1]' keyword=EXPOSURE operation=add value={on_exp} unit=s comment='updated exposure time'")
os.system(f"fthedit 'off_pulse.pha[1]' keyword=EXPOSURE operation=add value={off_exp} unit=s comment='updated exposure time'")

os.system(f"ftgrouppha infile=on_pulse.pha outfile=on_pulse_grp1.pha grouptype=min groupscale=1.0")

g = open("load_pulsedspec.xcm","w+")
g.write("data 1:1 on_pulse_grp1.pha \n")
g.write(f"resp {args.resp} \n")
g.write(f"arf {args.arf} \n")
g.write("back off_pulse.pha \n")
g.write("setplot energy \n")
g.write("ignore bad \n")
g.write("ignore **-0.2 10.0-** \n")
g.write("query yes \n")
g.write("abund wilm \n")
g.close()

g = open("load_pulsedspec.py","w+")
g.write("AllData('1:1 on_pulse_grp1.pha')\n")
g.write("spec1 = AllData(1)\n")
#g.write("spec1.ignore('bad')\n")
g.write("spec1.response = '"+args.resp+"' \n")
g.write("spec1.response.arf = '"+args.arf+"' \n")
g.write("spec1.background = 'off_pulse.pha' \n")
g.write("AllData.ignore('bad')\n")
g.write("AllData.ignore('**-0.2 10.0-**')\n")
g.close()




    


