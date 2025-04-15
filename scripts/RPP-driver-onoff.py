#!/usr/bin/env python
import argparse
import sys, os

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""Do RPP catalog spectral analysis.
    """,
)
parser.add_argument(
    "srcname",
    help="The source name (or any string) that you want to appear in plots and output file names.",
    type=str,
)
parser.add_argument(
    "--loadfile",
    help="Python load script made by RPP-pulsedspec.py.",
    default="load_pulsedspec.py",
)
args = parser.parse_args()

# --------------------------------------------

src = args.srcname
loadfile = args.loadfile

cmd = (
    'RPP-spec-onoff.py "tbabs*bbody" 0.01,0.044,1.0 --timestamp --loadfile '
    + loadfile
    + " --srcname "
    + src
    + " --outname rpp-onoff-bb-"
    + src
)
print(cmd)
os.system(cmd)

cmd = (
    'RPP-spec-onoff.py "tbabs*powerlaw" 0.01,1.9,0.01 --timestamp --loadfile '
    + loadfile
    + " --srcname "
    + src
    + " --outname rpp-onoff-pow-"
    + src
)
print(cmd)
os.system(cmd)

cmd = (
    'RPP-spec-onoff.py "tbabs*(bbody + bbody)" 0.01,0.044,1.0,0.044,1.0 --timestamp --loadfile '
    + loadfile
    + " --srcname "
    + src
    + " --outname rpp-onoff-2bb-"
    + src
)
print(cmd)
os.system(cmd)

cmd = (
    'RPP-spec-onoff.py "tbabs*(bbody + powerlaw)" 0.01,0.044,1.0,1.9,0.01 --timestamp --loadfile '
    + loadfile
    + " --srcname "
    + src
    + " --outname rpp-onoff-bbpow-"
    + src
)
print(cmd)
os.system(cmd)

cmd = (
    'RPP-spec-onoff.py "tbabs*(bbody + bbody + powerlaw)" 0.01,0.044,1.0,0.044,1.0,1.9,0.01 --timestamp --loadfile '
    + loadfile
    + " --srcname "
    + src
    + " --outname rpp-onoff-2bbpow-"
    + src
)
print(cmd)
os.system(cmd)
