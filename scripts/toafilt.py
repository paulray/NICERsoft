#! /usr/bin/env python
import  numpy as np
import argparse
import astropy
import astropy.units as u
import pint.toa
import sys
import os.path

import pint.logging
from loguru import logger as log
pint.logging.setup(level=pint.logging.script_level)


desc="Read TOAs from a .tim file and filter out TOAs based on simple conditions, then write new .tim file"
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("timname", help=".tim file to read TOAs from")
parser.add_argument("--outfile", help="Name of output file (default=STDOUT)", default=None)
parser.add_argument("--clobber", help="Overwrite outfile if it already exists.", action="store_true")
parser.add_argument(
    "--minhtest",
    help="Minimum value of H-test to keep (default=9)",
    type=float,
    default=9.0,
)
parser.add_argument(
    "--minexp",
    help="Minimum exposure keep (default=16 seconds)",
    type=float,
    default=16.0,
)
parser.add_argument(
    "--maxerr",
    help="Maximum TOA error to retain (unit=microsecons, default=25)",
    type=float,
    default=25.0,
)

## Parse arguments
args = parser.parse_args()


ts = pint.toa.get_TOAs(args.timname,include_bipm=False,planets=True)
ts.print_summary()

htest, htidx = ts.get_flag_value("htest",as_type=float)
htest = np.array(htest)


ncliphtest = (htest <= args.minhtest).sum()
ts_clip = ts[htest>args.minhtest]

ncliperr = (ts_clip.get_errors()>=args.maxerr*u.us).sum()
ts_clip = ts_clip[ts_clip.get_errors()<args.maxerr*u.us]

exp, expidx = ts_clip.get_flag_value("exposure",as_type=float)
exp = np.array(exp)

nclipexp = (exp <= args.minexp).sum()
ts_clip = ts_clip[exp>args.minexp]

print(f"Clipped {ncliphtest} TOAs for H-test; {ncliperr} for large errors and {nclipexp} for short exposures, out of {len(ts.get_mjds())}",file=sys.stderr)

if args.outfile:
    if os.path.isfile(args.outfile) and not args.clobber:
        log.error("Output file exists! Specify --clobber to overwrite")
    else:
        with open(args.outfile,"w") as ff:
            ts_clip.write_TOA_file(args.outfile)
