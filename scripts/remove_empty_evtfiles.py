#!/usr/bin/env python

import argparse

import astropy.io.fits as pyfits

parser = argparse.ArgumentParser(
    description="Read through list of EVENT fits files and remove any that have no events (0 rows in the EVENTS extension)"
)
parser.add_argument("infile", help="Text file containing list of event file names")
parser.add_argument(
    "outfile", help="Output text file containing filtered list of event file names"
)
args = parser.parse_args()

outf = open(args.outfile, "w")

print("{0:32s} {1:>12s} {2:>12s}".format("Filename", "#EVT", "#GTI"))
for fn in open(args.infile).readlines():
    fn = fn.strip()
    nevt = pyfits.getval(fn, "NAXIS2", extname="EVENTS")
    ngti = pyfits.getval(fn, "NAXIS2", extname="GTI")
    print("{0:32s} {1:12d} {2:12d}".format(fn, nevt, ngti))
    if nevt > 0:
        print(fn, file=outf)

outf.close()
