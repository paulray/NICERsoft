#!/usr/bin/env python
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import argparse
from astropy import log
from os import path
from glob import glob
from subprocess import check_call
import shutil
from astropy.table import Table

from nicer.values import *
from nicer.NicerFileSet import *


def get_parser():
    parser = argparse.ArgumentParser(
        description="Renormalize the profile in a template file (data after the META_STOP line) "
        "so that its minimum is zero and its mean is one, writing the result to fixed.tem."
    )
    parser.add_argument("template", help="Input template file")
    return parser


def main(argv=None):
    args = get_parser().parse_args(argv)

    fname = args.template
    log.info(fname)
    fin = open(fname, "r")
    fout = open("fixed.tem", "w")
    while True:
        line = fin.readline().strip()
        print(line, file=fout)
        if line.startswith("META_STOP"):
            break

    log.info("Broke")
    phases, lcSmooth = np.loadtxt(fin, unpack=True)

    lcSmooth -= lcSmooth.min()
    lcSmooth *= len(lcSmooth) / np.trapz(lcSmooth)
    print("trapz(lcSmooth)/len(lcSmooth)= ", np.trapz(lcSmooth) / len(lcSmooth))
    for x, b in zip(np.linspace(0.0, 1.0, len(lcSmooth), endpoint=False), lcSmooth):
        print("%f %f" % (x, b), file=fout)


if __name__ == "__main__":
    main()
