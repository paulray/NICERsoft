#!/usr/bin/env python

import numpy as np
from pint.models import get_model

import argparse


def get_parser():
    parser = argparse.ArgumentParser(
        description="Compare 2 par file"
    )
    parser.add_argument("par1", help="First par file", default=None)
    parser.add_argument("par2", help="Second par file", default=None)
    parser.add_argument("--verbose", help="Print more details.", action="store_true", default=False)
    return parser


def main(argv=None):
    args = get_parser().parse_args(argv)

    m1 = get_model(args.par1)
    m2 = get_model(args.par2)

    print(m1.compare(m2))


if __name__ == "__main__":
    main()
