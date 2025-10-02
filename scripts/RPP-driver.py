#!/usr/bin/env python
import argparse
import sys, os
import re

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
    help="Python load script made by nicerl3-spect.",
    default="merged_cutmpu7_loadRPP.py",
)
parser.add_argument(
    "--bbody1_kt",
    help="The starting value of kT (keV) for the first blackbody.",
    default="0.2",
)
parser.add_argument(
    "--bbody2_kt",
    help="The starting value of kT (keV) for the second blackbody (ignored in 1BB models).",
    default="0.5",
)
parser.add_argument(
    "--bbody_norm",
    help="The starting value of the blackbody normalization (will be used for both blackbodies in 2BB models).",
    default="1.e-6",
)
parser.add_argument(
    "--powerlaw_index",
    help="The starting value of the powerlaw photon index.",
    default="2.0",
)
parser.add_argument(
    "--powerlaw_norm",
    help="The starting value of the powerlaw normalization.",
    default="0.01",
)

parser.add_argument(
    "--add-gauss", dest="add_gauss", action="append", default=[],
    help=(
        "Add one Gaussian (repeatable). Format: LineE_keV,Sigma_keV,Norm. "
        "If at least one --add-gauss is provided, Gaussians will be enabled. "
        "Omitting the other Gaussian flags uses fit_rpp.py defaults."
    )
)

parser.add_argument(
    "--gauss-e-window", type=float, default=None,
    help=(
        "±ΔE (keV) window around initial LineE for each Gaussian. "
        "If omitted, RPP-spec.py uses its default: 0.05 keV."
    )
)

parser.add_argument(
    "--gauss-sigma-bounds", type=str, default=None,
    help=(
        "Sigma bounds for all Gaussians as 'smin,smax' in keV. "
        "If omitted, RPP-spec.py uses its default: 0.001,0.1 keV."
    )
)


args = parser.parse_args()

# --------------------------------------------

src = args.srcname
loadfile = args.loadfile

nh_startval = '0.1' #units of 1.e22 cm^-2
bbody1_startvals = args.bbody1_kt+','+args.bbody_norm
bbody2_startvals = args.bbody2_kt+','+args.bbody_norm
plaw_startvals = args.powerlaw_index+','+args.powerlaw_norm
bbody_startvals_list = [bbody1_startvals, bbody2_startvals]

models = [
    "tbabs*bbody",
    "tbabs*powerlaw",
    "tbabs*(bbody + bbody)",
    "tbabs*(bbody + powerlaw)",
    "tbabs*(bbody + bbody + powerlaw)"
]
models_short = ['bb','pow','2bb','bbpow','2bbpow']

for m, ms in zip(models, models_short):
    components = re.sub(r'[^a-zA-Z]+', ' ', m).split()

    ii_bb   = [i for i, x in enumerate(components) if x == "bbody"]
    ii_pow  = [i for i, x in enumerate(components) if x == "powerlaw"]
    ii_tbabs= [i for i, x in enumerate(components) if x == "tbabs"]

    for ii in range(len(ii_tbabs)):
        components[ii_tbabs[ii]] = nh_startval
    for ii in range(len(ii_bb)):
        components[ii_bb[ii]] = bbody_startvals_list[ii]
    for ii in range(len(ii_pow)):
        components[ii_pow[ii]] = plaw_startvals

    startvals = ",".join(components)

    cmd = (
        "RPP-spec.py" + ' "' + m + '" ' + startvals +
        " --timestamp" +
        " --loadfile " + loadfile +
        " --srcname " + src +
        " --outname " + src + "-scorp-" + ms
    )

    # Adds Gaussian components: only under user request
    for g in args.add_gauss:
        cmd += " --add-gauss " + g

    # These two are optional: if the user passes them, they override the default values.
    if args.gauss_e_window is not None:
        cmd += " --gauss-e-window " + str(args.gauss_e_window)
    if args.gauss_sigma_bounds is not None:
        cmd += " --gauss-sigma-bounds " + args.gauss_sigma_bounds

    print(cmd)
    os.system(cmd)

