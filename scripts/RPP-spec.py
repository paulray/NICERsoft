#!/usr/bin/env python
from xspec import *

# from loguru import logger as log
import matplotlib.pyplot as plt
import argparse
import os, sys
import time

# Parsing of CL -------------------------------------------------------


# Define a custom argument type for a list of floats
def list_of_floats(a):
    return list(map(float, a.split(",")))


# Define a custom argument type for a list of strings
def list_of_strings(a):
    return a.split(",")


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""Do RPP catalog spectral analysis.
    """,
)

parser.add_argument(
    "model_str", help="String describing the model, e.g. 'tbabs* (bb + pow)'", type=str
)
parser.add_argument(
    "model_startvals",
    help="Comma-separated list of starting parameter values; no spaces.",
    type=list_of_floats,
)
parser.add_argument(
    "--loadfile",
    help="Python load script made by nicerl3-spect.",
    default="merged_cutmpu7_loadRPP.py",
)
parser.add_argument(
    "--outname", help="Basename for output files.", default="myRPP-spec"
)
parser.add_argument(
    "--srcname",
    help="Source name, if you want it to appear on spectral plot.",
    default="",
)
parser.add_argument(
    "--freeze",
    help="Comma-separated list of parameters to freeze; no spaces.",
    type=list_of_strings,
    default="nxb(7)",
)
parser.add_argument(
    "--thaw",
    help="Comma-separated list of parameters to thaw; no spaces.",
    type=list_of_strings,
    default="sky(12)",
)
parser.add_argument(
    "--warn",
    help="Warn if an existing .xcm file will be overwritten (requires interactive user input). Default: overwrite without warning (non-interactive).",
    default=False,
    action="store_true",
)
parser.add_argument(
    "--timestamp",
    help="Add the date and time to output filenames.",
    default=False,
    action="store_true",
)
# Add option to read in a yaml file for starting values

args = parser.parse_args()

# Load models -------------------------------------------------------

Plot.splashPage = False
Xset.chatter = 10
# Xset.logChatter = 0
Xset.abund = "wilm"

exec(open(args.loadfile).read())
spec1.ignore("0.0-0.25")
m1 = Model(args.model_str, modName="m")
m1.setPars(args.model_startvals)
AllModels.show()
# AllData.show()

sky = AllModels(1, modName="sky")
nxb = AllModels(1, modName="nxb")

# Do the fitting -------------------------------------------------------

if args.timestamp:
    logfile = args.outname + "-" + time.strftime("%Y%m%d-%H%M%S") + ".log"
    xcmfile = args.outname + "-" + time.strftime("%Y%m%d-%H%M%S") + ".xcm"
else:
    logfile = args.outname + ".log"
    xcmfile = args.outname + ".xcm"

frozen_list = []

# Thaw params in SCORPEON, as requested
frozen_list.append("\n\nTHAWED on command line by user pre-fit:\n")
for pname in args.thaw:
    p = eval(pname)
    p.frozen = False
    tmp = " ".join(list(map(str, p.values)))
    tmp = "\n%s thawed with (val,delta,min,bot,top,max): %s" % (pname, tmp)
    frozen_list.append(tmp)

# Freeze as requested
frozen_list.append("\n\nFROZEN on command line by user pre-fit:\n")
for pname in args.freeze:
    p = eval(pname)
    p.frozen = True
    tmp = " ".join(list(map(str, p.values)))
    tmp = "\n%s frozen with (val,delta,min,bot,top,max): %s" % (pname, tmp)
    frozen_list.append(tmp)

Fit.statMethod = "pgstat"
Fit.nIterations = 200
Fit.query = "no"
Fit.perform()
Fit.show()

# Check if parameters are pegged at zero and freeze them. For sky components that can only vary within a certain range and can't be zero, freeze them at the current value. Each time the fit is redone.
frozen_list.append(
    "\n\nFROZEN automatically b/c of sigma = -1 or sigma > 10*value after initial fit:\n"
)
# frozen_list.append('\n\nFROZEN automatically b/c of sigma = -1 after initial fit:\n')
for m in [m1, sky, nxb]:
    for ii in range(m.startParIndex, m.nParameters + 1):
        if m(ii).sigma == -1 or m(ii).sigma > 10 * m(ii).values[0]:
            # if m(ii).sigma == -1:
            # If bot,top values allow zero, freeze it at zero
            if m(ii).values[3] <= 0 and m(ii).values[5] >= 0:
                m(ii).values = "0.0,0"
            else:  # freeze it at current value
                m(ii).frozen = True

            tmp = " ".join(list(map(str, m(ii).values)))
            tmp = "\n%s(%d) frozen with (val,delta,min,bot,top,max): %s" % (
                m.name,
                ii,
                tmp,
            )
            frozen_list.append(tmp)

            Fit.perform()

# print(Fit.statistic,Fit.testStatistic)
Xset.openLog(logfile)
AllData.show()
AllModels.show()
Fit.show()
AllModels.calcFlux("0.4 2.0")
flux_loint = AllData(1).flux[0]  # first element is the erg/cm^2/s energy flux
AllModels.calcFlux("2.0 10.0")
flux_hiint = AllData(1).flux[0]  # first element is the erg/cm^2/s energy flux
Xset.closeLog()

fd = open(logfile, "a+")
for f in frozen_list:
    fd.write(f)
fd.close()

if os.path.exists(xcmfile) and not args.warn:
    os.remove(xcmfile)
Xset.save(xcmfile)  # this asks for user input if xcmfile exists

# Save YAML file----------------------------------------------------

# print(m1.expression)
# print(m1.componentNames)
# print(m1.nParameters)
# print(m1.startParIndex)

if args.timestamp:
    fd = open(args.outname + "-" + time.strftime("%Y%m%d-%H%M%S") + ".yml", "w")
else:
    fd = open(args.outname + ".yml", "w")

bbody_count = 0
for cname in m1.componentNames:
    c = getattr(m1, cname)

    if "bbody" in cname:
        bbody_count = bbody_count + 1

    for pname in c.parameterNames:
        p = getattr(c, pname)
        pval = p.values[0]
        plo = p.values[0] - p.sigma
        phi = p.values[0] + p.sigma

        if cname == "TBabs":
            fd.write("nH: %g\n" % pval)
            fd.write("nH_high: %g\n" % phi)
            fd.write("nH_low: %g\n" % plo)
        elif "bbody" in cname:
            if bbody_count < 2:
                fd.write("%s_thermal: %g\n" % (pname, pval))
                fd.write("%s_thermal_high: %g\n" % (pname, phi))
                fd.write("%s_thermal_low: %g\n" % (pname, plo))
            else:
                fd.write("%s_thermal%d: %g\n" % (pname, bbody_count, pval))
                fd.write("%s_thermal%d_high: %g\n" % (pname, bbody_count, phi))
                fd.write("%s_thermal%d_low: %g\n" % (pname, bbody_count, plo))
        elif "powerlaw" in cname:
            if "PhoIndex" in pname:
                pname = "gamma"
            fd.write("%s_powerlaw: %g\n" % (pname, pval))
            fd.write("%s_powerlaw_high: %g\n" % (pname, phi))
            fd.write("%s_powerlaw_low: %g\n" % (pname, plo))

fd.write("flux_0.4-2: %g\n" % flux_loint)
fd.write("flux_2-10: %g\n" % flux_hiint)
fd.close()

# Make plots -------------------------------------------------------

# Don't generate XSPEC plots, just use this to get the plot data out
Plot.device = "/null"
# Plot.device = "/xs" # For debugging
Plot.xAxis = "keV"
Plot.add = True
Plot.xLog = True
Plot.yLog = True
Plot.setRebin()
Plot.addCommand("res y 0.01 20")
# Plot.addCommand("wenv output")
# Plot("ldata ratio resid")
# Plot("ldata resid")
Plot("ldata")
# Plot.show()
# sys.exit()
# Plot.delCommand(1)
en = Plot.x()
rates = Plot.y()
rates_err = Plot.yErr()
folded = Plot.model()

Plot.xAxis = "keV"
Plot.add = True
Plot.xLog = True
Plot.yLog = False
Plot.setRebin()
# Plot.addCommand("res y 0.01 20")
# Plot.addCommand("wenv output")
# Plot("ldata ratio resid")
# Plot("ldata resid")
Plot("ratio")
# Plot.show()
# Plot.delCommand(1)
en = Plot.x()
ratio = Plot.y()
ratio_err = Plot.yErr()
# folded = Plot.model()

fig, (ax1, ax2) = plt.subplots(
    2, 1, sharex=True, gridspec_kw={"hspace": 0, "height_ratios": [3, 1]}
)
fig.suptitle("%s   Model: %s" % (args.srcname, args.model_str))

ax1.errorbar(en, rates, yerr=rates_err, color="black", linestyle="")
ax1.step(en, rates, label="Data", color="black")
ax1.plot(en, folded, label="Full Model", color="orange")
ax1.plot(en, m1.folded(1), label="Source", color="red")
ax1.plot(en, sky.folded(1), label="Sky", color="blue")
ax1.plot(en, nxb.folded(1), label="NXB", color="green")
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_ylim(bottom=1.0e-9)
ax1.tick_params(
    axis="x", which="both", bottom=True, top=True, direction="in", labelbottom=False
)
ax1.legend(ncol=5, mode="expand")
# ax1.set_xlabel('Energy (keV)')
ax1.set_ylabel("counts / s / keV")
ax1.set_title(
    "Chi2: %.2f  DOF: %d  Chi2/DOF: %.2f  S[0.4-2]: %.2g  S[2-10]: %.2g"
    % (Fit.testStatistic, Fit.dof, Fit.testStatistic / Fit.dof, flux_loint, flux_hiint),
    fontsize="medium",
)

ax2.errorbar(en, ratio, yerr=ratio_err, color="black")
ax2.axhline(1, color="blue", linewidth=0.5)
ax2.set_xscale("log")
ax2.set_yscale("linear")
ax2.tick_params(axis="x", which="both", top=True, direction="in", labelbottom=True)
ax2.set_xlabel("Energy (keV)")
ax2.set_ylabel("Data/Full Model")

# plt.show()
if args.timestamp:
    plt.savefig(args.outname + "-" + time.strftime("%Y%m%d-%H%M%S") + ".png")
else:
    plt.savefig(args.outname + ".png")
