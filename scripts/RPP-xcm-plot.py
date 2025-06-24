#!/usr/bin/env python

from xspec import *

# from loguru import logger as log
import matplotlib.pyplot as plt
import argparse
import os, sys
import time
import numpy as np

# Parsing of CL -------------------------------------------------------

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""Load SCORPEON .xcm file saved with Xspec; make plot and .yml file for RPP catalog.
    """,
)
parser.add_argument(
    "xcmfile", help="Xspec .xcm file to load.", type=str
)
parser.add_argument(
    "--loadfile",
    help="Python load script made by nicerl3-spect.",
    default="merged_cutmpu7_loadRPP.py",
)
parser.add_argument(
    "--srcname",
    help="Source name, if you want it to appear on spectral plot.",
    default="",
)
parser.add_argument(
    "--timestamp",
    help="Add the date and time to output filenames.",
    default=False,
    action="store_true",
)
parser.add_argument(
    "--outname", help="Basename for output files.", default="rpp-scorp"
)

args = parser.parse_args()

# Load models -------------------------------------------------------

timestamp = time.strftime("%Y%m%d-%H%M%S")

if args.timestamp:
    logfile = args.outname + "-" + timestamp + ".log"
else:
    logfile = args.outname + ".log"

Plot.splashPage = False
Xset.openLog(logfile)

exec(open(args.loadfile).read())
Xset.restore(args.xcmfile)
#AllModels.show()

# Refit to get parameter uncertainties; the result will be basically the same as the starting fit from the .xcm file.
Fit.statMethod = "pgstat"
Fit.nIterations = 1000
Fit.query = "no"
Fit.perform()
Fit.show()
Xset.closeLog()

# Can't assume the names the user has given the models in the .xcm file; they may not have names at all. Figure out which is which by the component names:
m1 = None
sky = None
nxb = None

for mn in AllModels.sources:
    mname = AllModels.sources[mn]
    m = AllModels(1,mname)
    print('\nModel #:',mn,' Model name:',m.name)
    print('Component names:',m.componentNames)
    #print('Number of parameters:',m.nParameters)

    for cname in m.componentNames:
        if 'powerlaw' in cname or 'bbody' in cname:
            m1 = m
        elif 'sky' in cname:
            sky = m
        elif 'nxb' in cname:
            nxb = m

# Calculate fluxes -------------------------------------------------------
            
AllModels.calcFlux("0.4 2.0")
flux_loint = AllData(1).flux[0]  # first element is the erg/cm^2/s energy flux
AllModels.calcFlux("2.0 10.0")
flux_hiint = AllData(1).flux[0]  # first element is the erg/cm^2/s energy flux

# Save YAML file----------------------------------------------------

if args.timestamp:
    fd = open(args.outname + "-" + timestamp + ".yml", "w")
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
        psigma = p.sigma
        plo = max(p.values[0] - p.sigma, p.values[2])
        phi = min(p.values[0] + p.sigma, p.values[5])

        if cname == "TBabs":
            fd.write("nH: %g\n" % pval)
            fd.write("nH_sigma: %g\n" % psigma)
            fd.write("nH_high: %g\n" % phi)
            fd.write("nH_low: %g\n" % plo)
        elif "bbody" in cname:
            if bbody_count < 2:
                fd.write("%s_thermal: %g\n" % (pname, pval))
                fd.write("%s_thermal_sigma: %g\n" % (pname, psigma))
                fd.write("%s_thermal_high: %g\n" % (pname, phi))
                fd.write("%s_thermal_low: %g\n" % (pname, plo))
            else:
                fd.write("%s_thermal%d: %g\n" % (pname, bbody_count, pval))
                fd.write("%s_thermal%d_sigma: %g\n" % (pname, bbody_count, psigma))
                fd.write("%s_thermal%d_high: %g\n" % (pname, bbody_count, phi))
                fd.write("%s_thermal%d_low: %g\n" % (pname, bbody_count, plo))
        elif "powerlaw" in cname:
            if "PhoIndex" in pname:
                pname = "gamma"
            fd.write("%s_powerlaw: %g\n" % (pname, pval))
            fd.write("%s_powerlaw_sigma: %g\n" % (pname, psigma))
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

en_dif = np.array(en)[1:]-np.array(en)[0:-1]
en_dif = np.concatenate((en_dif,np.array([en_dif[-1]])))
#print('np.size(en_dif): ',np.size(en_dif))

fig, (ax1, ax2) = plt.subplots(
    2, 1, sharex=True, gridspec_kw={"hspace": 0, "height_ratios": [3, 1]}
)
fig.suptitle("%s   Model: %s" % (args.srcname, m1.expression))

ax1.errorbar(en, rates, yerr=rates_err, color="black", linestyle="")
ax1.step(en, rates, label="Data", color="black")
ax1.plot(en, folded, label="Full Model", color="orange")
ax1.plot(en, np.divide(m1.folded(1),en_dif), label="Source", color="red")
ax1.plot(en, np.divide(sky.folded(1),en_dif), label="Sky", color="blue")
ax1.plot(en, np.divide(nxb.folded(1),en_dif), label="NXB", color="green")
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_ylim(bottom=1.0e-6)
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
ax2.set_ylim(0.75, 1.25)
ax2.set_xlabel("Energy (keV)")
ax2.set_ylabel("Data/Full Model")

#plt.show()

if args.timestamp:
    plt.savefig(args.outname + "-" + timestamp + ".png")
else:
    plt.savefig(args.outname + ".png")
