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
    description="""Load On/Off .xcm file saved with Xspec; make plot and .yml file for RPP catalog.
    """,
)

parser.add_argument(
    "xcmfile", help="Xspec .xcm file to load.", type=str
)
parser.add_argument(
    "--loadfile",
    help="Python load script made by RPP-pulsedspec.py.",
    default="load_pulsedspec.py",
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
    "--outname", help="Basename for output files.", default="rpp-onoff"
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
#AllData.show()
#AllModels.show()
spec1 = AllData(1)

# Refit to get parameter uncertainties; the result will be basically the same as the starting fit from the .xcm file.
Fit.statMethod = "cstat"
Fit.nIterations = 1000
Fit.query = "no"
Fit.perform()
Fit.show()
Xset.closeLog()

# When we get the data points without background subtraction and reset the spec.background = '' before plotting, that resets the Fit.* variables, so let's save what we need to show on the plots:
test_statistic = Fit.testStatistic
dof = Fit.dof

# Can't assume the name the user has given the model in the .xcm file; it may not have a name at all. Find it by the component names:
m1 = None

for mn in AllModels.sources:
    mname = AllModels.sources[mn]
    m = AllModels(1,mname)
    print('\nModel #:',mn,' Model name:',m.name)
    print('Component names:',m.componentNames)
    #print('Number of parameters:',m.nParameters)

    for cname in m.componentNames:
        if 'powerlaw' in cname or 'bbody' in cname:
            m1 = m

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

minSig = 3
maxBins = 10
errType = 'sqrt'

Plot.device = "/null"
#Plot.device = "/xs" # For debugging
Plot.xAxis = "keV"
Plot.background = True
Plot.add = True
Plot.xLog = True
Plot.yLog = True
Plot.setRebin(minSig=minSig,maxBins=maxBins,errType=errType)
#Plot.addCommand("res y 0.01 20")
# Plot.addCommand("wenv output")
# Plot("ldata ratio")
Plot("ldata")
#Plot("ldata",'model m') #this doesn't work although it's from example online
Plot.show()
# Plot.delCommand(1)
en = Plot.x()
rates = Plot.y()
rates_err = Plot.yErr()
folded = Plot.model()
bkg = Plot.backgroundVals()

Plot.xAxis = "keV"
Plot.add = True
Plot.xLog = True
Plot.yLog = False
#Plot.setRebin(minSig=minSig,maxBins=maxBins,errType=errType)
# Plot.addCommand("res y 0.01 20")
# Plot.addCommand("wenv output")
# Plot("ldata ratio resid")
# Plot("ldata resid")
Plot("ratio")
# Plot.show()
# Plot.delCommand(1)
en = Plot.x()
ratio = Plot.y()
ratio_err = np.array(Plot.yErr())
ratio_err[ratio_err<0.0] = 0.0
# folded = Plot.model()

# Get data points without background subtraction
spec1.background = '' # NOTE: THIS RESETS THE FIT.* VARIABLES!
Plot.device = "/null"
#Plot.device = "/xs" # For debugging
Plot.xAxis = "keV"
Plot.background = False
Plot.add = True
Plot.xLog = True
Plot.yLog = True
Plot.setRebin(minSig=minSig,maxBins=maxBins,errType=errType)
#Plot.addCommand("res y 0.01 20")
# Plot.addCommand("wenv output")
# Plot("ldata ratio")
Plot("ldata")
#Plot("ldata",'model m') #this doesn't work although it's from example online
Plot.show()
# Plot.delCommand(1)
en_nobkg = Plot.x()
rates_nobkg = Plot.y()
rates_err_nobkg = Plot.yErr()

fig, (ax1, ax2) = plt.subplots(
    2, 1, sharex=True, gridspec_kw={"hspace": 0, "height_ratios": [3, 1]}
)
fig.suptitle("%s   Model: %s" % (args.srcname, m1.expression))

ax1.errorbar(en, rates, yerr=rates_err, color="black", marker='.',linestyle="",label='Data(Src)',zorder=1)
ax1.errorbar(en_nobkg, rates_nobkg, yerr=rates_err_nobkg, color="green", marker='.',linestyle="",label='Data(Src+Bkg)',zorder=1)
ax1.step(en, bkg,label='Bkg (Off pulse)')
ax1.plot(en, folded, label="Model(Src)", color="red",zorder=2)
#ax1.plot(en_norebin, m1.folded(1), label="Src", color="red",zorder=2)
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_ylim(bottom=1.0e-3)
ax1.tick_params(
    axis="x", which="both", bottom=True, top=True, direction="in", labelbottom=False
)
#ax1.legend(ncol=5, mode="expand")
ax1.legend()
# ax1.set_xlabel('Energy (keV)')
ax1.set_ylabel("counts / s / keV")
ax1.set_title(
    "Chi2: %.2f  DOF: %d  Chi2/DOF: %.2f  S[0.4-2]: %.2g  S[2-10]: %.2g"
    % (test_statistic, dof, test_statistic / dof, flux_loint, flux_hiint),
    fontsize="medium",
)

ax2.errorbar(en, ratio, yerr=ratio_err, color="black")
ax2.axhline(1, color="blue", linewidth=0.5)
ax2.set_xscale("log")
ax2.set_yscale("linear")
ax2.tick_params(axis="x", which="both", top=True, direction="in", labelbottom=True)
ax2.set_ylim(0.75, 1.25)
ax2.set_xlabel("Energy (keV)")
ax2.set_ylabel("Data/Model")

plt.show()

if args.timestamp:
    plt.savefig(args.outname + "-" + timestamp + ".png")
else:
    plt.savefig(args.outname + ".png")
