#!/usr/bin/env python
from xspec import *

import matplotlib.pyplot as plt

# import numpy as np
# import astropy.units as u


Plot.splashPage = False
exec(open("merged_cutmpu7_loadRPP.py").read())
m1 = Model("tbabs*nsatmos")
m1.setPars(0.01, 5.8, "1.4,0", "11.0,0", "0.420,0", 0.01)
m1.show()
AllData.show()

sky = AllModels(1, modName="sky")
nxb = AllModels(1, modName="nxb")

# Thaw param 12 (swcxovii_norm)
sky(12).frozen = False

Fit.statMethod = "pgstat"
Fit.nIterations = 200
Fit.query = "no"
Fit.perform()

# Now check if parameters are pegged at zero and freeze them

p = m1(1)  # n_H
if p.sigma == -1:
    # Freeze norm at 0.0 and refit
    p.values = "0.0,0"
    Fit.perform()


p = nxb(4)  # trel_norm
if p.sigma == -1:
    # Freeze norm at 0.0 and refit
    p.values = "0.0,0"
    Fit.perform()

AllModels.calcFlux("0.3 2.0 err")
fluxes = AllData(1).flux
# Eflux, Eflux,min, Eflux,max, Phflux, Phflux,min, Phflux,max, ... for each of the 3 sources


# Plot.device="xxx_fit.png/png"
Plot.device = "/null"
Plot.xAxis = "keV"
Plot.add = True
Plot.xLog = True
Plot.yLog = True
Plot.setRebin()
# Plot.addCommand("res y 0.001 2")
Plot("data")
# Plot.delCommand(1)
en = Plot.x()
rates = Plot.y()
rates_err = Plot.yErr()
folded = Plot.model()

fig, ax = plt.subplots(1, 1)
ax.errorbar(en, rates, yerr=rates_err)
ax.plot(en, folded)
ax.set_xscale("log")
ax.set_yscale("log")

plt.show()
#
# Plot.device="xxx_model.png/png"
# Plot.device = "/null"
# Plot("model")
