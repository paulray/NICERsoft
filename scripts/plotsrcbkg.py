#!/usr/bin/env python
from xspec import *
import argparse
import os
from os import path
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as pyfits

parser = argparse.ArgumentParser(description="Plot source and background spectra")
parser.add_argument("tot",help="PHA file name with src+bkg spectrum")
parser.add_argument("bkg",help="PHA file name with bkg estimate spectrum")
args = parser.parse_args()

tot = Spectrum(args.tot)
bkg = Spectrum(args.bkg)

objname = pyfits.getval(args.tot,'OBJECT')
exposure = pyfits.getval(args.tot,'EXPOSURE')

tot.response = path.join(os.environ['CALDB'],'data/nicer/xti/cpf/rmf/nixtiref20170601v001.rmf')
tot.response.arf = path.join(os.environ['CALDB'],'data/nicer/xti/cpf/arf/nixtiaveonaxis20170601v002.arf')

bkg.response = path.join(os.environ['CALDB'],'data/nicer/xti/cpf/rmf/nixtiref20170601v001.rmf')
bkg.response.arf = path.join(os.environ['CALDB'],'data/nicer/xti/cpf/arf/nixtiaveonaxis20170601v002.arf')


tot.ignore('**-0.3,8.0-**')
bkg.ignore('**-0.3,8.0-**')

#Plot.device = '/xs'
Plot.device = '/null'
Plot.xAxis = 'keV'
Plot('data')

energies = np.array(Plot.x(1))
p_tot = np.array(Plot.y(1))
p_toterr = np.array(Plot.yErr(1))
p_bkg = np.array(Plot.y(2))

# Plot using Matplotlib:
fig, axs = plt.subplots(nrows=2,ncols=1,sharey=True)
axs[0].errorbar(energies, p_tot,yerr=p_toterr,alpha=0.5)
axs[0].plot(energies, p_bkg,'r',alpha=0.5)
axs[0].set_xscale('log')
axs[0].grid(True)
axs[0].set_title('{0} Spectrum (Exposure = {1})'.format(objname,exposure))

axs[1].errorbar(energies, p_tot-p_bkg, yerr=p_toterr)
axs[1].set_xlabel('Energy (keV)')
axs[1].set_xscale('log')
axs[1].grid(True)
plt.show()
