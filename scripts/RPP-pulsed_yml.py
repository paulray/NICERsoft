#!/usr/bin/env python
import os, sys
import numpy as np
import glob as glob
import yaml

### Inputs ###

logfiles = glob.glob('*_pulsed_final.log')
if len(logfiles) > 1:
    print('Error: There are multiple files ending in "_pulsed_final.log". Please remove old versions of these logs before running RPP-pulsed_yml.py.')
    sys.exit(1)

logfile = glob.glob('*_pulsed_final.log')[0]
print(logfile)
psrname = logfile.split('_')[0]


### Open and read log file ###

log = open(logfile, 'r')
lines = log.readlines()
log.close()


### Find the final "show all" command ###
# There may have been multiple instances of "show all".
# This finds the final instance, which has the best-fit model parameters.

showall_inds = np.array([],dtype=int)
for i in range(0,len(lines)):
    if 'show all' in lines[i]:
        showall_inds = np.append(showall_inds, i)


### Find the final "current model list" line ###

currentmodel_inds = np.array([],dtype=int)
for i in range(showall_inds[-1],len(lines)):
    if 'Current model list' in lines[i]:
        currentmodel_inds = np.append(currentmodel_inds, i)


### Find the model parameter lines ###

parcomp_inds = np.array([], dtype=int)
for i in range(currentmodel_inds[-1],len(lines)):
    if lines[i].startswith('# par  comp'):
        parcomp_inds = np.append(parcomp_inds, i)
    if lines[i].startswith('#______'):
        parcomp_inds = np.append(parcomp_inds, i)

if len(parcomp_inds) > 2:
    print('WARNING: parcomp_inds should be a 2-element array, but this instance of parcomp_inds has more than 2 elements! Check that the correct model lines are being used for writing out the fit parameters!')

model_inds = np.arange(parcomp_inds[0]+1,parcomp_inds[1])
model_lines = np.array([], dtype=str)
for i in model_inds:
    model_lines = np.append(model_lines, lines[i].strip())
print(model_lines)


'''
# This is a sample of the model lines:
['#   1    1   TBabs      nH         10^22    2.46033E-02  +/-  8.54016E-03'
 '#   2    2   bbody      kT         keV      4.46025E-02  +/-  3.03212E-03'
 '#   3    2   bbody      norm                9.38307E-05  +/-  4.05075E-05'
 '#   4    3   powerlaw   PhoIndex            2.37633      +/-  0.344943'
 '#   5    3   powerlaw   norm                8.30041E-05  +/-  1.25409E-05']
'''


### Read error from log ###

for i in range(0,len(lines)):
    if lines[i].startswith('!XSPEC12>error'):
        error_start_ind = i
    elif lines[i].startswith('!XSPEC12>steppar 1'):
        error_end_ind = i

error_lines = np.array([], dtype=str)
for i in range(error_start_ind, error_end_ind):
    if lines[i].startswith('#     '):
        error_lines = np.append(error_lines, lines[i])
print(error_lines)
    
    
"""        
TO DO
# Read steppar from log and fit parabola

steppar_inds = np.array([], dtype=int)

for i in range(0,len(lines)):
    if 'PG-Statistic    Delta' in lines[i]:
        parnum = lines[i+1].strip().split()[-1]
        print(parnum) # 2D array?
        #for j in range(i+3,i+15):
            
        #lines[i+3]

"""



### Read fluxes ###

'''
Flux lines look like this:

#
!XSPEC12>flux 0.4 2.
# Model Flux 5.0665e-05 photons (5.9116e-14 ergs/cm^2/s) range (0.40000 - 2.0000 keV)
#
!XSPEC12>flux 2. 10.
# Model Flux 1.8176e-06 photons (9.6324e-15 ergs/cm^2/s) range (2.0000 - 10.000 keV)
'''


for i in range(0,len(lines)):
    if lines[i].startswith('!XSPEC12>flux 0.4 2.'):
        phflux_04_2 = lines[i+1].split(' ')[3]
        enflux_04_2 = lines[i+1].split(' ')[5].strip('(')
    elif lines[i].startswith('!XSPEC12>flux 2. 10.'):
        phflux_2_10 = lines[i+1].split(' ')[3]
        enflux_2_10 = lines[i+1].split(' ')[5].strip('(')
print(phflux_04_2, enflux_04_2, phflux_2_10, enflux_2_10)



### Write to .yml ###

# We have the following possibilities for model and parameter combos:
# TBabs nH + bbodyrad kT + bbodyrad norm -> bbcount = 2, plcount = 0
# TBabs nH + bbodyrad kT + bbodyrad norm + powerlaw PhoIndex + powerlaw norm -> bbcount=2, plcount=2
# TBabs nH + bbodyrad kT + bbodyrad norm + bbodyrad kT + bbodyrad norm + powerlaw PhoIndex + powerlaw norm -> bbcount=4, plcount=2
# TBabs nH + powerlaw PhoIndex + powerlaw norm -> bbcount=0, plcount=2

# Get spectral parameters
model_names = np.array([], dtype=str)
param_names = np.array([], dtype=str)
param_vals = np.array([], dtype=str)
bbcount = 0
plcount = 0
for m in model_lines:
    modname = m.split()[3].strip()
    if 'bbody' in modname:
        bbcount += 1
    if 'powerlaw' in modname:
        plcount += 1
    model_names = np.append(model_names, m.split()[3].strip())
    param_names = np.append(param_names, m.split()[4].strip())
    param_vals = np.append(param_vals, m.split()[-3].strip())

# Get error values
# Right now, using the "error" command in xspec for the errors
# There are the same number of error_lines as model_lines
# Will need to change this procedure when instead using steppar
confint_low = np.array([], dtype=str)
confint_high = np.array([], dtype=str)
for r in error_lines:
    confint_low = np.append(confint_low, r.split()[2].strip())
    confint_high = np.append(confint_high, r.split()[3].strip())
print(confint_low)
print(confint_high)

# Make dictionary and write to yaml
# For values coming out of a numpy array, you need to change their type from numpy.str or numpy.float to just str and then to float, e.g., 'nH': float(str(param_vals[0]))
if bbcount == 2 and plcount == 0:
    dict = {'name': psrname,
            'nH': float(str(param_vals[0])),
            'nH_low': float(str(confint_low[0])),
            'nH_high': float(str(confint_high[0])),
            'kT_thermal': float(str(param_vals[1])),
            'kT_thermal_low': float(str(confint_low[1])),
            'kT_thermal_high': float(str(confint_high[1])),
            'norm_thermal': float(str(param_vals[2])),
	    'norm_thermal_low': float(str(confint_low[2])),
            'norm_thermal_high': float(str(confint_high[2])),
            'flux_0.4-2': float(enflux_04_2),
            'flux_2_10': float(enflux_2_10)
    }
elif bbcount == 0 and plcount == 2:
    dict = {'name': psrname,
            'nH': float(str(param_vals[0])),
            'nH_low': float(str(confint_low[0])),
            'nH_high': float(str(confint_high[0])),
            'gamma_powerlaw': float(str(param_vals[1])),
            'gamma_powerlaw_low': float(str(confint_low[1])),
            'gamma_powerlaw_high': float(str(confint_high[1])),
            'norm_powerlaw': float(str(param_vals[2])),
            'norm_powerlaw_low': float(str(confint_low[2])),
            'norm_powerlaw_high': float(str(confint_high[2])),
            'flux_0.4-2': float(enflux_04_2),
            'flux_2_10': float(enflux_2_10)
    }
elif bbcount == 2 and plcount == 2:
    dict = {'name': psrname,
            'nH': float(str(param_vals[0])),
            'nH_low': float(str(confint_low[0])),
            'nH_high': float(str(confint_high[0])),
            'kT_thermal': float(str(param_vals[1])),
            'kT_thermal_low': float(str(confint_low[1])),
            'kT_thermal_high': float(str(confint_high[1])),
            'norm_thermal': float(str(param_vals[2])),
            'norm_thermal_low': float(str(confint_low[2])),
            'norm_thermal_high': float(str(confint_high[2])),
            'gamma_powerlaw': float(str(param_vals[3])),
            'gamma_powerlaw_low': float(str(confint_low[3])),
            'gamma_powerlaw_high': float(str(confint_high[3])),
            'norm_powerlaw': float(str(param_vals[4])),
            'norm_powerlaw_low': float(str(confint_low[4])),
            'norm_powerlaw_high': float(str(confint_high[4])),
            'flux_0.4-2': float(enflux_04_2),
            'flux_2_10': float(enflux_2_10)
    }
elif bbcount == 4 and plcount == 2:
    dict = {'name': psrname,
            'nH': float(str(param_vals[0])),
            'nH_low': float(str(confint_low[0])),
            'nH_high': float(str(confint_high[0])),
            'kT_thermal': float(str(param_vals[1])),
            'kT_thermal_low': float(str(confint_low[1])),
            'kT_thermal_high': float(str(confint_high[1])),
            'norm_thermal': float(str(param_vals[2])),
	    'norm_thermal_low': float(str(confint_low[2])),
            'norm_thermal_high': float(str(confint_high[2])),
            'kT_thermal2': float(str(param_vals[3])),
            'kT_thermal2_low': float(str(confint_low[3])),
            'kT_thermal2_high': float(str(confint_high[3])),
            'norm_thermal2': float(str(param_vals[4])),
	    'norm_thermal2_low': float(str(confint_low[4])),
            'norm_thermal2_high': float(str(confint_high[4])),
            'gamma_powerlaw': float(str(param_vals[5])),
            'gamma_powerlaw_low': float(str(confint_low[5])),
            'gamma_powerlaw_high': float(str(confint_high[5])),
            'norm_powerlaw': float(str(param_vals[6])),
            'norm_powerlaw_low': float(str(confint_low[6])),
            'norm_powerlaw_high': float(str(confint_high[6])),
            'flux_0.4-2': float(enflux_04_2),
            'flux_2_10': float(enflux_2_10)
    }


yaml_file = open('%s_pulsed.yml'%psrname,'w')
yaml.dump(dict, yaml_file)
yaml_file.close()
