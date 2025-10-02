import os, sys
import numpy as np
import glob as glob
import yaml

### Inputs ###

logfiles = glob.glob('*_scorpeon_final.log')
if len(logfiles) > 1:
    print('Error: There are multiple files ending in "_scorpeon_final.log". Please remove old versions of these logs before running RPP-scorpeon_yml.py.')
    sys.exit(1)

logfile = glob.glob('*_scorpeon_final.log')[0]
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
model_inds = np.array([], dtype=int)
for i in range(currentmodel_inds[-1],len(lines)):
    if lines[i].startswith('# par  comp'):
        parcomp_inds = np.append(parcomp_inds, i)
        
model_inds = np.arange(parcomp_inds[0]+1,parcomp_inds[1]-6)
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

# First need to find the lines containing the error values.
# error_start_ind_tmp is the earliest possible line containing an error value.
# error_end_ind is the line where steppar is run, and therefore is the latest possible error line.
for i in range(0,len(lines)):
    if lines[i].startswith('!XSPEC12>error'):
        error_start_ind_tmp = i
    elif lines[i].startswith('!XSPEC12>steppar 1'):
        error_end_ind = i

# Now we need to make sure that error_start_ind is truly the correct starting place.
# The reason is that if a better fit is found in the process of running 'error' in XSPEC,
#  then it will do another fit and will re-run error. So there could be multiple lines
#  containing an error value for the same parameter. We need the very last set of error lines
#  for our final error values.
# Thus we search for the final instance of an error line for parameter 1. Error lines begin with
#    '#     1 '
#  so when we find the last instance of a line starting in this way, then we know params 2, etc.
#  will follow in the next error lines.
for i in range(error_start_ind_tmp, error_end_ind):
    if lines[i].startswith('#     1 '):  # This is how any error lines for param 1 begin
        error_start_ind = i

error_lines = np.array([], dtype=str)
for i in range(error_start_ind, error_end_ind):
    if lines[i].startswith('#     ') and lines[i].endswith(')\n'):
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

for i in range(0,len(lines)):
    if lines[i].startswith('!XSPEC12>flux 0.4 2.'):
        phflux_04_2 = lines[i+2].split(' ')[3]
        enflux_04_2 = lines[i+2].split(' ')[5].strip('(')
    elif lines[i].startswith('!XSPEC12>flux 2. 10.'):
        phflux_2_10 = lines[i+2].split(' ')[3]
        enflux_2_10 = lines[i+2].split(' ')[5].strip('(')
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
            'gamma_powerlaw': float(str(param_vals[3])),
            'gamma_powerlaw_low': float(str(confint_low[3])),
            'gamma_powerlaw_high': float(str(confint_high[3])),
            'norm_powerlaw': float(str(param_vals[4])),
            'norm_powerlaw_low': float(str(confint_low[4])),
            'norm_powerlaw_high': float(str(confint_high[4])),
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


yaml_file = open('%s_scorpeon.yml'%psrname,'w')
yaml.dump(dict, yaml_file)
yaml_file.close()
