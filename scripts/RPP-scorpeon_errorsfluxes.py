#!/usr/bin/env python
import os, sys
import numpy as np
from optparse import OptionParser


# To run this script with minimal inputs:
#  > python readlog_writexcm.py -n <pulsar name> -x <input .xcm filename>

# It reads the log file to get the model parameters
# It only uses the input .xcm to write the new .xcm -- uses the input .xcm as the starting point.

### Inputs ###

parser = OptionParser()
parser.add_option("-n", "--psrname", dest="psrname", type="string", help="Pulsar name")
#parser.add_option("-l", "--input_log", dest="input_log", type="string", help="Name of input log file")
#parser.add_option("-x", "--input_xcm", dest="input_xcm", type="string", help="Name of input .xcm file")
parser.add_option("-f", "--paramerr_factor", dest="paramerr_factor", type=float, help="Factor by which to multiply the uncertainty values of the model parameters, in order to run steppar from (param_value-param_unc*factor) to (param_value+param_unc*factor)", default=1.5)
parser.add_option("-s", "--steppar_nsteps", dest="steppar_nsteps", type=int, help="Number of steps to use in steppar", default=12)

(options, args) = parser.parse_args()
psrname = options.psrname
#input_log = options.input_log
#input_xcm = options.input_xcm
paramerr_factor = options.paramerr_factor
steppar_nsteps = options.steppar_nsteps

input_log = '%s_scorpeon_bestfit.log' % psrname
input_xcm = '%s_scorpeon_bestfit.xcm' % psrname

### Open and read log file ###

log = open(input_log, 'r')
lines = log.readlines()
log.close()


### Find the final "show all" command and the final "current model list" line ###

showall_inds = np.array([],dtype=int)
for i in range(0,len(lines)):
    if 'show all' in lines[i]:
        showall_inds = np.append(showall_inds, i)

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
This is a sample of the model lines:
['#   1    1   TBabs      nH         10^22    2.46033E-02  +/-  8.54016E-03'
 '#   2    2   bbody      kT         keV      4.46025E-02  +/-  3.03212E-03'
 '#   3    2   bbody      norm                9.38307E-05  +/-  4.05075E-05'
 '#   4    3   powerlaw   PhoIndex            2.37633      +/-  0.344943'
 '#   5    3   powerlaw   norm                8.30041E-05  +/-  1.25409E-05']

However, if a parameter is frozen, then the line will appear as, e.g.:
  #   1    1   TBabs      nH         10^22    2.46033E-02  frozen
'''


### Write new .xcm file ###

newxcm = open('%s_scorpeon_errorsfluxes.xcm'%psrname, 'w')
newxcm.write('log %s_scorpeon_final.log\n'%psrname)
newxcm.write('cpd /xwin\n')
newxcm.write('statistic pgstat\n')
newxcm.write('@%s\n' % input_xcm)
newxcm.write('fit 100\n')
newxcm.write('error 1. 1-%d\n' % len(model_lines))
for m in model_lines:
    print(m)
    paramnum = m.strip('\n').split()[1]
    if m.split()[-1] == 'frozen':
        paramval = float(m.split()[-2])
        paramerr = 0.0
    else:
        paramval = float(m.split()[-3])
        paramerr = float(m.split()[-1])
    paramlo = str(paramval-paramerr*paramerr_factor)
    paramhi = str(paramval+paramerr*paramerr_factor)
    newxcm.write('steppar %s %s %s %d\n' % (paramnum,paramlo, paramhi, steppar_nsteps))
    print('param#  paramval  paramerr  lowval  highval')
    print(paramnum, paramval, paramerr, paramlo, paramhi)


newxcm.write('show all\n')
newxcm.write('flux 0.4 2.\n')
newxcm.write('flux 2. 10.\n')
newxcm.write('save all %s_scorpeon_final.xcm\n'%psrname)
newxcm.close()
