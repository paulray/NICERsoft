import os, sys
import numpy as np
from optparse import OptionParser


# To run this script with minimal inputs:
#  > python readlog_writexcm.py -n <pulsar name>

# It reads the correctly-named log file (PSRNAME_pulsed_bestfit.log) to get the model parameters.
# The uses the correctly-named .xcm file (PSRNAME_pulsed_bestfit.xcm) as the starting point for writing the new .xcm file.


### Inputs ###

parser = OptionParser()
parser.add_option("-n", "--psrname", dest="psrname", type="string", help="Pulsar name")
parser.add_option("-f", "--paramerr_factor", dest="paramerr_factor", type=float, help="Factor by which to multiply the uncertainty values of the model parameters, in order to run steppar from (param_value-param_unc*factor) to (param_value+param_unc*factor)", default=1.5)
parser.add_option("-s", "--steppar_nsteps", dest="steppar_nsteps", type=int, help="Number of steps to use in steppar", default=12)

(options, args) = parser.parse_args()
psrname = options.psrname
paramerr_factor = options.paramerr_factor
steppar_nsteps = options.steppar_nsteps

input_log = '%s_pulsed_bestfit.log' % psrname
input_xcm = '%s_pulsed_bestfit.xcm' % psrname

### Open and read log file ###

log = open(input_log, 'r')
lines = log.readlines()
log.close()


### Find the final "show all" command and the final "current model list" line ###

'''
This is what we're looking for in the .log file:

#Current model list:
#
#========================================================================
#Model TBabs<1>(bbodyrad<2> + powerlaw<3>) Source No.: 1   Active/On
#Model Model Component  Parameter  Unit     Value
# par  comp
#   1    1   TBabs      nH         10^22    6.87094E-02  +/-  4.73635E-03
#   2    2   bbodyrad   kT         keV      0.175399     +/-  5.56310E-03
#   3    2   bbodyrad   norm                3.79978      +/-  0.610930
#   4    3   powerlaw   PhoIndex            3.01897      +/-  9.75097E-02
#   5    3   powerlaw   norm                1.53023E-05  +/-  8.84851E-07
#________________________________________________________________________
'''

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

newxcm = open('%s_pulsed_errorsfluxes.xcm'%psrname, 'w')
newxcm.write('log %s_pulsed_final.log\n'%psrname)
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
newxcm.write('save all %s_pulsed_final.xcm\n'%psrname)
newxcm.close()
