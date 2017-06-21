#!/usr/bin/env python
from __future__ import division, print_function
from glob import glob
import argparse
import subprocess
from astropy import log
import sys
sys.path.append("/data/NICER/preliminary/DOY170PSR") 

parser = argparse.ArgumentParser(description="Process a set of NICER data into plots")
parser.add_argument("obsids",help="ObsIDs to process", nargs="+")
args = parser.parse_args()
log.info('Parsed some args')

for obsid in args.obsids:
    print(obsid)
    fnames = glob('/data/NICER/preliminary/DOY170PSR/{0}*.evt'.format(obsid))
    if len(fnames) == 0:
	log.info('didnt find anything i guess so bye bye')        
	continue
    log.info('found some stuff, now gonna go call THE MASTER plotter')
    #subprocess.call(['python', '/home/fhuynh/Documents/Practice_scripts/NICER_Plots_Creation/NICERsoft/scripts/master_plotter.py','--eng','--sci', '-s'] + fnames)
    subprocess.call(['python', '/home/fhuynh/Documents/Practice_scripts/NICER_Plots_Creation/NICERsoft/scripts/master_plotter.py','--eng','--sci', '--filtall', '--emin', '0.3', '--emax', '10.0'] + fnames)

