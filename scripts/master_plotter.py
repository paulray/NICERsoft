import matplotlib.pyplot as plot
import numpy as np
import pylab
import matplotlib.gridspec as gridspec
import argparse
from functionality import *
from sci_plots import sci_plots
from eng_plots import eng_plots
from astropy import log
plot.close('all')

#"What's efficiency?" Prof Darmofal re: code 

parser = argparse.ArgumentParser(description = "Plot the NICER data nicely.")
parser.add_argument("-s", "--save", help = "Save the output plots nicely.", action = "store_true")
parser.add_argument("-d", "--display", help = "Nicely display the nice NICER plots.", action = "store_true")
parser.add_argument("-c", "--sci", help = "Makes some nice science plots", action = "store_true")
parser.add_argument("-e", "--eng", help = "Makes some nice engineering plots", action = "store_true")
parser.add_argument("-p", "--pi", help = "Turn the pha_slow into nice PI values", action = "store_true")
args = parser.parse_args()

#TO CHANGE THE FILES, EDIT THE LIST BELOW:
filenames =['/export/home/zarzouma/pipeline/1706140000/xti/event_cl/ni1706140000_0mpu0_uf.evt',
'/export/home/zarzouma/pipeline/1706140000/xti/event_cl/ni1706140000_0mpu1_uf.evt',
'/export/home/zarzouma/pipeline/1706140000/xti/event_cl/ni1706140000_0mpu2_uf.evt',
'/export/home/zarzouma/pipeline/1706140000/xti/event_cl/ni1706140000_0mpu3_uf.evt',
'/export/home/zarzouma/pipeline/1706140000/xti/event_cl/ni1706140000_0mpu4_uf.evt',
'/export/home/zarzouma/pipeline/1706140000/xti/event_cl/ni1706140000_0mpu5_uf.evt',
'/export/home/zarzouma/pipeline/1706140000/xti/event_cl/ni1706140000_0mpu6_uf.evt'
]

log.info('Reading files')
data1, event_flags, info = smush(filenames)

IDS = np.array([0, 1, 2, 3, 4, 5,6, 7, 10, 11, 12, 13, 14, 15, 16, 17, 20, 21, 22, 23, 24, 25, 26, 27, 30, 31, 32, 33, 34, 35, 36, 37, 40, 41, 42, 43, 44, 45, 46, 47, 50, 51, 52, 53, 54, 55, 56, 57, 60, 61, 62, 63, 64, 65, 66, 67])

log.info('Computing reset rates')
avg_rate, ID_rates = reset_rate(data1, event_flags, IDS)

#Making the plots
if args.eng:
    figure1, num_events = eng_plots(data1, event_flags, info, avg_rate, ID_rates)
    plot.show()
elif args.sci:
    IDS, num_events, stdev, colors = hist_use(data1)
    figure2 = sci_plots(data1, event_flags, info, num_events, avg_rate, 0, ID_rates)
    plot.show()
elif args.pi:
    figure1, num_events = eng_plots(data1, event_flags, info, avg_rate, ID_rates)
    figure2 = sci_plots(data1, event_flags, info, num_events, avg_rate, 1, ID_rates)
    plot.show()
else:
    figure1, num_events = eng_plots(data1, event_flags, info, avg_rate, ID_rates)
    figure2 = sci_plots(data1, event_flags, info, num_events, avg_rate, 0, ID_rates)
    
if args.save:
    figure1.set_size_inches(16,12)
    figure2.set_size_inches(16,12)
    figure1.savefig(str(info["OBJECT"])+ str(info["DATE-OBS"]) +  "_ENG_.png", dpi = 100)
    figure2.savefig(str(info["OBJECT"])+ str(info["DATE-OBS"]) +  "_SCI_.png", dpi = 100)
elif args.display:
    plot.show()
