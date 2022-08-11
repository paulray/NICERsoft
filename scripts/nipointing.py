#!/usr/bin/env python
import argparse
import numpy as np
import matplotlib.pyplot as plt
from astroquery.simbad import Simbad

# from astropy.io import fits
from os import path
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
from astropy import log
from nicer.values import *

desc = """
Calculate the optimal (offset) pointing that maximizes the S/N for the target by minimizing the contamination from nearby targets.

Possibility to use SIMBAD coordinates for the source (with exact name). Information about nearby sources must be provided in a text file with 3 columns (RA, DEC, and estimated NICER count-rate).  Sources within 7 or 8 arcmin should be considered.

SIMBAD Queries require package 'astroquery'.
"""

parser = argparse.ArgumentParser(description=desc)
parser.add_argument(
    "name",
    help="SIMBAD Source name for query (must be exact). If --ra and --dec are provided, the query is ignored",
    type=str,
)
parser.add_argument(
    "countrate", help="Expected NICER count rate of your target", type=float
)
parser.add_argument(
    "near_srcs",
    help="Text file with RA (hms), DEC (dms), and estimated NICER count-rates of nearby sources",
    type=str,
)
parser.add_argument("--ra", help="RA in HH:MM:SS.s", type=str, default=None)
parser.add_argument("--dec", help="DEC in DD:MM:SS.s", type=str, default=None)
parser.add_argument(
    "--step",
    help="Number of grid points to search in RA/DEC",
    type=float,
    default=100.0,
)
parser.add_argument(
    "--saveplot", help="Saves plot instead of display", action="store_true"
)
args = parser.parse_args()


def GetCoordPSR(name):
    # Get the PSR coordinates from SIMBAD database
    return SkyCoord(
        ra=Simbad.query_object(name)["RA"][0],
        dec=Simbad.query_object(name)["DEC"][0],
        unit=(u.hourangle, u.deg),
    )


def AngSeparation(reference, obj):
    # Calculates the angular separation between reference position and obj
    return reference.separation(obj)


def ScaledCtRate(D, OptCtRate, effareaX, effareaY):
    # Scale a given count rate (OptCtRate) given an angular distance D
    # EffAreaX corresponds to Effective Area
    # EffAreaY corresponds to OffAxisAngle
    return OptCtRate * np.interp(D, effareaX, effareaY)


def SignaltoNoise(SrcCtsRate, BkgSrcRates, InstBkgd, ExpTime):
    # Calculate the S/N rate given
    #  1) Source count rate
    #  2) Sum the count rates due to background sources.
    #  3) the exposure time (ExpTime)
    #  4) the instrumental+particule background
    SNR = (SrcCtsRate * ExpTime) / np.sqrt(
        ExpTime * (SrcCtsRate + np.sum(BkgSrcRates) + InstBkgd)
    )
    return SNR


## load vigneting information:
## Effective area (EffArea) as function of OffAxisAngle
VignetingFile = path.join(datadir, "nicer_vignetting.dat")
if not path.isfile(VignetingFile):
    log.error(
        "This file nicer_vignetting.dat is missing. Check the {} directory".format(
            path.abspath(datadir)
        )
    )
    exit()
else:
    EffArea, OffAxisAngle = np.loadtxt(VignetingFile, unpack=True, usecols=(0, 1))


# Load info about nearby sources. Format: ra(deg), dec(deg), countrate(c/s)
NearbySources = np.loadtxt(args.near_srcs, unpack=True, dtype=str, usecols=(0, 1))
NearbySourcesCR = np.genfromtxt(args.near_srcs, unpack=True, dtype=float, usecols=(2))
SRCposition = SkyCoord(
    ra=NearbySources[0], dec=NearbySources[1], unit=(u.hourangle, u.deg)
)  ### Change SRCposition to NearbySourcesPOS

# Get pulsar info
PSRcountrate = args.countrate

if args.ra is not None and args.dec is not None:
    print("RA and DEC provided by user...Ignoring SIMBAD coordinates")
    PSRposition = SkyCoord(ra=args.ra, dec=args.dec, unit=(u.hourangle, u.deg))
else:
    PSRposition = GetCoordPSR(args.name)

print(PSRposition.ra.hms, PSRposition.dec.dms)


# Set other parameters
INSTbkgd = 0.2
EXPtime = 1000000


# Print info for Nominal pointing of NICER on the target PSR
SRCnominalDIST = AngSeparation(
    SRCposition, SkyCoord(ra=PSRposition.ra, dec=PSRposition.dec)
).arcmin
SRCscaleRates = ScaledCtRate(SRCnominalDIST, NearbySourcesCR, EffArea, OffAxisAngle)
print(
    "Target S/N at Nominal Pointing "
    + str(SignaltoNoise(PSRcountrate, SRCscaleRates, INSTbkgd, EXPtime))
)
print("Target count rate at Nominal pointing = " + str(PSRcountrate) + " cts/sec")
print(
    "Total count rate from nearby sources at Nominal pointing = "
    + str(np.sum(SRCscaleRates))
    + " cts/sec"
)
print("     Individual Nearby Sources rates: ")
print(str(SRCscaleRates))
print("     Nearby sources distance from Target (')")
print(SRCnominalDIST)
print("----------------------------------------------------------------------")


# Sampling positions:  (-3,3) arcmin in RA, (-3,3) arcmin in DEC,
DeltaRA = Angle(np.linspace(-3, 3.01, args.step), unit=u.deg) / 60
DeltaDEC = Angle(np.linspace(-3, 3.01, args.step), unit=u.deg) / 60
# Other arrays needed
sampleRA = np.zeros(len(DeltaRA) * len(DeltaDEC))
sampleDEC = np.zeros(len(DeltaRA) * len(DeltaDEC))
snr = np.zeros(len(DeltaRA) * len(DeltaDEC))
PSRrates = np.zeros(len(DeltaRA) * len(DeltaDEC))
SRCrates = np.zeros(len(DeltaRA) * len(DeltaDEC))

count = 0
for i in DeltaRA:
    for j in DeltaDEC:
        # Define NICER pointing given the offset
        NICERpointing = SkyCoord(ra=PSRposition.ra + i, dec=PSRposition.dec + j)

        # Separation between PSR and NICER
        PSRseparation = AngSeparation(PSRposition, NICERpointing)

        # Separation between Nearby sources and NICER pointing
        SRCseparation = AngSeparation(SRCposition, NICERpointing)

        # New count rates given the NICER pointing
        PSRcountrateScaled = ScaledCtRate(
            PSRseparation.arcmin, PSRcountrate, EffArea, OffAxisAngle
        )
        SRCcountrateScaled = ScaledCtRate(
            SRCseparation.arcmin, NearbySourcesCR, EffArea, OffAxisAngle
        )
        sampleRA[count] = NICERpointing.ra.deg
        sampleDEC[count] = NICERpointing.dec.deg
        PSRrates[count] = PSRcountrateScaled
        SRCrates[count] = np.sum(SRCcountrateScaled)

        # Get S/N ratio given the new count rates for PSR and nearby sources.
        snr[count] = SignaltoNoise(
            PSRcountrateScaled, SRCcountrateScaled, INSTbkgd, EXPtime
        )
        count = count + 1


# Find which offset maximizes the S/N ratio
OptimalPointingIdx = np.where(snr == max(snr))[0][0]
SRCoptimalSEPAR = AngSeparation(
    SRCposition,
    SkyCoord(
        ra=sampleRA[OptimalPointingIdx] * u.degree,
        dec=sampleDEC[OptimalPointingIdx] * u.degree,
    ),
).arcmin
SRCoptimalRATES = ScaledCtRate(SRCoptimalSEPAR, NearbySourcesCR, EffArea, OffAxisAngle)


# Print info for the optimal NICER pointing that maximizes the S/N ratio
print("Target S/N at Optimal Pointing " + str(snr[OptimalPointingIdx]))
print(
    "Target count rate at Optimal pointing = "
    + str(PSRrates[OptimalPointingIdx])
    + " cts/sec"
)
print(
    "Total count rate from nearby sources at Optimal pointing = "
    + str(SRCrates[OptimalPointingIdx])
    + " cts/sec"
)
print("     Individual Nearby SourcesSources: ")
print(str(SRCoptimalRATES))
print("----------------------------------------------------------------------")
print(
    "Optimal Pointing:  "
    + str(sampleRA[OptimalPointingIdx])
    + "  "
    + str(sampleDEC[OptimalPointingIdx])
)
# print("                :  " + str(sampleRA[OptimalPointingIdx].ra.hms) + "  " + str(sampleDEC[OptimalPointingIdx].dec.dms))
print(
    "corresponds to offset of "
    + str(
        AngSeparation(
            PSRposition,
            SkyCoord(
                ra=sampleRA[OptimalPointingIdx] * u.degree,
                dec=sampleDEC[OptimalPointingIdx] * u.degree,
            ),
        ).arcmin
    )
    + " arcmin"
)
print("----------------------------------------------------------------------")

# Plot the map of S/N ratio as function of NICER pointing
fig = plt.figure(figsize=(10, 6.5))
ax = fig.add_subplot(111)
# plt.gca().invert_xaxis()
plt.plot(SRCposition.ra, SRCposition.dec, marker=".", color="black", linestyle="")
plt.plot(PSRposition.ra, PSRposition.dec, marker="*", color="green", linestyle="")
plt.plot(
    sampleRA[OptimalPointingIdx],
    sampleDEC[OptimalPointingIdx],
    marker="o",
    color="red",
    linestyle="",
)


# label of the nearby sources
# n = np.arange(1,len(NearbySourcesCR)+1,1)
# for i, txt in enumerate(n):
#    ax.annotate(str(txt), (NearbySources[0,i],NearbySources[1,i]))

plt.scatter(sampleRA, sampleDEC, c=snr, s=8, edgecolor="")
plt.xlabel("RA", fontsize="large")
plt.ylabel("DEC", fontsize="large")
plt.title("S/N map for " + str(args.name))
cbar = plt.colorbar()
cbar.set_label("S/N")

if args.saveplot:
    plt.savefig("OptimalPoiting_{}.png".format(args.name.replace(" ", "")))
else:
    plt.show()
