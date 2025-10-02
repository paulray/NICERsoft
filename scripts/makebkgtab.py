#!/usr/bin/env python

import numpy as np
from astropy.table import Table
import argparse
from astropy import log
from astropy.table import Table, vstack
from os import path
from glob import glob
from nicer.values import *

parser = argparse.ArgumentParser(
    description="Construct BKG table from a pipe directory"
)
parser.add_argument(
    "pipedirs",
    help="Input _pipe directory process (should be a blank field BKG observation)",
    nargs="+",
)
parser.add_argument("--plot", action="store_true", help="Plot")
# parser.add_argument("--emin", help="Minimum energy to include (keV, default=0.25)", type=float, default=0.25)
args = parser.parse_args()

# Length of a chunk
chunklen = 30.0

if args.plot:
    import matplotlib.pyplot as plt

for pipedir in args.pipedirs:
    mkfname = glob(path.join(pipedir, "ni*.mkf"))[0]
    evtname = glob(path.join(pipedir, "cleanfilt.evt"))[0]
    outname = evtname.replace(".evt", ".bkgtab")
    #    outfile = open(evtname.replace('.evt','.bkgtab'),'w')
    #    print("{0:15s} {1:6s} {2:6s} {3:6s} {4:8s} {5:5s} {6:8s} {7:9s} {8:9s} {9:9s} {10:9s} {11:9s} {12:9s}".format(
    #    "# MET","Band1", "Band2", "Band3", "CORSAX", "SUN", "SUNANG", "OVERONLY","NOISE25","RATIOREJ","R1517", "IBG", "HREJ"),
    #    file=outfile)
    log.info("Processing {0} and {1}".format(evtname, mkfname))

    # Collect GTIs
    gtitable = Table.read(evtname, hdu=2)
    if "TIMEZERO" in gtitable.meta:
        tz = gtitable.meta["TIMEZERO"]
        # If there are multiple TIMEZERO entries in the header, just take the last
        if not np.isscalar(tz):
            tz = tz[-1]
        log.info("Applying TIMEZERO of {0} to gtitable in NicerFileSet".format(tz))
        gtitable["START"] += tz
        gtitable["STOP"] += tz
        gtitable.meta["TIMEZERO"] = 0.0
    log.info("Got the good times from GTI")
    gtitable["DURATION"] = gtitable["STOP"] - gtitable["START"]
    # Only keep GTIs longer than chunklen
    log.info("Discarding GTI shorter than {0} seconds".format(chunklen))
    idx = np.where(gtitable["DURATION"] > chunklen)[0]
    gtitable = gtitable[idx]

    if gtitable["DURATION"].sum() < chunklen:
        log.info("Skipping file with no good time.")
        continue
    # print(gtitable)

    # Read evt and mkf files into Tables
    etable = Table.read(evtname, hdu=1)
    if "TIMEZERO" in etable.meta:
        log.info("Applying TIMEZERO of {0} to etable".format(etable.meta["TIMEZERO"]))
        etable["TIME"] += etable.meta["TIMEZERO"]
        etable.meta["TIMEZERO"] = 0.0

    mktable = Table.read(mkfname, hdu=1)
    if "TIMEZERO" in mktable.meta:
        log.info("Applying TIMEZERO of {0} to mktable".format(mktable.meta["TIMEZERO"]))
        mktable["TIME"] += mktable.meta["TIMEZERO"]
        mktable.meta["TIMEZERO"] = 0.0

    # Make selectors for energy bands
    emin = 0.25
    emax = 0.4
    b1 = etable["PI"] > emin / PI_TO_KEV
    b2 = etable["PI"] < emax / PI_TO_KEV
    band1idx = np.where(b1 & b2)[0]

    emin = 0.4
    emax = 2.0
    b1 = etable["PI"] > emin / PI_TO_KEV
    b2 = etable["PI"] < emax / PI_TO_KEV
    band2idx = np.where(b1 & b2)[0]

    emin = 2.0
    emax = 8.0
    b1 = etable["PI"] > emin / PI_TO_KEV
    b2 = etable["PI"] < emax / PI_TO_KEV
    band3idx = np.where(b1 & b2)[0]

    del b1, b2

    firstgti = True
    # Loop over GTIs
    for curgti in gtitable:
        # Split each GTI into chunks
        nchunks = int(np.floor(curgti["DURATION"] / chunklen))
        print(curgti["DURATION"], nchunks, nchunks * chunklen)
        chunkmets = curgti["START"] + np.arange(nchunks) * chunklen
        # Use histogram to make light curves in energy bands 0.25-0.4, 0.4-2.0, 2.0-8.0
        chunkbins = np.append(chunkmets, chunkmets[-1] + chunklen)
        band1lc, edges = np.histogram(etable["TIME"][band1idx], bins=chunkbins)
        band2lc, edges = np.histogram(etable["TIME"][band2idx], bins=chunkbins)
        band3lc, edges = np.histogram(etable["TIME"][band3idx], bins=chunkbins)

        # Conver lcs to rates
        band1lc = np.array(band1lc, dtype=float) / chunklen
        band2lc = np.array(band2lc, dtype=float) / chunklen
        band3lc = np.array(band3lc, dtype=float) / chunklen

        # Extract MKF variables and take mean (or other appropriate combination) in each chunk

        # Add these here: ELV, BR_EARTH, SAA_TIME, SAA, FPM_UNDERONLY_COUNT, SAT_LAT, SAT_LON

        cor_sax = np.zeros(len(band1lc))
        sun_angle = np.zeros(len(band1lc))
        sunshine = np.zeros(len(band1lc))
        overonly = np.zeros(len(band1lc))
        underonly = np.zeros(len(band1lc))
        noise25 = np.zeros(len(band1lc))
        ratiorej = np.zeros(len(band1lc))
        rate1517 = np.zeros(len(band1lc))
        ibg = np.zeros(len(band1lc))
        hrej = np.zeros(len(band1lc))
        kp = np.zeros(len(band1lc))
        lat = np.zeros(len(band1lc))
        lon = np.zeros(len(band1lc))
        ra = np.zeros(len(band1lc))
        dec = np.zeros(len(band1lc))
        saa_time = np.zeros(len(band1lc))
        mag_angle = np.zeros(len(band1lc))
        i = 0
        for chmet in chunkmets:
            b1 = mktable["TIME"] > chmet
            b2 = mktable["TIME"] < chmet + chunklen
            chidx = np.where(b1 & b2)[0]
            # Now define the bkg predictor quantities from the MKF file
            cor_sax[i] = mktable["COR_SAX"][chidx].mean()
            sun_angle[i] = mktable["SUN_ANGLE"][chidx].mean()
            sunshine[i] = mktable["SUNSHINE"][chidx].mean()
            # For count rates, convert from perFPM to full NICER count rate
            overonly[i] = 52 * mktable["FPM_OVERONLY_COUNT"][chidx].mean()
            underonly[i] = 52 * mktable["FPM_UNDERONLY_COUNT"][chidx].mean()
            noise25[i] = 52 * mktable["FPM_NOISE25_COUNT"][chidx].mean()
            ratiorej[i] = 52 * mktable["FPM_RATIO_REJ_COUNT"][chidx].mean()
            rate1517[i] = 52 * mktable["FPM_XRAY_PI_1500_1700"][chidx].mean()
            kp[i] = mktable["KP"][chidx].mean()
            lat[i] = mktable["SAT_LAT"][chidx].mean()
            lon[i] = mktable["SAT_LON"][
                chidx
            ].min()  # Use min() here to prevent problems when LON wraps
            dec[i] = mktable["DEC"][chidx].mean()
            ra[i] = mktable["RA"][
                chidx
            ].min()  # Use min() here to prevent problems when RA wraps
            saa_time[i] = mktable["SAA_TIME"][chidx].min()  # Using min()
            mag_angle[i] = mktable["MAG_ANGLE"][chidx].mean()

            # FPM_TRUMP_SEL_1500_1800 (IBG)
            ibg[i] = mktable["FPM_TRUMP_SEL_1500_1800"][chidx].sum() / chunklen

            # FPM_RATIO_REJ_300_1800 (HREJ)
            hrej[i] = mktable["FPM_RATIO_REJ_300_1800"][chidx].sum() / chunklen
            i += 1

        if firstgti:
            f_chunkmets = chunkmets
            f_band1lc = band1lc
            f_band2lc = band2lc
            f_band3lc = band3lc
            f_cor_sax = cor_sax
            f_overonly = overonly
            f_underonly = underonly
            f_sun_angle = sun_angle
            f_sunshine = sunshine
            f_noise25 = noise25
            f_ratiorej = ratiorej
            f_rate1517 = rate1517
            f_ibg = ibg
            f_hrej = hrej
            f_kp = kp
            f_lat = lat
            f_lon = lon
            f_ra = ra
            f_dec = dec
            f_saa_time = saa_time
            f_mag_angle = mag_angle
            firstgti = False
        else:
            f_chunkmets = np.append(f_chunkmets, chunkmets)
            f_band1lc = np.append(f_band1lc, band1lc)
            f_band2lc = np.append(f_band2lc, band2lc)
            f_band3lc = np.append(f_band3lc, band3lc)
            f_cor_sax = np.append(f_cor_sax, cor_sax)
            f_overonly = np.append(f_overonly, overonly)
            f_underonly = np.append(f_underonly, underonly)
            f_sun_angle = np.append(f_sun_angle, sun_angle)
            f_sunshine = np.append(f_sunshine, sunshine)
            f_noise25 = np.append(f_noise25, noise25)
            f_ratiorej = np.append(f_ratiorej, ratiorej)
            f_rate1517 = np.append(f_rate1517, rate1517)
            f_ibg = np.append(f_ibg, ibg)
            f_hrej = np.append(f_hrej, hrej)
            f_kp = np.append(f_kp, kp)
            f_lat = np.append(f_lat, lat)
            f_lon = np.append(f_lon, lon)
            f_ra = np.append(f_ra, ra)
            f_dec = np.append(f_dec, dec)
            f_saa_time = np.append(f_saa_time, saa_time)
            f_mag_angle = np.append(f_mag_angle, mag_angle)

    f_obsid = np.array([etable.meta["OBS_ID"]] * len(f_band1lc))
    bkgtab = Table(
        [
            f_chunkmets,
            f_band1lc,
            f_band2lc,
            f_band3lc,
            f_cor_sax,
            f_sunshine,
            f_sun_angle,
            f_overonly,
            f_underonly,
            f_noise25,
            f_ratiorej,
            f_rate1517,
            f_ibg,
            f_hrej,
            f_kp,
            f_lat,
            f_lon,
            f_ra,
            f_dec,
            f_saa_time,
            f_obsid,
            f_mag_angle,
        ],
        names=(
            "met",
            "band1",
            "band2",
            "band3",
            "cor_sax",
            "sunshine",
            "sun_angle",
            "overonly",
            "underonly",
            "noise25",
            "ratiorej",
            "rate1517",
            "ibg",
            "hrej",
            "kp",
            "lat",
            "lon",
            "ra",
            "dec",
            "saa_time",
            "obsid",
            "magangle",
        ),
    )
    bkgtab.write(outname, format="fits", overwrite=True)

    if args.plot:
        fig, ax = plt.subplots()
        ax.plot(chunkmets - chunkmets[0], band1lc)
        plt.show()
