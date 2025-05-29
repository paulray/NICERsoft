#!/usr/bin/env python
"""
Module to download NICER data from the heasarc copy hosted at AWS S3

Can be run through the CLI via the command getniceraws. The usual getniceraws -h will provide a list of arguments
Example usage:

getniceraws -sn "psr b0540-69" -od ./psrb0540/nicer/ -s 2022-12-29 -e 2025-01-25

where -sn (--srcname) and -od (--outdir) are required arguments, and represent the name of the source and the output
directory where the observation ids will be downloaded. The -sn should match a known source name; it will be parsed by
astropy SkyCoord for RA and DEC coordinates, which will be passed to Heasarc().query_region. Otherwise, the user could
provide the RA and DEC directly through the argument -rd (--radec; e.g., -rd 85.046667 -69.331722). The -od is where
the observation ids will be downloaded. Note that the script does not create any directories by default
(e.g. srcname), except the observation ids themselves. If in the -od an observation_id directory exists and is not
empty, the script will skip the download of that observation. Finally, the -s (--start) and -e (--end) allow the user
to specify time range for the download, rather than the full nicer list of observations (default).
"""

import os
import subprocess
import argparse

from astroquery.heasarc import Heasarc
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time


def query_nicer_observations(srcname=None, ra=None, dec=None, radius_deg=0.2):
    heasarc = Heasarc()
    if srcname:
        coord = SkyCoord.from_name(srcname)
    elif ra is not None and dec is not None:
        coord = SkyCoord(ra, dec, unit="deg")
    else:
        raise ValueError("Provide either srcname or RA/DEC.")
    result = heasarc.query_region(coord, catalog='nicermastr', radius=radius_deg * u.deg)
    return result


def filter_by_date(obs_table, start_date, end_date):
    """Filter Astropy table between two ISO date strings."""
    t_start = Time(start_date).mjd
    t_end = Time(end_date).mjd
    return obs_table[(obs_table['time'] >= t_start) & (obs_table['time'] <= t_end)]


def get_year_month_from_mjd(mjd_array):
    """Convert an array of MJDs to list of 'YYYY_MM' strings."""
    times = Time(mjd_array, format="mjd").to_datetime()
    return [t.strftime("%Y_%m") for t in times]


def download_obs_from_aws(obsid, year_month, output_base):
    """
    Download NICER observation from AWS public bucket, unless it already exists or doesn't exist remotely.
    """
    obsid_str = str(obsid)
    output_dir = os.path.join(output_base, obsid_str)

    if os.path.isdir(output_dir) and any(os.scandir(output_dir)):
        print(f"[SKIP] {obsid_str} already exists at {output_dir}")
        return

    s3_path = f"s3://nasa-heasarc/nicer/data/obs/{year_month}/{obsid_str}/"

    # Check if OBSID exists on AWS
    check_cmd = [
        "aws", "s3", "ls", "--no-sign-request", s3_path
    ]
    check_result = subprocess.run(check_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if check_result.returncode != 0 or check_result.stdout.strip() == b"":
        print(f"[WARN] OBSID {obsid_str} not found on AWS site {s3_path} — skipping")
        return

    # Proceed to create output dir and download
    os.makedirs(output_dir, exist_ok=True)
    cmd = [
        "aws", "s3", "cp", "--no-sign-request", "--recursive",
        s3_path, output_dir
    ]

    print(f"[INFO] Downloading {obsid_str} from {s3_path}")
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if result.returncode != 0:
        print(f"[ERROR] Failed to download {obsid_str} from {s3_path}:\n{result.stderr.decode()}")

    else:
        print(f"[OK] Downloaded {obsid_str}")

    return


def main():
    parser = argparse.ArgumentParser(
        description="Download NICER observations from AWS for a given sourcename or coordinates and time range."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-sn", "--srcname", type=str, help="Target source name (e.g., 'PSR B0531+21')")
    group.add_argument("-rd", "--radec", nargs=2, type=float, metavar=('RA', 'DEC'),
                       help="RA and DEC in degrees (e.g., --radec 83.63 22.01)")
    parser.add_argument("-od", "--outdir", type=str, required=True,
                        help="Output directory for downloaded observation ids (e.g., 'psr_b0531+21')")
    parser.add_argument("-r", "--radius", type=float, default=0.05,
                        help="Search radius in degrees (default: 0.05)")
    parser.add_argument("-s", "--start", type=str, default="2017-01-01",
                        help="Start date (YYYY-MM-DD; default: 2017-01-01)")
    parser.add_argument("-e", "--end", type=str, default="2032-12-31",
                        help="End date (YYYY-MM-DD; default: 2032-12-31)")
    args = parser.parse_args()

    # Query NICER observations
    if args.srcname:
        result_table = query_nicer_observations(srcname=args.srcname, radius_deg=args.radius)
    else:
        ra, dec = args.radec
        result_table = query_nicer_observations(ra=ra, dec=dec, radius_deg=args.radius)

    if len(result_table) == 0:
        print("No observations found for the specified source.")
        return

    # Filter by time range
    timeflt_table = filter_by_date(result_table, args.start, args.end)
    if len(timeflt_table) == 0:
        print("No observations found in that time range.")
        return
    elif len(timeflt_table) > 1:
        print("Observations found for source in specified time range: \n")
        print(timeflt_table)
        print("\n")

    # Get YYYY_MM strings
    year_months = get_year_month_from_mjd(timeflt_table['time'].value)

    # Download each observation
    for row, ym in zip(timeflt_table, year_months):
        download_obs_from_aws(obsid=row['obsid'],
                              year_month=ym,
                              output_base=args.outdir)


if __name__ == "__main__":
    main()
