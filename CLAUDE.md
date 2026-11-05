# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

NICERsoft is a collection of user-contributed analysis tools for the NICER (Neutron star Interior Composition Explorer) X-ray telescope mission. This is separate from the official NICER team tools in HEASoft. The codebase focuses on pulsar timing analysis, spectral analysis, and various data processing pipelines for NICER X-ray observations.

## Repository Structure

- `scripts/` - Executable command-line tools for data processing and analysis
- `nicer/` - Importable Python modules providing core functionality
- `data/` - Reference data files (KP indices, etc.)

## Prerequisites

Required external dependencies:
- **HEASoft** with NICER Data Analysis Software (required for most functionality)
- **PyXspec** (requires HEASoft compiled from source for spectral analysis scripts)
- Python packages: numpy, matplotlib, scipy, astropy, pint-pulsar, loguru, tqdm

## Installation

The package can be installed using pip:

```bash
# For development (editable install):
pip install -e .

# For regular installation:
pip install .
```

This will:
- Install all Python dependencies automatically
- Make all scripts in `scripts/` available in your PATH (e.g., `psrpipe.py`, `photon_toa.py`)
- Make the `nicer` module importable from anywhere

**Legacy method** (still works but not recommended):
- Add `<basedir>/NICERsoft/scripts` to PATH
- Add `<basedir>/NICERsoft` to PYTHONPATH

## Core Architecture

### Data Processing Flow

1. **Input**: Raw NICER observation directories with standard structure:
   - `xti/event_cl/` - Cleaned event files (`*mpu7_cl.evt*`)
   - `xti/event_uf/` - Unfiltered event files (`*mpu*_uf.evt*`)
   - `auxil/` - Auxiliary files (orbit `.orb`, filter `.mkf`, housekeeping `.hk`)

2. **NicerFileSet Class** (`nicer/NicerFileSet.py`):
   - Central data management class for NICER observations
   - Automatically locates event files, orbit files, filter files, and housekeeping data
   - Handles TIMEZERO corrections and GTI (Good Time Interval) management
   - Provides latitude/longitude interpolation for spacecraft position

3. **Processing Pipelines**:
   - **psrpipe.py**: Main pulsar data processing pipeline
     - Applies GTI filtering based on elevation, pointing, and optionally darkness
     - Filters events by flags (SLOW/FAST), energy range, and detector masks
     - Can compute pulse phases using PINT if a `.par` timing model is provided
     - Outputs cleaned event files and extracted spectra/lightcurves
   - **bkgpipe.py**: Background estimation pipeline

### Key Modules

**nicer/values.py**: Defines fundamental constants and conversions
- MET/GPS time epoch definitions
- PI to energy conversions (PI_TO_KEV, KEV_TO_PI)
- Standard energy bands (SOFT: 0.4-2 keV, HARD: 2-10 keV)
- EVENT_FLAGS bit definitions
- Detector ID arrays and dead detector lists

**nicer/plotutils.py**: Diagnostic plotting utilities
- `find_hot_detectors()`: Automatically identifies noisy detectors using sigma-clipping
- Event count histograms per detector
- Creates diagnostic plots for data quality assessment

**nicer/fillgaps.py**: Gap-filling algorithm
- Fills data gaps while preserving noise characteristics
- Implements Howe & Schlossberger algorithm for frequency stability data

**nicer/fourier.py**: Fourier analysis utilities for periodic signals

**nicer/fitsutils.py**: FITS file handling utilities

### Pulsar Timing Analysis

**photon_toa.py**: Generate pulse Time of Arrivals (TOAs) from unbinned photon events
- Uses PINT pulsar timing library
- Supports NICER, RXTE, NuSTAR, XMM, Swift missions
- Can compute TOAs from template cross-correlation or Fourier harmonic phases
- Outputs TOAs in Tempo2 format
- Includes SNR calculation for generated TOAs

**niphaseogram.py**: Create pulse phase diagrams

**nicerfits2presto.py**: Convert NICER event files to PRESTO format for periodicity searches

### Spectral Analysis (RPP Catalog Suite)

The RPP (Rotating Periodic Pulsars) scripts form a coherent spectral analysis pipeline:

- **RPP-prep.py**: Prepare data for spectral analysis
- **RPP-driver.py**: Main driver for automated spectral fitting
  - Fits blackbody + powerlaw models (1BB, 2BB, etc.)
  - Supports adding Gaussian emission/absorption lines
  - Uses PyXspec for spectral fitting
- **RPP-profile.py**: Generate pulse profile plots (added off-pulse region markers)
- **RPP-spec.py**: Spectral extraction and fitting
- **RPP-pulsed*.py**: Pulsed flux and spectral analysis
- **RPP-xcm-plot*.py**: Plot results from XCM spectral fits

### Data Download and Selection

**ni_data_download.py**: Download NICER data from archives

**get_nicer_aws.py**: Fetch NICER data from AWS

**mkgti.py**: Create Good Time Intervals based on filtering criteria

**nioptcuts.py**: Apply optimal cuts to NICER data

### Utility Scripts

**photon_fit.py**: Fit models to photon arrival times

**cr_cut.py**: Apply count rate cuts to filter high-background periods

**issorb.py**, **issorball.py**: ISS orbit-related utilities

**add_kp.py**: Add Kp geomagnetic index data to observations

**calcdeadtime.py**: Calculate detector dead time

**InteractiveLC.py**: Interactive lightcurve analysis tool

## Common Development Patterns

### Working with Event Tables

Event tables (etables) are Astropy Table objects with these key columns:
- `TIME`: Event times in MET (Mission Elapsed Time)
- `PI`: Pulse Invariant energy channel (convert using PI_TO_KEV)
- `DET_ID`: Detector ID (0-67, with some dead detectors)
- `EVENT_FLAGS`: Bitmask for event classification
- `RAWX`, `RAWY`: Detector coordinates

Convert PI to energy: `energy_kev = etable['PI'] * PI_TO_KEV`

### Time Systems

- **MET**: NICER Mission Elapsed Time, seconds since 2014-01-01T00:00:00 UTC
- **GPS**: GPS time, seconds since 1980-01-06T00:00:00 UTC
- Use `MET0` and `GPS0` from `nicer.values` for conversions
- Always account for `TIMEZERO` in FITS headers when present

### Detector Masking

Hot/noisy detectors should be identified and masked:
```python
from nicer.plotutils import find_hot_detectors
bad_dets = find_hot_detectors(etable)
```

Dead detectors at launch: `{11, 20, 22, 60}` (from `DEAD_DETS`)

### GTI Filtering

Good Time Intervals are critical for clean data:
- Standard cuts: pointing within 0.015 degrees, elevation > 20 degrees
- Optional: filter on sunshine/darkness, overshoot rate, bad ratio events
- Filter files (`.mkf`) contain spacecraft housekeeping for GTI selection

## Important Conventions

- Scripts are designed to be run from command line with argparse interfaces. Each Python script builds its parser in `get_parser()` and runs in `main(argv=None)` behind an `if __name__ == "__main__":` guard; the docs are generated from `get_parser()`, so keep module-level code free of anything needing data, HEASoft or PyXspec
- New scripts must be added to a page in `docs/scripts/` (the docs build warns otherwise); `nicer/` code uses NumPy-style docstrings
- Most scripts use `loguru` logger (imported as `log`) instead of standard logging
- FITS file operations use `astropy.io.fits` (imported as `pyfits`)
- Phases are typically computed using the PINT pulsar timing package
- Energy ranges default to 0.25-12.0 keV unless otherwise specified

## Git Commit Style

Recent commits show a concise, imperative style:
- "Fix snr definition"
- "Add an SNR flag to generated TOAs"
- "Added ability for photon_toa.py to compute TOAs based on the phases of Fourier harmonics"

Focus on what was changed and why, using imperative mood for the first line.

## Documentation

Documentation is built with Sphinx (MyST Markdown, sphinx_rtd_theme) from `docs/` and published to GitHub Pages (https://paulray.github.io/NICERsoft/) by `.github/workflows/docs.yml` on pushes to master. Build locally with `pip install -e ".[docs]"` and `make html` (or `make livehtml`) in `docs/`.

Additional useful tools from George Younes: https://github.com/georgeyounes/NICERUTIL
