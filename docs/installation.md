# Installation

## Prerequisites

- **Python 3.8 or later.**
- **[HEASoft](https://heasarc.gsfc.nasa.gov/docs/software/heasoft/)** with the
  NICER Data Analysis Software, initialized (for example with `heainit`) in the
  shell you run the scripts from. Most of the pipelines call FTOOLS such as
  `nimaketime`, `niextract-events` and `extractor`.
- **PyXspec**, for the spectral-analysis (RPP) scripts. PyXspec needs HEASoft
  to be built from source.
- **[cartopy](https://scitools.org.uk/cartopy/)**, for the scripts that draw
  maps (`psrpipe.py`, `plothkmap.py`, `quickmap.py` and others). It can be
  hard to install with pip; on a Mac with MacPorts,
  `sudo port install py311-cartopy` (with your Python version) works, as does a
  `conda` install. `psrpipe.py --nomap` skips the maps if cartopy is not
  available.

A few scripts need other packages; they will tell you if one is missing.

## Installing NICERsoft

Clone the repository and install it with pip:

```bash
git clone https://github.com/paulray/NICERsoft.git
cd NICERsoft
pip install -e .
```

This installs the Python dependencies (numpy, scipy, matplotlib, astropy,
[PINT](https://nanograv-pint.readthedocs.io/), loguru, tqdm), puts the scripts
from `scripts/` on your `PATH`, and makes the `nicer` package importable. The
`-e` (editable) flag means changes to the repository take effect without
reinstalling.

### Without pip

You can also use NICERsoft without installing it: add `<dir>/NICERsoft/scripts`
to your `PATH` and `<dir>/NICERsoft` to your `PYTHONPATH`, where `<dir>` is
wherever you cloned it. You will need to install the dependencies yourself.

## Checking the installation

Every script prints its options with `--help`, for example:

```bash
psrpipe.py --help
```
