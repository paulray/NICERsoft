# NICER data

## Observation directories

NICER data are organized by observation ID (ObsID), a 10-digit number. Each
ObsID directory downloaded from the HEASARC has this structure:

```text
<obsid>/
├── auxil/
│   ├── ni<obsid>.orb          orbit file
│   ├── ni<obsid>.mkf          filter file (spacecraft and instrument housekeeping)
│   └── ni<obsid>_apid0260.hk  GPS (SPS) housekeeping
└── xti/
    ├── event_cl/
    │   ├── ni<obsid>_0mpu7_cl.evt    cleaned events (from nicerl2)
    │   └── ni<obsid>_0mpu7_ufa.evt   merged, unfiltered events (from nicerl2)
    ├── event_uf/
    │   └── ni<obsid>_0mpu?_uf.evt    raw unfiltered events, one per MPU
    └── hk/                           MPU housekeeping
```

Files may be gzipped (`.gz`); NICERsoft finds them either way.

The pipelines in NICERsoft start from the `ufa` (merged, unfiltered) event
files and the `.mkf` filter file, both made by the standard `nicerl2`
processing, and do their own screening (see {doc}`filtering`).

## Reading an observation in Python

{class}`nicer.NicerFileSet.NicerFileSet` finds the files in an ObsID directory
and loads the events, GTIs and filter file. It is what
{ref}`nicerql.py <nicerql>` uses. It takes an {class}`argparse.Namespace` of
options, normally from `nicerql.py`; see its documentation for the attributes
it uses.

## Event tables

Most of the `nicer` functions work on event tables: {class}`astropy.table.Table`
objects read from the `EVENTS` extension of an event file. The main columns are:

| Column | Meaning |
| --- | --- |
| `TIME` | Event time, in MET seconds (see {doc}`conventions`). Some functions expect it to be renamed `MET`. |
| `PI` | Energy channel, in units of 10 eV. |
| `PI_FAST` | Energy channel from the fast chain. |
| `DET_ID` | Detector that recorded the event (see {doc}`conventions`). |
| `EVENT_FLAGS` | Bit flags describing the event. |
| `RAWX`, `RAWY` | Detector position in the focal plane. |
| `PULSE_PHASE` | Pulse phase, added by `photonphase` (for example via `psrpipe.py --par`). |

For example, to read the events and convert their energies to keV:

```python
from astropy.table import Table
from nicer.values import PI_TO_KEV

etable = Table.read("cleanfilt.evt", hdu="EVENTS")
energy_kev = etable["PI"] * PI_TO_KEV
```
