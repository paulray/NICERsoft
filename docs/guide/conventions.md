# Times, energies and detectors

The constants mentioned here are defined in {mod}`nicer.values`.

## Time systems

**MET** (Mission Elapsed Time) is the time used in NICER files: seconds since
2014-01-01T00:00:00 UTC, available as {data}`nicer.values.MET0`. The exact
reference time and time system of a file are given by its `MJDREFI`,
`MJDREFF` and `TIMESYS` header keywords, which should be used for precise
conversions.

**GPS time**, used in some housekeeping, counts seconds since
1980-01-06T00:00:00 UTC ({data}`nicer.values.GPS0`).

### TIMEZERO

NICER files can have a `TIMEZERO` header keyword, which must be added to the
`TIME` column (and to GTI `START` and `STOP` times) to get the correct time.
The `nicer` functions that read files add `TIMEZERO` for you. If you read files
yourself, check for it:

```python
if "TIMEZERO" in etable.meta:
    etable["TIME"] += etable.meta["TIMEZERO"]
    etable.meta["TIMEZERO"] = 0.0
```

## Energies

`PI` channels are in units of 10 eV, so to convert to keV multiply by
{data}`nicer.values.PI_TO_KEV` (0.01); {data}`nicer.values.KEV_TO_PI` converts
the other way. Unless told otherwise, the pipelines keep events from 0.25 keV
up to the value of `--emax`.

The RPP catalog analysis uses a soft band of 0.4–2 keV
({data}`~nicer.values.SOFT_EMIN`, {data}`~nicer.values.SOFT_EMAX`) and a hard
band of 2–10 keV ({data}`~nicer.values.HARD_EMIN`,
{data}`~nicer.values.HARD_EMAX`).

## Event flags

The `EVENT_FLAGS` column is a set of bits. Astropy reads it as a boolean array
with one row per event, and {mod}`nicer.values` names the bits:

| Constant | Bit | Set when |
| --- | --- | --- |
| `FLAG_FIRST_MPU` | 2 | the "first MPU" flag is set |
| `FLAG_SLOW` | 3 | the slow chain triggered |
| `FLAG_FAST` | 4 | the fast chain triggered |
| `FLAG_SWTRIG` | 5 | software (forced) trigger |
| `FLAG_OVERSHOOT` | 6 | overshoot (large energy deposit, usually a particle) |
| `FLAG_UNDERSHOOT` | 7 | undershoot (detector reset, mostly from optical loading) |

For example, to select events that triggered both chains:

```python
from nicer.values import FLAG_SLOW, FLAG_FAST

both = etable["EVENT_FLAGS"][:, FLAG_SLOW] & etable["EVENT_FLAGS"][:, FLAG_FAST]
```

The pipelines select events with the FTOOLS expression
`EVENT_FLAGS=bx1x000`, which keeps slow-only and slow+fast events and removes
overshoots, undershoots and forced triggers.

## Detectors

NICER has 56 detectors (focal plane modules), read out by 7 MPUs. A detector's
`DET_ID` is 10 times its MPU number plus its number within the MPU, so the IDs
run 0–7, 10–17, ... 60–67 ({data}`nicer.values.IDS`). Detectors 11, 20, 22
and 60 were dead at launch ({data}`nicer.values.DEAD_DETS`).

Noisy detectors can be found with {func}`nicer.plotutils.find_hot_detectors`,
which flags detectors with far more events than the rest, and removed with the
`--mask` option of {ref}`psrpipe.py <psrpipe>`:

```python
from nicer.plotutils import find_hot_detectors

bad_dets = find_hot_detectors(etable)
```
