# Selecting good data

Good Time Intervals (GTIs) are the times when the data are usable. NICER's
background varies a lot around the orbit, so choosing GTIs carefully matters,
especially for faint sources.

## What psrpipe.py does

For each ObsID, {ref}`psrpipe.py <psrpipe>`:

1. Runs {ref}`nicerql.py <nicerql>` to make quick-look plots of the unscreened
   data (`*_prefilt_*.png`).
2. Makes GTIs from the filter (`.mkf`) file with `nimaketime`, keeping times
   when:
   - the pointing offset `ANG_DIST` is below 0.015° (`--angdist`);
   - the elevation above the Earth limb `ELV` is above 20° and above the
     bright Earth `BR_EARTH` is above 30°;
   - the star tracker is valid and the instrument is tracking the target;
   - the spacecraft is outside the NICER SAA region;
   - at least 7 detector modules are on (`--minfpm`);
   - the undershoot rate per module is below 500 (`--maxundershoot`, or no
     limit with `--nounderfilt`) and the overshoot rate per module is below
     1.5 (`--maxovershoot`).
3. Combines these with the GTIs already in the event file, and drops GTIs
   shorter than 16 s (`--mingti`).
4. Extracts the events in those GTIs with `PI` between `--emin` and `--emax`
   and `EVENT_FLAGS=bx1x000` (slow-only and slow+fast events), removing any
   detectors given with `--mask`. The result is `cleanfilt.evt`.
5. If `--par` is given, adds a `PULSE_PHASE` column with PINT's
   `photonphase`.
6. Extracts a spectrum and light curve with `extractor`.

## Optional cuts

These `psrpipe.py` options add further conditions to the GTIs:

| Option | Keeps times when |
| --- | --- |
| `--dark` | the ISS is in the Earth's shadow (`SUNSHINE=0`) |
| `--minsun X` | the Sun angle is above X degrees, or the ISS is in shadow |
| `--cormin X` | the cutoff rigidity `COR_SAX` is above X GeV (typical values 1.5–2) |
| `--kpmax X` | the geomagnetic Kp index is below X (needs a `KP` column in the `.mkf` file; see {ref}`add_kp.py <add-kp>`) |
| `--medianundershoot X` | the median undershoot rate is below X |
| `--overonlyexpr` | the overshoot rate is below `1.52*COR_SAX**(-0.633)` |
| `--filtpolar` | the ISS is outside the polar horn regions |

`--day` and `--daynight` use the event files that `nicerl2` makes with
`threshfilter=DAY` or `ALL`, which include data taken in daylight after the
2023 light leak.

## Other tools

- {ref}`cr_cut.py <cr-cut>` removes times when the count rate is above a
  threshold, which removes background flares. `psrpipe.py --merge --crcut`
  applies it to the merged event file.
- {ref}`ni_Htest_sortgti.py <ni-htest-sortgti>` sorts GTIs by their
  background rate and finds the set that maximizes the H-test detection
  significance of a pulsar, which can help a lot for faint pulsars.
- {ref}`nioptcuts.py <nioptcuts>` finds the energy range that maximizes the
  detection significance.
- {ref}`mkgti.py <mkgti>` makes a GTI file from a start and stop MET, for use
  with `niextract-events timefile=...`.
- {ref}`splitevt.py <splitevt>` splits an event file into chunks of a given
  length.

In Python, {func}`nicer.plotutils.apply_gti` selects the events in a GTI
table.
