# Quick start

This walks through a typical pulsar analysis: screen some observations, fold
them with a timing model, make a pulse template, and measure TOAs. It assumes
HEASoft is initialized and NICERsoft is installed (see {doc}`installation`).

## 1. Get the data

Download the observations (ObsIDs) you want from the HEASARC archive, for
example with {ref}`get_nicer_aws.py <get-nicer-aws>`:

```bash
get_nicer_aws.py --srcname "PSR J0030+0451" --outdir j0030
```

Each ObsID directory has the standard layout described in {doc}`guide/data`.
NICERsoft starts from the products of the standard `nicerl2` processing, so
run `nicerl2` on the ObsIDs first if they don't already have them.

## 2. Screen the data and compute pulse phases

{ref}`psrpipe.py <psrpipe>` selects good times and events for each ObsID, makes
quick-look plots, and, if you give it a timing model, adds a `PULSE_PHASE`
column:

```bash
psrpipe.py --par J0030+0451.par j0030/[0-9]*
```

For each ObsID this makes a directory `<obsid>_pipe/` containing:

- `cleanfilt.evt`: the screened events (with `PULSE_PHASE` if `--par` was given)
- `cleanfilt.pha` and `cleanfilt.lc`: a spectrum and light curve made with
  `extractor`
- quick-look plots (`*_prefilt_*.png`) from {ref}`nicerql.py <nicerql>`
- a copy of the orbit file

See {doc}`guide/filtering` for the screening criteria and the options that
control them. Add `--merge` to also merge all of the ObsIDs into a single event
file in `merged_pipe/`.

## 3. Look at the pulse profile

```bash
niphaseogram.py merged_pipe/cleanfilt.evt
```

This plots the phases against time, which shows quickly whether the timing
model is good.

## 4. Make a template

Fit a multi-Gaussian template to the profile interactively:

```bash
nitemplate.py --prof j0030 merged_pipe/cleanfilt.evt
```

This writes the template to `j0030.gauss`.

## 5. Measure TOAs

```bash
photon_toa.py --orbfile 1060020101_pipe/ni1060020101.orb.gz \
    --tint 50000 --outfile j0030.tim \
    1060020101_pipe/cleanfilt.evt j0030.gauss J0030+0451.par
```

This writes TOAs in Tempo2 format, one per 50 ks of exposure. If the event file
already has a `PULSE_PHASE` column, `--phasecol PULSE_PHASE` uses those phases
instead of recomputing them, which is much faster for large files. See
{doc}`guide/timing` for more.
