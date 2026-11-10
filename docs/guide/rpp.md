# RPP spectral analysis

The `RPP-*.py` scripts are the spectral analysis pipeline for the NICER
catalog of rotation-powered pulsars (RPPs). They fit the total and pulsed
spectra of a pulsar with blackbody and power-law models. They need HEASoft and
PyXspec.

## Preparing the data

{ref}`RPP-prep.py <rpp-prep>` does all of the preparation, starting from a
directory of ObsIDs and a timing model. Run it in the directory where you want
the output:

```bash
RPP-prep.py --srcname J0030 J0030+0451.par /data/j0030
```

It:

1. Runs {ref}`psrpipe.py <psrpipe>` on every ObsID, with the cuts used for the
   catalog (`--nomap --emin 0.22 --emax 15.01 --tidy --cormin 1.5 --kpmax 5 --mingti 100 --maxovershoot 1.5 --maxundershoot 600 --medianundershoot=100`), to make phase-tagged event files.
2. Merges the event files (`merged.evt`) and makes a 2–10 keV light curve.
3. Removes background flares with {ref}`cr_cut.py <cr-cut>`, cutting at a
   percentile of the count-rate distribution (`--cps_cut_percentile`, default 95%), to make `merged_cut.evt`.
4. Plots a phaseogram with {ref}`niphaseogram.py <niphaseogram>`.
5. Makes a spectrum and a SCORPEON background model with `nicerl3-spect`,
   writing XSPEC and PyXspec load scripts (`load_pulsedspec.xcm` and `load_pulsedspec.py`, respectively).
6. Finds the energy range with the most significant pulsations, with
   {ref}`nioptcuts.py <nioptcuts>`, unless `--emin` and `--emax` are given.
7. Analyzes the pulse profile with {ref}`RPP-profile.py <rpp-profile>`, which
   finds the off-pulse phase range with Bayesian blocks and writes the results
   to a YAML file (`<psr name>_profinfo.yml`). It also estimates the TOA
   accuracy for a 5 ks observation: the Cramér–Rao lower bound computed from
   the measured profile (including background and unpulsed emission), in the
   H-test band (`toa_*` keys in `optres`) and in the energy band that
   minimizes it (`toares`). Values are in phase units; multiply by the period
   for time. `--toa-mc` adds a Monte Carlo check of the bound, which fails
   for faint pulsars and double-peaked profiles.
8. Makes on-pulse and off-pulse spectra with
   {ref}`RPP-pulsedspec.py <rpp-pulsedspec>`, for fitting the pulsed spectrum.

## Fitting the total spectrum

{ref}`RPP-driver.py <rpp-driver>` runs {ref}`RPP-spec.py <rpp-spec>` for each
of the standard models (absorbed blackbody, power law, two blackbodies,
blackbody plus power law, and two blackbodies plus power law), with the
SCORPEON background model. `RPP-spec.py` can also add Gaussian emission or
absorption lines. Each fit makes a plot and a YAML file of results, including
fluxes.

If you would rather do a fit by hand in XSPEC, save it with `save all
mymodel.xcm` and use {ref}`RPP-spec-from-xcm.py <rpp-spec-from-xcm>` to refit the spectrum and plot or {ref}`RPP-xcm-plot.py <rpp-xcm-plot>` to plot only, to make the same products.

## Fitting the pulsed spectrum

The pulsed spectrum is the on-pulse spectrum with the off-pulse spectrum as
background. {ref}`RPP-driver-onoff.py <rpp-driver-onoff>` fits the standard
models to it with {ref}`RPP-spec-onoff.py <rpp-spec-onoff>`, and
{ref}`RPP-xcm-plot-onoff.py <rpp-xcm-plot-onoff>` handles a fit saved from
XSPEC.

## Errors and fluxes from XSPEC logs

For fits done by hand in XSPEC:

- {ref}`RPP-scorpeon_errorsfluxes.py <rpp-scorpeon-errorsfluxes>` and
  {ref}`RPP-pulsed_errorsfluxes.py <rpp-pulsed-errorsfluxes>` read the log of
  a best fit and write an XSPEC script that computes parameter errors,
  `steppar` scans and fluxes.
- {ref}`RPP-scorpeon_yml.py <rpp-scorpeon-yml>` and
  {ref}`RPP-pulsed_yml.py <rpp-pulsed-yml>` collect the results from the log
  of that script into a YAML file.
