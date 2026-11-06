# Spectral analysis (RPP)

Scripts for the rotation-powered pulsar catalog spectral analysis. They need HEASoft and PyXspec. See {doc}`../guide/rpp` for how they fit together.

(rpp-prep)=
## RPP-prep.py

```{argparse}
:filename: ../scripts/RPP-prep.py
:func: get_parser
:prog: RPP-prep.py
```

(rpp-profile)=
## RPP-profile.py

```{argparse}
:filename: ../scripts/RPP-profile.py
:func: get_parser
:prog: RPP-profile.py
```

(rpp-pulsedspec)=
## RPP-pulsedspec.py

Makes on-pulse and off-pulse spectra, using the off-pulse phase range found by `RPP-profile.py`.

```{argparse}
:filename: ../scripts/RPP-pulsedspec.py
:func: get_parser
:prog: RPP-pulsedspec.py
```

(rpp-driver)=
## RPP-driver.py

Fits the standard set of models to the total spectrum, with the SCORPEON background model, by running `RPP-spec.py` for each.

```{argparse}
:filename: ../scripts/RPP-driver.py
:func: get_parser
:prog: RPP-driver.py
```

(rpp-spec)=
## RPP-spec.py

Fits one model to the total spectrum, with the SCORPEON background model.

```{argparse}
:filename: ../scripts/RPP-spec.py
:func: get_parser
:prog: RPP-spec.py
```

(rpp-spec-from-xcm)=
## RPP-spec-from-xcm.py

```{argparse}
:filename: ../scripts/RPP-spec-from-xcm.py
:func: get_parser
:prog: RPP-spec-from-xcm.py
```

(rpp-xcm-plot)=
## RPP-xcm-plot.py

```{argparse}
:filename: ../scripts/RPP-xcm-plot.py
:func: get_parser
:prog: RPP-xcm-plot.py
```

(rpp-driver-onoff)=
## RPP-driver-onoff.py

Fits the standard set of models to the pulsed (on-pulse minus off-pulse) spectrum, by running `RPP-spec-onoff.py` for each.

```{argparse}
:filename: ../scripts/RPP-driver-onoff.py
:func: get_parser
:prog: RPP-driver-onoff.py
```

(rpp-spec-onoff)=
## RPP-spec-onoff.py

Fits one model to the pulsed (on-pulse minus off-pulse) spectrum.

```{argparse}
:filename: ../scripts/RPP-spec-onoff.py
:func: get_parser
:prog: RPP-spec-onoff.py
```

(rpp-xcm-plot-onoff)=
## RPP-xcm-plot-onoff.py

```{argparse}
:filename: ../scripts/RPP-xcm-plot-onoff.py
:func: get_parser
:prog: RPP-xcm-plot-onoff.py
```

(rpp-scorpeon-errorsfluxes)=
## RPP-scorpeon_errorsfluxes.py

```{argparse}
:filename: ../scripts/RPP-scorpeon_errorsfluxes.py
:func: get_parser
:prog: RPP-scorpeon_errorsfluxes.py
```

(rpp-scorpeon-yml)=
## RPP-scorpeon_yml.py

```{argparse}
:filename: ../scripts/RPP-scorpeon_yml.py
:func: get_parser
:prog: RPP-scorpeon_yml.py
```

(rpp-pulsed-errorsfluxes)=
## RPP-pulsed_errorsfluxes.py

```{argparse}
:filename: ../scripts/RPP-pulsed_errorsfluxes.py
:func: get_parser
:prog: RPP-pulsed_errorsfluxes.py
```

(rpp-pulsed-yml)=
## RPP-pulsed_yml.py

```{argparse}
:filename: ../scripts/RPP-pulsed_yml.py
:func: get_parser
:prog: RPP-pulsed_yml.py
```

(plotsrcbkg)=
## plotsrcbkg.py

```{argparse}
:filename: ../scripts/plotsrcbkg.py
:func: get_parser
:prog: plotsrcbkg.py
```
