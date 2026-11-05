# Pipelines and data selection

Scripts that screen NICER observations, make GTIs and quick-look plots, and prepare event files for analysis. See {doc}`../guide/filtering` for an overview.

(psrpipe)=
## psrpipe.py

The main processing pipeline. See {doc}`../guide/filtering` for what it does.

```{argparse}
:filename: ../scripts/psrpipe.py
:func: get_parser
:prog: psrpipe.py
```

(nicerql)=
## nicerql.py

Quick-look science, engineering, background and map plots for one ObsID. `psrpipe.py` runs it automatically.

```{argparse}
:filename: ../scripts/nicerql.py
:func: get_parser
:prog: nicerql.py
```

(bkgpipe)=
## bkgpipe.py

```{argparse}
:filename: ../scripts/bkgpipe.py
:func: get_parser
:prog: bkgpipe.py
```

(merge)=
## merge.py

```{argparse}
:filename: ../scripts/merge.py
:func: get_parser
:prog: merge.py
```

(cr-cut)=
## cr_cut.py

```{argparse}
:filename: ../scripts/cr_cut.py
:func: get_parser
:prog: cr_cut.py
```

(mkgti)=
## mkgti.py

```{argparse}
:filename: ../scripts/mkgti.py
:func: get_parser
:prog: mkgti.py
```

(splitevt)=
## splitevt.py

```{argparse}
:filename: ../scripts/splitevt.py
:func: get_parser
:prog: splitevt.py
```

(remove-empty-evtfiles)=
## remove_empty_evtfiles.py

```{argparse}
:filename: ../scripts/remove_empty_evtfiles.py
:func: get_parser
:prog: remove_empty_evtfiles.py
```

(add-kp)=
## add_kp.py

Adds a `KP` column to `.mkf` files, which `psrpipe.py --kpmax` uses.

```{argparse}
:filename: ../scripts/add_kp.py
:func: get_parser
:prog: add_kp.py
```

(calcdeadtime)=
## calcdeadtime.py

```{argparse}
:filename: ../scripts/calcdeadtime.py
:func: get_parser
:prog: calcdeadtime.py
```
