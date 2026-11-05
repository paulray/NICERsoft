# Contributing

Contributions of new tools and fixes are welcome. Please open a pull request at
<https://github.com/paulray/NICERsoft>.

- Executable scripts go in `scripts/`. Scripts in any language are fine, but
  Python scripts get documented automatically (see below).
- Importable Python modules go in `nicer/`.
- Most code uses [loguru](https://loguru.readthedocs.io/) for logging
  (`from loguru import logger as log`) and `astropy.io.fits` (imported as
  `pyfits`) for FITS files.

## Writing a script

Python scripts should build their command-line parser in a `get_parser()`
function and do their work in `main(argv=None)`:

```python
#!/usr/bin/env python
import argparse


def get_parser():
    parser = argparse.ArgumentParser(description="What the script does.")
    parser.add_argument("evfile", help="Input event file")
    parser.add_argument("--emin", help="Minimum energy (keV)", type=float, default=0.25)
    return parser


def main(argv=None):
    args = get_parser().parse_args(argv)
    ...


if __name__ == "__main__":
    main()
```

The documentation is built from `get_parser()`, so the description and help
strings are what users will see; make them complete. Keep the module level free
of anything that needs data, HEASoft or PyXspec: imports are fine, but checks
like `check_call("nicerversion")` belong in `main()`.

Then add the script to the right page in `docs/scripts/`, with a label and an
`argparse` directive like the others:

````markdown
(my-script)=
## my_script.py

```{argparse}
:filename: ../scripts/my_script.py
:func: get_parser
:prog: my_script.py
```
````

The docs build fails with a warning if a script in `scripts/` that has a
`get_parser()` function is missing from these pages.

## Docstrings

Functions and classes in `nicer/` should have
[NumPy-style docstrings](https://numpydoc.readthedocs.io/en/latest/format.html),
which are used to make the {doc}`api/index` pages:

```python
def light_curve(etable, startmet, stopmet, binsize):
    """Bin event times into a light curve.

    Parameters
    ----------
    etable : astropy.table.Table
        Event table with a ``MET`` column, sorted by time.
    binsize : float
        Bin size in seconds.

    Returns
    -------
    bins : numpy.ndarray
        Left edge of each bin, in seconds since ``startmet``.
    sums : numpy.ndarray
        Number of events in each bin.
    """
```

To add a new module to the API pages, add it to the list in
`docs/api/index.md`.

## Building the documentation

The documentation is written in Markdown ([MyST](https://myst-parser.readthedocs.io/))
and built with [Sphinx](https://www.sphinx-doc.org/). To build it:

```bash
pip install -e ".[docs]"
cd docs
make html
open _build/html/index.html
```

`make livehtml` starts a local server that rebuilds the pages and reloads the
browser whenever you save a file.

HEASoft and PyXspec are not needed to build the docs: the build replaces the
modules it can't import with stand-ins (see `autodoc_mock_imports` in
`docs/conf.py`). If you add a script that imports a package that isn't
installed by `pip install -e ".[docs]"`, add it to that list.

The documentation is built on every pull request, and published to GitHub
Pages when changes are merged to `master`.
