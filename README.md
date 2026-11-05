# NICERsoft
User-contributed analysis tools for the NICER mission. This is separate from the
NICER-team-developed tools that are part of the HEASoft distribution. Nothing
here should be construed as formally recommended by the NICER instrument team.

Please contribute tools you develop!

(Also note the George Younes has developed some useful tools, see here: https://github.com/georgeyounes/NICERUTIL )

## Disclaimer

All of this code is user contributed! Use at your own risk! No warranty is
expressed or implied. There may be bugs!  Do not use any of this code without
understanding what it is doing.

## Documentation

Documentation is at https://paulray.github.io/NICERsoft/

To build it locally, run `pip install -e ".[docs]"`, then `make html` in the
`docs` directory, and open `docs/_build/html/index.html`.

## Installation

Install the package using pip (recommended):

```bash
cd NICERsoft
pip install -e .
```

This will install all dependencies and make the scripts available in your PATH.

Alternatively, you can use the legacy method by adding `<basedir>/NICERsoft/scripts` to your PATH
and `<basedir>/NICERsoft` to your PYTHONPATH, where `<basedir>` is wherever you cloned NICERsoft.

## NOTES:

* Executable scripts go in the `scripts` subdirectory

* Importable python modules go in the `nicer` subdirectory

* Scripts in any language (Python, PERL, bash, tcsh, etc...) are OK to contribute!

## PREREQUISITES:

* pip-installable packages:
  * numpy, matplotlib, scipy, astropy, pint-pulsar, loguru, tqdm


* Heasoft (with NICER Data Analysis Software)


* Scripts using PyXspec require Heasoft to be compiled from source.
