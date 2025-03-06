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

See the Wiki for some documentation: https://github.com/paulray/NICERsoft/wiki

## Usage

To use these scripts add `<basedir>/NICERsoft/scripts` to your PATH
and add `<basedir>/NICERsoft` to your PYTHONPATH, where ``<basedir>`` is wherever you
cloned NICERsoft.

## NOTES:

* Executable scripts go in the `scripts` subdirectory

* Importable python modules go in the `nicer` subdirectory

* Scripts in any language (Python, PERL, bash, tcsh, etc...) are OK to contribute!

## PREREQUISITES:

* pip-installable packages:
  * numpy, matplotlib, scipy, astropy, pint-pulsar, loguru, tqdm


* Heasoft (with NICER Data Analysis Software)


* Scripts using PyXspec require Heasoft to be compiled from source.
