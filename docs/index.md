# NICERsoft

NICERsoft is a collection of user-contributed analysis tools for data from the
[NICER](https://heasarc.gsfc.nasa.gov/docs/nicer/) X-ray telescope on the
International Space Station. It is separate from the NICER team's software that
is distributed as part of [HEASoft](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/), and nothing here should be taken as formally
recommended by the NICER instrument team.

The tools are mostly aimed at pulsar timing and spectral analysis:

- **Data processing pipelines** that screen NICER observations, select good
  times and events, and make quick-look plots ({doc}`scripts/pipelines`).
- **Pulsar timing**: pulse phases, phaseograms, templates and photon TOAs
  ({doc}`scripts/timing`).
- **Spectral analysis** for the rotation-powered pulsar (RPP) catalog, using
  PyXspec ({doc}`scripts/spectral`).
- A Python package, `nicer`, with the functions the scripts are built on
  ({doc}`api/index`).

```{warning}
All of this code is user contributed. Use it at your own risk; there may be
bugs. Don't use any of it without understanding what it is doing.
```

Contributions are welcome. See {doc}`contributing`.

```{toctree}
:maxdepth: 2
:caption: Getting started

installation
quickstart
```

```{toctree}
:maxdepth: 2
:caption: User guide

guide/data
guide/conventions
guide/filtering
guide/timing
guide/rpp
```

```{toctree}
:maxdepth: 2
:caption: Reference

scripts/index
api/index
```

```{toctree}
:maxdepth: 1
:caption: Development

contributing
```

Other useful NICER tools by George Younes:
[NICERUTIL](https://github.com/georgeyounes/NICERUTIL).
