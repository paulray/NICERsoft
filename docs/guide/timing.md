# Pulsar timing

NICERsoft uses [PINT](https://nanograv-pint.readthedocs.io/) for pulse phases
and TOAs. A typical workflow is shown in {doc}`../quickstart`; this page gives
more detail on each step.

## Pulse phases

Phases are computed by PINT's `photonphase`, using a timing model (`.par` file)
and the NICER orbit file. {ref}`psrpipe.py <psrpipe>` does this automatically
when given `--par`, adding a `PULSE_PHASE` column to `cleanfilt.evt`. `--polycos`
uses polynomial approximations to the timing model, which is much faster for
large numbers of events. For a file that has already been processed, run
`photonphase` directly:

```bash
photonphase --orbfile ni1060020101.orb.gz --addphase cleanfilt.evt J0030+0451.par
```

## Looking at the pulsations

- {ref}`niphaseogram.py <niphaseogram>` plots phase against time, which shows
  whether the timing model holds over the whole data set.
- {ref}`plotphist.py <plotphist>` plots a binned pulse profile.
- {ref}`ni2dpulse.py <ni2dpulse>` plots phase against energy.
- {ref}`fitharms.py <fitharms>` fits the profile with harmonics and
  {ref}`calcrms.py <calcrms>` computes the fractional RMS amplitude.
- {ref}`nioptcuts.py <nioptcuts>` and
  {ref}`ni_Htest_sortgti.py <ni-htest-sortgti>` choose energy and GTI cuts that
  maximize the detection significance.

## Templates

TOAs are measured by fitting a template to the pulse profile.
{ref}`nitemplate.py <nitemplate>` fits a template of Gaussian components
interactively and writes it to a `.gauss` file. Templates in the format with a
`META_STOP` header line can be renormalized with
{ref}`fix_template.py <fix-template>`.

## TOAs

{ref}`photon_toa.py <photon-toa>` measures TOAs from the events, a template
and a timing model, and writes them in Tempo2 format. The main options are:

- `--tint` and `--maxint`: how much exposure goes into each TOA. With neither,
  one TOA is made from all of the events.
- `--grid`: make TOAs on a regular grid of times instead.
- `--phasecol PULSE_PHASE`: use phases already in the event file instead of
  recomputing them.
- `harmN` as the template name: measure the TOA from the phase of the Nth
  Fourier harmonic instead of fitting a template.
- `--topo`: make topocentric (spacecraft) TOAs, including the spacecraft
  position on each TOA line; otherwise the TOAs are barycentric.

Each TOA has flags giving the exposure, H-test value and signal-to-noise
ratio, among others. {ref}`toafilt.py <toafilt>` removes TOAs with a low
H-test value, short exposure or large uncertainty.

{ref}`photon_fit.py <photon-fit>` fits a template to a stored phase column and
reports the source and background properties without making TOAs, and
{ref}`pardiff.py <pardiff>` compares two timing models.

## Other formats

{ref}`nicerfits2presto.py <nicerfits2presto>` converts events to
[PRESTO](https://github.com/scottransom/presto) format for periodicity
searches.
