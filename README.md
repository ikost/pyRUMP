# [pyRUMP](https://ikost.github.io/pyRUMP/)

[![PyPI](https://img.shields.io/pypi/v/pyrump.svg)](https://pypi.org/project/pyrump/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/ikost/pyRUMP/blob/main/LICENSE)

![Measured RBS spectrum compared against a pyRUMP fit of a Si / SiO2 / Ru / Mn2.76Pt1 / Ru stack, with the Pt, Ru, Mn, Si and O surface edges marked](https://raw.githubusercontent.com/ikost/pyRUMP/main/docs/assets/mnpt-fit.png)
*Measured RBS data (grey) vs. a pyRUMP fit (red) of a Si / SiO<sub>2</sub> / Ru / Mn<sub>2.76</sub>Pt<sub>1</sub> / Ru stack, with the elements' surface edges marked.*

![The PIXE spectrum measured with the same beam compared against pyRUMP's PIXE simulation of the same stack, with the Si, Ru, Mn and Pt X-ray lines marked](https://raw.githubusercontent.com/ikost/pyRUMP/main/docs/assets/mnpt-pixe.png)
*The PIXE spectrum taken with the same beam (grey) vs. pyRUMP's PIXE simulation (red) of the sample fitted to the RBS above, with the X-ray lines marked.*

A clean Python reimplementation of **RUMP**, the Rutherford backscattering
spectrometry (RBS) simulation and analysis package originally written by
L. R. Doolittle and M. O. Thompson at Cornell.

The original is ~22k lines of unmaintained C from the late 1980s, with no
active support. This reimplementation is written using Claude Code and
reproduces RUMP's physics as a tested, importable Python library, with both
a batch CLI and RUMP's own interactive shell. It runs on Windows, macOS, and
Linux.

Rewriting from scratch with AI assistance is also a chance to move past a
straight port: a handful of bugs in the original C are fixed rather than
just reproduced (validated against the C oracle — see
[RUMP quirks and defects](https://ikost.github.io/pyRUMP/dev/rump-quirks/)),
new physics has been added (Andersen screening alongside RUMP's own
L'Ecuyer correction, via [`SCREENING`](https://ikost.github.io/pyRUMP/manual/config/#screening-new)),
and the workflow is gaining commands RUMP never had, like automatic
per-fit snapshots with [`AUTOSNAP`](https://ikost.github.io/pyRUMP/manual/pert/#autosnap-new),
and [`SNAPSHOT`](https://ikost.github.io/pyRUMP/manual/shell/#snapshot-snap-new),
which saves a session — fitted or tuned by hand — as a macro that brings it
back. SIMNRA's `.xnra` files can be read and written, spectrum, instrument
settings and sample together, with
[`GETNRA`/`WRITENRA`](https://ikost.github.io/pyRUMP/manual/buffers/#getnra-writenra-new);
NDF support is planned.

**New in pyRUMP 2.0: PIXE.** The beam that gives the RBS spectrum also
makes X-rays, and pyRUMP now simulates them: the K, L and M lines of every
element in the SIM sample, with ECPSSR ionisation cross sections, absorption
in the sample and the detector's filters, and the detector's resolution and
escape peaks, in a [PIXE](https://ikost.github.io/pyRUMP/manual/pixe/) window beside the RBS one. PERT
fits both spectra from the same sample, each element to the spectrum that
sees it best, so alloys that overlap in RBS, such as Fe–Ni or W–Ta, take
their composition from their X-ray lines and their thickness from RBS. The
[FeNi worked example](https://ikost.github.io/pyRUMP/getting-started/feni-rbs-pixe/) goes through one
measurement from start to finish.

## Quick start

```bash
pip install pyrump
```

Requires Python 3.10+; 3.14 is recommended, and is what pyRUMP is tested
on. numpy, scipy, and matplotlib are installed automatically, along with the
physics data tables pyRUMP needs at runtime.

```
pyrump                         # the interactive shell, from any directory
Your wish? cd examples
Your wish? xeq MnPt.RBS         /* load the measured spectrum, above    */
Your wish? sim get MnPt.lcm     /* load its 5-layer sample description  */
Your wish? compare              /* data vs. simulation, with residuals  */
```

Buffer 0 is always the simulation and recomputes itself when the sample or
the active buffer's parameters change — there is no "simulate" command,
exactly as in the original.

Or drive it as one-off batch commands:

```bash
pyrump simulate sample.lcm --energy 2.0 --beam 4He -o out.rbs
pyrump fit sample.lcm measured.rbs --vary thickness:0 --window 190 226
pyrump plot measured.rbs --compare out.rbs -o comparison.png
pyrump convert measured.rbs measured.dat
```

## Documentation

The full manual — interactive shell reference, CLI and Python API, worked
examples, the physics writeup, and the list of RUMP quirks and defects found
while porting — is at **https://ikost.github.io/pyRUMP/**.

## Contributing

```bash
pip install -e ".[dev]"
pytest              # unit tests, no external dependencies
ruff check .
```

Oracle-comparison tests (`pytest -m oracle`) need the RUMP C source, which
isn't redistributed here — see [Design and validation](https://ikost.github.io/pyRUMP/dev/validation/).
They skip cleanly when it's absent, so it's not needed for everyday development.

## Licensing and provenance

pyRUMP is MIT licensed, and an **independent reimplementation**: it is not
affiliated with, endorsed by, or derived from the RUMP source distribution.
See [Licensing and provenance](https://ikost.github.io/pyRUMP/dev/about/#licensing-and-provenance)
for the full story, including bundled data-table provenance.
