# [pyRUMP](https://ikost.github.io/pyRUMP/)

[![PyPI](https://img.shields.io/pypi/v/pyrump.svg)](https://pypi.org/project/pyrump/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

![Measured spectrum compared against a pyRUMP simulation, showing a Si / Ru / Mn2.74Pt1 / Ru stack](https://raw.githubusercontent.com/ikost/pyRUMP/main/docs/assets/mnpt-fit.png)
*Measured data (black) vs. a pyRUMP simulation (red) of a Si / Ru / Mn<sub>2.74</sub>Pt<sub>1</sub> / Ru stack.*

A clean Python reimplementation of **RUMP**, the Rutherford backscattering
spectrometry (RBS) simulation and analysis package originally written by
L. R. Doolittle and M. O. Thompson at Cornell.

The original is ~22k lines of unmaintained C from the late 1980s, with no
active support. This reimplementation is written using Claude Code and
reproduces RUMP's physics as a tested, importable Python library, with both
a batch CLI and RUMP's own interactive shell. It runs on Windows, macOS, and
Linux.

## Quick start

```bash
pip install pyrump
```

Requires Python 3.9+; numpy, scipy, and matplotlib are installed
automatically, along with the physics data tables pyRUMP needs at runtime.

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
isn't redistributed here — see [Design and validation](https://ikost.github.io/pyRUMP/validation/).
They skip cleanly when it's absent, so it's not needed for everyday development.

## Licensing and provenance

pyRUMP is MIT licensed, and an **independent reimplementation**: it is not
affiliated with, endorsed by, or derived from the RUMP source distribution.
See [Licensing and provenance](https://ikost.github.io/pyRUMP/about/#licensing-and-provenance)
for the full story, including bundled data-table provenance.
