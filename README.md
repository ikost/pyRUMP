# pyRUMP

[![PyPI](https://img.shields.io/pypi/v/pyrump.svg)](https://pypi.org/project/pyrump/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A clean Python reimplementation of **RUMP**, the Rutherford backscattering
spectrometry (RBS) simulation and analysis package originally written by
L. R. Doolittle and M. O. Thompson at Cornell. 

This implementations is written using Claude Code.

The original is ~22k lines of unmaintained C from the late 1980s, with a 1996-era
HTML manual and no active support. pyRUMP reproduces its physics as a tested,
importable library, with both a batch CLI and RUMP's own interactive shell.

## Install

```bash
pip install -e .
```

Python 3.9+, numpy, scipy. The four physics data tables pyRUMP needs at
runtime ship with the package, so nothing further is needed for simulation,
fitting, or the interactive shell.

## Quick start

```bash
pyrump                         # the interactive shell, from any directory
```

```
Your wish? get 2A.rbs           /* read a spectrum and its metadata  */
Your wish? sim                  /* edit the sample description       */
SIM Command: get ITO.lcm
SIM Command: return
Your wish? compare              /* data vs simulation, with residuals */
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

## Changelog

Release notes live on the [GitHub Releases page](https://github.com/ikost/pyRUMP/releases).

## Contributing

```bash
pip install -e ".[dev,plot]"
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
