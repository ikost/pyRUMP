# pyRUMP

A clean Python reimplementation of **RUMP**, the Rutherford backscattering
spectrometry (RBS) simulation and analysis package originally written by
L. R. Doolittle and M. O. Thompson at Cornell.

The original is ~22k lines of unmaintained C from the late 1980s, with a 1996-era
HTML manual and no active support. pyRUMP reproduces its physics as a tested,
importable library.

> **Status: the forward model is complete.** pyRUMP simulates RBS spectra —
> including straggling and detector resolution — that match the original to
> ~3e-6 in total counts and ~1e-5 of peak per channel, with depth-dependent
> composition profiles. Still to come: absorber/pileup, file I/O and fitting.
> See [Milestones](#milestones).

## Design

**Faithful first, corrected by choice.** The default reproduces the shipped C
bug-for-bug, because that is what every published RUMP result was produced with.
Known defects — and there are several — are reproduced exactly, with the
mathematically correct behaviour available behind explicit flags. See
[docs/rump-quirks.md](docs/rump-quirks.md).

**Validated against the original.** The legacy C is compiled into a shared
library and called directly from the test suite, so each stage is compared
function-by-function rather than by eyeballing a final spectrum.

## Installation

```bash
pip install -e ".[dev,plot]"
```

Python 3.11+, numpy, scipy.

## Validation

pyRUMP is tested against the original C at two levels.

**Unit oracle** — the physics translation units (`ziegler.c`, `stopping.c`,
`sigma.c`, …) are compiled into `libpyrump_oracle` and called via cffi. No TTY,
no graphics, no buffers. When a number disagrees, this isolates the cause to one
function.

```bash
python tests/oracle/build_oracle.py
pytest -m oracle
```

**End-to-end oracle** — the original `rump` binary is driven through a
pseudo-terminal to produce reference spectra.

Both require the legacy C tree, which is **not redistributed** (see
[Licensing](#licensing)). Point `PYRUMP_C_REFERENCE` at it, or place it at
`C-code/`. Tests skip cleanly when it is absent.

```bash
pytest              # unit tests
pytest -m oracle    # comparison against the C
```

### Current agreement

| Quantity | Agreement | Limited by |
|---|---|---|
| ZBL85 stopping, all 92 targets, H/He/Cu/Au beams, 10 keV–10 MeV | **6.1e-7** rel | float32 tables in the C |
| Fitted stopping polynomial (what the simulation consumes) | **1.3e-5** rel | float32 coefficient storage |
| Polynomial evaluation, given identical coefficients | **1.1e-14** rel | float64 round-off |
| Cross-sections and kinematics | **1e-10** rel | closed forms, nothing to fit |
| Bricks — 630 across 36 configurations | **5e-7** energies, **6e-6** heights | float32 coefficients |
| Full spectrum, total counts | **3e-6** rel | float32 coefficients |
| Full spectrum, per channel | **1e-5** of peak | float32 brick edges |
| With straggling and detector resolution | **3e-6** total, **1e-5** of peak | float32 brick edges |
| Depth profiles, all 11 evaluable forms | **2.6e-5** brick heights | float32 coefficients |

The oracle is the `float` build. RUMP cannot be built in double precision — its
table readers use `scanf("%f")` against `REAL` fields, so `-DREAL_IS_DOUBLE`
silently corrupts every table. That caps how tightly any float64 port can agree,
and the tolerances above are set by that floor rather than by choice.

## Milestones

| | Milestone | Status |
|---|---|---|
| M0 | Reference oracle (pty driver + cffi library) | done |
| M1 | Elements, isotopes, compound densities | done |
| M2 | Stopping powers: ZBL85, Konac, Mylar, priority chain | done |
| M3 | STOP_SQRT polynomial refit, Bragg summation, session cache | done |
| M4 | Kinematics, geometry, cross-sections | done |
| M5 | Slab march → bricks | done |
| M6 | Brick → channel fill | done |
| M7 | Straggling (closed-form erf) | done |
| M8 | Detector convolution | done |
| M9 | Depth profiles (13 EQUATION forms) | done |
| M10 | Absorber, foil, pileup, fuzz | next |
| M11 | File I/O (`.RBS`, ASCII, `.adt`/R33) | |
| M12 | Fitting (PERT) | |
| M13 | CLI, plotting, `.lcm` subset | |

## Licensing

pyRUMP is MIT licensed.

It is an **independent reimplementation**. The RUMP 2.0 C source is *not*
included and *not* redistributed: its licence permits modification but forbids
commercial exploitation, which is incompatible with an open-source distribution.
It is used here only as a local validation oracle.

RUMP and Genplot were trademarks of Computer Graphic Service, Ltd. CGS ceased
operating as a business in June 2012 and `genplot.com` no longer resolves; the
authors stated at the time that GENPLOT and RUMP remain free to download and
use. That removes the trademark concern, but **not** copyright in the original
source, which remains with its authors — hence the C tree is still not
redistributed here.

pyRUMP is not affiliated with, endorsed by, or derived from the RUMP 2.0 source
distribution.

Data tables (`atom4.dat`, `pscoef.dat`, `newstop.kal`, `density.tab`, `*.adt`)
are read from the legacy tree during development only, and are not
redistributed. Their provenance is third-party and independent of CGS —
`pscoef.dat` is the ZBL/TRIM `SCOEF` table, the `*.adt` files are IBANDL
evaluations — so they will be regenerated from primary sources (CIAAW/NIST,
published ZBL tables, IBANDL) before any public release. See `NOTICE`.

## References

- L. R. Doolittle, *Algorithms for the rapid simulation of Rutherford
  backscattering spectra*, Nucl. Instr. Meth. **B9** (1985) 344–351.
- L. R. Doolittle, *A new approach to Rutherford backscattering analysis*,
  Nucl. Instr. Meth. **B15** (1986) 227–231.
- J. F. Ziegler, J. P. Biersack, U. Littmark, *The Stopping and Range of Ions in
  Solids*, Pergamon (1985).
- G. Konac, S. Kalbitzer et al., Nucl. Instr. Meth. **B136–138** (1998) 159–165.
