# pyRUMP

A clean Python reimplementation of **RUMP**, the Rutherford backscattering
spectrometry (RBS) simulation and analysis package originally written by
L. R. Doolittle and M. O. Thompson at Cornell.

The original is ~22k lines of unmaintained C from the late 1980s, with a 1996-era
HTML manual and no active support. pyRUMP reproduces its physics as a tested,
importable library.

> **Status: complete.** All thirteen milestones are done. pyRUMP simulates RBS
> spectra matching the original to ~3e-6 in total counts, fits them with the
> same Poisson objective, reads and writes RUMP's native `.RBS` and `.lcm` files
> byte-identically, and ships a CLI. 503 tests.

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

## Use

```python
from pyrump.atomic.tables import PeriodicTable
from pyrump.sim.engine import Beam, UniformSample, simulate
from pyrump.model.geometry import Geometry
from pyrump.model.spectrum import Calibration

spectrum = simulate(
    UniformSample([1000.0], [14], [[1.0]]),   # 1000e15 at/cm^2 of silicon
    Beam(e0_MeV=2.0, z=2, mass=4.0026),       # 2 MeV alphas
    Geometry(theta=0.0, phi=10.0),            # phi is 180 - scattering angle
    registry, periodic_table, Calibration(kevch=5.0, npt=1024),
)
```

or from the shell:

```bash
pyrump simulate sample.lcm --energy 2.0 --beam 4He -o out.rbs
pyrump fit sample.lcm measured.rbs --vary thickness:0 --window 190 226
pyrump plot measured.rbs --compare out.rbs -o comparison.png
pyrump convert measured.rbs measured.dat
```

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
| Absorber, fuzz, multiple scattering | **3e-6** total, **4e-5** of peak | float32 brick edges |
| `.RBS` files read vs RUMP's own reader | **bit-identical** | — |
| Poisson objective vs `EvalChiPoisson` | **1e-5** reduced chi2 | float32 in the C |
| `.lcm` round-trip vs RUMP's own writer | **byte-identical** | — |

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
| M10 | Absorber, pile-up, fuzz, multiple scattering | done |
| M11 | File I/O (`.RBS` binary, ASCII) | done |
| M12 | Fitting (PERT) | done |
| M13 | CLI, plotting, `.lcm` subset | done |

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
