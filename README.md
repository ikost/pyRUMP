# pyRUMP

A clean Python reimplementation of **RUMP**, the Rutherford backscattering
spectrometry (RBS) simulation and analysis package originally written by
L. R. Doolittle and M. O. Thompson at Cornell.

The original is ~22k lines of unmaintained C from the late 1980s, with a 1996-era
HTML manual and no active support. pyRUMP reproduces its physics as a tested,
importable library.

> **Status: complete.** All fourteen milestones are done, including the
> interactive shell. pyRUMP simulates RBS spectra matching the original to
> ~3e-6 in total counts, fits them with the same Poisson objective, reads and
> writes RUMP's native `.RBS` and `.lcm` files byte-identically, and offers
> both a batch CLI and RUMP's interactive command environment. 668 tests.

## Try it

```bash
pip install -e ".[dev,plot]"
pyrump                         # the interactive shell, from any directory
```

No further setup needed — the physics data tables ship with the package.

## Documentation

- **[docs/usage.md](docs/usage.md)** — command reference, Python API, and worked
  examples: identifying elements from a spectrum, simulating a known structure,
  fitting a thickness.
- **[docs/algorithm.md](docs/algorithm.md)** — what the simulation computes and
  how, what it approximates, and where the shipped code diverges from the
  published algorithm.
- **[docs/rump-quirks.md](docs/rump-quirks.md)** — 20 documented defects and
  surprises in the original C, and what pyRUMP does about each.

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

Python 3.9+, numpy, scipy.

## Use

```python
from pathlib import Path
import pyrump
from pyrump.atomic.density import DensityTable
from pyrump.atomic.tables import PeriodicTable
from pyrump.io.kalbitzer import parse_kalbitzer
from pyrump.model.geometry import Geometry
from pyrump.model.spectrum import Calibration
from pyrump.sim.engine import Beam, UniformSample, simulate
from pyrump.stopping.kalbitzer import KalbitzerStopping
from pyrump.stopping.registry import StoppingRegistry
from pyrump.stopping.ziegler import ZieglerStopping

DATA = Path(pyrump.__file__).parent / "data"   # bundled with the package
table = PeriodicTable.load(DATA / "atom4.dat", DATA / "pscoef.dat")
registry = StoppingRegistry(
    table.elements,
    kalbitzer=KalbitzerStopping(parse_kalbitzer(DATA / "newstop.kal"), table.elements),
    ziegler=ZieglerStopping(table.elements),
)
densities = DensityTable.load(DATA / "density.tab")

spectrum = simulate(
    UniformSample([1000.0], [14], [[1.0]]),   # 1000e15 at/cm^2 of silicon
    Beam(e0_MeV=2.0, z=2, mass=4.0026),       # 2 MeV alphas
    Geometry(theta=0.0, phi=10.0),            # phi is 180 - scattering angle
    registry, table, Calibration(kevch=5.0, npt=1024),
)
```

or as one-off commands:

```bash
pyrump simulate sample.lcm --energy 2.0 --beam 4He -o out.rbs
pyrump fit sample.lcm measured.rbs --vary thickness:0 --window 190 226
pyrump plot measured.rbs --compare out.rbs -o comparison.png
pyrump convert measured.rbs measured.dat
```

Where an element's surface edge lands, at 2 MeV He and 170 degrees:

```python
from pyrump.physics.kinematics import kinematic_factor

gold = table.by_symbol("Au")
mass = max(gold.isotopes, key=lambda i: i.fraction).mass
print(kinematic_factor(4.0026, mass, 170.0) * 2000, "keV")   # surface-edge energy
```

Recovering a thickness by fitting simulated data against itself:

```python
from pyrump.fit.lm import fit
from pyrump.fit.parameters import FitInputs, thickness
from pyrump.fit.windows import Window, WindowSet
from pyrump.model.detector import Measurement

measurement = Measurement(omega_msr=1.0, charge_uC=10.0, fwhm_keV=15.0)
calibration = Calibration(kevch=5.0, npt=1024)
geometry = Geometry(theta=0.0, phi=10.0)

def run(inputs):
    return simulate(
        inputs.sample, inputs.beam, inputs.geometry, registry, table,
        inputs.calibration, inputs.measurement,
    ).counts

truth = FitInputs(UniformSample([1200.0], [14], [[1.0]]), Beam(e0_MeV=2.0, z=2, mass=4.0026),
                   geometry, calibration, measurement)
guess = FitInputs(UniformSample([900.0], [14], [[1.0]]), Beam(e0_MeV=2.0, z=2, mass=4.0026),
                   geometry, calibration, measurement)

result = fit(run, run(truth), guess, [thickness(0)], windows=WindowSet(error=[Window(200, 260)]))
print(result.parameters["thickness[0]"], "recovered, truth was 1200")
```

For more — identifying elements from a real spectrum, reading/writing
`.RBS` and `.lcm` files, depth profiles — see
[docs/usage.md](docs/usage.md), which this is distilled from.

### The interactive shell

Running `pyrump` bare starts RUMP's own command environment, with the original
command names and minimum abbreviations, from any directory:

```
Your wish? cd data              /* cd/ls/pwd, as the original had    */
Your wish? get 2A.rbs           /* read a spectrum and its metadata  */
Your wish? plot 1               /* persistent matplotlib window      */
Your wish? region 100 400       /* adjust it in place                */
Your wish? sim                  /* edit the sample description       */
SIM Command: get ITO.lcm
SIM Command: return
Your wish? compare              /* data vs simulation, with residuals */
Your wish? pert                 /* fit it                            */
PERT Command: window 355 375
PERT Command: thickness 1
PERT Command: go
Your wish? display              /* depth profile                     */
```

Spectra live in numbered buffers; buffer 0 is the simulation and recomputes
itself when the sample or the active buffer's parameters change — there is no
"simulate" command, exactly as in the original. `XEQ` replays a file of commands
and `SCRIPT` writes one, so an analysis can be checked in and reproduced.

See [docs/usage.md](docs/usage.md#interactive-shell) for the full command set.

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
| M14 | Interactive shell: buffers, SIM and PERT levels, macros | done |

**Known limitation**: non-Rutherford (tabulated-resonance) cross sections —
RUMP's `.adt`/R33 nuclear cross-section tables — have a complete, tested
reader (`pyrump.io.adt`) but aren't wired into simulation yet; pyRUMP
currently computes pure Rutherford + L'Ecuyer-screened scattering only. A
later milestone.

## Contributing

```bash
pip install -e ".[dev,plot]"
pytest              # unit tests, no external dependencies
ruff check .
```

Oracle-comparison tests (`pytest -m oracle`, and the wider set of tests that
compare against the legacy C for extra confidence) need the RUMP 2.0 C source,
which isn't redistributed here — see [Validation](#validation). They skip
cleanly when it's absent, so it's not needed for everyday development.

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

Four data tables (`atom4.dat`, `pscoef.dat`, `newstop.kal`, `density.tab`) are
bundled with pyRUMP (`src/pyrump/data/`). Their provenance is third-party and
independent of CGS — `pscoef.dat` is the ZBL/TRIM `SCOEF` table, `atom4.dat`
is elements and isotopes — and has been checked against current CIAAW/NIST
values and literature; see `src/pyrump/data/SOURCES.md` for full provenance
and `NOTICE` for citations. Non-Rutherford cross-section tables (`*.adt`) are
IBANDL evaluations and are **not** bundled — obtain them separately from
IBANDL if you need that data.

## References

- L. R. Doolittle, *Algorithms for the rapid simulation of Rutherford
  backscattering spectra*, Nucl. Instr. Meth. **B9** (1985) 344–351.
- L. R. Doolittle, *A new approach to Rutherford backscattering analysis*,
  Nucl. Instr. Meth. **B15** (1986) 227–231.
- J. F. Ziegler, J. P. Biersack, U. Littmark, *The Stopping and Range of Ions in
  Solids*, Pergamon (1985).
- G. Konac, S. Kalbitzer et al., Nucl. Instr. Meth. **B136–138** (1998) 159–165.
