# Using pyRUMP

Command reference, Python API, and worked examples. For what the simulation
actually computes, see [algorithm.md](algorithm.md).

Every number in the examples below was produced by running the code, not
written from memory.

---

## Contents

- [Setup](#setup)
- [Interactive shell](#interactive-shell)
- [Command line](#command-line)
- [Python API](#python-api)
- [Worked example: identifying what is in a sample](#worked-example-identifying-what-is-in-a-sample)
- [Worked example: simulating a known structure](#worked-example-simulating-a-known-structure)
- [Worked example: fitting a thickness](#worked-example-fitting-a-thickness)
- [Sample descriptions (`.lcm`)](#sample-descriptions-lcm)
- [Depth profiles](#depth-profiles)
- [Things that will catch you out](#things-that-will-catch-you-out)

---

## Setup

```bash
pip install -e ".[dev,plot]"
```

pyRUMP bundles the four data tables it needs at runtime — `atom4.dat`,
`pscoef.dat`, `newstop.kal`, `density.tab` — so nothing further is needed for
simulation, fitting, or the interactive shell.

The worked examples below that compare against RUMP's own shipped sample
spectrum and `.lcm` file (`Fixed/2A.rbs`, `Fixed/ITO.lcm`) do need the legacy
RUMP distribution, since those files aren't redistributed with pyRUMP (see
[Licensing](../README.md#licensing)). Point `--data DIR` or `PYRUMP_DATA` at
it if you have a copy:

```bash
export PYRUMP_DATA=/path/to/rump/data
```

---

## Interactive shell

Running `pyrump` with no arguments starts the interactive shell — RUMP's own
working style, from any directory:

```
$ pyrump
pyRUMP interactive shell -- tables from /path/to/rump/data
Type ? for commands, QUIT to leave.
Your wish? get 2A.rbs
Your wish? plot 1
Your wish? region 100 400
Your wish? sqrt
Your wish? sim
SIM Command: get ITO.lcm
SIM Command: show
SIM Command: return
Your wish? compare
Your wish? display
Your wish? quit
```

Command names and their **minimum abbreviations** follow the original
(`REGion`, `OVerlay`, `COMPare`), so `reg 100 400` and `region 100 400` are the
same command. `?` lists everything, with the required characters upper-cased.

Matching is *first-hit in table order*, exactly as `LexCmdl` (lexp2.c:639)
does — it never reports an ambiguous abbreviation, it takes the first entry that
fits. `pyrump shell` is the explicit form, and takes `--norc`, `--batch`, and a
macro file to run at startup.

### Getting around

The shell has RUMP's own filesystem commands (a port of the "General System
Commands" table at `lexp/system.c:175`), so you can move to your data rather than
restarting in the right directory:

| Command | Effect |
| --- | --- |
| `PWD` / `WHERE` | print the working directory |
| `CD <dir>` / `CHDIR` | change directory; **no argument goes home** |
| `PUSHDIR <dir>` / `POPDIR` | change directory remembering the old one, and come back |
| `LS [pattern]` / `DIRECTORY` | list files; `ls *.rbs` filters |
| `LL [pattern]` | long listing, with size and date |
| `TYPE <file>` / `CAT` / `MORE` | show a text file, paged when interactive |
| `CLS` | clear the screen |

Wildcards are expanded by the command itself, never by an OS shell, so `ls *.rbs`
behaves the same on Linux, macOS and Windows. Tab completion works on both
command names and paths.

These are reachable from `SIM` and `PERT` too — as in the original, a command the
sub-level does not know returns you to the RUMP level and runs there.

There is deliberately **no shell escape** (the original's `!` / `DOS` / `CSH`):
it would let any `.cmd` macro run arbitrary commands on your machine.

### Buffers

Spectra live in numbered buffers, one of which is ACTIVE and is what most
commands act on implicitly. Buffer **0 is the simulation**; data starts at 1.

| Command | Effect |
| --- | --- |
| `GET <file\|n>` | read a file into a buffer, or point at buffer *n* |
| `BUFFERS` | list the buffers, marking the active one |
| `ACTIVE` | print the active buffer's full parameter set |
| `COPY a b` / `MOVE a b` | copy / exchange |
| `RELEASE` / `NEWALL` | drop one / all |
| `WRITE f.rbs` / `WRASCII f.dat` | save the active buffer |

Unlike the original there is no ten-buffer limit and nothing is destroyed to
make room. Buffer 0 has **no simulate command**: it is recomputed whenever the
sample or the active buffer's parameters change, which is how RUMP behaved.
`RECALCULATE` forces it.

### Plotting

The plot is one persistent matplotlib window whose state survives between
commands. `PLOT` erases and draws, `OVERLAY` adds to it.

`REGION lo hi`, `COUNTS lo hi`, `EXPAND`, `BLOWUP`, `LINEAR`, `SQRT`, `LOG`,
`NORMALIZE`, `RAW`, `LABELS`, and `ENERGY` change how it is drawn and redraw
immediately. `PARMS` prints the current settings. `COMPARE` plots the active
buffer against the simulation with Poisson residuals; `DISPLAY` plots the
sample's composition against depth.

matplotlib is an optional dependency — install it with `pip install 'pyrump[plot]'`
or the shell will say so when you first plot.

### SIM and PERT

`SIM` edits the sample description; `PERT` fits it. Both are sub-levels with
their own prompt, and — as in the original — a command the sub-level does not
recognise is passed out to the RUMP level, which returns you there automatically.
`SIM <command>` also works as a one-shot from the top level.

```
Your wish? sim
SIM Command: thick 500 A
SIM Command: composition Si 1 /
SIM Command: next
SIM Command: thick 2000 A
SIM Command: composition Au 1 /
SIM Command: show
SIM Command: save mysample.lcm
SIM Command: return
```

SIM's sample-definition commands are the *same code* that parses `.lcm` files
(`SampleEditor` in `pyrump/script/lcm.py`), so what you type and what the file
holds cannot drift apart. `SIM SAVE` writes RUMP's own format.

PERT selects what may vary and over which channels, then `GO`:

```
Your wish? pert
PERT Command: window 355 375      /* compare only here            */
PERT Command: norm 140 200        /* rescale data to remove dose error */
PERT Command: thickness 1         /* vary layer 1's thickness     */
PERT Command: go
```

```
  reduced chi-square 1.2849 on 20 dof
  12 evaluations, Both `ftol` and `xtol` termination conditions are satisfied.
  data scaled by 0.99441 over the norm window
  thickness[0]                      299.198  +/- 0.3801   (was 200)
```

Fitted values are written back into the sample description, so `SIM SHOW` and
`SIM SAVE` reflect them. Two differences from the original: the data may be in
any buffer, not just buffer 1, and `MULTI` is the default because the solver is
a simultaneous least-squares fit (`SINGLE` loops one parameter at a time).

### Macros

`XEQ <file>` runs a file of commands through the same interpreter the prompt
uses, so an analysis can be checked in as a text file and replayed. `CALL` and
`EXECUTE` are synonyms.

`SCRIPT <file>` logs what you type into exactly such a file, and `SCRIPT OFF`
stops. `LOGFILE` and `RECORD` are synonyms. `~/.pyrumprc` is run at startup
unless you pass `--norc`.

> Note `SCRIPT`/`LOGFILE` need at least four characters, which is how the
> original kept them clear of `LOG` — the logarithmic yield axis. Typing `log`
> gets you the axis, `logf` the session log.

---

## Command line

```
pyrump [--data DIR] {shell,simulate,fit,convert,plot} ...
```

### `pyrump simulate`

Simulate a spectrum from a sample description.

```bash
pyrump simulate sample.lcm --energy 2.0 --beam 4He -o out.rbs
```

| Option | Default | Meaning |
|---|---|---|
| `--beam` | `4He` | Beam species: `4He`, `1H`, `He`, `28Si`… |
| `--energy` | `2.0` | Beam energy, MeV |
| `--theta` | `0.0` | Incidence angle from the sample normal, degrees |
| `--phi` | `10.0` | **180° minus the scattering angle** (10 → 170°) |
| `--psi` | `0.0` | Exit angle from the normal (GENERAL geometry only) |
| `--geometry` | `cornell` | `cornell`, `ibm` or `general` |
| `--kevch` | `5.0` | keV per channel |
| `--kev0` | `0.0` | keV at channel zero |
| `--channels` | `1024` | Number of channels |
| `--fwhm` | `15.0` | Detector resolution, keV |
| `--omega` | `1.0` | Detector solid angle, msr |
| `--charge` | `10.0` | Integrated charge, µC |
| `--current` | `0.0` | Beam current, nA — enables pile-up with `--tau` |
| `--tau` | `0.0` | Shaping time, µs |
| `-o, --output` | stdout | `.rbs` for binary, anything else for ASCII |
| `--two-column` | off | Write `channel value` instead of one column |

Output extension chooses the format. With no `-o`, energy/counts pairs go to
stdout and a summary to stderr.

### `pyrump fit`

Adjust sample parameters until the simulation matches a measurement.

```bash
pyrump fit start.lcm measured.rbs --vary thickness:0 --window 190 226
```

Takes every `simulate` option, plus:

| Option | Meaning |
|---|---|
| `--vary NAME` | Parameter to fit; repeat for several |
| `--window LOW HIGH` | Channel range to fit over; repeat for up to 10 |

`--vary` accepts `thickness:N` for layer *N*, or any of `mev`, `theta`, `phi`,
`psi`, `fwhm`, `tau`, `current`, `correction`, `kev/ch`, `kev(0)`, `straggle`,
`multiple_scatter`.

**Choose the window deliberately.** It should cover the part of the spectrum
that constrains what you are fitting, and no more — see
[Things that will catch you out](#things-that-will-catch-you-out).

### `pyrump convert`

```bash
pyrump convert measured.rbs measured.dat        # binary → ASCII
pyrump convert spectrum.dat spectrum.txt --two-column
```

Writing `.rbs` requires beam and geometry metadata, so it only works from a
`.rbs` source — ASCII files do not carry it.

### `pyrump plot`

```bash
pyrump plot measured.rbs                            # interactive window
pyrump plot measured.rbs -o spectrum.png            # save
pyrump plot measured.rbs --compare simulated.rbs -o comparison.png
```

With `--compare` you get the data, the simulation over it, and a residual panel
showing the **Poisson residuals the fit minimises** — not `data − model`.

---

## Python API

The CLI is a thin wrapper; the library is the primary interface.

```python
from pyrump.atomic.density import DensityTable
from pyrump.atomic.tables import PeriodicTable
from pyrump.io.kalbitzer import parse_kalbitzer
from pyrump.model.detector import Measurement
from pyrump.model.geometry import Geometry
from pyrump.model.spectrum import Calibration
from pyrump.sim.engine import Beam, UniformSample, simulate
from pyrump.stopping.kalbitzer import KalbitzerStopping
from pyrump.stopping.registry import StoppingRegistry
from pyrump.stopping.ziegler import ZieglerStopping

from pathlib import Path
import pyrump

DATA = Path(pyrump.__file__).parent / "data"   # bundled with the package

table = PeriodicTable.load(DATA / "atom4.dat", DATA / "pscoef.dat")
registry = StoppingRegistry(
    table.elements,
    kalbitzer=KalbitzerStopping(parse_kalbitzer(f"{DATA}/newstop.kal"), table.elements),
    ziegler=ZieglerStopping(table.elements),
)

spectrum = simulate(
    UniformSample(
        thicknesses=[1000.0],        # 1e15 atoms/cm^2
        element_z=[14],              # silicon
        compositions=[[1.0]],
    ),
    Beam(e0_MeV=2.0, z=2, mass=4.0026),
    Geometry(theta=0.0, phi=10.0),   # phi = 180 - scattering angle
    registry,
    table,
    Calibration(kevch=5.0, kev0=0.0, npt=1024),
    Measurement(omega_msr=1.0, charge_uC=10.0, fwhm_keV=15.0),
)

print(spectrum.total(), "counts")
```

Building the registry takes a moment; **build it once and reuse it**, especially
when fitting.

### Reading and writing files

```python
from pyrump.io.rbs import read_rbs, write_rbs
from pyrump.io.ascii import read_ascii, write_ascii

measured = read_rbs("data.rbs")
measured.counts          # np.ndarray
measured.calibration     # keV/channel, offset
measured.geometry        # angles, geometry convention
measured.e0_MeV, measured.zbeam, measured.mbeam
```

---

## Worked example: identifying what is in a sample

Given an unknown spectrum, the first question is which elements are present.
Each element's **surface edge** sits at *K*·*E*₀, so predicted edge positions
identify the peaks.

Using `2A.rbs`, one of the files shipped with RUMP:

```python
import numpy as np
from pyrump.io.rbs import read_rbs
from pyrump.physics.kinematics import kinematic_factor

s = read_rbs("C-code/rump/data/Fixed/2A.rbs")
print(s.identifier)
print(f"{s.e0_MeV} MeV, Z={s.zbeam}, scattering angle {s.geometry.scattering_angle}")
```

```
Binghampton_target_02A.RBS  RBS LT =  905.98 RT  962.42
3.0 MeV, Z=1, scattering angle 160.0
```

So: 3 MeV protons at 160°, 7.815 keV/channel with a 65.6 keV offset. Now
predict where each candidate element's edge would fall:

```python
E0 = s.e0_MeV * 1000
for symbol in ("C", "O", "Si", "Ti", "Fe", "In", "Sn", "Au"):
    element = table.by_symbol(symbol)
    mass = max(element.isotopes, key=lambda i: i.fraction).mass
    K = kinematic_factor(s.mbeam, mass, s.geometry.scattering_angle)
    print(f"{symbol:3s} K={K:.4f}  E={K*E0:7.1f} keV  channel {s.calibration.channel_of(K*E0):6.1f}")
```

```
C   K=0.7213  E= 2164.0 keV  channel  268.5
O   K=0.7829  E= 2348.6 keV  channel  292.1
Si  K=0.8695  E= 2608.5 keV  channel  325.4
Ti  K=0.9217  E= 2765.0 keV  channel  345.4
Fe  K=0.9325  E= 2797.4 keV  channel  349.6
In  K=0.9665  E= 2899.6 keV  channel  362.6
Sn  K=0.9679  E= 2903.7 keV  channel  363.2
Au  K=0.9803  E= 2941.0 keV  channel  367.9
```

The measured spectrum has falling edges at channels **267, 292 and 363**, which
match **carbon, oxygen, and indium/tin**. Indium and tin are 0.6 channels apart
here and cannot be separated — a general limitation for neighbouring heavy
elements, and the reason a fit constrains their *ratio* rather than resolving
them independently.

That composition — In, Sn, O over a C/O/H substrate — is indium tin oxide on a
polymer, which is exactly what `ITO.lcm` in the same directory describes.

---

## Worked example: simulating a known structure

RUMP ships both the measurement and a matching sample description, so we can
simulate one against the other:

```python
from pyrump.script.lcm import read_lcm, to_sample
from pyrump.sim.engine import Beam, simulate

densities = DensityTable.load(f"{DATA}/density.tab")
observed  = read_rbs(f"{DATA}/Fixed/2A.rbs")
sample    = to_sample(read_lcm(f"{DATA}/Fixed/ITO.lcm"), table, densities)

simulated = simulate(
    sample,
    Beam(e0_MeV=observed.e0_MeV, z=observed.zbeam, mass=observed.mbeam),
    observed.geometry, registry, table,
    observed.calibration, observed.measurement,
)
```

Comparing yields in each edge region:

| Region | Channels | Measured | Simulated |
|---|---|---|---|
| C edge | 255–270 | 266 902 | 27 205 |
| O edge | 280–295 | 63 160 | 8 150 |
| In/Sn edge | 350–366 | 35 261 | 2 273 |
| **total** | | **3 489 801** | **380 110** |

The **structure is right** — the edges land in the right channels and the
relative intensities are close — but the absolute yield is **9.2× low**.

That is not a simulation error. Running the same case through the original C
gives 381 120 counts, agreeing with pyRUMP to 2.6e-3. Both codes say the same
thing, so the discrepancy lives in the measurement's normalisation: the actual
collected charge, solid angle, or detector efficiency differs from the values
recorded in the file.

**This is the normal situation in RBS**, and it is why RUMP has both a `CORR`
factor and a normalisation window. Rather than trusting the charge integration,
you fit the scale:

```python
from pyrump.fit.windows import Window, WindowSet

windows = WindowSet(
    error=[Window(255, 370)],          # fit over the interesting region
    normalisation=Window(255, 370),    # and let the scale float
)
```

The normalisation window forces the total counts over that range to agree by
scaling the data, before χ² is evaluated — so a charge-integration error stops
biasing the fitted thicknesses.

---

## Worked example: fitting a thickness

Simulate a 2400 Å silicon layer, then recover it from a 2000 Å starting guess.

```bash
pyrump simulate truth.lcm  --energy 2.0 -o data.rbs
pyrump fit      start.lcm data.rbs --energy 2.0 --vary thickness:0 --window 190 226
```

```
reduced chi-square 0.0000 on 36 dof
10 evaluations, `xtol` termination condition is satisfied.
  thickness[0]                     1194.6  +/- 2.7494
```

2400 Å of silicon is 2400 × 0.4977 = **1194.5** in 10¹⁵ atoms/cm², so the fit
recovers it to better than 0.1%.

Note that thickness is reported in **areal density**, not Ångström. RBS measures
atoms per unit area; converting to a physical thickness needs an assumed
density, which is a separate and often less certain quantity.

The same fit from Python, with two parameters:

```python
from pyrump.fit.lm import fit
from pyrump.fit.parameters import FitInputs, thickness, parameter
from pyrump.fit.windows import Window, WindowSet

inputs = FitInputs(sample=sample, beam=beam, geometry=geometry,
                   calibration=calibration, measurement=measurement)

result = fit(
    lambda i: simulate(i.sample, i.beam, i.geometry, registry, table,
                       i.calibration, i.measurement).counts,
    measured_counts, inputs,
    [thickness(0), parameter("fwhm")],
    windows=WindowSet(error=[Window(200, 260)]),
)

print(result.parameters)      # {'thickness[0]': 1197.0, 'fwhm': 14.77}
print(result.uncertainties)
print(result.correlation)     # parameters are rarely independent
```

Always look at the correlation matrix. Thickness and resolution were −0.35
correlated here; strongly correlated parameters mean the data does not
constrain them separately, however tight the individual error bars look.

---

## Sample descriptions (`.lcm`)

Plain text, and the same format RUMP itself reads and writes — pyRUMP
round-trips RUMP's files byte-identically.

```
Sim Reset
Layer 1
 Thick 151 ITO
 Composition In 2 O 3 Sn 0.1 /
Next
 Thick 10 um
 Composition O 4 C 14 H 10 /
Maxpth 200
Foil disable
```

| Command | Meaning |
|---|---|
| `Sim Reset` | Start a new sample |
| `Layer 1` / `Next` | Begin a layer |
| `Thick <v> <unit>` | Thickness — see units below |
| `Composition <El> <n> … /` | Stoichiometry; the `/` terminates the list |
| `Sublayer <n>` | Force a sublayer count |
| `Sthickness <v> <unit>` | Or set sublayer thickness |
| `Equation <name> <params…>` | Depth profile |
| `Species <El> <n> … /` | What the profile blends toward |
| `Fuzzy <amount> <steps>` | Surface roughness |
| `Maxpth <v>` | Default sublayer thickness, 10¹⁵ at/cm² |
| `Absorber <n>` | First *n* layers are a dead layer/window, not sample |

**Thickness units** may be a length (`A`, `nm`, `um`), an explicit areal density
(`/CM2`, `M/CM2`), or a **compound name** from `density.tab`. A compound absent
from the table silently falls back to silicon's density — which is what makes
`Thick 151 ITO` come out at 75.5 × 10¹⁵ at/cm² in both RUMP and pyRUMP.

Commands pyRUMP does not implement (plotting, buffers, the `G_*` global-profile
subsystem) are collected in `script.ignored` rather than raising, so a real file
still loads.

---

## Depth profiles

A layer can vary with depth instead of being uniform:

```
 Equation Linear 0 0.2
 Species Au 1 /
```

blends from pure matrix at the surface to 20% gold at the back. The mixing rule
normalises **both** compositions first:

```
composition(x) = f(x)·species_normalised + (1 − f(x))·matrix_normalised
```

Available: `Constant`, `Linear`, `Erfc`/`Error`, `Exponential`, `Semi-infinite`,
`Thinfilm`, `BuriedThinFilm`, `Thickfilm`/`Thicfilm`, `Timedependent`,
`Gaussian`/`Implant`, `Edgeworth`.

`Gaussian`, `Thinfilm` and `BuriedThinFilm` are *integral* forms — they place an
exact dose per sublayer regardless of how coarse the grid is. The others sample
the sublayer centre and need enough sublayers to be accurate; each equation
carries a recommended count, which overrides `Maxpth`.

`Spline` and `Usereqn` are **not implemented** — they need GENPLOT's spline
fitter and expression evaluator. They raise rather than silently returning zero.

---

## Things that will catch you out

**`phi` is not the scattering angle.** It is 180° minus it. A detector at 170°
means `phi = 10`. Use `geometry.scattering_angle` for the physical value.

**Fitting windows must cover channels where the model has counts.** Poisson
likelihood is undefined where the model predicts zero, so those channels
contribute neither to χ² nor to the gradient. A window reaching past the
spectrum silently throws away most of its own evidence, and a parameter can sit
motionless while the fit reports success:

```
warning: 81 of 81 windowed channels had zero predicted counts
  thickness[0]                     995.54     ← did not move
```

pyRUMP reports the count and warns; RUMP's manual mentions it in one line.
Narrow the window and refit.

**Thickness is areal density**, 10¹⁵ atoms/cm². Converting to nanometres needs
an assumed atomic density.

**Straggling is off by default**, matching RUMP. Set `straggle=1.0` for the
Bohr value. Note the in/out combination is an approximation and there is no Chu
correction.

**Absolute yields depend on charge, solid angle and efficiency**, which are
rarely known to better than a few percent. Use a normalisation window rather
than trusting them — the worked example above is 9× out for exactly this reason.

**Neighbouring heavy elements are not separable.** In and Sn differ by 0.6
channels at 3 MeV; fit their ratio, not each independently.

**Build the stopping registry once.** It refits polynomials per beam energy, so
recreating it inside a fit loop is slow for no reason.
