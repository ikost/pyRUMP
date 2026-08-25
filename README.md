# pyRUMP

A clean Python reimplementation of **RUMP**, the Rutherford backscattering
spectrometry (RBS) simulation and analysis package originally written by
L. R. Doolittle and M. O. Thompson at Cornell. 

This implementations is written using Claude Code.

The original is ~22k lines of unmaintained C from the late 1980s, with a 1996-era
HTML manual and no active support. pyRUMP reproduces its physics as a tested,
importable library, with both a batch CLI and RUMP's own interactive shell.

## Changelog

### 1.0.1 (unreleased)

- `pyrump --version` / `-v` prints the installed version.
- The interactive shell prints its version in the startup banner.
- `QUIT`/`BYE` now ask for confirmation before leaving the interactive shell
  (skipped for piped input and macros, so scripting is unaffected).
- In SIM and PERT mode, `QUIT`/`q` now just returns to the RUMP main menu
  (like `RETURN`) instead of falling through to the exit prompt; the exit
  confirmation only triggers once you're already back at the main menu.
- `COMPARE` (data vs. simulation with residuals) is now usable directly from
  SIM and PERT, without leaving that level first. SIM already documented this
  as a synonym in the original RUMP; PERT is a new pyRUMP addition, handy for
  checking fit quality right after `GO`.

### 1.0.0

- Initial release.

## Install

```bash
pip install -e .
```

Python 3.9+, numpy, scipy. The four physics data tables pyRUMP needs at
runtime ship with the package, so nothing further is needed for simulation,
fitting, or the interactive shell.

## Contents

- [Changelog](#changelog)
- [Quick start](#quick-start)
- [Interactive shell](#interactive-shell)
- [CLI reference](#cli-reference)
- [Python API](#python-api)
- [Worked examples](#worked-examples)
- [Sample descriptions (`.lcm`)](#sample-descriptions-lcm)
- [Depth profiles](#depth-profiles)
- [Things that will catch you out](#things-that-will-catch-you-out)
- [Design and validation](#design-and-validation)
- [How the simulation works](#how-the-simulation-works)
- [RUMP quirks and defects found while porting](#rump-quirks-and-defects-found-while-porting)
- [Milestones](#milestones)
- [Contributing](#contributing)
- [Licensing and provenance](#licensing-and-provenance)
- [References](#references)

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
exactly as in the original. See [Interactive shell](#interactive-shell) for
the full session and command set.

Or drive it as one-off batch commands:

```bash
pyrump simulate sample.lcm --energy 2.0 --beam 4He -o out.rbs
pyrump fit sample.lcm measured.rbs --vary thickness:0 --window 190 226
pyrump plot measured.rbs --compare out.rbs -o comparison.png
pyrump convert measured.rbs measured.dat
```

See [CLI reference](#cli-reference) for every option, or [Python API](#python-api)
to call the library directly.


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

### Session and mode commands

| Command | Effect |
| --- | --- |
| `?` / `HELP` | list the commands available at the current level |
| `SIM` | enter the sample-description editor, its own prompt |
| `PERT` | enter the fitting sub-processor, its own prompt |
| `RETURN` | leave `SIM`/`PERT` back to the RUMP level |
| `DATA [dir]` | print, or reload the atomic tables from, a data directory |
| `QUIT` / `BYE` | leave pyRUMP (asks to confirm, when run interactively) |

`SIM <command>` also runs one SIM command without leaving the RUMP level, e.g.
`sim thick 500 A` — handy inside a one-line macro or when you only need to
tweak one thing.

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
| `GET <file\|n>` | read a file into a buffer, or point at buffer *n* (`READ` reads a file only) |
| `POINTAT <n>` | point at buffer *n*, by number only |
| `BUFFERS` | list the buffers, marking the active one |
| `ACTIVE` | print the active buffer's full parameter set |
| `EMPTY [n]` | reset buffer *n* (default: active, or the first free one) to blank |
| `COPY a b` / `MOVE a b` | copy / exchange |
| `RELEASE [n]` / `NEWALL` | drop one buffer (default: active) / drop all |
| `WRITE f.rbs` / `WRASCII f.dat` | save the active buffer, binary or text |
| `RECALCULATE` | force buffer 0 (the simulation) to recompute |

Unlike the original there is no ten-buffer limit and nothing is destroyed to
make room. Buffer 0 has **no simulate command**: it is recomputed whenever the
sample or the active buffer's parameters change, which is how RUMP behaved.

```
Your wish? get measured.rbs     /* reads the file into a new buffer   */
Your wish? get 0                /* point back at the simulation       */
Your wish? copy 0 2              /* snapshot the simulation into buffer 2 */
```

### Buffer and spectrum parameters

Each buffer carries its own beam, geometry, calibration and measurement
metadata. Every one of these **prints the current value with no argument, and
sets it (echoing the new value) with one** — and chains onto any further
command left on the line, so `Choff 0 FWHM 15` works in one go, exactly as
RUMP's own `WRASCII` output writes it back.

| Command | Sets |
| --- | --- |
| `BEAM 4He++` | beam species and charge state |
| `MEV <energy>` | beam energy, MeV |
| `THETA <deg>` | sample tilt |
| `PHI <deg>` | 180° minus the scattering angle |
| `PSI <deg>` | exit angle (GENERAL geometry only) |
| `GEOMETRY cornell\|ibm\|general` | detector geometry convention |
| `CONVERSION <keV/ch> [keV(0)]` | energy calibration |
| `CORRECTION <factor>` | normalization fudge factor |
| `CHARGE <uC>` | integrated beam dose |
| `CURRENT <nA>` | average beam current — enables pile-up with `TAU` |
| `CHOFF <n>` | channel number of the first data point |
| `FWHM <keV>` | detector resolution |
| `OMEGA <msr>` | detector solid angle |
| `TAU <us>` | MCA shaping time |
| `IDENTIFIER <text>` | free-text spectrum description |
| `DATE <text>` | when the spectrum was measured |
| `FILENAME <name>` | recorded source filename |

```
Your wish? beam 4He++
  beam Z=2 mass=4.0026 charge state 2
Your wish? mev 2.0
  MeV = 2
Your wish? conversion 5.0 0
  5 keV/channel, offset 0 keV
```

`SWALLOW [-twocolumn]`, used inside an `XEQ` macro, reads the macro file's
following lines straight into the active buffer as channel data (or
channel/value pairs), stopping at the first blank line — how a RUMP-written
`.cmd` file reconstructs a spectrum inline.

### Plotting

The plot is one persistent matplotlib window whose state survives between
commands.

| Command | Effect |
| --- | --- |
| `PLOT [buffer\|file]` | erase and plot a buffer (default: active) or file |
| `OVERLAY [buffer\|file]` | add another trace to the current plot |
| `REPLOT` | redraw the current plot, unchanged |
| `AXIS` | draw empty axes, with no data |
| `COMPARE` | active buffer vs. the simulation, with Poisson residuals |
| `DISPLAY` | sample composition vs. depth (from the SIM description) |
| `REGION lo hi` | channel range shown |
| `EXPAND lo hi` | narrow the current region and redraw |
| `COUNTS lo [hi]` | yield range shown |
| `BLOWUP <max>` | shorthand for `COUNTS 0 <max>` |
| `LINEAR` / `SQRT` / `LOG` | yield axis scale |
| `NORMALIZE` / `RAW` | normalized vs. raw yield units |
| `LABELS [off]` | axis labels on or off |
| `ENERGY [off]` | x axis in energy (keV) rather than channel |
| `PARMS` / `PARAMETERS` | print the current plot settings |

```
Your wish? plot 1               /* erase and plot buffer 1            */
Your wish? overlay 0            /* add the simulation on top          */
Your wish? region 100 400
Your wish? sqrt                 /* redraws immediately, sqrt yield    */
```

matplotlib installs by default with pyrump.

### Analysis

Element identification, calibration, and quantification, ported from RUMP's
`anlytc.c` command family. All of these act on the active buffer; region
arguments are plain 0-based channel indices, matching `INTEGRAL`'s existing
convention (not RUMP's own `first`-relative channel numbering).

| Command | Effect |
| --- | --- |
| `ELEMENT el [el ...]` | expected K, energy and channel of each element's surface edge |
| `MATRIX el` | expected energy, channel **and matrix height** for one element |
| `WHATISIT <channel>` | identify the elements whose surface edge is nearest a channel |
| `INFO el` | full report: density, K, cross section, stopping factors, isotopes |
| `INTEGRAL lo hi` | gross/net counts over a channel range (background-corrected net) |
| `THICKNESS lo hi el` | INTEGRAL plus conversion to atoms/cm² and Angstroms |
| `BACKGROUND lo1 hi1 lo2 hi2 order [-inplace] [-noplot]` | fit and strip a polynomial background |
| `SMOOTH [-sv\|-conv\|-fft] [-range lo hi] [n]` | smooth the active buffer |
| `FFT lo hi width` | same as `SMOOTH -fft -range lo hi width` |
| `WIDTH_THICK ch1 ch2 el` | thickness from a peak's half-height width |
| `CALIBRATE ch1 el1 ch2 el2 [energy channel]` | set keV/channel and keV(0) from two known peaks |
| `INTSET [Round\|Interp\|Surface\|Estimated\|Query\|?]` | INTEGRAL/THICKNESS rounding and alpha mode |
| `CURSOR` | not available in this shell — there is no interactive graphics device |
| `PROFILE` | not implemented — never was, even in the original |

`SMOOTH -conv`'s characteristic width uses RUMP's own (nonstandard)
`sigma = (FWHM/2)/sqrt(ln 2)/kevch` — not the usual `FWHM/(2*sqrt(2 ln 2))` —
reproduced deliberately, not corrected. `SMOOTH`'s default range is the whole
buffer; `-conv`'s iteration count and `-fft`'s width both come from a trailing
number, interpreted according to whichever mode is active.

`INTSET` picks two independent modes that both `INTEGRAL` and `THICKNESS`
honor: whether a region's boundaries are rounded to the nearest channel or
interpolated between them, and (for `THICKNESS` only) whether its second,
"compensated" pass uses an estimated alpha or asks you for one. Two-peak
calibration, then a thickness that uses it:

```
Your wish? calibrate 226 Si 369 Au
 Energy=2.0000 MeV    Conversion:4.9896 keV/ch   3.8329 keV(0)

Your wish? intset estimated
Your wish? thickness 180 220 Si
 Discrete integration on buffer 1
 Region:  180.0 to  220.0  Gross:     3498.31  Net:      -47.34  (#/uC/msr)
 Si surface approximation, density  2.32 g/cc
  (Gross)  2.2598e+18 Atoms/cm**2 ( 4539.9 Angstroms)
  ( Net ) -3.0582e+16 Atoms/cm**2 (  -61.4 Angstroms)
 Compensated calculation (Chu et al. page 65)
  (Gross)  2.0903e+18 Atoms/cm**2 ( 4199.4 Angstroms)
  ( Net ) -3.7411e+16 Atoms/cm**2 (  -75.2 Angstroms)
```

The negative "Net" values above aren't a bug: `180`-`220` sits on a flat part
of this (simulated, noiseless) spectrum, and the discrete net-background
correction assumes a sloped continuum either side of the region it's
integrating — pick regions either side of a real peak, not the middle of a
plateau, for a meaningful net figure.

Non-Rutherford (tabulated-resonance) cross sections aren't wired into any of
these — see the [known-limitations note](#milestones).

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

#### SIM commands

| Command | Effect |
| --- | --- |
| `?` / `HELP` | list the SIM commands |
| `RETURN` / `ABORT` / `QUIT` | return to the RUMP level |
| `LAYER <n>` | move to layer *n* |
| `NEXT` | move to the next layer, opening one if needed |
| `OPEN` | insert a blank layer above the current one |
| `RESET` | reset the sample to empty space |
| `SHOW` | display the sample description |
| `STATUS` | summarize layers, maxpth, straggle, multiple |
| `THICKNESS <v> <unit>` | this layer's thickness |
| `COMPOSITION El n [El n …] /` | this layer's stoichiometry |
| `SPECIES El n [El n …] /` | the impurity species an `EQUATION` blends toward |
| `EQUATION <name> <params…>` | depth-profile equation for this layer |
| `EQLIST` | list the known equation names |
| `FUZZ <amount> <steps>` | roughen the interface above this layer |
| `SUBLAYER <n>` | force a sublayer count |
| `STHICKNESS <v> <unit>` | or set the thickness of each sublayer |
| `MAXPTH <v>` | default sublayer thickness, 10¹⁵ at/cm² |
| `STRAGGLE <v>` | Bohr straggling multiplier |
| `ABSORBER <n>` | first *n* layers are a dead layer/window, not sample |
| `MULTIPLE <v>` | multiple-scattering amount |
| `GET <file>` | read a sample description from a `.lcm` file |
| `SAVE <file>` | write the sample description to a `.lcm` file |
| `DENSITY [pattern]` | list known thickness units, or matching `density.tab` compounds |
| `SPLOT` | overlay the simulation on the current plot |
| `COMPARE` | active buffer vs. the simulation, with residuals |

```
SIM Command: layer 1
SIM Command: thick 500 A
SIM Command: composition Si 1 /
SIM Command: density ito       /* look up how "ITO" resolves as a Thick unit */
```

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

#### PERT commands

| Command | Effect |
| --- | --- |
| `?` / `HELP` | list the PERT commands |
| `RETURN` / `QUIT` | return to the RUMP level |
| `GO` | run the fit |
| `PARMS` | display the current selection and windows |
| `CLEAR` | forget every selected parameter and window |
| `WINDOW lo hi` / `WINDOW clear` | add / clear an error window, in channels (up to 10) |
| `NORMALIZE lo hi` / `NORMALIZE off` | set / clear the normalisation window |
| `SINGLE` / `MULTI` | fit one parameter at a time / all together (default) |
| `VOLUME [off]` | verbose progress messages |
| `THICKNESS <layer>` | vary a layer's thickness |
| `COMPOSITION <layer> <El>` | vary one element's composition in a layer |
| `SPECIES <layer> <El>` | vary the `EQUATION` species composition |
| `EQUATION <layer> <n>` | vary equation parameter *n* |
| `MEV` / `FWHM` / `THETA` / `CORRECTION` / `STRAGGLE` | vary that beam, detector or sample parameter |
| `FUZZ` | not implemented — raises an error |
| `COMPARE` | active buffer vs. the simulation, with residuals |

```
PERT Command: window 355 375   /* compare only here                  */
PERT Command: thickness 1      /* vary layer 1's thickness           */
PERT Command: mev              /* also vary the beam energy          */
PERT Command: go
```

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

`FAITHFUL OFF`/`FAITHFUL ON` toggles the session between the shipped C's
bug-for-bug behaviour (the default) and the corrected physics available at
that point in the port — see
[Design and validation](#design-and-validation). It's a session setting, not
persisted on its own, so put it in `~/.pyrumprc` for a standing per-user
default:

```
$ cat ~/.pyrumprc
faithful off
```

`--faithful on`/`--faithful off` overrides that for one invocation, applied
after `~/.pyrumprc` runs but before any macro passed on the command line —
the macro can still set `FAITHFUL` itself if it needs to.

## CLI reference

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



## Python API

The CLI is a thin wrapper; the library is the primary interface.

```python
from pathlib import Path
import pyrump
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

## Worked examples

Every number below was produced by running the code, not written from memory.
For what the simulation actually computes, see
[How the simulation works](#how-the-simulation-works).

### Identifying what is in a sample

Given an unknown spectrum, the first question is which elements are present.
Each element's **surface edge** sits at *K*·*E*₀, so predicted edge positions
identify the peaks.

Using `2A.rbs`, one of the files shipped with RUMP (see
[Licensing and provenance](#licensing-and-provenance) — it isn't redistributed
with pyRUMP, so point `PYRUMP_DATA` at your own copy of the legacy `rump/data/`
tree to reproduce this):

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

### Simulating a known structure

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

That is not a simulation error: running the same case through the original C
gives 381 120 counts, agreeing with pyRUMP to 2.6e-3. Both codes say the same
thing, so the discrepancy lives in the measurement's normalisation — the
actual collected charge, solid angle, or detector efficiency differs from the
values recorded in the file. **This is the normal situation in RBS**, which is
why RUMP has both a `CORR` factor and a normalisation window. Rather than
trusting the charge integration, you fit the scale:

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

### Fitting a thickness

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
recovers it to better than 0.1%. Note that thickness is reported in **areal
density**, not Ångström — RBS measures atoms per unit area, and converting to
a physical thickness needs an assumed density, a separate and often less
certain quantity.

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

## Design and validation

**Faithful first, corrected by choice.** The default reproduces the shipped C
bug-for-bug, because that is what every published RUMP result was produced
with. Known defects — and there are several — are reproduced exactly, with
the mathematically correct behaviour available behind explicit flags. The
shell exposes this as a session setting, `session.settings.faithful`,
toggled with the `FAITHFUL` command and persisted through `~/.pyrumprc` (see
[Macros](#macros)) rather than a separate branch or fork — corrected and
faithful behaviour live in the same codebase so they stay comparable against
the C oracle side by side. See
[RUMP quirks and defects found while porting](#rump-quirks-and-defects-found-while-porting).

**Validated against the original**, at two levels. The legacy C is compiled
into a shared library and called directly from the test suite, so each stage
is compared function-by-function rather than by eyeballing a final spectrum.

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
[Licensing and provenance](#licensing-and-provenance)). Point `PYRUMP_C_REFERENCE`
at it, or place it at `C-code/`. Tests skip cleanly when it is absent.

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

## How the simulation works

A description of the forward model pyRUMP implements — what it computes, in what
order, and where it approximates. Written for someone who wants to know what the
numbers mean, not for someone reading the source.

The algorithm originates with L. R. Doolittle's 1985 paper. **pyRUMP follows the
shipped C, which diverges from that paper in places** — see
[Divergences from the published algorithm](#divergences-from-the-published-algorithm).

### The physical problem

A beam of light ions (usually 1–3 MeV ⁴He or ¹H) hits a sample. A few
backscatter from target nuclei and reach a detector at a fixed angle. The
detected energy encodes two things at once:

- **which nucleus** it bounced off — heavier targets take less energy from the
  projectile, so each element has a characteristic maximum energy
- **how deep** the collision happened — the ion loses energy travelling in and
  out again, so deeper collisions arrive lower in energy

A spectrum is therefore a depth profile smeared together across all elements
present. Simulating it forward, and adjusting the sample until the simulation
matches, is how the depth profile is recovered.

#### Kinematics

An elastic collision at scattering angle *φ* leaves the projectile with a
fixed fraction of its energy, the **kinematic factor**:

$$K = \left[\frac{\sqrt{1 - (x\sin\phi)^2} + x\cos\phi}{1+x}\right]^2,
\qquad x = \frac{m_1}{m_2}$$

For 2 MeV ⁴He at 170°, *K* is 0.566 on silicon and 0.922 on gold. Those factors
place the **surface edge** of each element: the highest energy at which it can
appear.

> **Angle convention.** RUMP stores `phi` as *180° minus* the scattering angle —
> a detector at 170° is entered as `phi = 10`. pyRUMP keeps that convention on
> `Geometry.phi` and exposes the physical angle as `Geometry.scattering_angle`.

#### Cross-section

How often a collision happens is Rutherford scattering, in the lab frame:

$$\sigma = \left(\frac{Z_1 Z_2 e^2}{4E}\right)^2 \frac{4}{\sin^4\phi}
  \frac{\left[\sqrt{1-(x\sin\phi)^2}+\cos\phi\right]^2}{\sqrt{1-(x\sin\phi)^2}}$$

The **Z₂²** dependence is why RBS is sensitive to heavy elements in a light
matrix and nearly blind the other way round. The **1/E²** means yield rises with
depth as the beam slows.

A screening correction (L'Ecuyer) reduces this slightly at low energy. RUMP has
no relativistic correction and does not implement Andersen screening.

### The brick

The central data structure, in the 1985 paper's own words:

> Each simulated spectrum is made up of the superimposed contributions from each
> isotope of each sublayer in the sample. Any such contribution will be referred
> to as a **brick**.

Every layer is cut into **sublayers**, and each sublayer contributes one brick
per isotope of every element present. A brick is a trapezoid in energy space:

```
        h_front  ___________
                /           \___
               /                \  h_back
              |                  |
        e_back                e_front     (energy increasing to the right)
```

- **e_front** — energy of particles scattered from the sublayer's *front* face
- **e_back** — from its *back* face, lower because the beam travelled further
- **h_front, h_back** — differential yield at each edge

A simulated spectrum is nothing more than every brick integrated onto the
detector's channels and summed. Natural silicon gives three bricks per sublayer
(²⁸Si, ²⁹Si, ³⁰Si), each with its own kinematic factor, so isotopes appear as
slightly displaced copies.

#### Why sublayers, and how many

Thicker sublayers mean fewer bricks and a faster simulation, but a coarser
approximation. RUMP's step size is **path-length based, not depth based**:

```
maxpath        = maxpth / max(|sec θ_in|, |sec θ_out|)
n_sublayers    = 1 + areal_thickness / maxpath
```

so tilting the sample automatically produces more, thinner sublayers. `maxpth`
defaults to 200 (in 10¹⁵ atoms/cm²) — roughly where the paper's error analysis
puts the energy-loss expansion at 10⁻⁵ fractional error.

A layer carrying a depth-profile equation ignores `maxpth` and uses a
per-equation recommended count instead (5 for `Constant`, 30 for `Thinfilm`).

### The pipeline

Stage order matters and is not arbitrary.

```
 stopping tables                fitted once per beam, then never recomputed
        │
        ▼
 slab discretization            layers → sublayers, composition per slab
        │
        ▼
 inbound march                  beam energy at every interface
        │
        ▼
 per-isotope depth loop         → bricks (the expensive part)
        │
        ▼
 brick → channel fill           trapezoid, or triangles if straggling
        │
        ▼
 detector convolution           Gaussian, applied to the whole spectrum
        │
        ▼
 yield normalisation            × Ω·Q/(charge state · CORR)
        │
        ▼
 pile-up  →  multiple scattering
```

#### 1. Stopping powers — the indirection that matters most

Energy loss per unit depth, *ε(E)*, comes from one of several models tried in
priority order:

1. **Konac/Kalbitzer** fits, for the specific ion/target pairs in `newstop.kal`
   (H, D, ³He, ⁴He on carbon and silicon)
2. **Ziegler ZBL85**, the general fallback for Z = 1…92
3. hard-coded polynomials for Mylar, the Z=93 pseudo-element

Note the order: **for ⁴He on silicon — the most common RBS measurement there is
— Ziegler is not used.** Konac wins, and the two differ by up to ~10%.

Then the step that surprises everyone:

> **RUMP never evaluates the stopping model during a simulation.** At the start
> of a run it samples the chosen model at 201 points and least-squares fits a
> degree-5 polynomial in √E. Every subsequent energy-loss calculation evaluates
> only that polynomial.

The fit window is tied to the beam energy — `[0.04·E₀, 1.15·E₀]` — so
coefficients change when the beam does, and the polynomial diverges rapidly
outside that range. Any reimplementation that calls the stopping model directly
disagrees with RUMP everywhere.

Per sublayer the elemental polynomials are combined by **Bragg's rule** —
linear additivity weighted by areal density, with no compound correction — which
also folds in the thickness, so the coefficients directly give eV through the
slab.

#### 2. Inbound march

The beam energy at each interface is computed once, before any element is
considered, because the incoming path does not depend on what it eventually
scatters from.

Energy loss across a sublayer uses a **third-order Taylor expansion** of
d*E*/d*a* = −*ε*(*E*):

$$E(a) = E_0 - a\varepsilon + \tfrac{1}{2}a^2\varepsilon\varepsilon'
        - \tfrac{1}{6}a^3\left(\varepsilon''\varepsilon^2 + \varepsilon'^2\varepsilon\right)$$

Truncating after the first term would be the familiar *surface approximation*.
The extra terms are what let RUMP use thick sublayers and stay fast — that is
the "rapid" in the 1985 paper's title. It is not a faster inner loop; it is
higher accuracy per sublayer, so that fewer are needed.

A **cutoff energy** (3% of the beam energy) ends the march: below it the
stopping fit is not trustworthy.

#### 3. Outbound path — the expensive half

Once scattered, a particle must be walked back out through every overlying
sublayer. Done for every sublayer and every isotope, this is the O(*N*²) term
that dominates the cost, and the reason RUMP works so hard to keep *N* small.

Two things come out of it: the exit energy, and `ratde` — the accumulated ratio
of stopping powers, which accounts for the scattered beam's energy spread
changing on the way out.

#### 4. Yield

Per sublayer and isotope:

$$h = \sigma(E)\, N_{\text{slab}}\, \frac{\sec\theta_{in}}{[\varepsilon]}
      \cdot \texttt{ratde} \cdot f_{\text{isotope}}$$

where **[ε]** is the stopping cross-section factor of Chu et al.,
$[\varepsilon] = K\varepsilon(E)\sec\theta_{in} + \varepsilon(KE)\sec\theta_{out}$,
which converts a depth interval into an energy interval.

#### 5. Filling channels

Without straggling, each trapezoid is integrated exactly onto the channels it
overlaps, including partial channels at both ends.

With straggling, the trapezoid is abandoned: it is split into **two triangles**,
each convolved with a Gaussian of its own width. Straggling is **Bohr only** and
**off by default** in RUMP.

#### 6. Detector resolution

A Gaussian convolution applied **once to the finished spectrum**, not per
sublayer — since 1994 the detector width was deliberately removed from the
straggling path. The kernel is channel-*integrated* (differences of the normal
CDF, not point samples) and truncated at 3σ.

Counts falling outside the channel range are discarded at both ends, so the
convolution is not count-conserving near the edges. Defensible — a real
multichannel analyser cannot record counts in channels it does not have — but it
matters when comparing totals.

#### 7. After the spectrum exists

- **Pile-up** — two events arriving within the detector's shaping time recorded
  as one of their combined energy. Needs the beam current and shaping time.
- **Multiple scattering** — an empirical low-energy tail with **no physical
  basis**; the C's own comment calls its scale factor "ad-hoc". Treat any number
  it produces as qualitative.
- **Fuzz** — surface roughness, as several simulations at Gaussian-weighted
  thicknesses. Every roughened layer multiplies the cost.

### What is approximated, and what is absent

Worth knowing before trusting a result.

| | |
|---|---|
| Straggling | Bohr only, **off by default**. The in/out paths are combined *linearly* rather than as `K²σ²_in + σ²_out` — an approximation, not the correct combination. No Chu correction. |
| Screening | L'Ecuyer only; no Andersen. |
| Relativity | No correction anywhere. |
| Multiple scattering | Empirical tail with no physical basis. |
| Channelling | Not modelled at all. |
| Nuclear reactions | Q ≠ 0 reactions are rejected outright. |
| Sample | Laterally uniform apart from the `Fuzz` roughness model. |
| Beam and detector | Treated as points — no finite spot or acceptance angle. |

The 1985 paper is explicit about most of this:

> The algorithms assume a laterally uniform sample and neglect the effects of
> channelling, core electron screening, nuclear reactions, nuclear resonances,
> and multiple scattering. They also ignore the effects brought about by a
> finite size beam spot and detector.

### Divergences from the published algorithm

The shipped C is not the program the papers describe.

**The parabolic brick is gone.** The 1985 paper's signature contribution — a
parabolic brick top with an exact analytic area from the integral ∫E(a)⁻²da —
sits behind `#if 0` in the C, labelled *"no longer using Doolittle qqq code"*.
The shipped path uses plain trapezoids. The integral is still computed on every
run and then discarded.

pyRUMP follows the C. It also carries the discarded integral on
`Bricks.area`, so the paper's algorithm could be offered later as an opt-in.

**Konac/Kalbitzer stopping** is not in the 1996 manual at all, yet it takes
priority over Ziegler for the most common measurements.

**Defects are reproduced deliberately.** pyRUMP's default is bug-for-bug
fidelity, because every published RUMP result was produced with those bugs
present. The catalogue is in
[RUMP quirks and defects found while porting](#rump-quirks-and-defects-found-while-porting)
— 20 entries. The one that affects results most directly is a wrong coefficient
index in the second-derivative macro (35–50% error in d²ε/dE², feeding the
third-order energy-loss term). `StoppingTable.derivative(..., faithful=False)`
gives the correct value.

The single deliberate departure is `SimStragf`, RUMP's seven-regime rational
approximation to the triangle⊗Gaussian integral. pyRUMP uses a closed `erf`
form instead: measured against numerical quadrature, the closed form is exact to
6e-14 where RUMP's fit carries 1.7e-6. The 1985 justification for the
approximation was explicitly about 1985 hardware.

### Accuracy

Stage-by-stage agreement with the original C is in
[Design and validation § Current agreement](#current-agreement) — every figure
there is limited by float32 storage in the C, not by pyRUMP. For a realistic
multi-layer sample with micron-thick polymer layers, agreement loosens to
~3e-3: thousands of sublayers accumulate single-precision differences, and part
of the beam falls below the stopping cutoff.

Full bibliography in [References](#references).

## RUMP quirks and defects found while porting

Running catalogue of behaviour in the legacy C that is surprising, wrong, or
simply undocumented. Each entry says what pyRUMP does about it.

Default policy (per the project plan): **reproduce faithfully, expose the fix
behind a flag.** A silent "correction" would make pyRUMP disagree with every
published RUMP result.

The mechanism for that flag is `session.settings.faithful` (toggled by the
`FAITHFUL` command, [Macros](#macros)), or a dedicated `Settings` field for a
correction that needs to be controlled independently. Not every entry below
is wired to it yet — each entry says whether it is.

### 1. The papers describe an algorithm the code no longer uses

The 1985 NIM-B paper's signature contribution — parabolic brick tops with an
exact analytic area from the "Rutherford integral" ∫E(a)⁻²da — is **dead code**.
`anlyz.c:496` sits inside `#if 0`, labelled *"no longer using Doolittle qqq
code"*. The shipped path is a plain trapezoid (`SimAnlyz4`, `anlyz.c:304`) or,
with straggling, two Gaussian-convolved triangles (`SimAnlyz3`, `anlyz.c:244`).

`layer[].qq` is still computed every run in `SimPrecal` and then discarded.

**pyRUMP:** trapezoid by default; the parabolic form is planned as an opt-in mode,
for which the dead C is the specification.

### 2. `SQRT_DDS_POWER` uses the wrong coefficient index

`stopping.h:47-49` computes d²S/dE² as

```c
(((((3.75*p[5])*e + 2*p[4])*e + 0.75*p[3])*e*e - 0.25*p[2]) / (e*e*e))
                                                       ^^^^ should be p[1]
```

The header's own comment two lines earlier gives the correct form,
`DDS = 1/(4e⁴) · Σ i(i-2)·aᵢ·eⁱ`, whose `i=1` term is `−p[1]·e`. There is no
`p[2]` term at all, since `i(i-2)` vanishes at `i=2`.

Verified both directions: the oracle reproduces the buggy macro exactly, and a
numerical second derivative reproduces the corrected form exactly. The error is
**35–50%**.

Impact is limited — d²S/dE² only enters the third-order term of the energy-loss
expansion (`creatr.c:1554`) — but it is real and it shifts the depth scale.

**pyRUMP:** `StoppingTable.derivative(..., order=2)` reproduces the bug;
`faithful=False` gives the correct value. Covered by
`test_second_derivative_reproduces_rump_bug`.

### 3. `-DREAL_IS_DOUBLE` silently corrupts every data table

RUMP's table readers hard-code `%f` in their `scanf` formats while writing into
`REAL` fields (`ziegler.c:100,114`; `atomio.c:162,163,173`). With `REAL` as
`double` these write four bytes into an eight-byte field, so the Ziegler and
atomic tables load as denormal garbage and `zstop` returns 0 without complaint.

`REAL_IS_DOUBLE` is never set by any shipped makefile, so the bug stayed latent.

**Consequence:** the plan's dual-precision oracle — build twice, diff to measure
RUMP's own float32 noise — is not possible without patching the readers, which
would make the oracle a *modified* RUMP. The `float` build is authoritative, and
tolerances are argued from single-precision reasoning instead.
`tests/oracle/oracle.py` refuses to load a corrupt build rather than returning
garbage.

### 4. Stopping tables are cached per (Z, mass) for the whole session

`RbsStpfind` (`stopping.c:274-279`) reuses an existing fitted table whenever the
new beam energy merely *fits inside* its window (`2·emin ≤ E ≤ emax`). Simulating
at 3 MeV and then at 2 MeV does **not** refit: the second run silently uses the
3 MeV window, and its coefficients differ from a fresh 2 MeV fit.

When only Z matches (different isotope), the table is reused with
`e_scale = table_mass / beam_mass` — the Amsel energy-scaling trick that lets one
table serve 3He/4He or H/D.

This is stateful behaviour that changes numbers, not an optimisation.

**pyRUMP:** `StoppingTableCache` reproduces the reuse rules, including
`e_scale`. The oracle gained `OracleResetStoppingTables()` so tests are
order-independent.

### 5. The fitted polynomial is only valid inside the beam-dependent window

`emin = 0.04·E_beam`, `emax = 1.15·E_beam` (`stopping.c:316-319`, STOP_SQRT).
Beyond `emax` the degree-5 fit diverges from the underlying model within a few
hundred keV — 28% at 3 MeV for a 2 MeV table — and RUMP returns the extrapolated
value without warning.

More importantly: **RUMP never evaluates Ziegler or Konac during a simulation.**
It fits once at startup and evaluates only the polynomial. Any port that calls
the stopping model directly will disagree with RUMP everywhere.

### 6. Konac/Kalbitzer outranks Ziegler

The priority chain (`stopping.c:479-515`) tries `newstop.kal` *before* ZBL. For
H/D/³He/⁴He on carbon or silicon — the most common RBS cases by far — RUMP is not
using Ziegler at all. The two models differ by up to ~10%.

`newstop.kal` is also not mentioned in the 1996 manual.

### 7. `ConvoluteDetector` loses counts at **both** edges

*(Corrected — an earlier version of this entry described only the high-energy
edge. The behaviour is symmetric.)*

Contributions that would land outside the channel range are discarded:

* the head loop drops the `i-k` half once `k > i` (`creatr.c:1233-1236`)
* the tail loop drops the `i+k` half once `k >= npt-i` (`creatr.c:1284-1291`)

Measured with a 30 keV FWHM at 5 keV/channel on a 200-channel spectrum, a delta
function of 1000 counts retains:

| Position | Retained |
|---|---|
| channel 100 (interior) | 999.15 — only the 3σ truncation is lost |
| channel 2 | 836.34 |
| channel 197 | 836.34 |

So this is not really a bug: a real MCA cannot record counts in channels it does
not have. But it does mean the integral of a convolved spectrum is below the
original whenever intensity sits within ~3σ of either end, which matters when
comparing totals.

**pyRUMP:** `convolve_edge='rump'` (default) reproduces it; `'renormalize'`
conserves counts by scaling each *source* channel's contribution by the kernel
weight that actually lands in range. Note the normalisation must be applied on
the source side — normalising the output instead rescales *received* weight and
inflates the edges rather than conserving.

### 8. The SIM `RECALCULATE` command is misspelled in the command table

`sim2.c:252` registers `{"recalculculate", 5, SIM_RECALCUL}` — note
`recalcul**cul**ate`. Typing `recalculate` does **not** match and silently falls
through to the top-level shell. Only the 5-character prefix `recal` works.

Relevant to anyone driving the legacy binary; `tests/oracle/driver.py` documents
and uses `recal`.

### 9. Straggling is Bohr-only, off by default, and combined incorrectly

`sample->straggle` defaults to 0 (straggling disabled). When enabled, the inbound
and outbound path variances are combined *linearly* as
`stragc = sec θ_in·K + sec θ_out` applied to the inbound Bohr variance
(`creatr.c:1661`), rather than as `K²σ²_in + σ²_out`. No Chu correction exists
anywhere in the simulation.

### 10. Physics that is absent

* no channelling
* no Q≠0 nuclear reactions (`reswork.c:355` rejects them explicitly)
* no Andersen screening — L'Ecuyer only
* no relativistic correction to the Rutherford cross-section
* multiple scattering is an ad-hoc exponential tail with no physical basis
  (`creatr.c:337-345`), described in the C itself as *"Ad-hoc scaling"*

### Dead or broken code to avoid porting

| Location | Status |
|---|---|
| `anlyz.c:496-642` | `#if 0` — the paper's parabolic algorithm; also syntactically incomplete at `:600-603` |
| `tables.c` | not in `OBJS`; does not compile (three typos). The live density-table loader is `sim2.c:1415` |
| `creatr.c:2100-2142` (`SimInloss`) | `#if 0`; useful as clean documentation of the ΔE expansion |
| `sigma.c:372-558` | `#if 0` — superseded Turos & Meyer cross-sections |
| `MyXsect/Xsect.c` | placeholder; writes to an undeclared variable, hook never called |
| `data/*.stp` | unreachable by default: stores `STOP_LINEAR` coefficients while the runtime is `STOP_SQRT` (`stopping.c:139` vs `:276`) |
| channelling hooks (`sigma_scale`, `dedx_scale`) | declared, never set; the only use is commented out at `creatr.c:1114` |

### 11. `reschk` is compiled out of `reswork.c` but referenced by `creatr.c`

`reswork.c:55` guards the resonance-table index with `#ifdef RESONANCE`, and the
`#define` lives in `xsect.h` — which `reswork.c` includes *after* some of its own
declarations. Depending on include order the array is not emitted, while
`creatr.c:1723` reads `reschk[z2]` unconditionally on the simulation path.

In the shipped build this happens to link because other translation units pull
`xsect.h` in first. It is fragile, and it bit the oracle build: a link that looks
successful can leave the resonance index pointing at nothing.

**pyRUMP:** `reswork.c` is linked into the oracle explicitly so the array exists
and is properly zero-initialised (no resonance tables loaded ⇒ pure Rutherford).

### 12. `SimStragf` rescales its own argument

`SimStragf(x, sig)` (anlyz.c:371) computes `newx = x * (1 + 3*sig)` as its first
real statement (anlyz.c:387). Its `x` is therefore in units of the **broadened**
width `|de| + 3*sqrt(2)*sigma`, not of the triangle base — which is what the
docstring's "triangle ... height 0 for x<0 or x>1" implies.

The rescaling exactly undoes the `fact = 1/(|de| + 3*sqrt(2)*sigma)`
normalisation the caller applies (anlyz.c:257), so the composition collapses to
`(E_j - E_peak)/|de|`.

Overlooking this makes the function appear wrong by up to **0.23** — comparable
to its entire range of 0.5 — while looking perfectly plausible in isolation.

**pyRUMP:** `stragf()` reproduces the C's contract including the rescaling;
`triangle_gaussian_integral()` is the underlying maths in unscaled coordinates.

### 13. Measured: `SimStragf`'s approximation error

The 1985 paper justifies the rational fit on 1985 hardware grounds:

> Analytic expressions for the functions f and g … are of limited utility.
> Direct evaluation is unnecessarily slow and often involves finding the small
> difference of large numbers. Single precision computation is inadequate …

Against numerical quadrature, over sigma in [0.05, 3] and x in [-1.5, 1.7]:

| | Worst absolute error |
|---|---|
| pyRUMP closed `erf` form | **6.4e-14** |
| RUMP `SimStragf` | **1.7e-6** |

So the approximation is good — its claim of "rapidly and accurately" holds — but
it is ~8 orders of magnitude looser than the closed form is in float64. pyRUMP
therefore does *not* port it: the closed form is more accurate, vectorises, and
avoids transcribing a table of hand-tuned constants.

This is the one place pyRUMP knowingly departs from bit-faithfulness. The
resulting spectrum still matches the C to ~3e-5 of peak.

### 14. `hfront` is recomputed per slab; only `efront` and `ratde` carry over

`SimCideal` guards the front-edge recomputation with `ok` (creatr.c:1791):

```c
if (! ok) {                       /* Recompute EFRONT & RATDE */
    efront = km2 * ein;
    SimFlyout(lay-1, &efront, &ratde);
}
rfront = strct[elno] * ratde * (secin/RbsEfact(...)) * fisot;   /* :1804 -- outside */
hfront = sigma * rfront;
```

Only the *geometry* is cached. The front **height** is recomputed every slab
from that slab's own areal density, because it is outside the guard.

With uniform composition the two are indistinguishable — `strct` is identical in
every slab — so a port can pass every uniform-layer test while getting this
wrong. It only shows up with a depth profile, where reusing the height shifts
the entire spectrum by one slab.

### 15. Layer density is an inverse-density average, and it sets the depth scale

`creatr.c:606-625` averages **cm³/atom**, not atoms/cm³ — the C's own comment
calls it "the idea of hard ball packing":

```
rho = ( sum_i (x_i / rho_i) / sum_i x_i ) ^ -1
```

The result converts areal thickness to physical thickness
(`cm_thick = cm2_thick/matrix_density/1E8`, creatr.c:692), so it sets the depth
scale of every depth-dependent profile equation. Getting it wrong — for example
using the 0.4997 silicon *fallback* constant instead of silicon's actual 0.49777
— leaves position-fraction forms (CONSTANT, LINEAR) matching perfectly while
every physical-depth form is off by a few percent.

The species composition gets its own density, used to convert Angstrom doses
(THINFILM, BURIEDTHINFILM) into 1e15 at/cm².

A pre-1997 `COMPATIBLE` mode averaged densities directly; pyRUMP implements only
`IMPROVED`, the shipped default.

### 16. `poly_e` reads one element past its array — and `-O2` makes it fatal

`poly_e(x, numer, iorder)` (gvcalc.c:4914) does

```c
numer += iorder;  tmp = *numer;
while (iorder--) tmp = tmp*x + *(--numer);
```

i.e. it touches `numer[0 .. iorder]` — **iorder+1** values. `NDTRI` calls it as
`poly_e(t*t, taylor, 10)` against a 10-element `taylor[]`, so it reads one past
the end.

At `-O0` the adjacent static happens to be tiny and Horner folds it to nothing,
so the bug is invisible. At `-O2` the compiler exploits the undefined behaviour
and `NDTRI` returns ±0.15 for **every** argument in roughly (0.15, 0.85) —
correct in the tails, badly wrong in the middle.

**The shipped RUMP passes no optimisation flag at all** (`makeosx.h`: `GOPTS`
has no `-O`, `CCOPT` is empty), so it runs the benign version. Anyone rebuilding
RUMP with optimisation enabled would silently corrupt every FUZZ profile.

**pyRUMP:** `tests/oracle/build_oracle.py` compiles `ndtri_probe.c` at `-O0`
specifically, with a test (`test_ndtri_matches_the_c`) that fails loudly if that
is ever lost.

### 17. The inbound march starts at `fsurf`, skipping the absorber

`SimPrecal` seeds `samm->layer[samm->fsurf].ehit = ee` with the **full** beam
energy (creatr.c:1530) and marches from there. Absorber layers are between the
sample and the detector, so the incoming beam never crosses them — only the
outgoing particle does, via `SimFlyout`.

```c
samm->layer[samm->fsurf].ehit = ee;                 /* creatr.c:1530 */
for (lay=samm->fsurf; lay<samm->num_layers; lay++)  /* creatr.c:1541 */
```

Marching through them on the way in double-counts their stopping and shifts
every edge — it costs ~20% of the total yield for a 200e15 at/cm² silicon
absorber, and grows from there. The absorber is also not tilted with the
sample: `SimFlyout` forces normal incidence through it (creatr.c:1971), since a
detector window does not rotate when the sample does.

### 18. `.RBS` data records use the generic type, not the explicit ones

The format defines five data-record types: `0011h` generic plus `0012h`-`0015h`
naming compression modes 0-3 explicitly. `read_data_records` (rbs_rdwr.c:774)
maps the generic form onto whichever mode the preceding init record declared.

**Every shipped fixture uses `0011h`.** A reader that implements only the
explicit types parses the whole file happily — headers, calibration, geometry,
identifier all correct — and returns a spectrum of zeros, because no data record
ever matches. There is no error and nothing looks wrong until you check the sum.

**pyRUMP:** handles both, and `test_generic_data_record_uses_the_declared_compression`
asserts the counts are non-zero rather than trusting the parse.

### 19. `MASSES:` and `ZEDS:` are ignored, and the shipped files disagree anyway

Both are `Q_IGNORE` in RUMP's header table (reswork.c:131-132); nuclide
identities come from `REACTION:` alone. That is just as well, because the
shipped files use incompatible conventions:

```
boron.adt    MASSES: 11 4 4 11                    <- target first
car_pp.adt   MASSES: 1.0078, 12, 1.0078, 12.0     <- projectile first, commas
```

A reader that trusts `MASSES:` gets boron as Z=4, m=11 — plausible-looking and
wrong. Parse `REACTION:` and ignore the rest.

Related: `QVALUE:` is often a comma-separated list (`0.00, 0.00, 0.00, ...`).
The C reads it with `atof()`, which consumes only the leading number.

### 20. RUMP cannot read its own bundled R33 example

`data/R33.Format` declares `Units: mb`. RUMP's accepted set is exactly
`b/sr`, `mb/sr`, `rtr`, `rr`, `relative` (reswork.c:325-334), so the file is
rejected outright.

It is a *format exemplar* from SigmaCalc, not a loadable table. pyRUMP refuses
it identically — matching the refusal is fidelity, not a gap.

### Notes on driving the engine from outside

`creatr.c`'s output stage is a **function pointer**, `SimFillSpectrum` (sample.h),
called once per brick. Redirecting it captures the engine's exact intermediate
results with the C completely unmodified — no `#ifdef` hooks, no patched copy.
`SimPileup`, `SimInitFillSpectrum` and `SimTermFillSpectrum` are pointers too.

Host state the engine needs, and where the real versions live:

| Symbol | Source | Note |
|---|---|---|
| `RbsNormK` | `bmanip.c:880` | sets the absolute yield scale — must be exact |
| `SimThickConvert` | `sim2.c:2349` | thickness-unit conversion |
| `RbsBuffers[0]` (`ALTBUF`) | `rumpdata.h:156` | the theory buffer |
| `RbsActiveBuf` (`ibuf`) | `rumpdata.h:151` | the data buffer, copied wholesale at `creatr.c:283` |
| `Rmp` | `rumpdata.h:201` | alias for `RbsDataBlock`; only `autsim` is read |
| `sigtab` / `coffe2` | `sim2.c:81-82` | **default `{-1,-1}`**; zero silently selects the manual-override branch and returns a cross-section of exactly zero |

Captured bricks carry the *scattered particle's* `z`/`mass` (for the stopper-foil
lookup at `anlyz.c:180`), **not** the target's — target identity is implicit in
the block ordering, one block per isotope, heaviest first.

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
compare against the legacy C for extra confidence) need the RUMP C source,
which isn't redistributed here — see [Design and validation](#design-and-validation).
They skip cleanly when it's absent, so it's not needed for everyday development.

## Licensing and provenance

pyRUMP is MIT licensed, and an **independent reimplementation**: it is not
affiliated with, endorsed by, or derived from the RUMP source distribution.

RUMP and Genplot were trademarks of Computer Graphic Service, Ltd. (CGS). CGS
ceased operating as a business in June 2012 and `genplot.com` no longer
resolves; the authors stated at the time that GENPLOT and RUMP remain free to
download and use, which removes the trademark concern but **not** copyright in
the original source, which remains with its authors — hence the C tree is
still not redistributed here.

Four data tables are bundled with pyRUMP (`src/pyrump/data/`), independent of
CGS and checked against current CIAAW/NIST values and literature:

- `pscoef.dat` — the ZBL/TRIM `SCOEF` stopping-coefficient table
- `newstop.kal` — Konac/Kalbitzer stopping-power fits
- `atom4.dat` — elements and isotopes
- `density.tab` — compound densities

See `src/pyrump/data/SOURCES.md` for full provenance, verification notes, and
the one correction made (a data-entry error in the GaP density); citations are
in [References](#references).

Non-Rutherford cross-section tables (`*.adt`) are IBANDL evaluations and are
**not** bundled — obtain them separately from IBANDL if you need that data.

## References

- L. R. Doolittle, *Algorithms for the rapid simulation of Rutherford
  backscattering spectra*, Nucl. Instr. Meth. **B9** (1985) 344–351.
- L. R. Doolittle, *A new approach to Rutherford backscattering analysis*,
  Nucl. Instr. Meth. **B15** (1986) 227–231.
- J. F. Ziegler, J. P. Biersack, U. Littmark, *The Stopping and Range of Ions in
  Solids*, Pergamon (1985) — source of `pscoef.dat`, the ZBL/TRIM SCOEF table.
- G. Konac, S. Kalbitzer, Ch. Klatt, D. Niemann, R. Stoll, Nucl. Instr. Meth.
  **B136–138** (1998) 159–165 — source of `newstop.kal`.
- W.-K. Chu, J. W. Mayer, M.-A. Nicolet, *Backscattering Spectrometry*,
  Academic Press (1978) — the [ε] stopping cross-section factor and kinematics.
- S. Baker, R. D. Cousins, Nucl. Instr. Meth. **221** (1984) 437 — the Poisson
  fitting objective.
- J. L'Ecuyer et al., Nucl. Instr. Meth. **160** (1979) 337 — screening
  correction.
- J. F. Ziegler, Nucl. Instr. Meth. **B136–138** (1998) 141 — screening and
  cross-section formulae.
- V. Quillet, F. Abel, M. Schott, Nucl. Instr. Meth. **B83** (1993) 47 —
  screening and cross-section formulae.
- A. F. Gurbich, Nucl. Instr. Meth. **B136–138** (1998) 60 — non-Rutherford
  cross-section evaluations, as distributed via IBANDL.
</content>
