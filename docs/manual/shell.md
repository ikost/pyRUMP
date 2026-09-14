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
(`REGion`, `OVerlay`), so `reg 100 400` and `region 100 400` are the same
command. `?` lists everything, with the required characters upper-cased. One
exception: `COMPARE` requires its full name at every level, since a partial
abbreviation collided with `COMPOSITION` in SIM and PERT; `CMP` works
everywhere instead as an explicit synonym.

Commands tagged `[new]` below have no original-RUMP counterpart — see the
[Changelog](https://github.com/ikost/pyRUMP/releases)'s versioning note for what that means for
`pip install`.

### Session and mode commands

| Command | Effect |
| --- | --- |
| `?` / `HELP` | list the commands available at the current level |
| `SIM` | enter the sample-description editor, its own prompt |
| `PERT` | enter the fitting sub-processor, its own prompt |
| `RETURN` | leave `SIM`/`PERT` back to the RUMP level |
| `DATA [dir]` | print, or reload the atomic tables from, a data directory |
| `QUIT` / `BYE` | leave pyRUMP (asks to confirm, when run interactively) |
| `FAITHFUL [on\|off]` `[new]` | toggle bug-for-bug vs. corrected physics (see [Design and validation](../physics/validation.md)) |

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
| `EMPTY [n]` | scroll a fresh blank buffer into buffer 1 (default), or reset buffer *n* in place |
| `COPY a b` / `MOVE a b` | copy / exchange |
| `RELEASE [n]` / `NEWALL` | drop one buffer (default: active) / drop all |
| `WRITE f.rbs` / `WRASCII f.dat` | save the active buffer, binary or text |
| `RECALCULATE` | force buffer 0 (the simulation) to recompute |

Reading a file you don't already have open always lands it in buffer 1 and
becomes ACTIVE, pushing every other data buffer up one slot -- matching the
original, where buffer 1 is "whatever was read most recently," not a fixed
slot. Unlike the original, though, nothing ever falls off the end and gets
destroyed to make room: the buffer list just keeps growing. Re-`GET`ting a
file already open in some buffer just re-selects it in place, without
scrolling anything (`cmds.htm`'s documented `PLOT` behaviour). Buffer 0 has
**no simulate command**: it is recomputed whenever the sample or the active
buffer's parameters change, which is how RUMP behaved.

```
Your wish? get measured.rbs     /* becomes buffer 1, ACTIVE            */
Your wish? get another.rbs      /* becomes the new buffer 1; measured.rbs is now buffer 2 */
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
| `OFFSET <keV(0)>` `[new]` | calibration offset alone, independent of `CONVERSION` |
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
| `CMP` | synonym for `COMPARE`, at the RUMP, SIM, and PERT levels alike |
| `DISPLAY` | sample composition vs. depth (from the SIM description) |
| `REGION lo hi` | channel range shown |
| `EXPAND lo hi` | narrow the current region and redraw |
| `COUNTS lo [hi]` | yield range shown |
| `BLOWUP <max>` | shorthand for `COUNTS 0 <max>` |
| `LINEAR` / `SQRT` / `LOG` | yield axis scale |
| `NORMALIZE` / `RAW` | normalized vs. raw yield units |
| `LABELS [off]` | axis labels on or off |
| `STRUCTLABEL [off]` | show the SIM sample's layer structure (substrate first) as the simulation's legend text, instead of "SIM" |
| `COMPFRAC [off]` `[new]` | within `STRUCTLABEL`, show each layer's composition as atomic fraction (sums to 1) instead of raw stoichiometry |
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
these — see the [known-limitations note](../dev/about.md#milestones).

### Macros

`XEQ <file>` runs a file of commands through the same interpreter the prompt
uses, so an analysis can be checked in as a text file and replayed. `CALL` and
`EXECUTE` are synonyms. A bare name with no extension is tried as-is, then as
`.cmd`, then as `.rbs`/`.RBS` -- the last two `[new]` because some RBS
acquisition software writes its output as a plain `EMPTY`/`SWALLOW` command
macro under that extension (see [Buffers](#buffers)), despite the name
suggesting real spectrum data. `XEQ`ing an actual binary `.rbs` file (that
belongs with `GET`) fails with a clear message rather than a decode error.

`SCRIPT <file>` logs what you type into exactly such a file, and `SCRIPT OFF`
stops. `LOGFILE` and `RECORD` are synonyms.

**`~/.pyrumprc`** is a plain macro file — no different from anything `XEQ`
runs — that's read once at startup, before you're dropped into the prompt,
unless you pass `--norc`. It lives at `Path.home() / ".pyrumprc"`: that's
`~/.pyrumprc` on Linux/macOS and `C:\Users\<you>\.pyrumprc` on Windows, same
filename either way, nothing further to configure. A minimal one is just
ordinary commands, one per line:

```
$ cat ~/.pyrumprc
faithful off
mev 3.5
theta 5
region 300 800
```

> Note `SCRIPT`/`LOGFILE` need at least four characters, which is how the
> original kept them clear of `LOG` — the logarithmic yield axis. Typing `log`
> gets you the axis, `logf` the session log.

`faithful off` toggles the session between the shipped C's bug-for-bug
behaviour (the default) and the corrected physics available at that point in
the port — see [Design and validation](../physics/validation.md). It's a
session setting, not persisted on its own, so `~/.pyrumprc` is how you make
it a standing per-user default. `--faithful on`/`--faithful off` overrides it
for one invocation, applied after `~/.pyrumprc` runs but before any macro
passed on the command line — the macro can still set `FAITHFUL` itself if it
needs to.

`mev 3.5`/`theta 5` are **default experiment settings**. `MEV`/`THETA`/`PHI`/
`PSI`/`OMEGA`/`CHARGE`/`CURRENT`/`FWHM`/`CORRECTION`/`CHOFF`/`CONVERSION`/
`OFFSET`/`GEOMETRY`/`BEAM` all normally act on the ACTIVE buffer — but before
any `GET`, there is no active buffer, so they fall back to a session-wide
default instead of erroring. That makes it possible to explore a `SIM`
sample's theoretical spectrum (`PLOT 0`) with no real data loaded at all.
The same defaults also fill in for a freshly-read ASCII spectrum, which
carries no beam/geometry/detector metadata of its own — so `GET`ting one
picks up your defaults instead of the code's hardcoded 2.0 MeV. A `.RBS`
file's own metadata always wins, and once any real buffer becomes ACTIVE,
these commands go back to editing it, exactly as before — the defaults are
only a fallback, never a silent override of real data.

`region 300 800` works from `~/.pyrumprc` for a different reason: `REGION`
(and `SCALE`/`LABELS`/`ENERGY`/`COUNTS`/`BLOWUP`, the other plot-state
commands) write straight to session-wide state that was never gated on an
active buffer in the first place, so they've always been usable before any
`GET` — no code changes were needed to support them here. It also now
shapes `COMPARE`, not just `PLOT`/`OVERLAY` (see the 1.1.0 changelog entry).
As a rule of thumb for anything not listed above: if a command already
writes to session-wide state rather than a specific buffer, it works from
`~/.pyrumprc` for free; only a command that hard-requires an active buffer
needs the fallback that `MEV`/`THETA`/etc. got in 1.1.0.

