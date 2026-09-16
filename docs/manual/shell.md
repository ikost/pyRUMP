## pyRUMP shell

Running `pyrump` with no arguments starts the pyRUMP shell — RUMP's own
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

The sections below follow `?`'s own grouping and order, most important
first — typing `?` at the RUMP prompt reproduces this same structure.

## General system commands

#### `?` / `HELP`

```
Your wish? help help
  HELP  list the commands
  usage: ? [name]
```

Lists every command at the current level, one section per table — e.g. at
the RUMP level, "Core workflow", "Settings", and so on, matching the
original's own grouping (`rump.c:351-360`). `HELP <name>` (or `? <name>`)
describes one command instead. A name this level doesn't know falls through
to whatever a sub-level would resolve it to, so `HELP THICK` works from the
RUMP prompt even though `THICK` is a `PERT` command.

#### `QUIT` / `BYE`

Leaves pyRUMP. Run interactively (not from an `XEQ` macro, which has no one
at the keyboard to answer), it asks "Really quit pyRUMP? [y/N]" first;
anything but `y`/`yes` cancels.

The rest of this section is a port of RUMP's own filesystem commands
(`lexp/system.c:175`), so you can move to your data rather than restarting
pyRUMP in the right directory. Reachable from every level — as in the
original, a command `SIM`/`PERT` doesn't know returns you to the RUMP level
and runs there.

Wildcards are expanded by the command itself, never by an OS shell, so `ls
*.rbs` behaves the same on Linux, macOS and Windows. Tab completion works on
both command names and paths.

There is deliberately **no shell escape** (the original's `!` / `DOS` /
`CSH`): it would let any `.cmd` macro run arbitrary commands on your
machine.

#### `LS` / `DIRECTORY` / `SL`

```
Your wish? help ls
  LS  list files, optionally matching a pattern
  usage: LS [pattern]
```

Lists files, optionally matching a glob (`ls *.rbs`). `SL` is a bare
synonym, not shown in `?`'s listing but still usable.

#### `LL`

```
Your wish? help ll
  LL  long listing with size and date
  usage: LL [pattern]
```

Long listing: size and modification time, columned. Deliberately not Unix
permission bits — they carry no meaning on Windows, and this listing looks
the same on every platform.

#### `CD` / `CHDIR`

```
Your wish? help cd
  CD  change directory (no argument: home)
  usage: CD [directory]
```

Changes directory; with no argument, goes home.

#### `PUSHDIR` / `POPDIR`

```
Your wish? help pushdir
  PUSHDir  change directory, remembering this one
  usage: PUSHDIR [dir]
```

`PUSHDIR <dir>` remembers the current directory before changing to `<dir>`
(or home, with no argument); `POPDIR` changes back to it.

#### `PWD` / `WHERE`

Prints the working directory.

#### `TYPE` / `CAT` / `MORE`

```
Your wish? help type
  TYpe  display a text file
  usage: TYPE <file>
```

Shows a text file, paged a screenful at a time when there's an actual
terminal to pause for — under `XEQ`, `--batch`, or in tests, it prints
straight through instead of blocking for a keypress.

#### `CLS`

Clears the screen.

#### `XEQ` / `CALL` / `EXECUTE`

```
Your wish? help xeq
  XEq  read RC43 .RBS file or execute a command file
  usage: XEQ <file>
```

Runs a file of commands through the same interpreter the prompt uses, so an
analysis can be checked in as a text file and replayed. A bare name with no
extension is tried as-is, then as `.cmd`, then as `.rbs`/`.RBS` -- the last
two `[new]` because some RBS acquisition software writes its output as a
plain `EMPTY`/`SWALLOW` command macro under that extension (see
[File formats](#file-formats)), despite the name suggesting real spectrum
data. `XEQ`ing an actual binary `.rbs` file (that belongs with `GET`) fails
with a clear message rather than a decode error.

```
Your wish? get measured.rbs     /* binary format -- reads records directly */
Your wish? xeq acquired.rbs     /* text macro -- replayed as commands      */
```

#### `ECHO` / `QUIET`

```
Your wish? help echo
  ECHO  echo commands as they run
  usage: ECHO [off]
```

`ECHO [off]` toggles whether commands are echoed as they run — most useful
inside an `XEQ` macro, to see what's actually executing. `QUIET` is the
original's own off-synonym for `ECHO`.

#### `SCRIPT` / `LOGFILE` / `RECORD`

```
Your wish? help script
  SCRIPT  record commands to a file for replay
  usage: SCRIPT [file|off]
```

Logs what you type into a file, for later replay with `XEQ`; `SCRIPT OFF`
(or `LOGFILE OFF`/`RECORD OFF`) stops. Needs at least four characters, which
is how the original kept it clear of `LOG` — the logarithmic yield axis.
Typing `log` gets you the axis, `logf` the session log.

Standing per-user defaults (the `FAITHFUL` toggle, default experiment
settings, plot state) are set once via `~/.pyrumprc` rather than every
session — see [Config](config.md).

### Core workflow

#### `GET` / `READ`

```
Your wish? help get
  GET  point at a buffer, or read a file into one
  usage: GET <file|n>
```

Reads a file into a buffer, or points at buffer *n* (`READ` reads a file
only). Reading a file you don't already have open always lands it in buffer
1 and becomes ACTIVE, pushing every other data buffer up one slot -- matching
the original, where buffer 1 is "whatever was read most recently," not a
fixed slot. Unlike the original, though, nothing ever falls off the end and
gets destroyed to make room: the buffer list just keeps growing. Re-`GET`ting
a file already open in some buffer just re-selects it in place, without
scrolling anything (`cmds.htm`'s documented `PLOT` behaviour). Buffer 0 has
**no simulate command**: it is recomputed whenever the sample or the active
buffer's parameters change, which is how RUMP behaved.

```
Your wish? get measured.rbs     /* becomes buffer 1, ACTIVE            */
Your wish? get another.rbs      /* becomes the new buffer 1; measured.rbs is now buffer 2 */
Your wish? get 0                /* point back at the simulation       */
Your wish? copy 0 2              /* snapshot the simulation into buffer 2 */
```

#### `XEQ`

Cross-referenced here from [General system commands](#xeq-call-execute):
`XEQ <file>` runs a command file, and loads data about as often as it runs a
macro, so it belongs alongside `GET` as much as alongside the filesystem
commands it's filed under.

#### `SIM`

Enters the sample-description editor: its own `SIM Command:` prompt, for
building or editing the layered target `PERT` fits against and `COMPARE`
plots. See [SIM and PERT](sim-pert.md).

`SIM <command>` runs one SIM command without leaving the RUMP level and
returning, e.g. `sim thick 500 A` — handy inside a one-line macro, or when
you only need to tweak one thing:

```
Your wish? sim thick 1 500 A
  layer 1 thickness = 500 A
```

#### `PERT`

Enters the fitting sub-processor: its own `PERT Command:` prompt, for
selecting what varies and running the least-squares search. Same one-shot
form as `SIM` — `PERT GO` re-runs the last selection without entering the
prompt. See [SIM and PERT](sim-pert.md).

#### `COMPARE` / `CMP`

Active buffer vs. the simulation, with Poisson residuals and a reduced
chi-square readout (over `PERT`'s error windows if set, else the visible
`REGION`). `CMP` is a synonym, usable at the RUMP, SIM, and PERT levels
alike — see the abbreviation note above for why `COMPARE` itself needs its
full name.

#### `PLOT`

```
Your wish? help plot
  PLot  erase and plot a buffer or file
  usage: PLOT [buffer|file]
```

Erases and plots a buffer (default: active) or file. The plot is one
persistent matplotlib window whose state survives between commands — see
[Plotting & display](#plotting-display) for everything that shapes it.

```
Your wish? plot 1               /* erase and plot buffer 1            */
```

#### `RECALCULATE`

Forces buffer 0 (the simulation) to recompute. Rarely needed by hand —
`SIM`/`PERT` changes already trigger it — but useful after something that
doesn't, e.g. reloading the atomic tables with `DATA`.

#### `RETURN`

Leaves `SIM` or `PERT` back to the RUMP level. Typing a command neither
sub-level recognizes does the same thing implicitly — it falls through to
RUMP and runs there — so `RETURN` is only needed to get back with nothing
else to run.

### Plotting & display

The plot is one persistent matplotlib window whose state survives between
commands. matplotlib installs by default with pyrump.

#### `OVERLAY`

```
Your wish? help overlay
  OVerlay  overlay a buffer or file on the current plot
  usage: OVERLAY [buffer|file]
```

Adds another trace to the current plot, without erasing it.

```
Your wish? overlay 0            /* add the simulation on top          */
```

#### `REPLOT`

Redraws the current plot, unchanged — useful after resizing the window or
after a setting that doesn't redraw on its own.

#### `FIGSAVE` / `HCOPY`

```
Your wish? help figsave
  FIGSave / HCOPY  save the current plot to an image file, e.g. FIGSAVE out.png
  usage: FIGSAVE <file>
```

Saves the current plot to an image file (`.png` by default; format follows
the extension). `HCOPY` is a synonym, the original's own name for the
command (a literal hard-copy to a plotter, in RUMP's day).

#### `REGION`

```
Your wish? help region
  REGion  set the channel range
  usage: REGION lo hi
```

Sets the channel range shown.

```
Your wish? region 100 400
```

#### `EXPAND`

```
Your wish? help expand
  EXpand  narrow the region and replot
  usage: EXPAND lo hi
```

Narrows the current region and redraws — a `REGION` that's relative to what's
already shown rather than absolute channel numbers.

#### `COUNTS`

```
Your wish? help counts
  COunts  set the yield range
  usage: COUNTS lo [hi]
```

Sets the yield range shown.

#### `LINEAR` / `SQRT` / `LOG`

Sets the yield axis scale. `sqrt` redraws immediately with the new scale:

```
Your wish? sqrt                 /* redraws immediately, sqrt yield    */
```

#### `NORMALIZE` / `RAW`

Toggles between normalized and raw yield units on the plot.

#### `LABELS`

```
Your wish? help labels
  LAbels  label the axes (LABELS OFF to suppress)
  usage: LABELS [off]
```

Turns axis labels on or off.

#### `STRUCTLABEL`

```
Your wish? help structlabel
  STRUctlabel  show the sample's layer structure as the simulation's legend text (STRUCTLABEL OFF for "SIM")
  usage: STRUCTLABEL [off]
```

Shows the `SIM` sample's layer structure (substrate first) as the
simulation's legend text, instead of "SIM".

#### `COMPFRAC` `[new]`

```
Your wish? help compfrac
  COMPfrac  show composition (STRUCTLABEL and SHOW) as atomic fraction, not raw stoichiometry
  usage: COMPFRAC [off]
```

Shows each layer's composition as atomic fraction (sums to 1) instead of raw
stoichiometry, in both `STRUCTLABEL` and SIM/PERT `SHOW`.

#### `ENERGY`

```
Your wish? help energy
  ENERgy  x axis in energy rather than channel
  usage: ENERGY [off]
```

Puts the x axis in energy (keV) rather than channel.

#### `AXIS`

Draws empty axes, with no data.

#### `BLOWUP`

```
Your wish? help blowup
  BLowup  expand the vertical scale
  usage: BLOWUP <max>
```

Shorthand for `COUNTS 0 <max>`.

#### `PARAMETERS` / `PARMS`

Prints the current plot settings.

#### `DISPLAY`

Plots the sample composition against depth, from the `SIM` description.

### Buffers

Spectra live in numbered buffers, one of which is ACTIVE and is what most
commands act on implicitly. Buffer **0 is the simulation**; data starts at
1. See [`GET`](#get-read) in Core workflow for how a freshly
read file gets slotted in.

#### `BUFFERS`

Lists the buffers, marking the active one.

#### `READ`

```
Your wish? help read
  REad  read a file into a buffer
  usage: GET <file|n>
```

Cross-referenced from [Core workflow](#get-read): `READ <file>`
reads a file into a buffer, same as `GET` except it never accepts a bare
buffer number.

#### `POINTAT`

```
Your wish? help pointat
  POintat  point at a buffer by number
  usage: POINTAT <n>
```

Points at buffer *n*, by number only — unlike `GET`, which also accepts a
filename.

#### `RELEASE` / `NEWALL`

```
Your wish? help release
  RELEASE  release the active buffer
  usage: RELEASE [n]
```

`RELEASE` drops one buffer (default: active); `NEWALL` drops all of them.

#### `EMPTY`

```
Your wish? help empty
  EMPty  reset a buffer to blank, or open a new one
  usage: EMPTY [n]
```

Scrolls a fresh blank buffer into buffer 1 (default), or resets buffer *n*
in place.

#### `COPY` / `MOVE`

```
Your wish? help copy
  COPY  copy one buffer to another
  usage: COPY <source> <target>

Your wish? help move
  MOVE  exchange two buffers
  usage: MOVE <a> <b>
```

`COPY` duplicates a buffer; `MOVE` exchanges two.

```
Your wish? copy 0 2              /* snapshot the simulation into buffer 2 */
```

#### `WRITE` / `WRASCII`

```
Your wish? help write
  WRITE  write the active buffer to a .rbs file
  usage: WRITE <file>

Your wish? help wrascii
  WRAscii  write the active buffer as text
  usage: WRASCII <file>
```

Saves the active buffer, binary or text.

### Sample & instrument parameters

Each buffer carries its own beam, geometry, calibration and measurement
metadata. Every one of these **prints the current value with no argument,
and sets it (echoing the new value) with one** — and chains onto any
further command left on the line, so `Choff 0 FWHM 15` works in one go,
exactly as RUMP's own `WRASCII` output writes it back.

```
Your wish? beam 4He++
  beam Z=2 mass=4.0026 charge state 2
Your wish? mev 2.0
  MeV = 2
Your wish? conversion 5.0 0
  5 keV/channel, offset 0 keV
```

#### `ACTIVE`

Prints the active buffer's full parameter set.

#### `BEAM`

```
Your wish? help beam
  BEAM  incident beam species, e.g. 4He++
  usage: BEAM 4He++
```

Sets the beam species and charge state.

#### `MEV`

```
Your wish? help mev
  MEV  beam energy
  usage: MEV <energy>
```

Sets the beam energy, MeV.

#### `THETA`

```
Your wish? help theta
  THEta  sample tilt
  usage: THETA <deg>
```

Sets the sample tilt.

#### `PHI`

```
Your wish? help phi
  PHI  supplement of the scattering angle
  usage: PHI <deg>
```

Sets 180° minus the scattering angle.

#### `PSI`

```
Your wish? help psi
  PSI  exit angle
  usage: PSI <deg>
```

Sets the exit angle (GENERAL geometry only).

#### `GEOMETRY`

```
Your wish? help geometry
  GEOMetry  cornell, ibm or general
  usage: GEOMETRY cornell|ibm|general
```

Sets the detector geometry convention.

#### `CONVERSION`

```
Your wish? help conversion
  CONVersion  keV per channel and offset
  usage: CONVERSION <keV/ch> [keV(0)]
```

Sets the energy calibration: slope and, optionally, offset together.

#### `SLOPE` `[new]`

```
Your wish? help slope
  SLOpe  keV/channel alone, independent of CONVERSION's offset
  usage: SLOPE <keV/ch>
```

Sets the calibration slope alone, independent of `CONVERSION`'s offset.

#### `OFFSET` `[new]`

```
Your wish? help offset
  OFFset  keV(0) alone, independent of CONVERSION's keV/ch
  usage: OFFSET <keV(0)>
```

Sets the calibration offset alone, independent of `CONVERSION`'s slope.

#### `CORRECTION`

```
Your wish? help correction
  CORrection  normalization fudge factor
  usage: CORRECTION <factor>
```

Sets the normalization fudge factor that absorbs charge-integration error
(RUMP's `CORR`). `PERT`'s `NORMALIZE` window can set this automatically from
a fit instead of you setting it by hand — see
[PERT commands](sim-pert.md#pert-commands).

#### `CHARGE`

```
Your wish? help charge
  CHarge  beam dose
  usage: CHARGE <uC>
```

Sets the integrated beam dose.

#### `CURRENT`

```
Your wish? help current
  CURRent  average beam current, for pileup
  usage: CURRENT <nA>
```

Sets the average beam current — enables pile-up modeling together with
`TAU`.

#### `CHOFF`

```
Your wish? help choff
  CHOff  channel number of the first data point
  usage: CHOFF <n>
```

Sets the channel number of the first data point.

#### `FWHM`

```
Your wish? help fwhm
  FWHM  detector resolution
  usage: FWHM <keV>
```

Sets the detector resolution.

#### `OMEGA`

```
Your wish? help omega
  OMEGA  detector solid angle
  usage: OMEGA <msr>
```

Sets the detector solid angle.

#### `TAU`

```
Your wish? help tau
  TAU  MCA shaping time
  usage: TAU <us>
```

Sets the MCA shaping time constant.

#### `IDENTIFIER`

```
Your wish? help identifier
  IDEntifier  description of the spectrum
  usage: IDENTIFIER <text>
```

Sets a free-text spectrum description.

#### `DATE`

```
Your wish? help date
  DATE  when the spectrum was measured
  usage: DATE <text>
```

Sets when the spectrum was measured.

#### `FILENAME`

```
Your wish? help filename
  FILEname  record the buffer's source filename
  usage: FILENAME <name>
```

Sets the recorded source filename.

#### `SWALLOW`

```
Your wish? help swallow
  SWALLOW  read the following macro lines as channel data
  usage: SWALLOW [-twocolumn]
```

Used inside an `XEQ` macro, reads the macro file's following lines straight
into the active buffer as channel data (or channel/value pairs), stopping at
the first blank line — how a RUMP-written `.cmd` file reconstructs a
spectrum inline.

### Analysis tools

Element identification, calibration, and quantification, ported from RUMP's
`anlytc.c` command family. All of these act on the active buffer; region
arguments are plain 0-based channel indices, matching `INTEGRAL`'s existing
convention (not RUMP's own `first`-relative channel numbering).

#### `CURSOR`

```
Your wish? help cursor
  CURsor  read channel/energy/yield at the nearest data point
  usage: CURSOR <channel>
```

Reads the channel, energy and yield at the nearest sampled channel in the
active buffer, printing which buffer that is up front (index and name) so
that's never a guess. RUMP's own `CURSOR` read whatever point a physical
graphics crosshair sat on; pyRUMP has no such device and no reason to fake
one with mouse clicks — matplotlib's own toolbar already gives a live x/y
readout for free while hovering a plot, so this takes the channel directly
instead and snaps to the nearest real sample, the same "read what's
actually there" the original's cursor served.

#### `ELEMENT`

```
Your wish? help element
  ELement  expected energy/channel of an element's surface peak
  usage: ELEMENT el [el ...]
```

Prints the expected K, energy and channel of each element's surface edge.

#### `MATRIX`

```
Your wish? help matrix
  MATrix  expected energy, channel and matrix height
  usage: MATRIX el
```

Prints the expected energy, channel **and matrix height** for one element.

#### `WHATISIT`

```
Your wish? help whatisit
  WHATisit  identify elements near a channel
  usage: WHATISIT <channel>
```

Identifies the elements whose surface edge is nearest a channel.

#### `INFO`

```
Your wish? help info
  INFo  detailed report on an element
  usage: INFO el
```

Full report: density, K, cross section, stopping factors, isotopes.

#### `INTEGRAL`

```
Your wish? help integral
  INTegral  sum counts over a channel range
  usage: INTEGRAL lo hi
```

Gross/net counts over a channel range (background-corrected net).

#### `THICKNESS`

```
Your wish? help thickness
  THICkness  integral plus thickness conversion
  usage: THICKNESS lo hi element
```

`INTEGRAL` plus conversion to atoms/cm² and Angstroms.

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

#### `BACKGROUND`

```
Your wish? help background
  BACKground  fit and subtract a polynomial background
  usage: BACKGROUND lo1 hi1 lo2 hi2 order [-inplace] [-noplot]
```

Fits and strips a polynomial background.

#### `SMOOTH`

```
Your wish? help smooth
  SMOoth  smooth the active buffer (-sv, -conv, -fft)
  usage: SMOOTH [-sv|-conv|-fft] [-range lo hi] [n]
```

Smooths the active buffer. `-conv`'s characteristic width uses RUMP's own
(nonstandard) `sigma = (FWHM/2)/sqrt(ln 2)/kevch` — not the usual
`FWHM/(2*sqrt(2 ln 2))` — reproduced deliberately, not corrected. The
default range is the whole buffer; `-conv`'s iteration count and `-fft`'s
width both come from a trailing number, interpreted according to whichever
mode is active.

#### `FFT`

```
Your wish? help fft
  FFT  FFT smooth (same as SMOOTH -FFT -RANGE)
  usage: FFT lo hi width
```

Same as `SMOOTH -fft -range lo hi width`.

#### `WIDTH_THICK`

```
Your wish? help width_thick
  WIDth_thick  thickness from a peak's half-height width
  usage: WIDTH_THICK ch1 ch2 element
```

Computes thickness from a peak's half-height width.

#### `CALIBRATE`

```
Your wish? help calibrate
  CALibrate  energy-calibrate from two known peaks
  usage: CALIBRATE ch1 el1 ch2 el2 [energy_eV marker_channel]
```

Sets keV/channel and keV(0) from two known peaks — see the `THICKNESS`
example above, which chains a `CALIBRATE` into it.

#### `INTSET`

```
Your wish? help intset
  INTSET  change INTEGRAL/THICKNESS rounding and alpha mode
  usage: INTSET [Round|Interp|Surface|Estimated|Query|?]
```

Picks two independent modes that both `INTEGRAL` and `THICKNESS` honor:
whether a region's boundaries are rounded to the nearest channel or
interpolated between them, and (for `THICKNESS` only) whether its second,
"compensated" pass uses an estimated alpha or asks you for one.

### Settings

#### `DATA`

```
Your wish? help data
  DATA  show or change the atomic data directory
  usage: DATA [dir]
```

With no argument, prints the directory the atomic tables (`atom4.dat`,
`pscoef.dat`, stopping-power tables, …) were loaded from. With one, reloads
every table from that directory instead and forces the simulation buffer to
recompute — for comparing two table sets without restarting pyRUMP.

#### `FAITHFUL` `[new]`

```
Your wish? help faithful
  FAIThful  toggle faithful (bug-for-bug) vs corrected physics (FAITHFUL OFF to correct)
  usage: FAITHFUL [on|off]
```

Toggles the session between the shipped C's bug-for-bug behaviour (the
default) and pyRUMP's corrected physics where the two diverge — see
[Design and validation](../physics/validation.md) for what "corrected"
covers. With no argument, reports the current state instead of changing it.
Has no original-RUMP counterpart, hence `[new]`.

#### `SCREENING` `[new]`

```
Your wish? help screening
  SCREening  select the Rutherford screening correction: NONE, LECUYER (default) or ANDERSEN
  usage: SCREENING [none|lecuyer|andersen]
```

Selects the Rutherford screening correction. `LECUYER` is RUMP's own and
the default. `ANDERSEN` is a pyRUMP addition (Andersen et al., Phys. Rev. A
21 (1980) 1891) — more accurate at forward angles, but with no RUMP oracle
to validate it against, and per its own literature, may be inaccurate below
a few hundred keV or at small scattering angles.

#### `MODE` `[new]`

```
Your wish? help mode
  MODE  SIM/PERT thickness convention, COMP or ATOMS -- see MODE with no argument
  usage: MODE [Comp|Atoms]
```

Sets how `SIM`/`PERT` describe a layer's thickness. `Comp` is RUMP's own
convention: a physical thickness (normally Angstroms) split across elements
by stoichiometric ratio — `SIM`/`PERT`'s `THICKNESS` and `COMPOSITION`.
`Atoms` instead holds each element's own areal density directly in
`COMPOSITION`, with `THICKNESS` just their sum in `/CM2` — `SIM`/`PERT`'s
`ATOMS`. Only one set of commands is usable at a time, gated by this
setting, so a fit can never mix the two conventions on the same layer.
Switching recalculates every layer between the two conventions through each
layer's own atomic density — not a relabelling — so the simulated spectrum
is unchanged either way. With no argument, prints the current mode.

#### `PROFILE`

Not implemented — never was, even in the original (`RbsNewprf`, dead code in
the shipped C). Reproduced verbatim: it prints "OOPS: Didn't think anyone
used this routine anymore - sorry not implemented" and does nothing.

### File formats

Two unrelated things can carry a `.rbs`/`.RBS` extension. `GET` and `XEQ`
each expect one of them, and loading the wrong one with the wrong command
fails with a message naming the mistake, rather than a raw decode error or
silently-wrong data.

**RUMP binary format** (`GET`/`WRITE`) is a TIFF-inspired stream of
self-describing, big-endian 32-bit-word records:

```
[length_in_words] [record_type] [data...] [checksum]
```

`length` counts the header and checksum words themselves, so a record holds
`length - 3` data words, at most 1024. Every word in a record, checksum
included, sums to zero modulo 2^32 -- the format's own integrity check,
verified on every read. Records carry the program identifier, comments, the
sample identifier, date and live/real time, the dead-time correction,
accelerator parameters (E0, Zbeam, Mbeam, charge state, charge, current), MCA
calibration (keV/channel, keV at channel 0, first channel, FWHM), RBS
geometry (kind, theta, phi, psi, solid angle) and finally the counts, in one
of four encodings: raw floats, raw 32-bit integers, differential integers, or
differential integers plus zero-run compression. The spec
(`html/RUMP/rbs_inf.htm`, Aug 1994) is explicitly open and ships a reference
C API, unlike the rest of the legacy distribution -- see
[`pyrump/io/rbs.py`](https://github.com/ikost/pyRUMP/blob/main/src/pyrump/io/rbs.py)
for the clean-room implementation that reads and writes it.

**RC43 acquisition macro** (`XEQ`) is what some RBS acquisition software (NEC's
RC43 among them) writes instead, saved under a `.RBS` extension despite the
name suggesting real binary spectrum data. It's an ordinary text macro --
`IDENTIFIER`, `DATE`, `CONVERSION`, `MEV`, `BEAM`, `GEOMETRY`,
`THETA`/`PHI`/`PSI`, `OMEGA`, `CHOFF`, `FWHM`, `CURRENT`, `CHARGE` -- ending
in `SWALLOW`, which reads everything after it straight into the buffer as
channel data (see [Sample & instrument parameters](#sample-instrument-parameters)
for `SWALLOW` itself). It's the same mechanism a RUMP-written `.cmd` file
uses to reconstruct a spectrum inline; RC43 just also uses it under a `.RBS`
name.

```
Your wish? get measured.rbs     /* binary format -- reads records directly */
Your wish? xeq acquired.rbs     /* text macro -- replayed as commands      */
```

See [Quick start](../getting-started/index.md#data-loading) for both formats
worked end to end, including what each of `GET`/`XEQ` prints when pointed at
the other one's file by mistake.
