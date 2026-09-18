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
both command names and paths, backed by the standard library `readline`
module. Windows has no built-in `readline`, so pyRUMP pulls in
[pyreadline3](https://pypi.org/project/pyreadline3/) there automatically as
a dependency; without it (e.g. an old install predating this), tab silently
does nothing.

There is deliberately **no shell escape** (the original's `!` / `DOS` /
`CSH`): it would let any `.cmd` macro run arbitrary commands on your
machine.

#### `LS` / `DIRECTORY` / `SL`

```
usage: LS [pattern]
```

Lists files, optionally matching a glob (`ls *.rbs`). `SL` is a bare
synonym, not shown in `?`'s listing but still usable.

#### `LL`

```
usage: LL [pattern]
```

Long listing: size and modification time, columned. Deliberately not Unix
permission bits — they carry no meaning on Windows, and this listing looks
the same on every platform.

#### `CD` / `CHDIR`

```
usage: CD [directory]
```

Changes directory; with no argument, goes home.

#### `PUSHDIR` / `POPDIR`

```
usage: PUSHDIR [dir]
```

`PUSHDIR <dir>` remembers the current directory before changing to `<dir>`
(or home, with no argument); `POPDIR` changes back to it.

#### `PWD` / `WHERE`

Prints the working directory.

#### `TYPE` / `CAT` / `MORE`

```
usage: TYPE <file>
```

Shows a text file, paged a screenful at a time when there's an actual
terminal to pause for — under `XEQ`, `--batch`, or in tests, it prints
straight through instead of blocking for a keypress.

#### `CLS`

Clears the screen.

#### `XEQ` / `CALL` / `EXECUTE`

```
usage: XEQ <file>
```

Runs a file of commands through the same interpreter the prompt uses, so an
analysis can be checked in as a text file and replayed. A bare name with no
extension is tried as-is, then as `.cmd`, then as `.rbs`/`.RBS` -- the last
two `[new]` because some RBS acquisition software writes its output as a
plain `EMPTY`/`SWALLOW` command macro under that extension (see
[File formats](file-formats.md)), despite the name suggesting real spectrum
data. `XEQ`ing an actual binary `.rbs` file (that belongs with `GET`) fails
with a clear message rather than a decode error.

```
Your wish? get measured.rbs     /* binary format -- reads records directly */
Your wish? xeq acquired.rbs     /* text macro -- replayed as commands      */
```

#### `ECHO` / `QUIET`

```
usage: ECHO [off]
```

`ECHO [off]` toggles whether commands are echoed as they run — most useful
inside an `XEQ` macro, to see what's actually executing. `QUIET` is the
original's own off-synonym for `ECHO`.

#### `SCRIPT` / `LOGFILE` / `RECORD`

```
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
usage: PLOT [buffer|file]
```

Erases and plots a buffer (default: active) or file. The plot is one
persistent matplotlib window whose state survives between commands — see
[Plotting & display](plotting.md) for everything that shapes it.

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

Plotting & display commands (`OVERLAY`, `REPLOT`, `FIGSAVE`, `REGION`,
`EXPAND`, `COUNTS`, `LINEAR`/`SQRT`/`LOG`, `NORMALIZE`/`RAW`, `LABELS`,
`STRUCTLABEL`, `COMPFRAC`, `ENERGY`, `AXIS`, `BLOWUP`, `PARAMETERS`,
`DISPLAY`) have moved to [Plotting & display](plotting.md).

Buffers (`BUFFERS`, `READ`, `POINTAT`, `RELEASE`/`NEWALL`, `EMPTY`, `COPY`/
`MOVE`, `WRITE`/`WRASCII`) and sample & instrument parameters (`ACTIVE`,
`BEAM`, `MEV`, `THETA`, `PHI`, `PSI`, `GEOMETRY`, `CONVERSION`, `SLOPE`,
`OFFSET`, `CORRECTION`, `CHARGE`, `CURRENT`, `CHOFF`, `FWHM`, `OMEGA`, `TAU`,
`IDENTIFIER`, `DATE`, `FILENAME`, `SWALLOW`) have moved to
[Buffers & instrument parameters](buffers.md).

Analysis tools (`CURSOR`, `ELEMENT`, `MATRIX`, `WHATISIT`, `INFO`,
`INTEGRAL`, `THICKNESS`, `BACKGROUND`, `SMOOTH`, `FFT`, `WIDTH_THICK`,
`CALIBRATE`, `INTSET`) have moved to [Analysis tools](analysis.md).

Settings (`DATA`, `FAITHFUL`, `SCREENING`, `MODE`, `PROFILE`) have moved to
[Config](config.md#settings).
