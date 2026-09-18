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

Command names and their **minimum abbreviations** as well as some useful synonyms for backward compatibility with original RUMP follow the original
(`REGion`, `OVerlay`), so `reg 100 400` and `region 100 400` are the same
command. `help` or `?` lists everything, with the required characters upper-cased. One
exception: `compare` requires its full name at every level, since a partial
abbreviation collided with `COMPOSITION` in SIM and PERT; `cmp` works
everywhere instead as an explicit synonym.

Commands tagged `[new]` below have no original-RUMP counterpart — see the
[Changelog](https://github.com/ikost/pyRUMP/releases)'s versioning note for what that means for.

## General system commands

####  `HELP` / `?` 

```
help [name]
```

Lists every command at the current level . `HELP <name>` (or `? <name>`)
describes one command instead. 

#### `QUIT` / `BYE` / 'q'

Leaves pyRUMP or goes one level up from `SIM` and `PERT`. It asks "Really quit pyRUMP? [y/N]" first; anything but `y`/`yes` cancels.


#### `LS` / `DIRECTORY` / `SL`

```
LS [pattern]
```

Lists files, optionally matching a glob (`ls *.rbs`). `SL` is a bare
synonym, not shown in `?`'s listing but still usable. Tab autocompletion and wildcards are implemented.

#### `LL`

```
LL [pattern]
```

Long listing: size and modification time, columned. Deliberately not Unix
permission bits — they carry no meaning on Windows, and this listing looks
the same on every platform. 

#### `CD` / `CHDIR`

```
CD [directory]
```

Changes directory; with no argument, goes home.
Tab autocompletion and wildcards are implemented.


#### `PUSHDIR` / `POPDIR`

```
PUSHDIR [dir]
```

`PUSHDIR <dir>` remembers the current directory before changing to `<dir>`
(or home, with no argument); `POPDIR` changes back to it.

#### `PWD` / `WHERE`

Prints the working directory.

####  `CAT` / `MORE` / `TYPE` 

```
CAT <file>
```

Shows a text file, paged a screenful at a time when there's an actual
terminal to pause for — under `XEQ`, `--batch`, or in tests, it prints
straight through instead of blocking for a keypress.

#### `CLS`

Clears the screen.

#### `XEQ` / `CALL` / `EXECUTE`

```
XEQ <file>
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
ECHO [off]
```

`ECHO [off]` toggles whether commands are echoed as they run — most useful
inside an `XEQ` macro, to see what's actually executing. `QUIET` is the
original's own off-synonym for `ECHO`.

#### `SCRIPT` / `LOGFILE` / `RECORD`

```
SCRIPT [file|off]
```

Logs what you type into a file, for later replay with `XEQ`; `SCRIPT OFF`
(or `LOGFILE OFF`/`RECORD OFF`) stops. Needs at least four characters, which
is how the original kept it clear of `LOG` — the logarithmic yield axis.
Typing `log` gets you the axis, `logf` the session log.

Standing per-user defaults (the `FAITHFUL` toggle, default experiment
settings, plot state) are set once via `~/.pyrumprc` rather than every
session — see [Config](config.md).

## Core workflow

### `GET` / `READ`

```
GET <file|n>
```

Reads a file into a buffer, or points at buffer *n*. `READ` is a synonym,
usable wherever `GET` is. Reading a file you don't already have open always lands it in buffer
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

### `XEQ`

Cross-referenced here from [General system commands](#xeq-call-execute):
`XEQ <file>` runs a command file, and loads data about as often as it runs a
macro, so it belongs alongside `GET` as much as alongside the filesystem
commands it's filed under.

### `SIM`

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

### `PERT`

Enters the fitting sub-processor: its own `PERT Command:` prompt, for
selecting what varies and running the least-squares search. Same one-shot
form as `SIM` — `PERT GO` re-runs the last selection without entering the
prompt. See [SIM and PERT](sim-pert.md).

### `COMPARE` / `CMP`

Active buffer vs. the simulation, with Poisson residuals and a reduced
chi-square readout (over `PERT`'s error windows if set, else the visible
`REGION`). `CMP` is a synonym, usable at the RUMP, SIM, and PERT levels
alike — see the abbreviation note above for why `COMPARE` itself needs its
full name.

### `PLOT`

```
 PLOT [buffer|file]
```

Erases and plots a buffer (default: active) or file. The plot is one
persistent matplotlib window whose state survives between commands — see
[Plotting & display](plotting.md) for everything that shapes it.

```
Your wish? plot 1               /* erase and plot buffer 1            */
```

### `RECALCULATE`

Forces buffer 0 (the simulation) to recompute. Rarely needed by hand —
`SIM`/`PERT` changes already trigger it — but useful after something that
doesn't, e.g. reloading the atomic tables with `DATA`.

### `RETURN` / 'Q'

Leaves `SIM` or `PERT` back to the RUMP level. Typing a command neither
sub-level recognizes does the same thing implicitly — it falls through to
RUMP and runs there — so `RETURN` is only needed to get back with nothing
else to run.
