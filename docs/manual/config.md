## Config

### `~/.pyrumprc`

A plain macro file — no different from anything `XEQ` runs — that's read
once at startup, before you're dropped into the prompt, unless you pass
`--norc`. It lives at `Path.home() / ".pyrumprc"`: that's `~/.pyrumprc` on
Linux/macOS and `C:\Users\<you>\.pyrumprc` on Windows, same filename either
way, nothing further to configure. A minimal one is just ordinary commands,
one per line:

```
$ cat ~/.pyrumprc
faithful off
mev 3.5
theta 5
region 300 800
```

`faithful off` toggles the session between the shipped C's bug-for-bug
behaviour (the default) and the corrected physics available at that point in
the port — see [Design and validation](../dev/validation.md). It's a
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
[Design and validation](../dev/validation.md) for what "corrected"
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
