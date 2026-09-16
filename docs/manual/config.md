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
