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

