## SIM — the sample editor

`SIM` edits the sample description: the stack of layers that buffer 0 (the
simulation) is computed from, and that [`PERT`](pert.md) fits. It is not a
"simulate" command — buffer 0 recomputes itself whenever the sample changes.

`SIM` opens its own `SIM Command:` prompt. A command SIM does not recognise
is passed to the RUMP level, which also takes you back there. From the RUMP
level, `SIM <command>` runs a single SIM command without entering the prompt:

```
Your wish? sim get MnPt.lcm
Your wish? sim show
```

Layer 1 is the surface; higher numbers are deeper. Every command acts on the
current layer, which `SHOW` marks with `>`.

```
Your wish? sim
SIM Command: thick 500 A              /* layer 1, at the surface */
SIM Command: composition Si 1 /
SIM Command: next                     /* layer 2, below it       */
SIM Command: thick 2000 A
SIM Command: composition Au 1 /
SIM Command: show
SIM Command: save mysample.lcm
SIM Command: return
```

SIM's commands are the same code that reads `.lcm` files (see
[Sample descriptions](samples.md)), so what you type and what the file holds
cannot drift apart.

### Getting around

#### `HELP` / `?`

```
usage: HELP [<name>]
```

Lists the SIM commands by group. `HELP <name>` describes one command, with
its usage; a name SIM does not know is looked up at the RUMP level.

```
SIM Command: help thick
```

#### `RETURN` / `QUIT` / `Q` / `ABORT`

Goes back to the RUMP level. Inside SIM, `QUIT` does not exit pyRUMP.

#### `SHOW`

Lists the layers, surface first, then the sample-wide settings. The brackets
show each layer in the other [`MODE`](config.md#mode-new) — exactly what
switching `MODE` would turn it into. In MODE COMP that is each element's
areal density in 10¹⁵ at/cm²:

```
SIM Command: show
 >  1            40 A        Ru 1           [29 /CM2  Ru 29.03]
    2           331 A        Mn 2.73 Pt 1   [254 /CM2  Mn 186.11 Pt 68.05]
    3            40 A        Ru 1           [29 /CM2  Ru 29.03]
    4           330 A        Si 1 O 2       [149 /CM2  Si 49.57 O 99.13]
    5          1000 nm       Si 1           [4978 /CM2  Si 4977.69]
  maxpth 1000   straggle 0   multiple 0   absorber 0
```

In MODE ATOMS it is the thickness in Angstroms and the composition as
fractions of 1:

```
SIM Command: show
 >  1            29 /CM2     Ru 29.03             [40 A  Ru 1]
    2           254 /CM2     Mn 186.11 Pt 68.05   [331 A  Mn 0.732 Pt 0.268]
    ...
```

With `COMPFRAC` on, MODE COMP compositions show as atomic fractions (see
[Plotting & display](plotting.md)). MODE ATOMS always shows each element's
10¹⁵ at/cm², since the fractions are already in the brackets.

#### `STATUS`

Summarises the sample: layer count, elements, sample-wide settings, and
whether the simulation is up to date.

```
SIM Command: status
  2 layers, maxpth 1000, straggle 0, multiple 0
  elements: Si Au
  simulation is stale
```

### Layers

#### `LAYER`

```
usage: LAYER [<n>]
```

Moves to layer *n*. One past the last layer gives a fresh blank layer to
fill in. With no argument, prints which layer you are on.

```
SIM Command: layer 2
  you are now working on layer # 2 of 5
```

#### `NEXT`

Moves one layer deeper. Past the last layer you get a fresh blank layer. A
layer left without a thickness disappears when you move away from it.

```
SIM Command: next
  you are now working on a fresh layer # 6 (of 5)
```

#### `OPEN` / `INSERT`

Inserts a blank layer above the current one and moves to it. Layers below it
are renumbered, and SIM warns if a `PERT` selection refers to one of them.

```
SIM Command: layer 1
SIM Command: open                     /* new layer 1; the old one is now 2 */
SIM Command: thick 20 A
SIM Command: comp C 1 /
```

#### `CLOSE` / `DELETE` / `CLEAR`

```
usage: CLOSE [<n>]
```

Removes the current layer, or layer *n*. Layers below it are renumbered,
with the same `PERT` warning as `OPEN`.

```
SIM Command: close 3
  layer 3 deleted
```

#### `RESET`

Empties the sample.

```
SIM Command: reset
sample reset to empty space
```

### Layer contents — MODE COMP

In MODE COMP (the default, and RUMP's own convention), a layer is a thickness
plus a stoichiometry. [`MODE`](config.md#mode-new) switches to MODE ATOMS,
right here in SIM; switching converts every layer, and the spectrum stays
the same.

#### `THICKNESS`

```
usage: THICKNESS <value> [<unit>]
```

Sets the layer's thickness. The unit is `A` unless given; `nm`, `um`,
`/CM2` (10¹⁵ at/cm²), `M/CM2` or a compound name from `density.tab` also
work (see [`DENSITY`](#density)).

```
SIM Command: thick 500 A
SIM Command: thick 250 /CM2
```

#### `COMPOSITION`

```
usage: COMPOSITION <element> <n> [<element> <n> ...] /
```

Sets the layer's stoichiometry. Only the ratios matter; `/` ends the list.
Element symbols can be typed in any case (`si`, `SI`); they are stored and
shown as `Si`, here and in `SPECIES`, `ATOMS` and `.lcm` files.

```
SIM Command: comp In 2 O 3 /
```

### Layer contents — MODE ATOMS

#### `ATOMS` `[new]`

```
usage: ATOMS <element> <amount> [<element> <amount> ...] /
```

Sets each element's own areal density, in 10¹⁵ at/cm². The layer's
thickness becomes their sum, in `/CM2`. Needs [`MODE ATOMS`](config.md#mode-new).

```
SIM Command: atoms Mn 185 Pt 65.5 /   /* 250.5 /CM2 in total */
```

### Profiles and interfaces

#### `SPECIES`

```
usage: SPECIES <element> <n> [<element> <n> ...] /
```

Sets what an `EQUATION` blends the layer toward.

```
SIM Command: species Au 1 /
```

#### `EQUATION`

```
usage: EQUATION <name> <p1> [<p2> ...]
```

Makes the layer's composition vary with depth: a mix of its `COMPOSITION`
and its `SPECIES`, weighted by the equation. See
[Depth profiles](samples.md#depth-profiles) for the equations and their
parameters.

```
SIM Command: species Au 1 /
SIM Command: equation linear 0 0.2    /* 0% Au at the front, 20% at the back */
```

#### `EQLIST` `[new]`

Lists the equation names `EQUATION` accepts.

```
SIM Command: eqlist
```

#### `FUZZ`

```
usage: FUZZ <amount> [<steps>]
```

Roughens the layer: simulates `<steps>` variants of its thickness, spread by
`<amount>` (in the layer's own unit), and sums them with Gaussian weights.
Each fuzzed layer multiplies the simulation time by `<steps>`. `FUZZ 0`
turns it off. The original took the layer number as its first argument;
here it acts on the current layer.

```
SIM Command: fuzz 30 5                /* +/- about 30 A, 5 variants */
```

#### `SUBLAYER`

```
usage: SUBLAYER <n>
```

Splits the layer into *n* internal slabs for the calculation. Takes
precedence over `STHICKNESS` and `MAXPTH`.

```
SIM Command: sublayer 10
```

#### `STHICKNESS`

```
usage: STHICKNESS <value> [<unit>]
```

Sets the thickness of each internal slab in this layer instead.

```
SIM Command: sthick 50 A
```

### Sample-wide

#### `MAXPTH`

```
usage: MAXPTH [<value>]
```

Maximum internal slab thickness, in 10¹⁵ at/cm², for layers without
`SUBLAYER` or `STHICKNESS`. Smaller is more accurate and slower. With no
argument, prints the current value.

```
SIM Command: maxpth 200
  maxpth 200
```

#### `STRAGGLE`

```
usage: STRAGGLE <factor>
```

Energy-loss straggling, as a multiple of Bohr's value. `0` (the default)
turns it off.

```
SIM Command: straggle 1
```

#### `ABSORBER`

```
usage: ABSORBER <n>
```

Treats the first *n* layers as a foil between the sample and the detector (a
dead layer, window or stopper foil): particles cross them only on the way
out, and they do not tilt with the sample. `0` (the default) means none.

```
SIM Command: absorber 1               /* layer 1 is the detector window */
```

#### `MULTIPLE`

```
usage: MULTIPLE <strength>
```

Adds an ad-hoc low-energy tail; it does not model real multiple scattering
(see [Physics](physics.md)). `0` (the default) turns it off. The original
called it `MULTIPLE_SCATTER`.

```
SIM Command: multiple 1
```

### Files and plotting

#### `GET`

```
usage: GET <file>
```

Reads a `.lcm` sample description, replacing the current sample. `.lcm` is
added when the name has no extension. Commands pyRUMP does not know are
skipped with a warning.

```
SIM Command: get MnPt
read MnPt.lcm: 5 layers
```

#### `SAVE`

```
usage: SAVE <file>
```

Writes the sample in RUMP's own `.lcm` format; `.lcm` is added when the name
has no extension.

```
SIM Command: save mysample
wrote mysample.lcm
```

#### `DENSITY`

```
usage: DENSITY [<pattern>]
```

Lists the thickness units, and the compounds in `density.tab` that also work
as a `THICKNESS` unit, with their atomic density in 10²³ at/cm³. A pattern
narrows the compound list. A compound that is not in the table is not an
error: it silently gets silicon's density.

```
SIM Command: density SiO
```

#### `SPLOT`

```
usage: SPLOT [<element> | <layer>]
```

Overlays the simulation on the current plot. With an element, only that
element's contribution; with a layer number, only that layer's.

```
SIM Command: splot Mn                 /* manganese only */
SIM Command: splot 2                  /* layer 2 only   */
```

#### `COMPARE` / `CMP`

Plots the active buffer against the simulation, with residuals, the same as
[`COMPARE`](shell.md#compare-cmp) at the RUMP level.

```
SIM Command: cmp
```
