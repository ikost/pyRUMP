## Plotting & display

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
