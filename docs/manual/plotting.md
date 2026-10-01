## Plotting & display

The plot is one persistent matplotlib window whose state survives between
commands. matplotlib installs by default with pyrump.

#### `OVERLAY`

```
usage: OVERLAY [buffer|file]
```

Adds another trace to the current plot, without erasing it — on a
`COMPARE`, to its top panel.

```
Your wish? overlay 0            /* add the simulation on top          */
```

#### `REPLOT`

Redraws the current plot, unchanged — useful after resizing the window or
after a setting that doesn't redraw on its own.

#### `FIGSAVE` / `HCOPY`

```
usage: FIGSAVE <file>
```

Saves the current plot to an image file (`.png` by default; format follows
the extension). Raster formats such as `.png` are written at 300 dpi;
`.pdf` and `.svg` are vector and scale freely. `HCOPY` is a synonym, the
original's own name for the command (a literal hard-copy to a plotter, in
RUMP's day).

#### `REGION`

```
usage: REGION lo hi
```

Sets the channel range shown.

```
Your wish? region 100 400
```

#### `EXPAND`

```
usage: EXPAND lo hi
```

Narrows the current region and redraws — a `REGION` that's relative to what's
already shown rather than absolute channel numbers.

#### `COUNTS`

```
usage: COUNTS lo [hi]
```

Sets the yield range shown.

#### `LINEAR` / `SQRT` / `LOG`

Sets the yield axis scale. `sqrt` redraws immediately with the new scale:

```
Your wish? sqrt                 /* redraws immediately, sqrt yield    */
```

On a `COMPARE` the scale applies to the top panel; the residuals stay linear.
`LOG` ignores a `COUNTS` floor of 0 or below, which a log axis can't show.

#### `NORMALIZE` / `RAW`

Toggles between normalized and raw yield units on the plot. It applies to
`COMPARE` too, whose residuals and chi-square are then taken on the
normalized yield rather than raw counts — the readout says so.

#### `LABELS`

```
usage: LABELS [off]
```

Turns axis labels on or off.

#### `STRUCTLABEL`

```
usage: STRUCTLABEL [off]
```

Shows the `SIM` sample's layer structure (substrate first) as the
simulation's legend text, instead of "SIM".

#### `COMPFRAC` `[new]`

```
usage: COMPFRAC [off]
```

Shows each layer's composition as atomic fraction (sums to 1) instead of raw
stoichiometry, in both `STRUCTLABEL` and SIM/PERT `SHOW`. `SHOW` in MODE
ATOMS is the exception: it keeps each element's 10¹⁵ at/cm², with the
fractions in brackets.

#### `ENERGY`

```
usage: ENERGY [off]
```

Puts the x axis in energy (keV) rather than channel.

#### `AXIS`

Draws empty axes, with no data.

#### `BLOWUP`

```
usage: BLOWUP <max>
```

Shorthand for `COUNTS 0 <max>`.

#### `PARAMETERS` / `PARMS`

Prints the current plot settings.

#### `DISPLAY`

Plots the sample composition against depth, from the `SIM` description.
