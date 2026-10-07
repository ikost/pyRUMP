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

Redraws the current plot — useful after resizing the window or after a
setting that doesn't redraw on its own. Any simulation on it is brought up
to date with the sample first, even with `LIVE OFF`.

#### `LIVE` `[new]`

```
usage: LIVE [off]
```

Keeps the simulation on the plot current. After every command that changes
the sample (any `SIM` edit, `SIM GET`, a `PERT` fit) or something the
simulation depends on (`MEV`, `THETA`, `FAITHFUL`, `SCREENING`, ...), the
simulation curves already on the plot are recomputed and the window
redrawn. This covers `PLOT 0`, `OVERLAY 0`, the simulation in a `COMPARE`
(with its residuals and chi-square) and each `SIM SPLOT` element or layer.
A `SPLOT` whose element or layer is no longer in the sample is dropped,
with a note. On by default.

* It **never opens a window**: with nothing plotted, or the window closed,
  it does nothing.
* An `XEQ` macro counts as one command, so a macro of twenty `SIM` edits
  recomputes once, at the end.
* `PERT`'s `GO` is not slowed down: the plot is redrawn once, with the
  fitted sample, after the fit.
* A sample that can't be simulated (all its layers deleted, say) is
  reported once, and the plot is left as it was.
* With [PIXE](pixe.md) enabled, the PIXE window follows too — it is one
  switch for both windows, and an open PIXE window follows even with the
  RBS window closed. PIXE-only settings (detector, `FILTER`, `CALIB`, `H`,
  ...) redraw the PIXE window themselves, with `LIVE` on or off, and never
  touch the RBS plot.

Each redraw costs a fresh simulation (one more for each `SPLOT` curve, and
the PIXE one), so for a slow sample turn it off with `LIVE OFF` and redraw
with `REPLOT` when you want to see it. `LIVE` works at any prompt without
leaving it. It is a standing preference, like `AUTOCMP`: put `LIVE OFF` in
`~/.pyrumprc` to keep it off; `SAVE` does not write it.

```
SIM Command: live off
live off
SIM Command: thick 60 A         /* the plot stays as it is            */
SIM Command: live               /* ...and catches up now              */
live on
```

#### `FIGSAVE` / `HCOPY`

```
usage: FIGSAVE <file>
```

Saves the current plot to an image file (`.png` by default; format follows
the extension). Raster formats such as `.png` are written at 300 dpi;
`.pdf` and `.svg` are vector and scale freely. `HCOPY` is a synonym, the
original's own name for the command (a literal hard-copy to a plotter, in
RUMP's day).

`[new]` The RBS plot is saved as `<file>_rbs`: `FIGSAVE fit` writes
`fit_rbs.png`. A name that already ends in `_rbs` isn't doubled. With
`PAIR ON` and the PIXE window open, the PIXE plot is saved next to it as
`fit_pixe.png` (see [PIXE](pixe.md)).

```
Your wish? figsave fit
wrote fit_rbs.png
```

#### `REGION`

```
usage: REGION lo hi
```

Sets the channel range shown. Any range is accepted, so `REGION` zooms in
and back out. With no plot open, the range is stored and used by the next
`PLOT`. With no arguments it prints the plot settings, the same as
`PARAMETERS`.

```
Your wish? region 100 400
```

#### `EXPAND`

```
usage: EXPAND lo hi
```

Zooms in on part of the current plot and redraws it. The channel numbers
are absolute, the same as `REGION`'s, but they must lie inside the range
already shown: `EXPAND` can only narrow the view. Use `REGION` to widen it
again. Limits given high-to-low are swapped. `EXPAND` needs a plot to be
open. If it is rejected, the plot keeps its current range.

`REGION` does everything `EXPAND` does and can also zoom out, so it is the
one to use in most cases. `EXPAND` is kept for compatibility with RUMP
macros. In the original it was an interactive tool: with no numbers typed,
you picked the two limits with the cursor on the plot. pyRUMP has no cursor
picking, so the channels must always be typed.

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

#### `ENERGY` `[new]`

```
usage: ENERGY [off]
```

Puts the x axis in energy (keV) rather than channel, as a single axis along
the bottom.

By default (and after `ENERGY OFF`) the plot uses RUMP's frame: channels along
the bottom and energy in MeV along the top, both on the scale of the first
buffer plotted. `COMPARE` shows the energy axis above its upper panel. RUMP
itself has no `ENERGY` command; it always draws both axes.

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

```
usage: DISPLAY [<depth> [<unit>]]
```

Plots the sample composition against depth, from the `SIM` description:
each element's atomic fraction, layer by layer, in the RBS window. Depth
is in 10¹⁵ at/cm², the unit RBS measures.

With no depth, the whole sample is shown, substrate included. A thick
substrate then takes up nearly all the axis, so give a depth to see the
films: `DISPLAY` stops there. The depth takes `SIM THICKNESS`'s units,
`A` unless given (`nm`, `um`, `/CM2`). A length is converted layer by
layer with each layer's own density, the one its thickness was given
with (a compound unit's, such as `SIO2`, or its elements'), and the
plot says what it came to:

```
Your wish? display 600
  depth 0-600 A = 0-422.2 /CM2 of the sample's 25256.4
```

RUMP has `DISPLAY` too, but there it does nothing: its `SimDrawSample`
is an empty stub.
