## PERT — fitting

`PERT` fits the [`SIM`](sim.md) sample to the active data buffer by least
squares. You choose which parameters may vary and which channels count, then
run `GO`. Fitted values are written back into the sample, so `SIM SHOW` and
`SIM SAVE` show them.

`PERT` opens its own `PERT Command:` prompt. A command PERT does not
recognise is passed to the RUMP level, which also takes you back there. From
the RUMP level, `PERT <command>` runs a single PERT command without entering
the prompt:

```
Your wish? pert get MnPt.pert go      /* load a saved setup and fit */
```

```
Your wish? pert
PERT Command: window 800 1200         /* fit these channels only        */
PERT Command: thickness 2             /* vary layer 2's thickness       */
PERT Command: composition 2 Mn        /* and its Mn content             */
PERT Command: go
```

Two differences from the original: the data may be in any buffer, not just
buffer 1, and `MULTI` is the default.

### Getting around

#### `HELP` / `?`

```
usage: HELP [<name>]
```

Lists the PERT commands. `HELP <name>` describes one command, with its
usage; a name PERT does not know is looked up at the RUMP level.

```
PERT Command: help thickness
```

#### `RETURN` / `QUIT` / `Q`

Goes back to the RUMP level. Inside PERT, `QUIT` does not exit pyRUMP.

#### `PARMS`

Shows the current setup: fit mode, windows, and what is varying.

```
PERT Command: parms
  mode        multiple variable
  autocmp     off
  autosnap    off
  error win   [1] 800-1200
  norm win    (none)
  varying:
    [1] kev(0)
    [2] layer 2 thickness
    [3] layer 2 composition Mn
    [4] fwhm  bounds 18-25
```

#### `SHOW`

Lists the sample, the same as [`SIM SHOW`](sim.md#show). Useful to check
layer numbers before selecting parameters.

```
PERT Command: show
```

### Running a fit

#### `GO`

Runs the fit and prints the result (see [Reading the fit report](#reading-the-fit-report)).
Needs a sample, a data buffer, and at least one varying parameter. Refuses to
run if a varying parameter does not match the current
[`MODE`](config.md#mode-new) (switching `MODE` normally drops those).

```
PERT Command: go
```

#### `SINGLE` / `MULTI`

`MULTI` (the default) fits all varying parameters together. `SINGLE` fits
them one at a time, in the order they were selected.

```
PERT Command: single
  single-variable mode
```

#### `VOLUME`

```
usage: VOLUME [off]
```

Prints a line for every simulation `GO` runs, so a long fit visibly makes
progress. Off by default.

```
PERT Command: volume
  messages on
```

#### `AUTOCMP` `[new]`

```
usage: AUTOCMP [off]
```

Runs `COMPARE` at the end of every `GO`. Off by default; `GET` and `CLEAR`
keep the setting.

```
PERT Command: autocmp
  autocmp on
```

#### `AUTOSNAP` `[new]`

```
usage: AUTOSNAP [off]
```

After every `GO`, takes a [`SNAPSHOT`](shell.md#snapshot-snap-new) of the
fit, writing the same files `SNAP` does, named after the sample: the first word of the data buffer's
`IDENTIFIER`, or its file name. If another dataset already has that name,
the data file's own name is used instead, with a warning (see
[`SNAPSHOT`](shell.md#snapshot-snap-new)). Off by default; `GET` and `CLEAR`
keep the setting.

* `<sample>.report` — the fit result, appended, so refits keep a history
* `<sample>.pert` — the PERT setup, as `SAVE` writes it
* `<sample>.lcm` — the fitted sample, as `SIM SAVE` writes it
* `<sample>_rbs.png` — the `COMPARE` plot, as `FIGSAVE` writes it (with
  `<sample>_pixe.png` when it saves the PIXE window too)
* `<sample>_fit.xeq` — the restore macro: `pyrump <sample>_fit.xeq` is this
  session again, with the fitted spectrum parameters too (`CORRECTION`,
  calibration …)

Tuning by hand after the fit isn't covered: type `SNAP` again to save it. Its
report block is marked *set by hand*.

```
PERT Command: autosnap
  autosnap on
```

`REPORT`, its old name, still works, so older macros and `~/.pyrumprc`
files keep working, but it is no longer listed by `HELP`.

#### `HIGHLIGHT` `[new]`

```
usage: HIGHLIGHT [off]
```

Shades the fit's windows on the plots, so it's plain which channels `GO`
fits. In the RBS window, the [`WINDOW`](#window) channels are light blue and
the [`NORMALIZE`](#normalize) window light orange. In the PIXE window, the
[`PIXWIN`](#pixwin-new) channels are light blue. It works in `PLOT` and
`COMPARE` alike, in both panels of a comparison, and on the channel or the
energy axis. Adding or clearing a window redraws the plot at once. On by
default; `HIGHLIGHT OFF` gives clean plots for figures. `GET` and `CLEAR`
keep the setting, and `SAVE` doesn't write it. A snapshot keeps it off if
you turned it off.

```
PERT Command: highlight off
  highlight off
```

#### `COMPARE` / `CMP`

Plots the active buffer against the simulation, with residuals, the same as
[`COMPARE`](shell.md#compare-cmp) at the RUMP level.

```
PERT Command: cmp
```

[`EXPORTCMP` / `EC`](shell.md#exportcmp-ec-new) writes the same comparison
to a column file from here too, without leaving PERT, and
[`EXPORT`](shell.md#export-new) writes the active buffer on its own.

### Saving and clearing

#### `GET`

```
usage: GET <file> [GO]
```

Replays a saved `.pert` setup, replacing the current one. A trailing `GO`
runs the fit right after. `.pert` is added when the name has no extension.

```
PERT Command: get MnPt go
```

#### `SAVE`

```
usage: SAVE <file>
```

Writes the current setup (windows, fit mode, varying parameters and their
bounds) to a `.pert` file, as the PERT commands that recreate it.

```
PERT Command: save usual
wrote usual.pert
```

#### `CLEAR` `[new]`

```
usage: CLEAR [<n>]
```

Forgets the whole setup, or only varying parameter *n*, numbered as `PARMS`
lists them.

```
PERT Command: clear 3                 /* stop varying parameter [3] */
```

### Windows

#### `WINDOW`

```
usage: WINDOW <first> <last>
usage: WINDOW CLEAR [<n>]
```

Adds an error window: only these channels count in the fit. Up to 10
windows; with none, the whole spectrum counts. The plots shade them in
light blue ([`HIGHLIGHT`](#highlight-new)). `WINDOW CLEAR` removes all
windows, `WINDOW CLEAR <n>` only window *n*. With no argument, shows the
setup like `PARMS`.

```
PERT Command: window 800 1200
  error windows [1] 800-1200
```

#### `NORMALIZE`

```
usage: NORMALIZE <first> <last>
usage: NORMALIZE CLEAR
```

Scales the data to the simulation over these channels, to absorb a charge or
dose error. After `GO`, the scale is stored in the buffer's
[`CORRECTION`](buffers.md#correction). Cannot be combined with varying
`CORRECTION`. The plots shade it in light orange
([`HIGHLIGHT`](#highlight-new)).

```
PERT Command: normalize 300 400
  normalisation window 300-400
```

The scale stands for the dose, so with [`PIXWIN`](#pixwin-new) windows it is
applied to the PIXE spectrum too.

#### `PIXWIN` `[new]`

```
usage: PIXWIN <first> <last>
usage: PIXWIN CLEAR [<n>]
```

Sets the PIXE channels that the parameters fitted to the PIXE spectrum are
fitted over. Channels are numbered as the `.PIX` file numbers them, as
PIXE `REGION` takes them. This is for alloys whose elements overlap in RBS
but not in their X-ray lines: Ta–W, Fe–Ni, Ni–Co.

**You choose, element by element, which spectrum fits it.** `COMPOSITION`,
`ATOMS` and `SPECIES` take `PIXE` after the element
([`COMPOSITION 2 Fe PIXE`](#composition)). Without it, or with `RBS`, the
element is fitted to the RBS spectrum, as in RUMP. Everything else —
thicknesses, calibration, beam — is always fitted to RBS. GO then
alternates:

1. the `PIXE` parameters, and H with [`PIXH`](#pixh-new), fitted to the PIXE
   spectrum over the `PIXWIN` channels, everything else held;
2. the rest fitted to the RBS spectrum over the [`WINDOW`](#window)
   channels, the PIXE parameters held;

and repeats both until no value moves by more than 0.1 % (at most 10
rounds), since the two depend on each other a little: the PIXE yields on
the thicknesses, the RBS spectrum on the composition. Each parameter is
fitted to one spectrum only, so the two never pull against each other.
GO reports each spectrum's chi-square, and how many rounds it took.

* **Put the windows on clear peaks.** The PIXE continuum isn't simulated, so
  a window over background would be fitted to nothing. One window may span a
  whole group of lines. Overlaps such as Fe Kβ on Ni Kα, or Ta Lβ near W Lα,
  are fine, because the simulation has every line.
* **Leave an element unvaried as the reference.** In MODE COMP only the
  ratios count: fit `Ta` to PIXE against a fixed `W`, and the thickness to
  RBS. GO refuses every element of a layer fitted to PIXE in MODE COMP.
* **A normalisation window** ([`NORMALIZE`](#normalize)) sets the dose for
  the PIXE data too.

GO refuses `PIXE` parameters without `PIXWIN` windows or without a PIXE
spectrum in the buffer (`PIXE GET`, or `PAIR ON`), and an element with no
simulated counts in the windows — O in an oxide, say — since the PIXE
spectrum can say nothing about it.

Up to 10 windows. `PIXWIN CLEAR` removes all of them, `PIXWIN CLEAR <n>`
only window *n*, and `PIXWIN` alone lists them. They're saved with `SAVE`,
so snapshots keep them too. The PIXE window shades them
([`HIGHLIGHT`](#highlight-new)), and its `COMPARE` scores the chi-square
over them. `PARMS` shows the spectrum each parameter is fitted to, and H:

```
PERT Command: window 380 500
PERT Command: pixwin 786 994        /* Ta L and W L lines, 7.9-10 keV */
PERT Command: composition 1 Ta pixe
PERT Command: thickness 1
PERT Command: pixh
PERT Command: parms
  mode        multiple variable
  autocmp     off
  autosnap    off
  highlight   on
  error win   [1] 380-500
  norm win    (none)
  PIXE win    [1] 786-994
  PIXE H      K 1   L 1.25   M 1   L fitted (Ta L lines in the PIXE windows)
  varying:
    [1] layer 1 composition Ta       PIXE
    [2] layer 1 thickness            RBS
    [3] PIXH                         PIXE
PERT Command: go
  Fitting film: Si [20000/cm2] - WTa2 [1000/cm2]
  PIXE, channels 786-994: layer 1 composition Ta, PIXH L
  RBS, channels 380-500:  layer 1 thickness

  fit took 0.33 s, 2 rounds PIXE -> RBS

  PIXE: reduced chi-square 1.1153 on 207 dof   (was 10.7757)
  RBS:  reduced chi-square 1.5234 on 120 dof   (was 1.5276)
  40 evaluations, converged
  layer 1 composition Ta            1.01684  +/- 0.02061   (was 2)
  layer 1 thickness               1000 /cm2  +/- 0.2016 /cm2   (was 1000 /cm2)
  PIXH L                           0.990994  +/- 0.008716   (was 1.25)

  film: Si [20000/cm2] - WTa1.02 [1000/cm2]
```

That is a simulated W–Ta film (truth: Ta 1, H 1) started at Ta 2 with H_L
25% off. From RBS alone the same fit gives Ta to ±0.057; the PIXE lines
narrow it to ±0.021. The whole analysis of a measured Fe–Ni film, from
loading the data to saving the session, is the
[FeNi worked example](../getting-started/feni-rbs-pixe.md).

`PIXWIN` and `PIXH` need all four letters, so that `PIXE` (or `PIX`) still
opens the PIXE prompt from PERT, as it does at every other prompt.

### Parameters

Each command below adds one parameter to vary; selecting the same one again
replaces it. Layer numbers are the ones `SHOW` lists.

Every parameter takes an optional `<min> <max>` bound at the end: give both
or neither. The fit stays inside the bound. `PARMS` shows it, and
`SAVE`/`GET` keep it. A `THICKNESS` bound is in the layer's own unit (the one
`SHOW` lists, usually Å); an `ATOMS` bound is in 10¹⁵ at/cm².

#### `THICKNESS`

```
usage: THICKNESS <layer> [<min> <max>]
```

Varies a layer's thickness. Needs [MODE COMP](config.md#mode-new). The
bound is in the layer's own unit, e.g. Å for a layer `SHOW` lists in `A`.

```
PERT Command: thickness 2 300 360
```

#### `COMPOSITION`

```
usage: COMPOSITION <layer> <element> [PIXE|RBS] [<min> <max>]
```

Varies one element's stoichiometry in a layer. The element must already be
in that layer. Needs [MODE COMP](config.md#mode-new). `PIXE` fits it to the
PIXE spectrum over the [`PIXWIN`](#pixwin-new) channels instead of the RBS
one.

```
PERT Command: composition 2 Mn 2 4
```

!!! warning "Always keep one element fixed"
    Only the ratios within a layer matter, so varying every element of a
    layer together has no effect on the spectrum and the fit fails. In an
    *N*-element layer, vary at most *N*−1 of them.

#### `ATOMS` `[new]`

```
usage: ATOMS <layer> <element> [PIXE|RBS] [<min> <max>]
```

Varies one element's own areal density in a layer, keeping the other
elements' amounts fixed, in 10¹⁵ at/cm². Needs [MODE ATOMS](config.md#mode-new).
`PIXE` fits it to the PIXE spectrum over the [`PIXWIN`](#pixwin-new)
channels. With H calibrated and held (no [`PIXH`](#pixh-new)), every
element of a layer can be fitted to PIXE this way, each amount absolute.

```
PERT Command: atoms 2 Mn
```

#### `SPECIES`

```
usage: SPECIES <layer> <element> [PIXE|RBS] [<min> <max>]
```

Varies one element of a layer's [`SPECIES`](sim.md#species). The element
must already be declared there. `PIXE` fits it to the PIXE spectrum.

```
PERT Command: species 2 Au
```

#### `EQUATION`

```
usage: EQUATION <layer> <n> [<min> <max>]
```

Varies parameter *n* of a layer's [`EQUATION`](sim.md#equation).

```
PERT Command: equation 2 2            /* the second equation parameter */
```

#### `STRAGGLE`

```
usage: STRAGGLE [<min> <max>]
```

Varies the sample's [straggling factor](sim.md#straggle).

```
PERT Command: straggle 0.5 2
```

#### `FUZZ`

Not implemented; it prints an error.

#### `MEV`

```
usage: MEV [<min> <max>]
```

Varies the beam energy, in MeV.

```
PERT Command: mev
```

#### `FWHM`

```
usage: FWHM [<min> <max>]
```

Varies the detector resolution, in keV.

```
PERT Command: fwhm 18 25
```

#### `THETA`

```
usage: THETA [<min> <max>]
```

Varies the sample tilt, in degrees. With `GEOMETRY IBM` or `CORNELL` the exit
angle follows it; with `GENERAL` it stays at `PSI` (see [Experimental
geometry](geometry.md)).

```
PERT Command: theta -12 -6
```

#### `CORRECTION`

```
usage: CORRECTION [<min> <max>]
```

Varies the buffer's normalisation factor ([`CORRECTION`](buffers.md#correction)).
Cannot be combined with a `NORMALIZE` window.

```
PERT Command: correction
```

#### `SLOPE` / `KEV/CH`

```
usage: SLOPE [<min> <max>]
```

Varies the energy calibration's slope, in keV per channel.

```
PERT Command: slope
```

#### `OFFSET` / `KEV(0)`

```
usage: OFFSET [<min> <max>]
```

Varies the energy calibration's offset, in keV, e.g. to follow a
sample-charging shift.

```
PERT Command: offset
```

#### `PIXH` `[new]`

```
usage: PIXH [<min> <max>]
```

Varies the PIXE instrumental constant H together with the parameters
fitted to the PIXE spectrum. GO varies H of each shell (K, L or M) whose
lines of those elements fall in the [`PIXWIN`](#pixwin-new) channels: Fe
and Ni K lines in the windows, H_K. The PIXE spectrum then gives only the
ratios of the elements, which is all it needs to. The fitted H is written
back to the PIXE setup (PIXE `H`).

Without `PIXH`, H stays at the PIXE prompt's values, as calibrated on a
standard: the PIXE spectrum then gives absolute amounts. `PARMS` says
which H is fitted, and why:

```
PERT Command: pixh
  varying PIXH
PERT Command: parms
  ...
  PIXE H      K 1   L 1.25   M 1   L fitted (Ta L lines in the PIXE windows)
```

H scales the counts of every element in its shell, so some element with
lines there must not be fitted to PIXE: its amount, from RBS, sets the
scale. GO refuses `PIXH` when every element with lines of a fitted shell
in the windows is fitted to PIXE (`ATOMS 1 W PIXE` and `ATOMS 1 Ta PIXE`
together, say).

### Reading the fit report

```
Your wish? pert get MnPt.pert go
  ...
  Fitting MA8408: Si [5um] - SiO2 [300A] - Ru [38A] - Mn2.83Pt [291A] - Ru [40A]

  fit took 0.60 s

  reduced chi-square 4.4259 on 792 dof   (was 165.1066)
  87 evaluations, converged
  fwhm                              19.8701  +/- 0.1221   (was 15)
  correction                       0.859777  +/- 0.001692   (was 1)
  kev(0)                            42.5962  +/- 0.1588   (was 48)
  layer 1 thickness                    40 A  +/- 0.3466 A   (was 40 A)   [29 /CM2]
  layer 2 thickness                   317 A  +/- 0.7872 A   (was 291 A)   [244 /CM2]
  layer 3 thickness                    41 A  +/- 0.3302 A   (was 38 A)   [30 /CM2]
  layer 2 composition Mn             2.7605  +/- 0.01587   (was 2.82931)
  kev/ch                            1.60215  +/- 0.0001141   (was 1.59)
  layer 4 thickness                   532 A  +/- 10.05 A   (was 300 A)   [240 /CM2]

  MA8408: Si [5um] - SiO2 [532A] - Ru [41A] - Mn2.76Pt [317A] - Ru [40A]
```

* **`Fitting MA8408: …`** — the sample name and the layer structure before the
  fit, substrate first. The last line shows the structure after it.
* **`reduced chi-square … on … dof`** — chi-square per degree of freedom,
  over the error windows; `(was …)` is the value before the fit. Near 1, the
  model explains the data as well as counting statistics allow. Much more
  than 1 means a systematic mismatch: a missing layer, a wrong calibration or
  geometry, or non-Rutherford scattering in the window. Much less than 1
  usually means the window is too narrow.
* **`… evaluations, converged`** — how many simulations the fit ran. `did not
  converge` means it stopped at its limit first: check the starting values
  or narrow the window before trusting the numbers.
* **One line per parameter** — the fitted value, its uncertainty, and the
  starting value. As in [`SIM SHOW`](sim.md#show), the other
  [`MODE`](config.md#mode-new)'s view follows in brackets. A `THICKNESS` is in
  the layer's own unit, with its 10¹⁵ at/cm² (`/CM2`) in brackets. An `ATOMS`
  amount is in `/CM2`, with the layer's resulting thickness in Å in
  brackets.
* **`data scaled by …`** — only with a `NORMALIZE` window: the scale applied,
  and the `CORRECTION` it was stored as.
