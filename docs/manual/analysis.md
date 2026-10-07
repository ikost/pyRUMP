## Analysis tools

Element identification, calibration, and quantification, ported from RUMP's
`anlytc.c` command family. All of these act on the active buffer; region
arguments are plain 0-based channel indices, matching `INTEGRAL`'s existing
convention (not RUMP's own `first`-relative channel numbering).

#### `CURSOR`

```
usage: CURSOR <channel>
```

Reads the channel, energy and yield at the nearest sampled channel in the
active buffer, printing which buffer that is up front (index and name) so
that's never a guess. RUMP's own `CURSOR` read whatever point a physical
graphics crosshair sat on; pyRUMP has no such device and no reason to fake
one with mouse clicks — matplotlib's own toolbar already gives a live x/y
readout for free while hovering a plot, so this takes the channel directly
instead and snaps to the nearest real sample, the same "read what's
actually there" the original's cursor served.

#### `ELEMENT`

```
usage: ELEMENT el [el ...]
```

Prints the expected K, energy and channel of each element's surface edge.
If a plot is showing, each edge is also ticked along the bottom of the plot
and labelled with its symbol (`^{A}X` for an isotope such as `28Si`), the
same kind of marker WHATISIT draws, so peaks can be labelled on the figure.
Edges outside the current region are skipped.

With [PIXE](pixe.md) on, each element's main X-ray lines follow its RBS
line: their energy, and their channel as the PIXE plot's channel axis
reads it (the active buffer's PIXE calibration). Lines outside the PIXE
spectrum are left out. If a PIXE plot is showing, the lines are labelled
on it as well, the way `MARKERS` labels the sample's, even with
`MARKERS OFF`:

```
Your wish? element Pt Fe
  Pt  Z=78  Mass=195.090  K(ion)=0.9219  Energy=  1751.6 keV  Channel=1076.478
      X-rays  Lα 9.435 keV ch 941.5   Lβ1 11.071 keV ch 1104.5   Lβ2 11.248 keV ch 1122.1
              Lγ1 12.942 keV ch 1290.7   Mα 2.050 keV ch 206.4   Mβ 2.127 keV ch 214.1
  Fe  Z=26  Mass= 55.847  K(ion)=0.7524  Energy=  1429.5 keV  Channel= 873.888
      X-rays  Kα 6.400 keV ch 639.4   Kβ 7.058 keV ch 705.0   Lα 0.705 keV ch 72.5
              Lβ1 0.718 keV ch 73.7
```

Both plots' marks last until the plot is next drawn.

#### `MATRIX`

```
usage: MATRIX el
```

Prints the expected energy, channel **and matrix height** for one element.

#### `WHATISIT`

```
usage: WHATISIT <channel>
```

Identifies the elements whose surface edge is nearest a channel.

#### `INFO`

```
usage: INFO el
```

Full report: density, K, cross section, stopping factors, isotopes.

#### `INTEGRAL`

```
usage: INTEGRAL lo hi
```

Gross/net counts over a channel range (background-corrected net).

#### `THICKNESS`

```
usage: THICKNESS lo hi element
```

`INTEGRAL` plus conversion to atoms/cm² and Angstroms.

```
Your wish? calibrate 226 Si 369 Au
 Energy=2.0000 MeV    Conversion:4.9896 keV/ch   3.8329 keV(0)

Your wish? intset estimated
Your wish? thickness 180 220 Si
 Discrete integration on buffer 1
 Region:  180.0 to  220.0  Gross:     3498.31  Net:      -47.34  (#/uC/msr)
 Si surface approximation, density  2.32 g/cc
  (Gross)  2.2598e+18 Atoms/cm**2 ( 4539.9 Angstroms)
  ( Net ) -3.0582e+16 Atoms/cm**2 (  -61.4 Angstroms)
 Compensated calculation (Chu et al. page 65)
  (Gross)  2.0903e+18 Atoms/cm**2 ( 4199.4 Angstroms)
  ( Net ) -3.7411e+16 Atoms/cm**2 (  -75.2 Angstroms)
```

The negative "Net" values above aren't a bug: `180`-`220` sits on a flat part
of this (simulated, noiseless) spectrum, and the discrete net-background
correction assumes a sloped continuum either side of the region it's
integrating — pick regions either side of a real peak, not the middle of a
plateau, for a meaningful net figure.

Non-Rutherford (tabulated-resonance) cross sections aren't wired into any of
these — see the [known-limitations note](../dev/about.md#milestones).

#### `BACKGROUND`

```
usage: BACKGROUND lo1 hi1 lo2 hi2 order [-inplace] [-noplot]
```

Fits and strips a polynomial background.

#### `SMOOTH`

```
usage: SMOOTH [-sv|-conv|-fft] [-range lo hi] [n]
```

Smooths the active buffer. `-conv`'s characteristic width uses RUMP's own
(nonstandard) `sigma = (FWHM/2)/sqrt(ln 2)/kevch` — not the usual
`FWHM/(2*sqrt(2 ln 2))` — reproduced deliberately, not corrected. The
default range is the whole buffer; `-conv`'s iteration count and `-fft`'s
width both come from a trailing number, interpreted according to whichever
mode is active.

#### `FFT`

```
usage: FFT lo hi width
```

Same as `SMOOTH -fft -range lo hi width`.

#### `WIDTH_THICK`

```
usage: WIDTH_THICK ch1 ch2 element
```

Computes thickness from a peak's half-height width.

#### `CALIBRATE`

```
usage: CALIBRATE ch1 el1 ch2 el2 [energy_eV marker_channel]
```

Sets keV/channel and keV(0) from two known peaks — see the `THICKNESS`
example above, which chains a `CALIBRATE` into it.

#### `INTSET`

```
usage: INTSET [Round|Interp|Surface|Estimated|Query|?]
```

Picks two independent modes that both `INTEGRAL` and `THICKNESS` honor:
whether a region's boundaries are rounded to the nearest channel or
interpolated between them, and (for `THICKNESS` only) whether its second,
"compensated" pass uses an estimated alpha or asks you for one.
