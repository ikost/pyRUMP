## Analysis tools

Element identification, calibration, and quantification, ported from RUMP's
`anlytc.c` command family. All of these act on the active buffer; region
arguments are plain 0-based channel indices, matching `INTEGRAL`'s existing
convention (not RUMP's own `first`-relative channel numbering).

#### `CURSOR`

```
Your wish? help cursor
  CURsor  read channel/energy/yield at the nearest data point
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
Your wish? help element
  ELement  expected energy/channel of an element's surface peak
  usage: ELEMENT el [el ...]
```

Prints the expected K, energy and channel of each element's surface edge.

#### `MATRIX`

```
Your wish? help matrix
  MATrix  expected energy, channel and matrix height
  usage: MATRIX el
```

Prints the expected energy, channel **and matrix height** for one element.

#### `WHATISIT`

```
Your wish? help whatisit
  WHATisit  identify elements near a channel
  usage: WHATISIT <channel>
```

Identifies the elements whose surface edge is nearest a channel.

#### `INFO`

```
Your wish? help info
  INFo  detailed report on an element
  usage: INFO el
```

Full report: density, K, cross section, stopping factors, isotopes.

#### `INTEGRAL`

```
Your wish? help integral
  INTegral  sum counts over a channel range
  usage: INTEGRAL lo hi
```

Gross/net counts over a channel range (background-corrected net).

#### `THICKNESS`

```
Your wish? help thickness
  THICkness  integral plus thickness conversion
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
Your wish? help background
  BACKground  fit and subtract a polynomial background
  usage: BACKGROUND lo1 hi1 lo2 hi2 order [-inplace] [-noplot]
```

Fits and strips a polynomial background.

#### `SMOOTH`

```
Your wish? help smooth
  SMOoth  smooth the active buffer (-sv, -conv, -fft)
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
Your wish? help fft
  FFT  FFT smooth (same as SMOOTH -FFT -RANGE)
  usage: FFT lo hi width
```

Same as `SMOOTH -fft -range lo hi width`.

#### `WIDTH_THICK`

```
Your wish? help width_thick
  WIDth_thick  thickness from a peak's half-height width
  usage: WIDTH_THICK ch1 ch2 element
```

Computes thickness from a peak's half-height width.

#### `CALIBRATE`

```
Your wish? help calibrate
  CALibrate  energy-calibrate from two known peaks
  usage: CALIBRATE ch1 el1 ch2 el2 [energy_eV marker_channel]
```

Sets keV/channel and keV(0) from two known peaks — see the `THICKNESS`
example above, which chains a `CALIBRATE` into it.

#### `INTSET`

```
Your wish? help intset
  INTSET  change INTEGRAL/THICKNESS rounding and alpha mode
  usage: INTSET [Round|Interp|Surface|Estimated|Query|?]
```

Picks two independent modes that both `INTEGRAL` and `THICKNESS` honor:
whether a region's boundaries are rounded to the nearest channel or
interpolated between them, and (for `THICKNESS` only) whether its second,
"compensated" pass uses an estimated alpha or asks you for one.
