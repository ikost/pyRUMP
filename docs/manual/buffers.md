## Buffers & instrument parameters

### Buffers

Spectra live in numbered buffers, one of which is ACTIVE and is what most
commands act on implicitly. Buffer **0 is the simulation**; data starts at
1. See [`GET`](shell.md#get-read) in Core workflow for how a freshly
read file gets slotted in.

#### `BUFFERS`

Lists the buffers, marking the active one.

#### `POINTAT`

```
usage: POINTAT <n>
```

Points at buffer *n*, by number only — unlike `GET`, which also accepts a
filename.

#### `RELEASE` / `NEWALL`

```
usage: RELEASE [n]
```

`RELEASE` drops one buffer (default: active); `NEWALL` drops all of them.

#### `EMPTY`

```
usage: EMPTY [n]
```

Scrolls a fresh blank buffer into buffer 1 (default), or resets buffer *n*
in place.

#### `COPY` / `MOVE`

```
usage: COPY <source> <target>
usage: MOVE <a> <b>
```

`COPY` duplicates a buffer; `MOVE` exchanges two.

```
Your wish? copy 0 2              /* snapshot the simulation into buffer 2 */
```

#### `WRITE` / `WRASCII`

```
usage: WRITE <file>
usage: WRASCII <file>
```

Saves the active buffer. `WRITE` writes the full RUMP binary `.rbs`
format -- identifier, comments, date, dead-time correction, accelerator
parameters, calibration, geometry, and counts -- the same self-describing
format `GET` reads (see [File formats](file-formats.md)). It always writes
compression mode 0 (unpacked floats) at format version 1.0, RUMP's own
default; raising it to 1.1 is gated behind `CONFIG WRITE_LEVEL` in the
original, which pyRUMP doesn't currently expose as a command.

`WRASCII` writes RUMP's own plain-text dialect instead (bmanip.c:595-650): a
keyword header -- identifier, date, charge/energy, calibration, angles,
detector solid angle/correction, channel offset/FWHM, current, geometry and
beam code -- then the literal line `Swallow`, then one count per line. It's
an ordinary text file, editable like the RC43 macros `XEQ` reads (see
[File formats](file-formats.md)), and matches the real RUMP binary's output
byte-for-byte for the same buffer. That said, the header is only for a human
reading the file: RUMP's own plain-ASCII reader (what `GET` falls back to for
anything that isn't the binary format) only ever picks up the first
non-numeric line as the identifier and ignores the rest, so `GET`ting a
`WRASCII` file back -- in pyRUMP or the original -- still falls back to your
`~/.pyrumprc` defaults for beam/geometry/calibration, not the header's own
values. `[new]` A bare filename with no extension gets `.dat` added
automatically; `.rbs`/`.RBS` is a footgun rather than a hard error --
`WRASCII` still writes the file but warns, since `GET` would otherwise
mistake it for an RC43 macro and refuse it (see [File formats](file-formats.md)).

`[new]` To open a spectrum in Origin, Excel or gnuplot instead, use
[`EXPORT`](shell.md#export-new). It writes named channel, energy, counts and
error columns, but `GET` can't read the file back.

#### `GETNRA` / `WRITENRA` `[new]`

```
usage: GETNRA <file> [-data] [-simulation]      (short: GN)
usage: WRITENRA <file>                          (short: WN)
```

Read and write SIMNRA's `.xnra` files (see
[File formats](file-formats.md#simnra-xnra) for the field mapping).

`GETNRA` loads the measured spectrum into buffer 1, the same as `GET`, and the
file's sample into SIM. If SIM already holds a sample, it asks before
replacing it. Inside an `XEQ` macro there is nobody to answer, so it replaces
the sample and says so. It prints every setting it had to ignore or convert,
the geometry in RUMP's terms, and the SIMNRA physics settings that pyRUMP
doesn't share, so differences between the two simulations are expected
rather than surprising.

* `-data` loads only the spectrum and leaves SIM alone.
* `-simulation` (or `-sim`) also loads SIMNRA's own simulated spectrum into a
  buffer of its own, after the measured data, to `OVERLAY` against pyRUMP's
  buffer 0. A file that holds only a simulation needs this flag.

`GET` also reads `.xnra` files, but only the spectrum; it never touches SIM.

`WRITENRA` writes the active buffer, the SIM sample and the simulation
(buffer 0). A bare filename gets `.xnra` added. A buffer read from an `.xnra`
file is written over its original document, so SIMNRA settings pyRUMP
doesn't model (cross-section choices, layer roughness, plot scaling) are
kept. A graded layer is written as the uniform sublayers pyRUMP simulates it
with, since SIMNRA has no profiles.

```
Your wish? gn MnPt.xnra       /* spectrum -> buffer 1, sample -> SIM */
Your wish? compare            /* pyRUMP's simulation against the data */
Your wish? wn MnPt-pyrump     /* writes MnPt-pyrump.xnra for SIMNRA   */
```

### Sample & instrument parameters

Each buffer carries its own beam, geometry, calibration and measurement
metadata. Every one of these **prints the current value with no argument,
and sets it (echoing the new value) with one** — and chains onto any
further command left on the line, so `Choff 0 FWHM 15` works in one go,
exactly as RUMP's own `WRASCII` output writes it back.

```
Your wish? beam 4He++
  beam Z=2 mass=4.0026 charge state 2
Your wish? mev 2.0
  MeV = 2
Your wish? conversion 5.0 0
  5 keV/channel, offset 0 keV
```

#### `ACTIVE`

Prints the active buffer's full parameter set.

#### `BEAM`

```
usage: BEAM 4He++
```

Sets the beam species and charge state.

#### `MEV`

```
usage: MEV <energy>
```

Sets the beam energy, MeV.

#### `THETA`

```
usage: THETA <deg>
```

Sets the sample tilt: the angle between the beam and the sample normal. See
[Experimental geometry](geometry.md).

#### `PHI`

```
usage: PHI <deg>
```

Sets 180° minus the scattering angle: the angle between the beam (looking
back at the source) and the detector. A detector at 170° is `PHI 10`. See
[Experimental geometry](geometry.md).

#### `PSI`

```
usage: PSI <deg>
```

Sets the exit angle, between the detector and the sample normal. Used with
`GEOMETRY GENERAL` only: `IBM` and `CORNELL` compute it from `THETA` and
`PHI`. See [Experimental geometry](geometry.md).

#### `GEOMETRY`

```
usage: GEOMETRY cornell|ibm|general
```

Sets how the exit angle `PSI` is obtained: `IBM` (the sample normal stays in
the plane of beam and detector), `CORNELL` (the tilt axis lies in that plane)
or `GENERAL` (`PSI` as typed). See [Experimental geometry](geometry.md), with
figures of the IBM and Cornell geometries.

#### `CONVERSION`

```
usage: CONVERSION <keV/ch> [keV(0)]
```

Sets the energy calibration: slope and, optionally, offset together.

#### `SLOPE` `[new]`

```
usage: SLOPE <keV/ch>
```

Sets the calibration slope alone, independent of `CONVERSION`'s offset.

#### `OFFSET` `[new]`

```
usage: OFFSET <keV(0)>
```

Sets the calibration offset alone, independent of `CONVERSION`'s slope.

#### `CORRECTION`

```
usage: CORRECTION <factor>
```

Sets the normalization fudge factor that absorbs charge-integration error
(RUMP's `CORR`). `PERT`'s `NORMALIZE` window can set this automatically from
a fit instead of you setting it by hand — see
[`PERT NORMALIZE`](pert.md#normalize).

#### `CHARGE`

```
usage: CHARGE <uC>
```

Sets the integrated beam dose.

#### `CURRENT`

```
usage: CURRENT <nA>
```

Sets the average beam current — enables pile-up modeling together with
`TAU`.

#### `CHOFF`

```
usage: CHOFF <n>
```

Sets the channel number of the first data point.

#### `FWHM`

```
usage: FWHM <keV>
```

Sets the detector resolution.

#### `OMEGA`

```
usage: OMEGA <msr>
```

Sets the detector solid angle.

#### `TAU`

```
usage: TAU <us>
```

Sets the MCA shaping time constant.

#### `IDENTIFIER`

```
usage: IDENTIFIER <text>
```

Sets a free-text spectrum description.

#### `DATE`

```
usage: DATE <text>
```

Sets when the spectrum was measured.

#### `FILENAME`

```
usage: FILENAME <name>
```

Sets the recorded source filename.

#### `SWALLOW`

```
usage: SWALLOW [-twocolumn]
```

Used inside an `XEQ` macro, reads the macro file's following lines straight
into the active buffer as channel data (or channel/value pairs), stopping at
the first blank line — how a RUMP-written `.cmd` file reconstructs a
spectrum inline.
