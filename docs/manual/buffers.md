## Buffers & instrument parameters

### Buffers

Spectra live in numbered buffers, one of which is ACTIVE and is what most
commands act on implicitly. Buffer **0 is the simulation**; data starts at
1. See [`GET`](shell.md#get-read) in Core workflow for how a freshly
read file gets slotted in.

#### `BUFFERS`

Lists the buffers, marking the active one.

#### `READ`

```
usage: GET <file|n>
```

Cross-referenced from [Core workflow](shell.md#get-read): `READ <file>`
reads a file into a buffer, same as `GET` except it never accepts a bare
buffer number.

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

Saves the active buffer, binary or text.

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

Sets the sample tilt.

#### `PHI`

```
usage: PHI <deg>
```

Sets 180° minus the scattering angle.

#### `PSI`

```
usage: PSI <deg>
```

Sets the exit angle (GENERAL geometry only).

#### `GEOMETRY`

```
usage: GEOMETRY cornell|ibm|general
```

Sets the detector geometry convention.

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
[PERT commands](sim-pert.md#pert-commands).

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
