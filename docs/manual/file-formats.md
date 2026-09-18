## File formats

Two unrelated things can carry a `.rbs`/`.RBS` extension. `GET` and `XEQ`
each expect one of them, and loading the wrong one with the wrong command
fails with a message naming the mistake, rather than a raw decode error or
silently-wrong data.

**RUMP binary format** (`GET`/`WRITE`) is a TIFF-inspired stream of
self-describing, big-endian 32-bit-word records:

```
[length_in_words] [record_type] [data...] [checksum]
```

`length` counts the header and checksum words themselves, so a record holds
`length - 3` data words, at most 1024. Every word in a record, checksum
included, sums to zero modulo 2^32 -- the format's own integrity check,
verified on every read. Records carry the program identifier, comments, the
sample identifier, date and live/real time, the dead-time correction,
accelerator parameters (E0, Zbeam, Mbeam, charge state, charge, current), MCA
calibration (keV/channel, keV at channel 0, first channel, FWHM), RBS
geometry (kind, theta, phi, psi, solid angle) and finally the counts, in one
of four encodings: raw floats, raw 32-bit integers, differential integers, or
differential integers plus zero-run compression. The spec
(`html/RUMP/rbs_inf.htm`, Aug 1994) is explicitly open and ships a reference
C API, unlike the rest of the legacy distribution -- see
[`pyrump/io/rbs.py`](https://github.com/ikost/pyRUMP/blob/main/src/pyrump/io/rbs.py)
for the clean-room implementation that reads and writes it.

**RC43 acquisition macro** (`XEQ`) is what some RBS acquisition software (NEC's
RC43 among them) writes instead, saved under a `.RBS` extension despite the
name suggesting real binary spectrum data. It's an ordinary text macro --
`IDENTIFIER`, `DATE`, `CONVERSION`, `MEV`, `BEAM`, `GEOMETRY`,
`THETA`/`PHI`/`PSI`, `OMEGA`, `CHOFF`, `FWHM`, `CURRENT`, `CHARGE` -- ending
in `SWALLOW`, which reads everything after it straight into the buffer as
channel data (see [Sample & instrument parameters](buffers.md#sample-instrument-parameters)
for `SWALLOW` itself). It's the same mechanism a RUMP-written `.cmd` file
uses to reconstruct a spectrum inline; RC43 just also uses it under a `.RBS`
name.

```
Your wish? get measured.rbs     /* binary format -- reads records directly */
Your wish? xeq acquired.rbs     /* text macro -- replayed as commands      */
```

**Plain ASCII spectrum** (`GET`) is anything `GET` is pointed at that isn't
one of the binary extensions above (`.rbs`, `.rump`, `.frs`, `.fres`,
`.pixe`) -- the format is sniffed from the file's content, not its name, so
`.dat`, `.asc`, `.ascii`, `.txt`, or no extension at all all work the same
way. Three dialects are recognized automatically:

* **one column** -- just the counts, one per line, in channel order
* **two column** -- `channel value` pairs
* **tab-delimited** -- as exported by a spreadsheet

Any leading lines that aren't numbers become the identifier; RUMP's own
[`WRASCII`](buffers.md#write-wrascii) output (a keyword header block ending
in the literal line `Swallow`, then one count per line) is also recognized
this way, matching what RUMP's own plain-ASCII reader does with it: only the
first non-numeric line becomes the identifier, the rest is read but not
otherwise acted on. So `GET`ting a plain ASCII file -- `WRASCII`'s own output
included -- always fills in your `~/.pyrumprc` defaults for
beam/geometry/detector metadata, never the header's own values
-- see [Config](config.md#pyrumprc).

```
Your wish? get counts.dat       /* one count per line, channel order       */
```

It's tempting to read `WRASCII`'s header (`Ident`/`Date`/`Charge`/`Conversion`/
`Theta`/... ending in `Swallow`) as the same kind of thing as an RC43 macro,
since the vocabulary and the `Swallow` terminator look alike -- but it isn't
one, and `XEQ`ing it fails immediately:

```
Your wish? xeq buffer.dat
ERROR: buffer.dat:2: unrecognized command: Spectrum
```

An RC43 macro's lines are real RUMP commands (`IDENTIFIER`, `CONVERSION`,
`MEV`, ...), meant to be replayed by the interpreter. `WRASCII`'s header
lines (`Spectrum`, `Ident`, ...) come from one hardcoded `fprintf` in the
original (`bmanip.c:619-641`) and were only ever meant to be read by a human
-- not one of them is a real command, so only `GET`'s lenient plain-ASCII
reader (which just treats any non-numeric line other than the first as
ignorable junk) can do anything with the file at all. `WRITE`/`GET` (binary)
and `WRASCII`/`GET` (text) are the two real pairs; there's no `WRASCII`/`XEQ`
pair to match RC43's `XEQ`-only dialect.

One consequence worth watching for: because both RC43 macros and `WRASCII`
output are plain text, and `GET` only tells them apart by checking whether a
file with a `.rbs`/`.RBS` extension is printable text at all (not by its
content), giving `WRASCII` output that same extension makes `GET` reject it
as a suspected RC43 macro -- `use XEQ, not GET, to load it` -- even though it
isn't one. `WRASCII` defaults to a `.dat` extension when you don't give one
of your own, and warns if you explicitly ask for `.rbs`/`.RBS`, but the
warning doesn't stop the write -- rename it before `GET`ting it back if you
want that to work.

See [Quick start](../getting-started/index.md#data-loading) for both binary
formats worked end to end, including what each of `GET`/`XEQ` prints when
pointed at the other one's file by mistake.
