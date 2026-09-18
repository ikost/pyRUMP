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

See [Quick start](../getting-started/index.md#data-loading) for both binary
formats worked end to end, including what each of `GET`/`XEQ` prints when
pointed at the other one's file by mistake.
