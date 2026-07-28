# pyRUMP — state of work

Written to let a cold session resume without re-deriving anything. Read this
plus [rump-quirks.md](rump-quirks.md) and you have the context.

**Status:** M0–M10 complete — the entire forward model. 382 tests pass in ~50 s.
Under git; `LICENSE` (MIT) and `NOTICE` in place.

---

## The one-paragraph version

RUMP is a 1980s C program for simulating Rutherford backscattering spectra,
unmaintained since its authors wound the business down in 2012. pyRUMP
reimplements it in Python. The legacy C compiles cleanly under modern clang, so
it is built into a shared library and called directly from the test suite —
every stage is validated function-by-function against the original rather than
by eyeballing a spectrum. Full spectra agree to ~3e-6 in total counts.

---

## Non-obvious things that took real work to find

Do not re-litigate these. Full catalogue with citations:
[rump-quirks.md](rump-quirks.md), 17 entries.

1. **The engine is `creatr.c`, not `sim2.c`.** `sim2.c` is the SIM command
   interpreter. Natural wrong turn; costs a day.
2. **The papers describe an algorithm the code abandoned.** The 1985 paper's
   parabolic bricks + exact Rutherford integral are `#if 0` at `anlyz.c:496`.
   Porting "from the papers" produces something that does not match RUMP.
3. **RUMP never evaluates Ziegler during a simulation.** It refits a degree-5
   polynomial in √E once at startup and evaluates only that.
4. **`SimFillSpectrum` is a function pointer.** Redirecting it captures every
   brick from a completely unmodified `creatr.c` — this is what made M5–M10
   tractable. See `tests/oracle/csrc/sim_probe.c`.
5. **The `double` oracle build is impossible.** RUMP's readers use
   `scanf("%f")` against `REAL` fields, so `-DREAL_IS_DOUBLE` silently corrupts
   every table. The `float` build is authoritative; all tolerances follow.
6. **`SQRT_DDS_POWER` has an index bug** (`p[2]` where the maths needs `p[1]`),
   35–50% error. Reproduced by default; `faithful=False` corrects it.
7. **Stopping tables are session-cached per (Z, mass)** and reused whenever the
   new beam energy fits the old window. Order-dependent; the oracle has
   `reset_stopping_tables()`.
8. **`SimStragf` rescales its own argument** (`x*(1+3*sig)`). Missing this makes
   it look wrong by 0.23 out of a range of 0.5.
9. **`hfront` is recomputed every slab**; only `efront`/`ratde` carry over.
   Invisible with uniform composition; shifts the spectrum one slab with a
   depth profile.
10. **Layer density is an inverse-density average** ("hard ball packing"), and
    it sets the depth scale of every profile equation. The 0.4997 constant in
    the C is a *fallback*, not silicon's actual 0.49777.
11. **The absorber is not traversed inbound.** It sits between sample and
    detector; `SimPrecal` seeds the march at `fsurf`. Getting this wrong costs
    ~20% of the yield.
12. **`NDTRI` reads one element past its coefficient array.** Harmless in
    RUMP's own build, catastrophic when transplanted. pyRUMP uses scipy.
13. **The SIM command is misspelled** `recalculculate` in the C's table; only
    the prefix `recal` works when driving the binary.

---

## Layout

```
src/pyrump/
  model/      element, geometry, spectrum (Calibration), detector (Measurement)
  atomic/     tables (PeriodicTable), density (incl. inverse-density average)
  io/         atom4, scoef, kalbitzer          <- legacy table parsers
  stopping/   ziegler, kalbitzer, mylar, registry, polyfit, table, cache, bragg
  physics/    kinematics, xsec/rutherford
  profiles/   equations (13 EQUATION forms + mixing rule)
  sim/        slabs, precal, ideal, outbound, bricks, engine,
              fill/{trapezoid,straggled}, convolve,
              absorber, fuzz, pileup, multiscatter
tests/
  oracle/     build_oracle.py, oracle.py (cffi), driver.py (pty), csrc/*.c
  unit/       one file per milestone
```

Entry point:

```python
pyrump.sim.engine.simulate(
    sample, beam, geometry, registry, periodic_table, calibration, measurement
) -> Spectrum
```

Pipeline order matters and follows `SimCreateDetails` (creatr.c:307-345):
fuzz replicas → fill (trapezoid or straggled) → detector convolution → yield
normalisation → pile-up → multiple-scattering tail.

---

## The oracle

Two mechanisms, both requiring the legacy C tree (`PYRUMP_C_REFERENCE`, or
`C-code/` in place). Tests skip cleanly without it.

```bash
python tests/oracle/build_oracle.py    # builds libpyrump_oracle_{float,double}
pytest                                 # everything
pytest -m oracle                       # only the C comparisons
```

* **cffi library** — physics translation units compiled into a dylib. No TTY,
  no graphics. Exposes `zstop`, fitted polynomials, cross-sections,
  `SimStragf`, `NDTRI`, and full `simulate_bricks()` / `simulate_spectrum()`
  with sample, beam, equation, absorber, fuzz and pile-up setters.
* **pty driver** — drives the real `bin/rump` binary interactively.

`build_oracle.py` discovers missing symbols from the linker and generates
aborting stubs, so adding a source file needs no hand-maintained stub list.
Anything genuinely reachable must be linked for real — `reschk`, `ArrayMinMax`,
`FitPolynomial`, `g_sppfa/g_sppsl` all are.

Two traps in the host state (`sim_probe.c`):
* `sigtab` defaults to `{-1,-1}`; zero silently selects a manual-override branch
  and returns a cross-section of exactly zero.
* `gvcalc.c` will not compile under clang (declares `srand48`/`drand48`
  non-static, then defines them static), so `NDTRI` is extracted verbatim into
  `ndtri_probe.c` with its array padded — see quirk 16.

---

## Agreement achieved

| Quantity | Agreement |
|---|---|
| ZBL85 stopping, 92 targets × 4 beams, 10 keV–10 MeV | 6.1e-7 |
| Fitted stopping polynomial | 1.3e-5 |
| Polynomial evaluation given identical coefficients | 1.1e-14 |
| Cross-sections, kinematics | 1e-10 |
| Bricks, 630 across 36 configurations | 5e-7 energies, 6e-6 heights |
| Depth profiles, 11 evaluable forms | 2.6e-5 brick heights |
| Full spectrum, total counts | 3e-6 |
| Full spectrum, per channel | 1e-5 of peak |
| Absorber / fuzz / pile-up / multiple scattering | 3e-6 total, 4e-5 of peak |

Everything is capped by float32 storage in the C. Tolerances are argued from
that floor, not chosen for convenience.

---

## Testing conventions worth keeping

Established by repeated failure — five times a "bug" was a bad test, not code:

* **Derive expected values from the same source the engine uses.** Literature
  tables use modern atomic weights; `atom4.dat` is from 1993. Isotopic vs
  average mass shifts an edge enough to change channel.
* **When a comparison fails, check both sides are configured identically before
  touching the physics.** One absorber "failure" was the setting applied to the
  oracle but not the sample.
* **A partial channel can hold 30% of the peak.** Do not mask on magnitude; use
  relative *or* small-vs-peak.
* **Test the function, not the fitted coefficients.** The monomial basis in √E
  is near-degenerate: coefficients drift 2e-4 while the curve moves 5e-7.
* **Separate formula error from fit error** by feeding the oracle's own
  coefficients through our evaluator.
* **Assert that an effect actually changes the spectrum**, or a comparison can
  pass because both sides did nothing.
* Where pyRUMP deliberately departs (only `SimStragf`), measure the difference
  against an independent reference rather than asserting it away.

---

## Next: M11 — file I/O

**`.RBS` binary** is fully specified in `C-code/html/RUMP/rbs_inf.htm`, and the
spec explicitly permits reimplementation. Reader/writer: `rump/rbs_rdwr.c`.

* 32-bit words, **big-endian on disk**
* record = `[length][type][data...][checksum]`, 3–1027 words; **all words in a
  record sum to 0 mod 2^32**
* first record must be type `0000h`; program id `10211210h` = RUMP
* types: `0000h` program, `0001h/0002h` comments, `0010h` field initiator,
  `0020h` array initiator, `0011h` data, `0101h` id, `0102h` livetime,
  `0103h` date, `0110h` correction, `0111h` accelerator, `0112h` MCA,
  `0120h` RBS geometry, `0121h` FRES geometry
* four compression modes: unpacked real, unpacked int, differential int,
  differential + zero-run. Differential: first element 4-byte absolute, then
  1-byte deltas, escape `80h` → 2-byte, escape `8000h` → 4-byte absolute
* geometry codes `0` Cornell, `1` IBM, `-1` general; **`Phi` is the supplement**
  of the scattering angle

Acceptance: the three `C-code/rump/data/Fixed/*.rbs` fixtures round-trip
byte-identically, and all four compression modes decode. The C's writer can
generate fixtures for modes the shipped files do not cover.

Also in M11: ASCII 1col/2col/`.xls` (`rdwr.c`), and `.adt`/R33 cross-section
tables in three dialects (`reswork.c:194`).

Then: M12 fitting (Baker-Cousins Poisson χ², error + normalisation windows, LM
via scipy — test the χ² surface, never the optimiser trajectory), M13 CLI.

---

## Licensing position (settled 2026-07-28)

**CGS ceased operating as a business in June 2012** (Cornell conflict-of-interest
policy), stating GENPLOT and RUMP remain free to download and use.
`genplot.com` no longer resolves. Established from web-search summaries; the
primary page could not be fetched.

* **Trademark: no longer a concern.** `pyRUMP` is fine.
* **Copyright: unchanged.** A dissolved company does not extinguish the authors'
  copyright, so the C tree is not redistributed and `.gitignore` excludes it.
* **`pscoef.dat`: unchanged, and independent of CGS.** It is the ZBL/TRIM SCOEF
  table (Ziegler/Pergamon). Regenerate from published tables before release, as
  recorded in `NOTICE`.
