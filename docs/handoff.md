# pyRUMP — state of work

Written to let a cold session resume without re-deriving anything. Read this
plus [rump-quirks.md](rump-quirks.md) and you have the context.

**Status:** M0–M9 complete. The forward model works end to end, including
depth-dependent composition profiles. 315 tests pass in ~50 s. Under git since
2026-07-28.

---

## The one-paragraph version

RUMP is a 1980s C program for simulating Rutherford backscattering spectra,
unmaintained and abandoned by its authors. pyRUMP reimplements it in Python.
The legacy C compiles cleanly under modern clang, so it is built into a shared
library and called directly from the test suite — every stage is validated
function-by-function against the original rather than by eyeballing a spectrum.
Full spectra now agree to ~3e-6 in total counts.

---

## Non-obvious things that took real work to find

Do not re-litigate these.

1. **The engine is `creatr.c`, not `sim2.c`.** `sim2.c` is the SIM command
   interpreter. Natural wrong turn; costs a day.

2. **The papers describe an algorithm the code abandoned.** The 1985 paper's
   parabolic bricks + exact Rutherford integral are `#if 0` at `anlyz.c:496`
   ("no longer using Doolittle qqq code"). The shipped C uses plain trapezoids.
   Porting "from the papers" produces something that does not match RUMP.

3. **RUMP never evaluates Ziegler during a simulation.** It refits a degree-5
   polynomial in √E once at startup and evaluates only that. Any port that
   calls the stopping model directly disagrees everywhere.

4. **`SimFillSpectrum` is a function pointer.** Redirecting it captures every
   brick from a completely unmodified `creatr.c`. This is what made M5–M8
   tractable; see `tests/oracle/csrc/sim_probe.c`.

5. **The `double` oracle build is impossible.** RUMP's readers use
   `scanf("%f")` against `REAL` fields, so `-DREAL_IS_DOUBLE` silently corrupts
   every table. The `float` build is authoritative, and all tolerances follow
   from that.

6. **`SQRT_DDS_POWER` has an index bug** (`p[2]` where the maths needs `p[1]`),
   35–50% error. Reproduced by default, `faithful=False` corrects it.

7. **Stopping tables are session-cached per (Z, mass)** and reused whenever the
   new beam energy fits the old window. Results are order-dependent; the oracle
   has `reset_stopping_tables()` for determinism.

8. **`SimStragf` rescales its own argument** (`x*(1+3*sig)`). Missing this makes
   it look wrong by 0.23 out of a range of 0.5.

9. **The SIM command is misspelled** `recalculculate` in the C's table; only
   the prefix `recal` works when driving the binary.

Full catalogue with citations: [rump-quirks.md](rump-quirks.md), 15 entries.

10. **`hfront` is recomputed every slab**, only `efront`/`ratde` carry over.
    Invisible with uniform composition; shifts the spectrum by one slab with a
    depth profile.

11. **Layer density is an inverse-density average** ("hard ball packing"), and
    it sets the depth scale of every profile equation. The 0.4997 constant in
    the C is a *fallback*, not silicon's actual 0.49777.

---

## Layout

```
src/pyrump/
  model/      element, geometry, spectrum (Calibration), detector (Measurement)
  atomic/     tables (PeriodicTable), density
  io/         atom4, scoef, kalbitzer          <- legacy table parsers
  stopping/   ziegler, kalbitzer, mylar, registry, polyfit, table, cache, bragg
  physics/    kinematics, xsec/rutherford
  sim/        slabs, precal, ideal, outbound, bricks, engine,
              fill/{trapezoid,straggled}, convolve
tests/
  oracle/     build_oracle.py, oracle.py (cffi), driver.py (pty),
              csrc/{stubs,ui_stubs,sim_probe,stragf_probe,oracle_api}.c
  unit/       one file per milestone
```

Entry point: `pyrump.sim.engine.simulate(sample, beam, geometry, registry,
periodic_table, calibration, measurement) -> Spectrum`.

---

## The oracle

Two independent mechanisms, both requiring the legacy C tree
(`PYRUMP_C_REFERENCE`, or `C-code/` in place). Tests skip cleanly without it.

```bash
python tests/oracle/build_oracle.py    # builds libpyrump_oracle_{float,double}
pytest                                 # everything
pytest -m oracle                       # only the C comparisons
```

* **cffi library** — physics translation units compiled into a dylib. No TTY,
  no graphics. Gives `zstop`, the fitted polynomials, cross-sections,
  `SimStragf`, and full `simulate_bricks()` / `simulate_spectrum()`.
* **pty driver** — drives the real `bin/rump` binary interactively.

`build_oracle.py` discovers missing symbols from the linker and generates
aborting stubs, so adding a source file does not require hand-maintaining a
stub list. Anything genuinely reachable must be linked for real, not stubbed —
`reschk`, `ArrayMinMax`, `FitPolynomial`, `g_sppfa/g_sppsl` all are.

Host state whose *values* matter (see `sim_probe.c`): `sigtab` defaults to
`{-1,-1}` — zero would silently select a manual-override branch and return a
cross-section of exactly zero.

---

## Agreement achieved

| Quantity | Agreement |
|---|---|
| ZBL85 stopping, 92 targets × 4 beams, 10 keV–10 MeV | 6.1e-7 |
| Fitted stopping polynomial | 1.3e-5 |
| Polynomial evaluation given identical coefficients | 1.1e-14 |
| Cross-sections, kinematics | 1e-10 |
| Bricks, 630 across 36 configurations | 5e-7 energies, 6e-6 heights |
| Full spectrum, total counts | 3e-6 |
| Full spectrum, per channel | 1e-5 of peak |

Everything is capped by float32 storage in the C. Tolerances are argued from
that floor, not chosen for convenience.

---

## Testing conventions worth keeping

Established by repeated failure — four times a "bug" turned out to be a bad
test, never the code:

* **Derive expected values from the same source the engine uses.** Literature
  tables use modern atomic weights; `atom4.dat` is from 1993. Isotopic vs
  average mass shifts an edge by enough to change channel.
* **A partial channel can hold 30% of the peak.** Do not mask on magnitude.
  Use relative *or* small-vs-peak.
* **Test the function, not the coefficients.** The monomial basis in √E is
  near-degenerate: coefficients drift 2e-4 while the curve moves 5e-7.
* **Separate formula error from fit error.** Feed the oracle's own coefficients
  through our evaluator to test the formula in isolation.
* Where pyRUMP deliberately departs (only `SimStragf`), measure the difference
  against an independent reference rather than asserting it away.

---

## Next: M10 — absorber, foil, pileup, fuzz

* **Absorber layers** (`sample->absorber_layers`, `samm->fsurf`) sit in front of
  the sample: traversed on the exit path only, and forced to normal incidence
  (`creatr.c:1971` — they are not tilted with the sample).
* **Fuzz** — lateral thickness non-uniformity. N Gaussian-weighted
  re-simulations via `ndtri` quantiles (`creatr.c:665-686`), summed with weight
  `samm->ampl`. Best implemented as a batched leading axis, not N reruns.
* **Pileup** — `SimNewPileup` (`creatr.c:1331`, Custer thesis) needs `tau` and
  `current`; plus a legacy `SimOldPileup`.
* **Multiple scattering** — an ad-hoc exponential tail with no physical basis
  (`creatr.c:337-345`). Mark `experimental`.
* **Stopper foil** — `stopfoil.c` is a pure lookup table and **no data file
  ships**; the format is documented only inside `#ifdef TESTING`. Probably
  document as unsupported rather than implement.

`sim_probe.c` will need absorber/fuzz/pileup fields on `OracleSetSample` to
validate these.

Then: M11 file I/O (`.RBS` binary is fully specified in `html/RUMP/rbs_inf.htm`),
M12 fitting, M13 CLI.

---

## Licensing position (settled 2026-07-28)

**CGS ceased operating as a business in June 2012** (Cornell conflict-of-interest
policy), stating that GENPLOT and RUMP remain free to download and use.
`genplot.com` no longer resolves. Established from web search summaries; the
primary page could not be fetched directly.

Consequences:

* **Trademark: no longer a concern.** `pyRUMP` is fine.
* **Copyright: unchanged.** A dissolved company does not extinguish the authors'
  copyright, so the C tree is still not redistributed and `.gitignore` still
  excludes it.
* **`pscoef.dat`: unchanged, and independent of CGS.** It is the ZBL/TRIM SCOEF
  table (Ziegler/Pergamon). Regenerate from published tables before release, as
  recorded in `NOTICE`.

Under git since 2026-07-28; `LICENSE` (MIT) and `NOTICE` are in place.
