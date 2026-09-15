## Design and validation

**Faithful first, corrected by choice.** The default reproduces the shipped C
bug-for-bug, because that is what every published RUMP result was produced
with. Known defects — and there are several — are reproduced exactly, with
the mathematically correct behaviour available behind explicit flags. The
shell exposes this as a session setting, `session.settings.faithful`,
toggled with the `FAITHFUL` command and persisted through `~/.pyrumprc` (see
[Macros](../manual/shell.md#macros)) rather than a separate branch or fork — corrected and
faithful behaviour live in the same codebase so they stay comparable against
the C oracle side by side. `FAITHFUL` governs physics only — for pyRUMP's own
command-surface additions beyond stock RUMP, see the
[Changelog](https://github.com/ikost/pyRUMP/releases)'s versioning note. See
[RUMP quirks and defects found while porting](../dev/rump-quirks.md).

**Validated against the original**, at two levels. The legacy C is compiled
into a shared library and called directly from the test suite, so each stage
is compared function-by-function rather than by eyeballing a final spectrum.

**Unit oracle** — the physics translation units (`ziegler.c`, `stopping.c`,
`sigma.c`, …) are compiled into `libpyrump_oracle` and called via cffi. No TTY,
no graphics, no buffers. When a number disagrees, this isolates the cause to one
function.

```bash
python tests/oracle/build_oracle.py
pytest -m oracle
```

**End-to-end oracle** — the original `rump` binary is driven through a
pseudo-terminal to produce reference spectra.

Both require the legacy C tree, which is **not redistributed** (see
[Licensing and provenance](../dev/about.md#licensing-and-provenance)). Point `PYRUMP_C_REFERENCE`
at it, or place it at `C-code/`. Tests skip cleanly when it is absent.

```bash
pytest              # unit tests
pytest -m oracle    # comparison against the C
```

### Current agreement

| Quantity | Agreement | Limited by |
|---|---|---|
| ZBL85 stopping, all 92 targets, H/He/Cu/Au beams, 10 keV–10 MeV | **6.1e-7** rel | float32 tables in the C |
| Fitted stopping polynomial (what the simulation consumes) | **1.3e-5** rel | float32 coefficient storage |
| Polynomial evaluation, given identical coefficients | **1.1e-14** rel | float64 round-off |
| Cross-sections and kinematics | **1e-10** rel | closed forms, nothing to fit |
| Bricks — 630 across 36 configurations | **5e-7** energies, **6e-6** heights | float32 coefficients |
| Full spectrum, total counts | **3e-6** rel | float32 coefficients |
| Full spectrum, per channel | **1e-5** of peak | float32 brick edges |
| With straggling and detector resolution | **3e-6** total, **1e-5** of peak | float32 brick edges |
| Depth profiles, all 11 evaluable forms | **2.6e-5** brick heights | float32 coefficients |
| Absorber, fuzz, multiple scattering | **3e-6** total, **4e-5** of peak | float32 brick edges |
| `.RBS` files read vs RUMP's own reader | **bit-identical** | — |
| Poisson objective vs `EvalChiPoisson` | **1e-5** reduced chi2 | float32 in the C |
| `.lcm` round-trip vs RUMP's own writer | **byte-identical** | — |

The oracle is the `float` build. RUMP cannot be built in double precision — its
table readers use `scanf("%f")` against `REAL` fields, so `-DREAL_IS_DOUBLE`
silently corrupts every table. That caps how tightly any float64 port can agree,
and the tolerances above are set by that floor rather than by choice.


### Accuracy

Stage-by-stage agreement with the original C is in
[Design and validation § Current agreement](#current-agreement) — every figure
there is limited by float32 storage in the C, not by pyRUMP. For a realistic
multi-layer sample with micron-thick polymer layers, agreement loosens to
~3e-3: thousands of sublayers accumulate single-precision differences, and part
of the beam falls below the stopping cutoff.

Full bibliography in [References](../dev/about.md#references).

### Atomic data provenance

The stage-by-stage agreement above is *internal* — pyRUMP against the shipped C, both
consuming the same bundled tables. It says nothing about whether those tables' underlying
physical constants are still correct. That's checked separately, against external reference
sources, and re-verified periodically rather than assumed from the tables' 1990s origin.

**Atomic masses and isotope abundances** (`atom4.dat`) — checked against CIAAW (Commission on
Isotopic Abundances and Atomic Weights) current standard atomic weights, Z=1-92. Every sampled
mass matches within <0.1% relative. Two elements have isotope-abundance intervals CIAAW has
revised since the table's origin — magnesium (post-2011 remeasurement) and zirconium (revised
2024) — but the effect on each element's average mass is negligible (<0.03% and <0.1%
respectively) and was left as-is rather than patched, to avoid disturbing the isotope table's
internal normalization.

**Elemental densities** (`atomic_density`) — checked against NIST's X-Ray Mass Attenuation
Coefficients reference densities (Hubbell & Seltzer, NISTIR 5632). Median relative difference
**0.3%** across 75 comparable elements. Three of the four elements differing by more than 5% are
allotrope-convention mismatches rather than errors (carbon: graphite vs. a lower bulk value;
phosphorus: white vs. red; selenium: grey vs. an amorphous form) — both values are individually
defensible, they just describe different physical forms of the element. The fourth, **calcium**,
had no such explanation and was corrected to NIST's reference value.

Full detail, including per-element figures, in [`SOURCES.md`](https://github.com/ikost/pyRUMP/blob/main/src/pyrump/data/SOURCES.md)
alongside the tables themselves.

