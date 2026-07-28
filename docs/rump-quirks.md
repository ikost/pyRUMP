# RUMP quirks and defects found while porting

Running catalogue of behaviour in the legacy C that is surprising, wrong, or
simply undocumented. Each entry says what pyRUMP does about it.

Default policy (per the project plan): **reproduce faithfully, expose the fix
behind a flag.** A silent "correction" would make pyRUMP disagree with every
published RUMP result.

---

## 1. The papers describe an algorithm the code no longer uses

The 1985 NIM-B paper's signature contribution — parabolic brick tops with an
exact analytic area from the "Rutherford integral" ∫E(a)⁻²da — is **dead code**.
`anlyz.c:496` sits inside `#if 0`, labelled *"no longer using Doolittle qqq
code"*. The shipped path is a plain trapezoid (`SimAnlyz4`, `anlyz.c:304`) or,
with straggling, two Gaussian-convolved triangles (`SimAnlyz3`, `anlyz.c:244`).

`layer[].qq` is still computed every run in `SimPrecal` and then discarded.

**pyRUMP:** trapezoid by default; the parabolic form is planned as an opt-in mode,
for which the dead C is the specification.

---

## 2. `SQRT_DDS_POWER` uses the wrong coefficient index

`stopping.h:47-49` computes d²S/dE² as

```c
(((((3.75*p[5])*e + 2*p[4])*e + 0.75*p[3])*e*e - 0.25*p[2]) / (e*e*e))
                                                       ^^^^ should be p[1]
```

The header's own comment two lines earlier gives the correct form,
`DDS = 1/(4e⁴) · Σ i(i-2)·aᵢ·eⁱ`, whose `i=1` term is `−p[1]·e`. There is no
`p[2]` term at all, since `i(i-2)` vanishes at `i=2`.

Verified both directions: the oracle reproduces the buggy macro exactly, and a
numerical second derivative reproduces the corrected form exactly. The error is
**35–50%**.

Impact is limited — d²S/dE² only enters the third-order term of the energy-loss
expansion (`creatr.c:1554`) — but it is real and it shifts the depth scale.

**pyRUMP:** `StoppingTable.derivative(..., order=2)` reproduces the bug;
`faithful=False` gives the correct value. Covered by
`test_second_derivative_reproduces_rump_bug`.

---

## 3. `-DREAL_IS_DOUBLE` silently corrupts every data table

RUMP's table readers hard-code `%f` in their `scanf` formats while writing into
`REAL` fields (`ziegler.c:100,114`; `atomio.c:162,163,173`). With `REAL` as
`double` these write four bytes into an eight-byte field, so the Ziegler and
atomic tables load as denormal garbage and `zstop` returns 0 without complaint.

`REAL_IS_DOUBLE` is never set by any shipped makefile, so the bug stayed latent.

**Consequence:** the plan's dual-precision oracle — build twice, diff to measure
RUMP's own float32 noise — is not possible without patching the readers, which
would make the oracle a *modified* RUMP. The `float` build is authoritative, and
tolerances are argued from single-precision reasoning instead.
`tests/oracle/oracle.py` refuses to load a corrupt build rather than returning
garbage.

---

## 4. Stopping tables are cached per (Z, mass) for the whole session

`RbsStpfind` (`stopping.c:274-279`) reuses an existing fitted table whenever the
new beam energy merely *fits inside* its window (`2·emin ≤ E ≤ emax`). Simulating
at 3 MeV and then at 2 MeV does **not** refit: the second run silently uses the
3 MeV window, and its coefficients differ from a fresh 2 MeV fit.

When only Z matches (different isotope), the table is reused with
`e_scale = table_mass / beam_mass` — the Amsel energy-scaling trick that lets one
table serve 3He/4He or H/D.

This is stateful behaviour that changes numbers, not an optimisation.

**pyRUMP:** `StoppingTableCache` reproduces the reuse rules, including
`e_scale`. The oracle gained `OracleResetStoppingTables()` so tests are
order-independent.

---

## 5. The fitted polynomial is only valid inside the beam-dependent window

`emin = 0.04·E_beam`, `emax = 1.15·E_beam` (`stopping.c:316-319`, STOP_SQRT).
Beyond `emax` the degree-5 fit diverges from the underlying model within a few
hundred keV — 28% at 3 MeV for a 2 MeV table — and RUMP returns the extrapolated
value without warning.

More importantly: **RUMP never evaluates Ziegler or Konac during a simulation.**
It fits once at startup and evaluates only the polynomial. Any port that calls
the stopping model directly will disagree with RUMP everywhere.

---

## 6. Konac/Kalbitzer outranks Ziegler

The priority chain (`stopping.c:479-515`) tries `newstop.kal` *before* ZBL. For
H/D/³He/⁴He on carbon or silicon — the most common RBS cases by far — RUMP is not
using Ziegler at all. The two models differ by up to ~10%.

`newstop.kal` is also not mentioned in the 1996 manual.

---

## 7. `ConvoluteDetector` loses counts at **both** edges

*(Corrected — an earlier version of this entry described only the high-energy
edge. The behaviour is symmetric.)*

Contributions that would land outside the channel range are discarded:

* the head loop drops the `i-k` half once `k > i` (`creatr.c:1233-1236`)
* the tail loop drops the `i+k` half once `k >= npt-i` (`creatr.c:1284-1291`)

Measured with a 30 keV FWHM at 5 keV/channel on a 200-channel spectrum, a delta
function of 1000 counts retains:

| Position | Retained |
|---|---|
| channel 100 (interior) | 999.15 — only the 3σ truncation is lost |
| channel 2 | 836.34 |
| channel 197 | 836.34 |

So this is not really a bug: a real MCA cannot record counts in channels it does
not have. But it does mean the integral of a convolved spectrum is below the
original whenever intensity sits within ~3σ of either end, which matters when
comparing totals.

**pyRUMP:** `convolve_edge='rump'` (default) reproduces it; `'renormalize'`
conserves counts by scaling each *source* channel's contribution by the kernel
weight that actually lands in range. Note the normalisation must be applied on
the source side — normalising the output instead rescales *received* weight and
inflates the edges rather than conserving.

---

## 8. The SIM `RECALCULATE` command is misspelled in the command table

`sim2.c:252` registers `{"recalculculate", 5, SIM_RECALCUL}` — note
`recalcul**cul**ate`. Typing `recalculate` does **not** match and silently falls
through to the top-level shell. Only the 5-character prefix `recal` works.

Relevant to anyone driving the legacy binary; `tests/oracle/driver.py` documents
and uses `recal`.

---

## 9. Straggling is Bohr-only, off by default, and combined incorrectly

`sample->straggle` defaults to 0 (straggling disabled). When enabled, the inbound
and outbound path variances are combined *linearly* as
`stragc = sec θ_in·K + sec θ_out` applied to the inbound Bohr variance
(`creatr.c:1661`), rather than as `K²σ²_in + σ²_out`. No Chu correction exists
anywhere in the simulation.

---

## 10. Physics that is absent

* no channelling
* no Q≠0 nuclear reactions (`reswork.c:355` rejects them explicitly)
* no Andersen screening — L'Ecuyer only
* no relativistic correction to the Rutherford cross-section
* multiple scattering is an ad-hoc exponential tail with no physical basis
  (`creatr.c:337-345`), described in the C itself as *"Ad-hoc scaling"*

---

## Dead or broken code to avoid porting

| Location | Status |
|---|---|
| `anlyz.c:496-642` | `#if 0` — the paper's parabolic algorithm; also syntactically incomplete at `:600-603` |
| `tables.c` | not in `OBJS`; does not compile (three typos). The live density-table loader is `sim2.c:1415` |
| `creatr.c:2100-2142` (`SimInloss`) | `#if 0`; useful as clean documentation of the ΔE expansion |
| `sigma.c:372-558` | `#if 0` — superseded Turos & Meyer cross-sections |
| `MyXsect/Xsect.c` | placeholder; writes to an undeclared variable, hook never called |
| `data/*.stp` | unreachable by default: stores `STOP_LINEAR` coefficients while the runtime is `STOP_SQRT` (`stopping.c:139` vs `:276`) |
| channelling hooks (`sigma_scale`, `dedx_scale`) | declared, never set; the only use is commented out at `creatr.c:1114` |

---

## 11. `reschk` is compiled out of `reswork.c` but referenced by `creatr.c`

`reswork.c:55` guards the resonance-table index with `#ifdef RESONANCE`, and the
`#define` lives in `xsect.h` — which `reswork.c` includes *after* some of its own
declarations. Depending on include order the array is not emitted, while
`creatr.c:1723` reads `reschk[z2]` unconditionally on the simulation path.

In the shipped build this happens to link because other translation units pull
`xsect.h` in first. It is fragile, and it bit the oracle build: a link that looks
successful can leave the resonance index pointing at nothing.

**pyRUMP:** `reswork.c` is linked into the oracle explicitly so the array exists
and is properly zero-initialised (no resonance tables loaded ⇒ pure Rutherford).

---

## 12. `SimStragf` rescales its own argument

`SimStragf(x, sig)` (anlyz.c:371) computes `newx = x * (1 + 3*sig)` as its first
real statement (anlyz.c:387). Its `x` is therefore in units of the **broadened**
width `|de| + 3*sqrt(2)*sigma`, not of the triangle base — which is what the
docstring's "triangle ... height 0 for x<0 or x>1" implies.

The rescaling exactly undoes the `fact = 1/(|de| + 3*sqrt(2)*sigma)`
normalisation the caller applies (anlyz.c:257), so the composition collapses to
`(E_j - E_peak)/|de|`.

Overlooking this makes the function appear wrong by up to **0.23** — comparable
to its entire range of 0.5 — while looking perfectly plausible in isolation.

**pyRUMP:** `stragf()` reproduces the C's contract including the rescaling;
`triangle_gaussian_integral()` is the underlying maths in unscaled coordinates.

---

## 13. Measured: `SimStragf`'s approximation error

The 1985 paper justifies the rational fit on 1985 hardware grounds:

> Analytic expressions for the functions f and g … are of limited utility.
> Direct evaluation is unnecessarily slow and often involves finding the small
> difference of large numbers. Single precision computation is inadequate …

Against numerical quadrature, over sigma in [0.05, 3] and x in [-1.5, 1.7]:

| | Worst absolute error |
|---|---|
| pyRUMP closed `erf` form | **6.4e-14** |
| RUMP `SimStragf` | **1.7e-6** |

So the approximation is good — its claim of "rapidly and accurately" holds — but
it is ~8 orders of magnitude looser than the closed form is in float64. pyRUMP
therefore does *not* port it: the closed form is more accurate, vectorises, and
avoids transcribing a table of hand-tuned constants.

This is the one place pyRUMP knowingly departs from bit-faithfulness. The
resulting spectrum still matches the C to ~3e-5 of peak.

---

## 14. `hfront` is recomputed per slab; only `efront` and `ratde` carry over

`SimCideal` guards the front-edge recomputation with `ok` (creatr.c:1791):

```c
if (! ok) {                       /* Recompute EFRONT & RATDE */
    efront = km2 * ein;
    SimFlyout(lay-1, &efront, &ratde);
}
rfront = strct[elno] * ratde * (secin/RbsEfact(...)) * fisot;   /* :1804 -- outside */
hfront = sigma * rfront;
```

Only the *geometry* is cached. The front **height** is recomputed every slab
from that slab's own areal density, because it is outside the guard.

With uniform composition the two are indistinguishable — `strct` is identical in
every slab — so a port can pass every uniform-layer test while getting this
wrong. It only shows up with a depth profile, where reusing the height shifts
the entire spectrum by one slab.

---

## 15. Layer density is an inverse-density average, and it sets the depth scale

`creatr.c:606-625` averages **cm³/atom**, not atoms/cm³ — the C's own comment
calls it "the idea of hard ball packing":

```
rho = ( sum_i (x_i / rho_i) / sum_i x_i ) ^ -1
```

The result converts areal thickness to physical thickness
(`cm_thick = cm2_thick/matrix_density/1E8`, creatr.c:692), so it sets the depth
scale of every depth-dependent profile equation. Getting it wrong — for example
using the 0.4997 silicon *fallback* constant instead of silicon's actual 0.49777
— leaves position-fraction forms (CONSTANT, LINEAR) matching perfectly while
every physical-depth form is off by a few percent.

The species composition gets its own density, used to convert Angstrom doses
(THINFILM, BURIEDTHINFILM) into 1e15 at/cm².

A pre-1997 `COMPATIBLE` mode averaged densities directly; pyRUMP implements only
`IMPROVED`, the shipped default.

---

## 16. `poly_e` reads one element past its array — and `-O2` makes it fatal

`poly_e(x, numer, iorder)` (gvcalc.c:4914) does

```c
numer += iorder;  tmp = *numer;
while (iorder--) tmp = tmp*x + *(--numer);
```

i.e. it touches `numer[0 .. iorder]` — **iorder+1** values. `NDTRI` calls it as
`poly_e(t*t, taylor, 10)` against a 10-element `taylor[]`, so it reads one past
the end.

At `-O0` the adjacent static happens to be tiny and Horner folds it to nothing,
so the bug is invisible. At `-O2` the compiler exploits the undefined behaviour
and `NDTRI` returns ±0.15 for **every** argument in roughly (0.15, 0.85) —
correct in the tails, badly wrong in the middle.

**The shipped RUMP passes no optimisation flag at all** (`makeosx.h`: `GOPTS`
has no `-O`, `CCOPT` is empty), so it runs the benign version. Anyone rebuilding
RUMP with optimisation enabled would silently corrupt every FUZZ profile.

**pyRUMP:** `tests/oracle/build_oracle.py` compiles `ndtri_probe.c` at `-O0`
specifically, with a test (`test_ndtri_matches_the_c`) that fails loudly if that
is ever lost.

---

## 17. The inbound march starts at `fsurf`, skipping the absorber

`SimPrecal` seeds `samm->layer[samm->fsurf].ehit = ee` with the **full** beam
energy (creatr.c:1530) and marches from there. Absorber layers are between the
sample and the detector, so the incoming beam never crosses them — only the
outgoing particle does, via `SimFlyout`.

Marching through them on the way in double-counts their stopping and shifts
every edge. The absorber is also not tilted with the sample: `SimFlyout` forces
normal incidence through it (creatr.c:1971), since a detector window does not
rotate when the sample does.

---

## 16. `NDTRI` reads one element past its coefficient array

`gvcalc.c:4764` declares the Taylor coefficients with **ten** elements:

```c
static double taylor[] = { -4.59e-10, 0.1667, 0.0583, ..., -0.000338 };   /* 10 */
...
t += t*poly_e(t*t, taylor, 10);
```

but `poly_e(x, numer, iorder)` starts at `numer[iorder]` and works down, so it
reads `taylor[10]` — one past the end.

In RUMP's own build the next object in memory is `ctay[0]` (~6e-9), whose
contribution at x^10 is negligible, so the bug never shows. Under any other
layout it reads whatever happens to be there: when the function was extracted
into the oracle it returned a constant ±0.15 across the entire central region
(|p−0.5| < 0.35), against a true range of ±∞.

`NDTRI` feeds the FUZZ feature (creatr.c:679), so a port that transplanted the
routine without noticing would get plausible-looking but wrong roughness
quantiles.

**pyRUMP:** uses `scipy.special.ndtri`. With the array padded to eleven elements
(the eleventh being 0.0, which is what the series actually calls for), RUMP's
routine agrees with scipy to **8.4e-10** across p ∈ [1e-6, 0.999] — so the
substitution is exact within the C's own precision.

---

## 17. The absorber is *not* traversed on the way in

`SimPrecal` seeds the inbound march at `samm->fsurf`, not at slab 0:

```c
samm->layer[samm->fsurf].ehit = ee;                 /* creatr.c:1533 */
for (lay=samm->fsurf; lay<samm->num_layers; lay++)  /* creatr.c:1541 */
```

Absorber layers sit between the sample and the **detector** — a dead layer or a
Mylar window — so the beam reaches the sample undegraded and the absorber
attenuates only the outgoing particle. They are also traversed at normal
incidence, because a detector window does not tilt when the sample does
(`creatr.c:1971`).

Marching the beam through the absorber on the way in as well costs ~20% of the
total yield for a 200e15 at/cm² Si absorber, and grows from there.

---

## Notes on driving the engine from outside

`creatr.c`'s output stage is a **function pointer**, `SimFillSpectrum` (sample.h),
called once per brick. Redirecting it captures the engine's exact intermediate
results with the C completely unmodified — no `#ifdef` hooks, no patched copy.
`SimPileup`, `SimInitFillSpectrum` and `SimTermFillSpectrum` are pointers too.

Host state the engine needs, and where the real versions live:

| Symbol | Source | Note |
|---|---|---|
| `RbsNormK` | `bmanip.c:880` | sets the absolute yield scale — must be exact |
| `SimThickConvert` | `sim2.c:2349` | thickness-unit conversion |
| `RbsBuffers[0]` (`ALTBUF`) | `rumpdata.h:156` | the theory buffer |
| `RbsActiveBuf` (`ibuf`) | `rumpdata.h:151` | the data buffer, copied wholesale at `creatr.c:283` |
| `Rmp` | `rumpdata.h:201` | alias for `RbsDataBlock`; only `autsim` is read |
| `sigtab` / `coffe2` | `sim2.c:81-82` | **default `{-1,-1}`**; zero silently selects the manual-override branch and returns a cross-section of exactly zero |

Captured bricks carry the *scattered particle's* `z`/`mass` (for the stopper-foil
lookup at `anlyz.c:180`), **not** the target's — target identity is implicit in
the block ordering, one block per isotope, heaviest first.
