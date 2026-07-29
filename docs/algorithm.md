# How the simulation works

A description of the forward model pyRUMP implements — what it computes, in what
order, and where it approximates. Written for someone who wants to know what the
numbers mean, not for someone reading the source.

The algorithm originates with L. R. Doolittle's 1985 paper. **pyRUMP follows the
shipped C, which diverges from that paper in places** — see
[Divergences from the paper](#divergences-from-the-published-algorithm).

---

## The physical problem

A beam of light ions (usually 1–3 MeV ⁴He or ¹H) hits a sample. A few
backscatter from target nuclei and reach a detector at a fixed angle. The
detected energy encodes two things at once:

- **which nucleus** it bounced off — heavier targets take less energy from the
  projectile, so each element has a characteristic maximum energy
- **how deep** the collision happened — the ion loses energy travelling in and
  out again, so deeper collisions arrive lower in energy

A spectrum is therefore a depth profile smeared together across all elements
present. Simulating it forward, and adjusting the sample until the simulation
matches, is how the depth profile is recovered.

### Kinematics

An elastic collision at scattering angle *φ* leaves the projectile with a
fixed fraction of its energy, the **kinematic factor**:

$$K = \left[\frac{\sqrt{1 - (x\sin\phi)^2} + x\cos\phi}{1+x}\right]^2,
\qquad x = \frac{m_1}{m_2}$$

For 2 MeV ⁴He at 170°, *K* is 0.566 on silicon and 0.922 on gold. Those factors
place the **surface edge** of each element: the highest energy at which it can
appear.

> **Angle convention.** RUMP stores `phi` as *180° minus* the scattering angle —
> a detector at 170° is entered as `phi = 10`. pyRUMP keeps that convention on
> `Geometry.phi` and exposes the physical angle as `Geometry.scattering_angle`.

### Cross-section

How often a collision happens is Rutherford scattering, in the lab frame:

$$\sigma = \left(\frac{Z_1 Z_2 e^2}{4E}\right)^2 \frac{4}{\sin^4\phi}
  \frac{\left[\sqrt{1-(x\sin\phi)^2}+\cos\phi\right]^2}{\sqrt{1-(x\sin\phi)^2}}$$

The **Z₂²** dependence is why RBS is sensitive to heavy elements in a light
matrix and nearly blind the other way round. The **1/E²** means yield rises with
depth as the beam slows.

A screening correction (L'Ecuyer) reduces this slightly at low energy. RUMP has
no relativistic correction and does not implement Andersen screening.

---

## The brick

The central data structure, in the 1985 paper's own words:

> Each simulated spectrum is made up of the superimposed contributions from each
> isotope of each sublayer in the sample. Any such contribution will be referred
> to as a **brick**.

Every layer is cut into **sublayers**, and each sublayer contributes one brick
per isotope of every element present. A brick is a trapezoid in energy space:

```
        h_front  ___________
                /           \___
               /                \  h_back
              |                  |
        e_back                e_front     (energy increasing to the right)
```

- **e_front** — energy of particles scattered from the sublayer's *front* face
- **e_back** — from its *back* face, lower because the beam travelled further
- **h_front, h_back** — differential yield at each edge

A simulated spectrum is nothing more than every brick integrated onto the
detector's channels and summed. Natural silicon gives three bricks per sublayer
(²⁸Si, ²⁹Si, ³⁰Si), each with its own kinematic factor, so isotopes appear as
slightly displaced copies.

### Why sublayers, and how many

Thicker sublayers mean fewer bricks and a faster simulation, but a coarser
approximation. RUMP's step size is **path-length based, not depth based**:

```
maxpath        = maxpth / max(|sec θ_in|, |sec θ_out|)
n_sublayers    = 1 + areal_thickness / maxpath
```

so tilting the sample automatically produces more, thinner sublayers. `maxpth`
defaults to 200 (in 10¹⁵ atoms/cm²) — roughly where the paper's error analysis
puts the energy-loss expansion at 10⁻⁵ fractional error.

A layer carrying a depth-profile equation ignores `maxpth` and uses a
per-equation recommended count instead (5 for `Constant`, 30 for `Thinfilm`).

---

## The pipeline

Stage order matters and is not arbitrary.

```
 stopping tables                fitted once per beam, then never recomputed
        │
        ▼
 slab discretization            layers → sublayers, composition per slab
        │
        ▼
 inbound march                  beam energy at every interface
        │
        ▼
 per-isotope depth loop         → bricks (the expensive part)
        │
        ▼
 brick → channel fill           trapezoid, or triangles if straggling
        │
        ▼
 detector convolution           Gaussian, applied to the whole spectrum
        │
        ▼
 yield normalisation            × Ω·Q/(charge state · CORR)
        │
        ▼
 pile-up  →  multiple scattering
```

### 1. Stopping powers — the indirection that matters most

Energy loss per unit depth, *ε(E)*, comes from one of several models tried in
priority order:

1. **Konac/Kalbitzer** fits, for the specific ion/target pairs in `newstop.kal`
   (H, D, ³He, ⁴He on carbon and silicon)
2. **Ziegler ZBL85**, the general fallback for Z = 1…92
3. hard-coded polynomials for Mylar, the Z=93 pseudo-element

Note the order: **for ⁴He on silicon — the most common RBS measurement there is
— Ziegler is not used.** Konac wins, and the two differ by up to ~10%.

Then the step that surprises everyone:

> **RUMP never evaluates the stopping model during a simulation.** At the start
> of a run it samples the chosen model at 201 points and least-squares fits a
> degree-5 polynomial in √E. Every subsequent energy-loss calculation evaluates
> only that polynomial.

The fit window is tied to the beam energy — `[0.04·E₀, 1.15·E₀]` — so
coefficients change when the beam does, and the polynomial diverges rapidly
outside that range. Any reimplementation that calls the stopping model directly
disagrees with RUMP everywhere.

Per sublayer the elemental polynomials are combined by **Bragg's rule** —
linear additivity weighted by areal density, with no compound correction — which
also folds in the thickness, so the coefficients directly give eV through the
slab.

### 2. Inbound march

The beam energy at each interface is computed once, before any element is
considered, because the incoming path does not depend on what it eventually
scatters from.

Energy loss across a sublayer uses a **third-order Taylor expansion** of
d*E*/d*a* = −*ε*(*E*):

$$E(a) = E_0 - a\varepsilon + \tfrac{1}{2}a^2\varepsilon\varepsilon'
        - \tfrac{1}{6}a^3\left(\varepsilon''\varepsilon^2 + \varepsilon'^2\varepsilon\right)$$

Truncating after the first term would be the familiar *surface approximation*.
The extra terms are what let RUMP use thick sublayers and stay fast — that is
the "rapid" in the 1985 paper's title. It is not a faster inner loop; it is
higher accuracy per sublayer, so that fewer are needed.

A **cutoff energy** (3% of the beam energy) ends the march: below it the
stopping fit is not trustworthy.

### 3. Outbound path — the expensive half

Once scattered, a particle must be walked back out through every overlying
sublayer. Done for every sublayer and every isotope, this is the O(*N*²) term
that dominates the cost, and the reason RUMP works so hard to keep *N* small.

Two things come out of it: the exit energy, and `ratde` — the accumulated ratio
of stopping powers, which accounts for the scattered beam's energy spread
changing on the way out.

### 4. Yield

Per sublayer and isotope:

$$h = \sigma(E)\, N_{\text{slab}}\, \frac{\sec\theta_{in}}{[\varepsilon]}
      \cdot \texttt{ratde} \cdot f_{\text{isotope}}$$

where **[ε]** is the stopping cross-section factor of Chu et al.,
$[\varepsilon] = K\varepsilon(E)\sec\theta_{in} + \varepsilon(KE)\sec\theta_{out}$,
which converts a depth interval into an energy interval.

### 5. Filling channels

Without straggling, each trapezoid is integrated exactly onto the channels it
overlaps, including partial channels at both ends.

With straggling, the trapezoid is abandoned: it is split into **two triangles**,
each convolved with a Gaussian of its own width. Straggling is **Bohr only** and
**off by default** in RUMP.

### 6. Detector resolution

A Gaussian convolution applied **once to the finished spectrum**, not per
sublayer — since 1994 the detector width was deliberately removed from the
straggling path. The kernel is channel-*integrated* (differences of the normal
CDF, not point samples) and truncated at 3σ.

Counts falling outside the channel range are discarded at both ends, so the
convolution is not count-conserving near the edges. Defensible — a real
multichannel analyser cannot record counts in channels it does not have — but it
matters when comparing totals.

### 7. After the spectrum exists

- **Pile-up** — two events arriving within the detector's shaping time recorded
  as one of their combined energy. Needs the beam current and shaping time.
- **Multiple scattering** — an empirical low-energy tail with **no physical
  basis**; the C's own comment calls its scale factor "ad-hoc". Treat any number
  it produces as qualitative.
- **Fuzz** — surface roughness, as several simulations at Gaussian-weighted
  thicknesses. Every roughened layer multiplies the cost.

---

## What is approximated, and what is absent

Worth knowing before trusting a result.

| | |
|---|---|
| Straggling | Bohr only, **off by default**. The in/out paths are combined *linearly* rather than as `K²σ²_in + σ²_out` — an approximation, not the correct combination. No Chu correction. |
| Screening | L'Ecuyer only; no Andersen. |
| Relativity | No correction anywhere. |
| Multiple scattering | Empirical tail with no physical basis. |
| Channelling | Not modelled at all. |
| Nuclear reactions | Q ≠ 0 reactions are rejected outright. |
| Sample | Laterally uniform apart from the `Fuzz` roughness model. |
| Beam and detector | Treated as points — no finite spot or acceptance angle. |

The 1985 paper is explicit about most of this:

> The algorithms assume a laterally uniform sample and neglect the effects of
> channelling, core electron screening, nuclear reactions, nuclear resonances,
> and multiple scattering. They also ignore the effects brought about by a
> finite size beam spot and detector.

---

## Divergences from the published algorithm

The shipped C is not the program the papers describe.

**The parabolic brick is gone.** The 1985 paper's signature contribution — a
parabolic brick top with an exact analytic area from the integral ∫E(a)⁻²da —
sits behind `#if 0` in the C, labelled *"no longer using Doolittle qqq code"*.
The shipped path uses plain trapezoids. The integral is still computed on every
run and then discarded.

pyRUMP follows the C. It also carries the discarded integral on
`Bricks.area`, so the paper's algorithm could be offered later as an opt-in.

**Konac/Kalbitzer stopping** is not in the 1996 manual at all, yet it takes
priority over Ziegler for the most common measurements.

**Defects are reproduced deliberately.** pyRUMP's default is bug-for-bug
fidelity, because every published RUMP result was produced with those bugs
present. The catalogue is in [rump-quirks.md](rump-quirks.md) — 20 entries.
The one that affects results most directly is a wrong coefficient index in the
second-derivative macro (35–50% error in d²ε/dE², feeding the third-order
energy-loss term). `StoppingTable.derivative(..., faithful=False)` gives the
correct value.

The single deliberate departure is `SimStragf`, RUMP's seven-regime rational
approximation to the triangle⊗Gaussian integral. pyRUMP uses a closed `erf`
form instead: measured against numerical quadrature, the closed form is exact to
6e-14 where RUMP's fit carries 1.7e-6. The 1985 justification for the
approximation was explicitly about 1985 hardware.

---

## Accuracy

pyRUMP is validated against the original C stage by stage. The legacy code is
compiled into a shared library and called directly from the test suite, so each
quantity is compared where it is produced rather than only at the end.

| Quantity | Agreement |
|---|---|
| ZBL85 stopping, 92 targets × 4 beam species | 6.1e-7 |
| Fitted stopping polynomial | 1.3e-5 |
| Cross-sections, kinematics | 1e-10 |
| Bricks, 630 across 36 configurations | 5e-7 energies, 6e-6 heights |
| Full spectrum, total counts | 3e-6 |
| Poisson fitting objective | 1e-5 |
| `.RBS` files, read and written | bit-identical |

Every one of these is limited by **float32 storage in the C**, not by pyRUMP.
The tolerances are argued from that floor rather than chosen for convenience.

For a realistic multi-layer sample with micron-thick polymer layers, agreement
loosens to ~3e-3 — thousands of sublayers accumulate single-precision
differences, and part of the beam falls below the stopping cutoff.

---

## References

- L. R. Doolittle, *Algorithms for the rapid simulation of Rutherford
  backscattering spectra*, Nucl. Instr. Meth. **B9** (1985) 344–351.
- L. R. Doolittle, *A new approach to Rutherford backscattering analysis*,
  Nucl. Instr. Meth. **B15** (1986) 227–231.
- W.-K. Chu, J. W. Mayer, M.-A. Nicolet, *Backscattering Spectrometry*,
  Academic Press (1978) — the source of the [ε] factor and the kinematics.
- J. F. Ziegler, J. P. Biersack, U. Littmark, *The Stopping and Range of Ions in
  Solids*, Pergamon (1985).
- S. Baker, R. D. Cousins, Nucl. Instr. Meth. **221** (1984) 437 — the fitting
  objective.
