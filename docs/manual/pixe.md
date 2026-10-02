## PIXE — particle-induced X-ray emission

!!! warning "In development for pyRUMP 2.0"
    This page is the working specification for PIXE support. It records the
    design as agreed so far, and will become the user manual once PIXE
    ships. Command names and details may still change.

    **Working so far:** the PIXE prompt with `GET`, `PAIR`, `DISABLE` and
    `SHOW`; the detector and calibration settings with their built-in
    defaults; reading Oxford `.PIX` files; and the PIXE window showing the
    measured spectrum. **Not yet:** the simulation, line markers, the
    bare-substrate background and fitting.

PIXE adds a second measurement to the one pyRUMP already models. The beam
that produces the RBS spectrum also ionises inner shells in the sample, and
an X-ray detector records the characteristic lines. pyRUMP simulates those
lines from the **same SIM sample** that the RBS simulation uses, so one
sample description is checked against both spectra.

Scope of the first version:

* **Beams:** protons and ⁴He at a few MeV.
* **Samples:** thin films (one or more layers) on a substrate, as described
  in [SIM](sim.md). The film peaks are simulated; the substrate and the
  continuum come from a measured **bare-substrate spectrum**.
* **Lines:** K, L and M lines of the film elements.
* **Absolute scale:** an instrumental constant **H** for each shell (K, L,
  M), as in GUPIX.

### Turning PIXE on and off

`PIXE` opens its own `PIXE Command:` prompt, like [SIM](sim.md) and
[PERT](pert.md). Entering it the first time **enables** PIXE: from then on
the PIXE spectrum and its simulation are kept up to date alongside the RBS
ones, and the PIXE plot window appears. `RETURN` goes back to the RUMP
level and leaves PIXE enabled.

`DISABLE` turns PIXE off and leaves the PIXE prompt in one step. The PIXE
window closes, and the settings are kept for next time. Entering `PIXE`
again enables it again.

From the RUMP level, `PIXE <command>` runs a single PIXE command without
entering the prompt.

```
Your wish? pixe
PIXE Command: get MnPt.PIX
PIXE Command: show
PIXE Command: return
Your wish? compare                     /* RBS window and PIXE window */
...
Your wish? pixe disable                /* back to RBS only           */
```

### Two plot windows

PIXE draws in **its own matplotlib window**, separate from the RBS window,
so the two can be sized and placed independently (for example one on each
screen). The RBS window behaves exactly as without PIXE.

The PIXE window shows:

* the measured PIXE spectrum, against energy in keV, on a logarithmic yield
  axis by default;
* the simulation: film peaks plus the background;
* optionally the separate contributions: each element's peaks, the
  background, the Si escape peaks;
* labelled markers for the lines (Mn Kα, Pt Mα, Ru Lα, …), like `ELEMENT`'s
  edge ticks;
* with `COMPARE`, a residuals panel below the spectrum.

Inside the PIXE prompt, `REGION`, `COUNTS`, `LOG`/`LINEAR`/`SQRT` and
`FIGSAVE` act on the PIXE window (with `REGION` in keV). At the RUMP level
they act on the RBS window, as before.

### PIXE data and buffers

A PIXE spectrum lives **in a buffer**, next to that buffer's RBS spectrum.
A buffer may hold an RBS spectrum, a PIXE spectrum, or both.

* **`GET <file>`** in the PIXE prompt reads a PIXE spectrum into the active
  buffer. With no buffer active, it creates a new one that holds only the
  PIXE spectrum.
* The beam (ion, energy), the beam tilt and the **charge** come from the
  buffer's RBS part. RC43 records both spectra in one run with one charge
  integrator, so the charge in the `.RBS` file is the charge for the `.PIX`
  file too. A buffer with no RBS part takes these from the defaults —
  `MEV`, `BEAM`, `THETA`, `CHARGE` as set in `~/.pyrumprc`, or pyRUMP's own
  hard-coded defaults.
* The PIXE spectrum keeps its own **live and real time**. The dead-time
  factor LT/RT is applied to the PIXE yield.
* **`PAIR ON|OFF`** (off by default): with `PAIR ON`, reading `x.RBS` with
  `GET` or `XEQ` also looks for `x.PIX` in the same folder (in any upper or
  lower case) and loads it into the same buffer.

**The bare-substrate spectrum is just another buffer.** Measure the bare
substrate under the same conditions, load it like any other measurement
(with `PAIR ON`, `XEQ bare.RBS` brings its `.PIX` along and so its charge),
then name it with `BARE <buffer>`. Its PIXE spectrum, scaled by charge and
dead time, becomes the background of every PIXE simulation.

### File format: Oxford Instruments `.PIX`

The PIXE reader recognises the file by its first line, `Oxford Instruments,
Inc`, not by its extension. A typical header:

```
Oxford Instruments, Inc
ID:           MnPt.PIX  PIXE X-Ray  LT =  951.132 RT  953.434 Gain  61.6
Acquisition Date:           08-19-2026 16:43:34
Elapsed Real Time:           953.434
Elapsed Live Time:           951.132
Conversion Gain:             2048 Calibration
High Voltage:   -120        A =            1.009699E-02
Coarse Gain:   44           B =           -3.864563E-02
Fine Gain:   1.10           C = -1.263E-004

Channel   Contents
 1             0
 2             0
 ...
```

* **Channels** are read with the numbers the file gives them. RC43 numbers
  them from 1; other Oxford exports start at 0.
* **Calibration:** in RC43 files, `A` is the gain in keV per channel and `B`
  the offset in keV. Other exports use `A` and `B` differently, so the
  header calibration is used only when it looks plausible (a gain of about
  1–50 eV per channel), and only as the **starting value**: gain and offset
  are refined in the fit. `C` is ignored — it is the same constant in every
  file seen so far, including a 1998 export, so it is not a calibration
  term.
* **Live and real time** give the dead-time correction.
* The file has **no charge**; see above.
* Header labels vary between exports (`Coarse`/`Course`, tabs or spaces,
  CRLF line ends), so the reader keys on the label text, not on positions.

### Detector setup

The X-ray detector is described in the PIXE prompt. Like any other setting,
it can be made permanent in `~/.pyrumprc`. End the block with `DISABLE`
to store the setup without turning PIXE on at startup, or with `RETURN` to
start every session with PIXE enabled:

```
pixe
 angle 45                  ! detector axis to sample normal, degrees
 solid 1.2                 ! solid angle, msr
 window Be 8               ! window material and thickness, µm
 crystal Si 450            ! crystal thickness, µm
 fwhm 131                  ! resolution at Mn Kα, eV
 pair on
disable
```

(The values above are placeholders, not a real detector.) One-shot
`PIXE <command>` lines work in `~/.pyrumprc` too, and never turn PIXE on.

`SHOW` prints the current setup as commands in this form, ready to paste
into `~/.pyrumprc`.

**Built-in defaults.** Without any setup, pyRUMP uses these — like the RBS
defaults, something sensible to start from. Where a value could be taken
from a real measurement, it comes from the MnPt example's SDD:

| Setting | Default | Source |
|---|---|---|
| `ANGLE` | 45° | typical |
| `SOLID` | 1 msr | typical |
| `WINDOW` | Be 8 µm | typical |
| `CRYSTAL` | Si 450 µm | typical |
| `FWHM` | 131 eV at Mn Kα | fitted to `MnPt.PIX` |
| `FANO` | 0.13 | fitted to `MnPt.PIX` (with Si Kα's 79 eV) |
| `CALIB` | 0.01009699 keV/ch, −0.03864563 keV | `MnPt.PIX`'s RC43 header |
| `PAIR` | off | |

`CALIB` is the fallback for a spectrum whose header has no usable
calibration. A spectrum's own header calibration, when it has one, is used
for that spectrum; `CALIB` given after `GET` changes the active buffer's
PIXE calibration too.

### Commands (draft)

| Group | Commands |
|---|---|
Commands marked † are not implemented yet.

| Group | Commands |
|---|---|
| Prompt | `RETURN`, `DISABLE`, `SHOW`, `HELP` |
| Data | `GET <file>`, `PAIR ON\|OFF`, `BARE <buffer>` †, `LIVETIME <live> <real>` † |
| Geometry | `ANGLE <deg>` (−90 to 90), `SOLID <msr>` |
| Detector | `WINDOW <element> <µm>`, `CRYSTAL <element> <µm>`, `FWHM <eV>`, `FANO <F>`, `DEADLAYER <µm>` †, `ESCAPE ON\|OFF` †, `TAIL ON\|OFF` † |
| Absorbers | `FILTER <element> <µm> [<hole %>]`, `FILTER CLEAR` |
| Calibration | `CALIB <gain keV/ch> <offset keV>`, `H K\|L\|M <value>` † |
| Elements | `LINES` † (table of lines, energies, cross sections and yields), `EXCLUDE <element>` †, `INCLUDE <element>` † |
| Background | `BGSCALE <s>` †, `SMOOTH <channels>` † |
| Plot | `PLOT`, `REGION <keV> <keV>` / `REGION ALL`, `COUNTS <low> <high>` / `COUNTS ALL`, `LOG`/`LINEAR`/`SQRT`, `FIGSAVE <file>`, `MARKERS ON\|OFF` †, `COMPONENTS ON\|OFF` † |
| Output | `EXPORT <file>` † |

Absorbers are single elements for now; compounds (Mylar, Kapton) will
follow.

### Physics

#### Ionisation cross sections

Inner-shell ionisation cross sections σ_i(E) come from **ECPSSR** (Brandt
& Lapicki): the plane-wave Born approximation corrected for energy loss,
Coulomb deflection, perturbed stationary states (binding and polarisation)
and relativistic effects. pyRUMP ships tables computed by its own
independent implementation — no ISICS or GUPIX code or data — for
protons (0.1–5 MeV) and ⁴He (0.2–12 MeV), K shell for Z = 6–92 and L1–L3
for Z = 18–92. **M shells are still to be added**; they are essential for
heavy elements under ⁴He, where the M lines are by far the strongest (Pt
Mα in the MnPt example).

The correction terms were checked term by term against Geant4's
independent ECPSSR code. Against ISICS values, the tables agree within ±5 %
for most cases from 0.5 MeV (H) and 1 MeV (He) up; see the limitations
below for where they don't. Values between grid points are interpolated in
ln σ versus ln E, never extrapolated.

#### From vacancies to X-rays

Fluorescence yields ω, Coster–Kronig probabilities f, radiative rates F,
line energies and absorption edges come from
[xraylib](https://github.com/tschoonj/xraylib). Vacancies are carried down
the shell by Coster–Kronig transitions before they decay. For the L shell:

$$n_1 = \sigma_{L1},\qquad n_2 = \sigma_{L2} + f_{12}\,n_1,\qquad
n_3 = \sigma_{L3} + f_{13}\,n_1 + f_{23}\,n_2$$

and likewise for M1–M5. A line ℓ that fills a vacancy in subshell i has
the production cross section σ_ℓ = n_i ω_i F_iℓ. (The direct f₁₃ is used,
not xraylib's combined f′₁₃, which would count the transfer twice.)

#### Yield

The film is cut into the same kind of sublayers as the RBS simulation, and
the beam's energy E_k in each one comes from the same stopping powers. The
yield of line ℓ of element Z is

$$Y_\ell = H\,N_\text{ion}\,\frac{\Omega}{4\pi}\,\varepsilon(E_\ell)\,
T_\text{abs}(E_\ell)\sum_k \frac{c_{Z,k}\,\Delta m_k}{\cos\alpha}\,
\sigma_\ell(E_k)\,A_{k,\ell}$$

with the absorption of the X-rays on their way out

$$A_{k,\ell} = \exp\!\Big(-\sum_{j<k}\frac{\mu_j\,\rho t_j}{\cos\theta_d}\Big)
\cdot\frac{1-e^{-x_k}}{x_k},\qquad x_k = \frac{\mu_k\,\rho t_k}{\cos\theta_d}$$

where α is the beam's angle to the sample normal, θ_d the detector's,
Δm_k the sublayer's areal density, c_Z,k the element's atomic fraction, μ
the mass attenuation coefficient and ε the detector efficiency. For truly
thin films A → 1 and E_k → E₀, and this is the plain thin-film formula;
summing over sublayers costs nothing and stays right when a cap layer
absorbs. H is the instrumental constant for the line's shell (K, L or M).
It absorbs the solid angle, charge calibration and database errors, and is
best determined from standards — including films whose amounts RBS has
already measured.

The substrate (the last SIM layer) is not part of this sum: its peaks are
in the bare-substrate spectrum.

#### Detector

* **Resolution:** σ² = σ²_noise + ε_Si F E with ε_Si = 3.64 eV, set by the
  FWHM at Mn Kα and the Fano factor F. In the MnPt example, Si Kα is 79 eV
  wide and the 5.9 keV line 131 eV, which this model reproduces with F ≈
  0.13.
* **Line shape:** a Gaussian per line, with an optional low-energy tail.
* **Si escape peaks** 1.740 keV below each line above the Si K edge.
* **Efficiency:** window, contact and dead layer transmission times the
  crystal's absorption, or an imported efficiency curve.

#### Background

$$M(E) = \sum_\ell Y_\ell\,P_\ell(E) \;+\; s\,\frac{Q}{Q_\text{bare}}\,
B_\text{bare}(E)\,T_\text{film}(E)$$

The bare-substrate spectrum B_bare is scaled by the ratio of charges (and
dead-time factors), and by the transmission T_film of the film above the
substrate. The free scale s, close to 1, absorbs the slightly lower beam
energy reaching the substrate and small charge errors. The bare spectrum
must be measured with the same beam, geometry, detector and filters.

### Fitting

PERT will vary PIXE quantities alongside the RBS ones: element amounts
(through the shared sample), the PIXE energy calibration, the resolution,
the background scale s and the H values. A fit can use the RBS spectrum,
the PIXE spectrum, or both at once.

### Known limitations

* **Light elements (C, N, O):** ECPSSR is least reliable here. For carbon
  the tables run from about 23 % above to 16 % below ISICS across the He
  energy range. Mass attenuation coefficients and detector efficiency below
  1 keV are also uncertain, so light-element results depend on H from
  standards.
* **M shells** carry larger cross-section uncertainties than K and L, again
  absorbed into H.
* **He beams on light elements** cause multiple ionisation, which shifts
  and broadens the K lines; this is not modelled.
* Films are laterally uniform and flat with sharp interfaces — the same
  assumptions as the RBS simulation.

### References

* W. Brandt, G. Lapicki, Phys. Rev. A 20 (1979) 465; Phys. Rev. A 23 (1981)
  1717 — ECPSSR theory.
* J.A. Maxwell, J.L. Campbell, W.J. Teesdale, Nucl. Instr. Meth. B 43 (1989)
  218; J.A. Maxwell, W.J. Teesdale, J.L. Campbell, Nucl. Instr. Meth. B 95
  (1995) 407 — the GUPIX spectrum model and layered-target yields.
* J.L. Campbell et al., Nucl. Instr. Meth. B 499 (2021) 77 — SDD line
  shapes, H(E) calibration.
* T. Schoonjans et al., Spectrochim. Acta B 66 (2011) 776 — xraylib.
