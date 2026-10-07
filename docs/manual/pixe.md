## PIXE — particle-induced X-ray emission

!!! note "New in pyRUMP 2.0"
    pyRUMP simulates the characteristic K, L and M lines of the SIM sample
    in a PIXE window beside the RBS one, and [PERT](pert.md#pixwin-new) fits
    compositions to the PIXE spectrum while the rest is fitted to RBS (see
    [Fitting](#fitting)). Not modelled yet: the continuum background, peak
    tails and pile-up, so fit windows go on clear peaks. Commands marked †
    in the [command summary](#commands) are planned for a later version.

PIXE adds a second measurement to the one pyRUMP already models. The beam
that produces the RBS spectrum also ionises inner shells in the sample, and
an X-ray detector records the characteristic lines. pyRUMP simulates those
lines from the **same SIM sample** that the RBS simulation uses, so one
sample description is checked against both spectra.

Scope of the first version:

* **Beams:** protons and ⁴He at a few MeV.
* **Samples:** thin films (one or more layers) on a substrate, as described
  in [SIM](sim.md). The peaks of the films and the substrate are
  simulated; the continuum is still to come — fitted, or taken from a
  measured reference spectrum.
* **Lines:** K, L and M lines of the film elements.
* **Absolute scale:** an instrumental constant **H** for each shell (K, L,
  M), as in GUPIX.

### Turning PIXE on and off

`PIXE` opens its own `PIXE Command:` prompt, like [SIM](sim.md) and
[PERT](pert.md). Entering it the first time **enables** PIXE: from then on
every RBS plot command shows PIXE too (see the next section). Entering
opens no window by itself. `RETURN` goes back to the RUMP level and leaves
PIXE enabled.

`DISABLE` turns PIXE off and leaves the PIXE prompt in one step. The PIXE
window closes, and the settings and the buffers' PIXE spectra are kept for
next time. Entering `PIXE` again enables it again.

From the RUMP level, `PIXE <command>` runs a single PIXE command without
entering the prompt, and never turns PIXE on.

```
Your wish? pixe pair on                /* x.RBS brings x.PIX along   */
Your wish? xeq MnPt.RBS                /* both spectra, no window yet */
Your wish? sim get MnPt.lcm
Your wish? pixe                        /* PIXE on                    */
PIXE Command: return
Your wish? compare                     /* RBS and PIXE comparisons   */
Your wish? figsave fit                 /* fit_rbs.png, fit_pixe.png  */
...
Your wish? pixe disable                /* back to RBS only           */
```

### Two plot windows

PIXE draws in **its own matplotlib window**, separate from the RBS window,
so the two can be sized and placed independently (for example one on each
screen). The RBS window behaves exactly as without PIXE.

**While PIXE is enabled, the PIXE window mirrors the RBS window** — the same
buffers, the PIXE spectrum where the RBS window shows the RBS one:

| At the RUMP level | RBS window | PIXE window |
|---|---|---|
| `PLOT 1` | buffer 1 | buffer 1's PIXE spectrum |
| `OVERLAY 2` | adds buffer 2 | adds buffer 2's PIXE spectrum |
| `PLOT 0` | the simulation | the PIXE simulation |
| `SIM SPLOT Mn`, `SIM SPLOT 2` | one element's or layer's part | the same part of the PIXE simulation |
| `COMPARE` / `CMP` | data, simulation, residuals | the same, for PIXE |
| `REPLOT` | redraw | redraw |

A buffer with no PIXE spectrum is named in the PIXE window ("buffer 2: no
PIXE spectrum") rather than left out silently. Reading data (`GET`, `XEQ`
with `PAIR ON`, `PIXE GET`), changing a PIXE setting and — with
[`LIVE`](plotting.md#live-new) on — changing the sample only update a PIXE
window that is already open; they never open one.

Inside the PIXE prompt, `PLOT [buffer]` (`PLOT 0`: the simulation) and
`COMPARE` draw in the PIXE window alone and leave the RBS window as it is.

The PIXE window shows:

* the spectra against **energy** (bottom axis) with the **channel** numbers
  on top, on a logarithmic yield axis by default;
* the simulation of the SIM sample's film lines; the background is still
  to come;
* markers for the lines of every element in the SIM sample (Mn Kα, Pt Mα,
  Ru Lα, …): a tick at each line's energy hanging from the top edge, led to
  a vertical label. Where lines crowd, the labels are spread apart and the
  ticks stay at the true energies; lines too close to tell apart share one
  label ("Pt Mα/Mβ"). The yield axis leaves room above the tallest peak for
  them. `MARKERS ON` (the default) marks Kα, Kβ, Lα, Lβ1, Lβ2, Lγ1, Mα and
  Mβ; `MARKERS ALL` adds Ll, Mζ and Mγ; `MARKERS OFF` hides them. RBS
  absorber layers are left out — the beam never reaches them. To mark an
  element that isn't in the sample, use
  [`ELEMENT`](analysis.md#element) at the RUMP prompt: it lists the
  element's lines and labels them here;
* with `COMPARE`, the Poisson residuals below the spectrum and the reduced
  chi-square over the channels shown, as in the RBS comparison.

Inside the PIXE prompt, `REGION`, `COUNTS`, `LOG`/`LINEAR`/`SQRT` and
`FIGSAVE` act on the PIXE window; at the RUMP level they act on the RBS
window, as before. `REGION` is in **channels**, as the spectrum file
numbers them, like the RBS window's: `REGION 220 1190`.

**FIGSAVE.** At the RUMP level, `FIGSAVE fit` saves the RBS window as
`fit_rbs.png` — and with `PAIR ON` and the PIXE window open, the PIXE window
next to it as `fit_pixe.png`. Inside the PIXE prompt, `FIGSAVE` saves the
PIXE window alone, under the name given.

### PIXE data and buffers

A PIXE spectrum lives **in a buffer**, next to that buffer's RBS spectrum.
A buffer may hold an RBS spectrum, a PIXE spectrum, or both.

* **`GET <file>`** in the PIXE prompt reads a PIXE spectrum into the active
  buffer. With no buffer active, it creates a new one that holds only the
  PIXE spectrum.
* The beam (ion, energy) and the **charge** come from the buffer's RBS
  part, and so does the sample tilt with `PAIR ON` or `THETA RBS` (see
  Geometry below; otherwise PIXE has its own `THETA`, 0 by default). RC43
  records both spectra in one run with one charge integrator, so the
  charge in the `.RBS` file is the charge for the `.PIX` file too. A buffer with no RBS part takes these from the defaults —
  `MEV`, `BEAM`, `THETA`, `CHARGE` as set in `~/.pyrumprc`, or pyRUMP's own
  hard-coded defaults.
* **`CORRECTION` applies to PIXE too.** It divides the dose for both
  spectra: N_ion = `CHARGE` / (e × charge state × `CORRECTION`). Changing
  CORR rescales the PIXE simulation as well as the RBS one. That includes a
  CORR typed by hand, CORR varied in a PERT fit, and the CORR a PERT
  normalisation window writes back. This is right when CORR corrects the
  charge (integrator error, secondary electrons, beam on the frame): the
  same ions reach both detectors. It is wrong when CORR stands in for
  something only RBS has, such as RBS dead time (pyRUMP's RBS simulation
  has no live-time factor of its own), an `OMEGA` error, or a cross-section
  or stopping-power error. PIXE would then be scaled by it as well. Keep
  CORR for the dose. PIXE `SHOW` prints the dose in force.
* The PIXE spectrum keeps its own **live and real time**. The dead-time
  factor LT/RT is applied to the PIXE yield.
* **`PAIR ON|OFF`** (off by default): with `PAIR ON`, reading `x.RBS` with
  `GET` or `XEQ` also looks for `x.PIX` in the same folder (in any upper or
  lower case) and loads it into the same buffer. The two spectra are then
  one run on one sample, so the PIXE tilt is always the RBS buffer's
  `THETA`, as under `THETA RBS`, and a `THETA` typed in the PIXE prompt is
  refused until `PAIR OFF`.

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

### Exports

`EXPORT` and `EXPORTCMP` (`EC`) write the PIXE spectrum, or the PIXE
comparison, as plain columns for Origin, Excel, gnuplot or pandas — always
to `<file>_pixe`, next to the RBS file of the same name:

| Where | `PAIR ON`, PIXE spectrum in the buffer | otherwise |
|---|---|---|
| RUMP level `EXPORT f`, `EC f` | `f.txt` (RBS) **and** `f_pixe.txt` (PIXE) | `f.txt` (RBS) |
| PIXE prompt `EXPORT f`, `EC f` | `f_pixe.txt` (PIXE only) | `f_pixe.txt` (PIXE only) |

`EXPORT`'s columns are `channel` (as the spectrum file numbers them),
`energy_keV` (lower edge of the channel), `counts` and `error` (√counts).
`EXPORTCMP` has `channel`, `energy_keV`, `counts`, `simulation`, `diff`,
`residual` (Poisson, in σ), then one column per element of the sample with
its simulated lines (`sim_Mn`, `sim_Pt`, …).

The `#` header records what is needed to read them: the PIXE file, live and
real time, calibration, beam, the beam and X-ray angles in force, the dose,
the whole PIXE setup (as `SHOW` prints it), and for `EXPORTCMP` the physics
settings, the sample and the reduced chi-square over the PIXE `REGION`.
A `.csv` file is comma-separated, anything
else tab-separated, and a bare filename gets `.txt`.

### Detector setup

The X-ray detector is described in the PIXE prompt. Like any other setting,
it can be made permanent in `~/.pyrumprc`. End the block with `DISABLE`
to store the setup without turning PIXE on at startup, or with `RETURN` to
start every session with PIXE enabled:

```
pixe
 geometry ibm              ! where the detector is: IBM, CORNELL or GENERAL
 theta 0                   ! sample tilt, degrees (0: normal incidence), or RBS
 phi 45                    ! detector to the beam (looking back), degrees
 psi 0                     ! detector to the sample normal: GENERAL only
 solid 25 7.125 in         ! 25 mm^2 active area, 7.125 in from the sample
 window Be 12.5            ! window material and thickness, µm
 crystal Si 500            ! crystal thickness, µm
 fwhm 122                  ! resolution at Mn Kα, eV
 fano 0.104                ! how the peak width grows with energy
 filter clear
 filter 1 mylar 62         ! absorber against bremsstrahlung, µm (effective)
 pair on
disable
```

One-shot `PIXE <command>` lines work in `~/.pyrumprc` too, and never turn
PIXE on. `SHOW` prints the current setup as commands in this form, ready to
paste into `~/.pyrumprc`, followed by the beam and X-ray angles in force,
and the dose, shared with RBS:

```
  ! dose, shared with RBS: CHARGE 10 uC / CORRECTION 0.8 = 12.5 uC, 7.802e+13 ions (charge state 1)
```

**Geometry.** The PIXE detector is placed the way the RBS detector is,
with `GEOMETRY`, `THETA`, `PHI` and `PSI` meaning the same (the
[Experimental geometry](geometry.md) page draws all three cases), but set
in the PIXE prompt, apart from the RBS buffer's: the two detectors stand in
different places.

* `PHI` is the detector's angle to the beam, measured as for RBS from the
  beam looking back at the source. An untilted sample sees the X-rays
  leave at `PHI` to its normal.
* `THETA` is the sample tilt: the beam comes in at |THETA| to the normal.
  `THETA RBS` takes the RBS buffer's `THETA` instead — one sample, tilted
  once for both detectors — and follows it while PERT fits `THETA`. A
  number unlinks it again. `PAIR ON` always takes the RBS buffer's, and
  refuses a number until `PAIR OFF`.
* `GEOMETRY` says where the X-rays leave, at PSI to the normal: `IBM`, the
  detector beside the beam in the plane the sample turns in, PSI =
  \|THETA + PHI\| (a negative `THETA` turns the sample towards the
  detector); `CORNELL`, cos PSI = cos THETA × cos PHI; `GENERAL`, `PSI` as
  typed.

The default, `GEOMETRY IBM`, `PHI 45`, `THETA 0`, has the beam at normal
incidence and the X-rays leaving at 45°. On the RC43 endstation the sample
is turned towards the PIXE detector by the RBS file's negative `THETA`:
with `PAIR ON` or `THETA RBS`, the MnPt example's `THETA -9` gives 9° in and
\|−9 + 45\| = 36° out. `SHOW`, and every one of these commands, prints
the two angles in force.

In the PIXE prompt, `THETA`, `PHI`, `PSI` and `GEOMETRY` are the PIXE
detector's; the RBS buffer's are set at the RUMP level, after `RETURN`.

**Resolution: `FWHM` and `FANO`.** A line's width has two parts:
electronic noise, the same for every line, and the statistics of the
charge the X-ray frees in the crystal, which grow with its energy E:

  FWHM(E)² = noise² + 2.355² × 3.64 eV × F × E

`FWHM` is the total width at Mn Kα (5.9 keV), in eV: the number a
detector's specification quotes. `FANO` is F, which splits that width
between the two parts. So a larger F makes the lines below 5.9 keV (Si K)
narrower and those above it (Pt L) wider. Si detectors have F ≈ 0.1. The
RC43 values, 122 eV and 0.104, come from the Si Kα (78.5 eV) and Mn Kα
(122 eV) widths of real spectra; F rarely needs changing. `HELP FWHM` and
`HELP FANO` in the PIXE prompt say the same.

**Solid angle.** `SOLID <msr>` sets it directly. `SOLID <area mm²>
<distance> [MM|IN]` computes it for a round detector seen on axis: 25 mm²
at 7.125 in (181 mm) is 0.763 msr.

**Absorbers.** `WINDOW` and `FILTER` take an element or a compound by
name: `MYLAR` (C₁₀H₈O₄, 1.40 g/cm³) or `KAPTON` (C₂₂H₁₀N₂O₅, 1.42 g/cm³),
NIST's compositions and densities, mixed by mass fraction. `CRYSTAL` takes
an element.

**Filters** are numbered from the sample outwards — filter 1 faces the
sample — and the X-rays cross them all:

```
PIXE Command: filter 1 mylar 50 hole 60%     /* sets (or replaces) filter 1 */
PIXE Command: filter 2 al 10 45%             /* a second one, once 1 exists */
PIXE Command: filter                         /* lists them, with the total
                                                transmission at 1.5-10 keV */
PIXE Command: filter clear 2                 /* removes filter 2            */
PIXE Command: filter clear                   /* removes them all            */
```

A filter can only be set next to the existing ones, so the numbering has no
gaps; removing one moves the ones after it up. A hole percentage (`HOLE`
and `%` are optional) makes a "funny filter": that fraction of the X-rays
passes unattenuated.

**Built-in defaults.** Without any setup, pyRUMP uses these — like the RBS
defaults, something sensible to start from. They describe a real setup:
the NEC RC43 endstation the MnPt example was measured on, with an Amptek
silicon drift detector.

| Setting | Default | Source |
|---|---|---|
| `GEOMETRY` | `IBM`, `PHI 45`, `THETA 0` (normal incidence) | RC43 endstation: the detector 45° from the beam |
| `SOLID` | 0.763 msr | 25 mm² at 7.125 in |
| `WINDOW` | Be 12.5 µm | Amptek SDD |
| `CRYSTAL` | Si 500 µm | Amptek FAST SDD (typical) |
| `FWHM` | 122 eV at Mn Kα | fitted: Mn Kα doublet in `MnPt.PIX` (specification: 130 eV) |
| `FANO` | 0.104 | fitted: Si Kα 78.5 eV and Mn Kα 122 eV |
| `FILTER` | Mylar 62 µm | effective thickness, fitted; see below (specified: 125 µm) |
| `CALIB` | 0.01009699 keV/ch, −0.03864563 keV | `MnPt.PIX`'s RC43 header |
| `PAIR` | off | |
| `MARKERS` | on | |

`CALIB` is the fallback for a spectrum whose header has no usable
calibration. A spectrum's own header calibration, when it has one, is used
for that spectrum; `CALIB` given after `GET` changes the active buffer's
PIXE calibration too.

**How the filter, resolution and calibration defaults were measured.** Two
spectra from this setup, four years apart, pin them down:

* *A 258 nm SiO₂-on-Si reference* (2022, 1.9 MeV ⁴He, 40 µC) has nothing
  but Si and O, so its thick-target Si K (1.74 keV) depends only on the
  filter. Through the specified 125 µm of Mylar, only 1.7×10⁻⁵ of it would
  arrive and the simulation would be 200 times too weak; the measured yield
  needs **63.0 µm**. `examples/MnPt` (2026) gives **60.8 µm** the same way.
  Both took the sample tilt from the RBS file (`THETA RBS`). The default
  is their mean, 62 µm — an *effective* thickness, which also
  takes up any error in the Si K cross section (H_K = 1 assumed). A 5 %
  change in it changes Si K by about 30 %, the higher-energy lines hardly.
* Si Kα is 78.5 eV wide in both spectra and Mn Kα, fitted as its Kα₁/Kα₂
  doublet, 122 ± 3 eV: a Fano factor of 0.104 with about 50 eV of
  electronic noise.
* The RC43 header calibration puts Si Kα within 0.4 eV of its true energy
  but Mn Kα 19 eV high — 0.46 % too much gain. Gain and offset need fitting;
  the header is a good start.

Its spectrum has nothing at 2.05 keV (0.6 % of Si K, from Si K's own tail),
so the 2.05 keV peak in MnPt (16.6 % of its Si K) comes from the sample —
the Pt Mα/Mβ pair.

### Commands

Commands marked † are planned for a later version.

| Group | Commands |
|---|---|
| Prompt | `RETURN`, `DISABLE`, `SHOW`, `HELP` |
| Data | `GET <file>`, `PAIR ON\|OFF`, `BARE <buffer>` †, `LIVETIME <live> <real>` † |
| Geometry | `GEOMETRY IBM\|CORNELL\|GENERAL`, `THETA <deg>\|RBS` (−90 to 90), `PHI <deg>` (0 to 180), `PSI <deg>` (`GENERAL` only), `SOLID <msr>` or `SOLID <area mm²> <distance> [MM\|IN]` |
| Detector | `WINDOW <element\|MYLAR\|KAPTON> <µm>`, `CRYSTAL <element> <µm>`, `FWHM <eV>`, `FANO <F>`, `DEADLAYER <µm>` †, `ESCAPE ON\|OFF`, `TAIL ON\|OFF` † |
| Absorbers | `FILTER` (list), `FILTER <n> <element\|MYLAR\|KAPTON> <µm> [[HOLE] <%>]`, `FILTER CLEAR [<n>]` |
| Calibration | `CALIB <gain keV/ch> <offset keV>`, `H <K> <L> <M>` or `H K\|L\|M <value>` |
| Simulation | `LINES [ALL]` (typed in full — `LIN` is `LINEAR`; table of lines, energies, cross sections, efficiency and counts), `EXCLUDE <element>` †, `INCLUDE <element>` † |
| Background | `BGSCALE <s>` †, `SMOOTH <channels>` † |
| Plot | `PLOT [buffer]`, `COMPARE` / `CMP`, `REGION <channel> <channel>` / `REGION ALL`, `COUNTS <low> <high>` / `COUNTS ALL`, `LOG`/`LINEAR`/`SQRT`, `MARKERS ON\|ALL\|OFF`, `FIGSAVE <file>`, `COMPONENTS ON\|OFF` † |
| Output | `EXPORT <file>`, `EXPORTCMP <file>` / `EC` (to `<file>_pixe`; at the RUMP level with `PAIR ON` as well) |

Other compounds (by formula and density) can be added to the list as
needed.

### Physics

#### Ionisation cross sections

Inner-shell ionisation cross sections σ_i(E) come from **ECPSSR** (Brandt
& Lapicki): the plane-wave Born approximation corrected for energy loss,
Coulomb deflection, perturbed stationary states (binding and polarisation)
and relativistic effects. pyRUMP ships tables computed by its own
independent implementation — no ISICS or GUPIX code or data — for
protons (0.1–5 MeV) and ⁴He (0.2–12 MeV), K shell for Z = 6–92 and L1–L3
for Z = 18–92 and M1–M5 for Z = 50–92, at 50 energies per decade — dense
enough that interpolation is good to 0.3 %. The generator is in
`tools/ecpssr/`, with its theory and validation.

The **M shells** matter most for heavy elements under ⁴He, where the M lines
are by far the strongest (Pt Mα in the MnPt example: about 20,000 b for M5
at 1.9 MeV, against 4.5 b for L3). They follow the M-shell treatment of
**ISICS2011** (ISICS version 5.1; Liu & Cipolla 1996, Cipolla 2007, 2009):
Slater-screened hydrogenic PWBA for 3s/3p/3d, the united-atom binding
correction ζ = U(Z₂+Z₁)/U(Z₂) — the subshell's binding energy in the
element Z₂+Z₁ over its own — and the same energy-loss and
Coulomb-deflection corrections, but no M-specific polarisation or
relativistic terms.

Compared with ISICS itself (versions 1.0 and 5.1, compiled locally, used for
validation only and not distributed) and given the same binding energies,
pyRUMP's tables agree within 1.1 % for M1–M5 and 0.1–1.5 % for K and L
(except Au K at the lowest velocities, where σ < 0.005 b).

The correction terms were checked term by term against Geant4's
independent ECPSSR code. Against ISICS values, the tables agree within ±5 %
for most cases from 0.5 MeV (H) and 1 MeV (He) up; see the limitations
below for where they don't. Values between grid points are interpolated in
ln σ versus ln E, never extrapolated.

#### From vacancies to X-rays

Fluorescence yields ω, Coster–Kronig probabilities f, radiative rates F,
line energies, absorption edges and mass attenuation coefficients come from
[xraylib](https://github.com/tschoonj/xraylib) 4.3.0. pyRUMP does not need
xraylib installed: `tools/pixe_atomic.py` generated the tables once, and
they ship with pyRUMP (`pyrump/data/pixe/`, with their provenance and
xraylib's licence in `data/SOURCES.md`), so results only change when the
tables are deliberately regenerated. The attenuation table is sampled so
that interpolating it reproduces xraylib within about 0.1 %, absorption
edges included. Vacancies are carried down
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

where α is the beam's and θ_d the X-rays' angle to the sample normal (see Geometry above),
Δm_k the sublayer's areal density, c_Z,k the element's atomic fraction, μ
the mass attenuation coefficient and ε the detector efficiency. For truly
thin films A → 1 and E_k → E₀, and this is the plain thin-film formula;
summing over sublayers costs nothing and stays right when a cap layer
absorbs. H is the instrumental constant for the line's shell (K, L or M).
It absorbs the solid angle, charge calibration and database errors, and is
best determined from standards — including films whose amounts RBS has
already measured. Fit the standard's `CORRECTION` on its RBS spectrum
**first**, then set H from its PIXE spectrum. That way the standard's
charge error stays in its CORR rather than in H. H then describes the PIXE
detector alone, and each sample's CORR carries only its own dose.

**The substrate** — the last SIM layer — is part of this sum like any
other layer, with the thickness SIM gives it, as far as the beam gets: where the beam has slowed below the cross-section tables (0.1 MeV
for protons, 0.2 MeV for ⁴He) the rest is dropped, its cross sections being
orders of magnitude down by then. For 1.9 MeV ⁴He in Si that is about 4 µm,
so a SIM substrate of 3–5 µm, as is usual for RBS, gives the full
thick-target yield; a thinner one gives less (MnPt.lcm's 1000 nm Si: about
half). To see the film's part alone, use `SIM SPLOT <layer>` or
`SIM SPLOT <element>`. When a measured bare-substrate background is in use
(still to come), it supplies the substrate's peaks and the simulation will
leave the substrate out by itself.

#### Detector

* **Resolution:** σ² = σ²_noise + ε_Si F E with ε_Si = 3.64 eV, set by the
  FWHM at Mn Kα and the Fano factor F (see Resolution under Detector setup).
  In the MnPt example, Si Kα is 78.5 eV wide and Mn Kα 122 eV, which this
  model reproduces with F = 0.104.
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

In PERT you choose, element by element, which spectrum fits it:
`COMPOSITION`, `ATOMS` and `SPECIES` take `PIXE` after the element
([`COMPOSITION 2 Fe PIXE`](pert.md#pixwin-new)) to be fitted to the PIXE
spectrum over the windows set with [`PIXWIN <first> <last>`](pert.md#pixwin-new),
while thicknesses and everything else are fitted to RBS. GO alternates the
two fits until neither moves. H stays as set here unless
[`PIXH`](pert.md#pixh-new) varies it, for the shells of the PIXE elements'
lines in the windows. That's how alloys that overlap in RBS (Ta–W, Fe–Ni,
Ni–Co) get their ratio from their X-ray lines; the
[FeNi worked example](../getting-started/feni-rbs-pixe.md) goes through
one.

Still to come: the PIXE energy calibration and resolution as fit
parameters, the background (with its scale s), and fitting the PIXE spectrum
on its own, without RBS.

### Known limitations

* **Light elements (C, N, O):** ECPSSR is least reliable here. For carbon
  the tables run from about 23 % above to 16 % below ISICS across the He
  energy range. Mass attenuation coefficients and detector efficiency below
  1 keV are also uncertain, so light-element results depend on H from
  standards.
* **M shells** use ISICS2011's simplified model and carry larger
  cross-section uncertainties than K and L — expect tens of percent for
  slow ⁴He — absorbed into H_M from a standard (e.g. an Au film whose
  thickness RBS has measured).
* **He beams on light elements** cause multiple ionisation, which shifts
  and broadens the K lines; this is not modelled.
* Films are laterally uniform and flat with sharp interfaces — the same
  assumptions as the RBS simulation.

### References

* W. Brandt, G. Lapicki, Phys. Rev. A 20 (1979) 465; Phys. Rev. A 23 (1981)
  1717 — ECPSSR theory.
* Z. Liu, S.J. Cipolla, Comput. Phys. Commun. 97 (1996) 315; S.J. Cipolla,
  Comput. Phys. Commun. 176 (2007) 157 and 180 (2009) 1716 — ISICS, whose
  2011 version (ISICS2011, CPC Program Library ADDS v5.1,
  [doi:10.17632/dmjpvt86cx.1](https://doi.org/10.17632/dmjpvt86cx.1)) gives
  the M-shell treatment used here and served to validate all the tables.
* J.A. Maxwell, J.L. Campbell, W.J. Teesdale, Nucl. Instr. Meth. B 43 (1989)
  218; J.A. Maxwell, W.J. Teesdale, J.L. Campbell, Nucl. Instr. Meth. B 95
  (1995) 407 — the GUPIX spectrum model and layered-target yields.
* J.L. Campbell et al., Nucl. Instr. Meth. B 499 (2021) 77 — SDD line
  shapes, H(E) calibration.
* T. Schoonjans et al., Spectrochim. Acta B 66 (2011) 776 — xraylib.
