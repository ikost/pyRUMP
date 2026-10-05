# ECPSSR ionization cross-section tables for pyrump (H⁺ and ⁴He)

These tables give ion-induced inner-shell **ionization** cross sections for PIXE simulation. They are calculated with the ECPSSR theory of Brandt & Lapicki, using an independent Python implementation written for this purpose. No ISICS or GUPIX code or data was used to produce them. ISICS-derived values were used only for validation (section 5).

| File | Content |
|---|---|
| `src/pyrump/data/pixe/ecpssr_H1.csv.gz` | protons: K (Z = 6–92), L1, L2, L3 (Z = 18–92), M1–M5 (Z = 50–92), E = 0.1–5 MeV (86 energies) |
| `src/pyrump/data/pixe/ecpssr_He4.csv.gz` | ⁴He: same shells and Z ranges, E = 0.2–12 MeV (90 energies) |
| `tools/ecpssr/validation_vs_ISICS.png` | ratio of these tables to ISICS-derived values |
| `tools/ecpssr/` | the generator (`ecpssr.py`, `pwba.py`, `gos.py`, `make_tables.py`); its cached L-shell form factors (`gos_20.npz`, `gos_21.npz`) are kept out of git |

pyRUMP reads the tables with `pyrump.pixe.xsect`; nothing here is imported at run time.

M shells (M1–M5, Z = 50–92) follow the simplified M-shell ECPSSR of ISICS2011
(ISICS v5.1, section 2b).

---

## 1. CSV format

One row per (ion, element, subshell, energy):

| column | meaning |
|---|---|
| `ion` | `H1` or `He4` |
| `Z1` | projectile charge (1 or 2) |
| `A1_u` | projectile (nuclear) mass in u: 1.007276467 (p), 4.001506179 (α) |
| `Z`, `element` | target atomic number and symbol |
| `shell` | `K`, `L1`, `L2`, `L3`, `M1` … `M5` |
| `U_keV` | subshell binding energy used (xraylib `EdgeEnergy`, v4.3.0) |
| `E_MeV` | **total** lab kinetic energy of the ion (not per nucleon) |
| `sigma_barn` | ionization cross section in barn (1 b = 10⁻²⁴ cm²) |

Energy grids: log-spaced, 50 points per decade — H 0.1–5 MeV, He 0.2–12 MeV. At that
density, linear interpolation in (ln E, ln σ) reproduces direct ECPSSR calls at the grid
midpoints within 0.3 % (worst case: L1 of Gd and Pt at the lowest energies). The first
tables, on a hand-picked grid of 31–34 energies, were off by up to 5 % there.

**Interpolation:** interpolate linearly in (ln E, ln σ), and do not extrapolate. A handful of rows are exactly 0: K shell of Z ≥ 85 at ≤ 0.12 MeV protons, where the energy-loss correction closes the channel. Treat those as 0.

**These are ionization, not X-ray production, cross sections.** pyRUMP converts them to line production cross sections (`pyrump.pixe.production`) with the fluorescence yields, direct Coster–Kronig probabilities (`FL12/FL13/FL23`, not `FLP13`) and radiative rates of its xraylib tables:

- σ_x(Kℓ) = σ_K ω_K F_Kℓ
- L3 vacancies: n₃ = σ_L3 + f₁₃ σ_L1 + f₂₃ (σ_L2 + f₁₂ σ_L1). Then σ_x(L3ℓ) = n₃ ω₃ F₃ℓ, and similarly for L1 and L2.

---

## 2. Theory implemented

The model is the ECPSSR theory of Brandt & Lapicki: perturbed-stationary-state binding and polarization, energy loss, Coulomb deflection and relativistic corrections applied to the hydrogenic PWBA.

$$\sigma^{ECPSSR}_s = C_s(x_C)\; f_s(z)\; \sigma^{PWBA}_s\!\left(m^R_s\,\eta_s,\ \zeta_s\theta_s\right)$$

All quantities are in atomic units, for subshell s with principal quantum number n:

| quantity | expression |
|---|---|
| screened charge | Z_s = Z₂ − 0.3 (K), Z₂ − 4.15 (L) |
| reduced binding | θ_s = n² U_s / (Z_s² Ry) |
| reduced energy | η_s = (v₁ / Z_s v₀)² |
| velocity parameter | ξ_s = 2n √η_s / θ_s |
| binding + polarization | ζ_s = 1 + 2Z₁/(Z_s θ_s) · [g_s(ξ_s) − h_s(ξ_s)] |
| polarization | h_s = 2n I(c_s n/ξ_s)/(θ_s ξ_s³), with c = 1.5 (K, L1) and 1.25 (L2, L3); I(x) is Brandt–Lapicki's piecewise function |
| binding functions | g_K, g_L1, g_L2,3: Brandt–Lapicki rational polynomials |
| relativistic | m^R = √(1 + 1.1y²) + y, with y = 0.40 (Z_s/c)² / (n ξ/ζ) (K, L1) and y = 0.15 (Z_s/c)² / (ξ/ζ) (L2, L3) |
| energy loss | z = √(1 − 4ζ/(Mθξ²)); f(z) = 2⁻ᵖ(p−1)⁻¹[(pz−1)(1+z)ᵖ + (pz+1)(1−z)ᵖ], with p = 9 (K, L1) and 11 (L2, L3) |
| Coulomb deflection | x_C = 2π d q₀ ζ / (z(1+z)), d q₀ = Z₁Z₂U/(M v₁³); C = p·E_{p+1}(x_C) |

Here M is the reduced mass of the projectile and target atom.

**PWBA:** σ^PWBA = 8π a₀² Z₁²/(η Z_s⁴) · I_s(η, θ), where

$$I_s = N_e\int_{\theta/2n^2}^{\infty} d\varepsilon \int_{\varepsilon/\sqrt\eta}^{\infty}\frac{dq}{q^3}\,G_{nl}(k^2,q),\qquad \tfrac{k^2}{2}=\varepsilon-\tfrac{1}{2n^2}$$

- **K shell:** G is the closed-form Walske / Merzbacher–Lewis 1s form factor.
- **L1 (2s) and L2/L3 (2p):** G is computed numerically. Energy-normalized Coulomb waves (power series near the origin, Numerov beyond) are summed over partial waves l′ ≤ 45, and the overlaps are integrated with spherical Bessel functions.

**Outer screening:** the ejected electron's energy is the hydrogenic one, k²/2 = ε − 1/(2n²). The observed binding energy θ only sets the lower limit of the energy-transfer integral. When θ < 1, energy transfers just below the hydrogenic threshold are included by analytic continuation of the form factor to k² < 0, with the factor 1/(1 − e^{−2π/k}) replaced by 1. This is the Khandelwal–Choi–Merzbacher convention used in the standard PWBA tables and in ISICS. It changes σ by up to a factor ~3 for light elements, so it matters.

The construction was checked in three ways:

- The high-velocity Bethe coefficient C₁(θ) of the K shell reproduces the tabulated values exactly (to 4 digits, θ = 0.4–1).
- The reduced PWBA functions F = θI/η agree with the tabulated K, L1 and L2,3 universal functions within ±1 % (K) and ±2 % (L) for η/θ² ≥ 0.01, over θ = 0.3–3.
- The low-velocity K limit reproduces Brandt–Lapicki's asymptotic formula F → (2⁹/45) ξ⁸ (1 + 1.72ξ²)⁻⁴.

---

## 3. What is *not* in the model (use with care)

- **No united-atom (UA) or Hartree–Slater (HS) corrections**, and no electron capture: this is the plain ECPSSR. Electron capture is negligible for H and He on the elements here.
- **Light elements at high velocity (C–Ne K shell; L shells of Z ≲ 30):** these tables are 10–20 % **below** ISICS (section 5), probably because of additional corrections in ISICS. The deviation grows with ion energy.
- **Heavy elements at very low velocity** (e.g. Au K below 1 MeV protons): the Coulomb-deflection factor gives values well below ISICS. These cross sections are ≤ 10⁻³ b and irrelevant for PIXE.
- **L1 subshell:** agreement with ISICS is only ±10 %. Its effect on Lα is small (it enters through CK transfer), but it matters for Lβ₃,₄ and Lγ₂,₃.
- **He on light elements:** multiple ionization shifts ω and line energies. That is outside the ionization cross section and should be treated in the line model.
- **Binding energies** are xraylib edge energies (solid-state values).

In practice, pyrump's H(E) calibration with standards (separate K and L curves) absorbs systematic errors of this size.

---

## 4. Typical values (for sanity checks)

| ion, E | Si K | Fe K | Pt L3 |
|---|---|---|---|
| H⁺, 2 MeV | 1.85 × 10⁴ b | 250 b | 80.5 b |
| ⁴He, 2 MeV | 1.07 × 10⁴ b | 16.5 b | 5.53 b |

---

## 2b. M shells: ISICS2011's model

For M1–M5 the Brandt–Lapicki binding, polarization and relativistic
functions are not available here, so the M shells follow the treatment of
ISICS2011, version 5.1 of ISICS (Z. Liu, S.J. Cipolla, Comput. Phys. Commun.
97 (1996) 315; S.J. Cipolla, Comput. Phys. Commun. 176 (2007) 157 and 180 (2009)
1716; CPC Program Library ADDS v5.1):

| quantity | M shells |
|---|---|
| screened charge | Z_s = Z₂ − 11.25 (3s, 3p), Z₂ − 21.15 (3d) |
| electrons | M1 2, M2 2, M3 4 (3p split ⅓ : ⅔), M4 4, M5 6 (3d split 0.4 : 0.6) |
| binding + polarization | ζ = U_s(Z₂+Z₁)/U_s(Z₂), the united-atom binding (the subshell's binding energy in element Z₂+Z₁ over its own); (1 + Z₁/Z_s)² when Z₁+Z₂ > 103, as in ISICS v1.0; no g(ξ) or h(ξ) |
| relativistic | none, m^R = 1 |
| energy loss, Coulomb deflection | as for K/L, with p = 9, 11, 13 for 3s, 3p, 3d |

ISICS does not apply the energy-loss factor f(z); it is kept here, as for K and L
(it is ≈ 1 for these heavy targets). The PWBA form factors for 3s, 3p, 3d are computed
numerically like 2s/2p (`gos_30/31/32.npz`, about 13 minutes in all), not from ISICS's
analytic expressions.

**Validation against ISICS itself** (v1.0 and v5.1 compiled locally from the CPC
distribution, used only for this comparison), with ISICS's own binding energies and
masses so that only the models are compared:

| shells | cases | ours / ISICS |
|---|---|---|
| M1–M5 | Sn, W, Pt, Au, Pb, U; H 1–3 MeV, ⁴He 1–3 MeV | 0.989–1.000 against v5.1 (worst 1.1 %, mostly < 0.5 %) |
| K, L1–L3 | Fe, Ag, Au; H and ⁴He 1–3 MeV | 0.97–1.00 for L; K the same, except Au K at low velocity (0.83–0.97), where only our f(z) differs and σ < 0.005 b. v1.0 and v5.1 give identical K and L |

The united-atom binding of v5.1 lowers the M cross sections by 1–10 % from v1.0's
(1.9 MeV ⁴He on Pt: M5 −3 %, M2 −10 %).

## 5. Validation against ISICS-derived values

The reference is the ISICS ECPSSR cross sections as reproduced by the polynomial fits of Taborda et al. (X-Ray Spectrom. 40 (2011) 127), which Geant4 distributes. That data carries a non-commercial, Geant4-only licence, so it was used only for this comparison and is **not** included here.

The comparison covers 17 elements per shell (Z = 6–92 for K, Z = 18–92 for L), for E ≥ 0.5 MeV (H) and E ≥ 1 MeV (He), and only rows with σ > 0.01 b.

| ion | shell | within ±5 % | within ±10 % | median ratio | range |
|---|---|---|---|---|---|
| H⁺ | K | 82 % | 91 % | 1.014 | 0.81 – 1.05 |
| H⁺ | L1 | 50 % | 89 % | 1.019 | 0.76 – 1.11 |
| H⁺ | L2 | 89 % | 96 % | 1.021 | 0.81 – 1.07 |
| H⁺ | L3 | 89 % | 96 % | 1.019 | 0.81 – 1.06 |
| ⁴He | K | 79 % | 96 % | 1.028 | 0.85 – 1.10 |
| ⁴He | L1 | 56 % | 87 % | 1.021 | 0.78 – 1.12 |
| ⁴He | L2 | 76 % | 98 % | 1.039 | 0.83 – 1.10 |
| ⁴He | L3 | 81 % | 97 % | 1.037 | 0.83 – 1.10 |

The largest deviations are light-element K shells and low-Z L shells at the highest energies (see `validation_vs_ISICS.png`). For Z ≳ 12 (K) and Z ≳ 36 (L2, L3) at 0.5–5 MeV protons, agreement is typically within ±4 %. He values run systematically about 3–5 % above ISICS.

---

## 6. Regenerating / extending

```bash
pip install numpy scipy xraylib
python tools/ecpssr/make_tables.py   # writes src/pyrump/data/pixe/; ~10 min with the cached gos_20/gos_21.npz
```

- Deleting the `.npz` files forces the numerical L-shell form factors to be recomputed (~11 min).
- `ecpssr.ecpssr(shell, Z1, A1, Z2, A2, E_MeV, U_keV)` returns a single cross section. It works for any bare light ion (d, ³He) and any energy with θζ ≥ 0.3; below that the PWBA form-factor grid ends and it raises rather than clamping.
- **M shells:** a fuller M-shell ECPSSR (binding, polarization and relativistic corrections for 3s/3p/3d, as in later ISICS versions) would replace ISICS v1.0's simplified ζ in `ecpssr.py`; the form factors are already there.

---

## 7. References

- W. Brandt, G. Lapicki, Phys. Rev. A 20 (1979) 465; Phys. Rev. A 23 (1981) 1717. ECPSSR theory.
- E. Merzbacher, H.W. Lewis, Handbuch der Physik 34 (1958) 166. PWBA, outer screening.
- G.S. Khandelwal, B.H. Choi, E. Merzbacher, At. Data 1 (1969) 103; B.H. Choi, E. Merzbacher, G.S. Khandelwal, At. Data 5 (1973) 291. PWBA K, L.
- D.W. Rice, G. Basbas, F.D. McDaniel, At. Data Nucl. Data Tables 20 (1977) 503. Extended PWBA K tables.
- G. Lapicki, J. Phys. Chem. Ref. Data 18 (1989) 111.
- Z. Liu, S.J. Cipolla, Comput. Phys. Commun. 97 (1996) 315; S.J. Cipolla, CPC 180 (2009) 1716. ISICS (validation reference only).
- A. Taborda et al., X-Ray Spectrom. 40 (2011) 127. Polynomial fits to ISICS (validation reference only).
- T. Schoonjans et al., Spectrochim. Acta B 66 (2011) 776. xraylib (binding energies).
