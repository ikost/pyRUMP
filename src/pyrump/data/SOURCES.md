# Data table provenance

These four tables are bundled so `pip install pyrump` works with no configuration. They
originate with the legacy RUMP 2.0 C distribution's `rump/data/` tree, cross-checked against
current primary/authoritative sources before shipping. See `README.md` at the repo root
("Licensing and provenance") for the non-redistributability of the RUMP C source itself,
which is unrelated to these tables.

## `pscoef.dat` — Ziegler/TRIM stopping coefficients

J. F. Ziegler, J. P. Biersack, U. Littmark, *The Stopping and Range of Ions in Solids*,
Pergamon Press (1985) — the ZBL/TRIM `SCOEF` table.

## `newstop.kal` — Konac/Kalbitzer stopping-power fits

G. Konac, S. Kalbitzer, Ch. Klatt, D. Niemann, R. Stoll, Nucl. Instr. Meth. B136-138 (1998)
159-165.

## `atom4.dat` — elements, atomic masses, isotopic abundances

Checked (August 2026) against current CIAAW (Commission on Isotopic Abundances and Atomic
Weights) standard atomic weights and isotope-abundance intervals for a broad element sample
spanning Z=1-92. Every sampled atomic mass matches CIAAW's current values well within
RBS-relevant precision (<0.1% relative).

Two elements have isotope abundances CIAAW has revised since this table's 1990s origin:
- **Magnesium** — revised post-2011 (isotope-ratio remeasurements); atom4.dat's Mg-24/25/26
  fractions differ from the current CIAAW interval by <1% relative per isotope, <0.03% on the
  average atomic mass.
- **Zirconium** — revised by CIAAW in 2024; atom4.dat's Zr-92/94/96 fractions differ from the
  current values by up to ~0.7% relative per isotope, negligible effect on the average mass.

Neither revision is large enough to affect RBS simulation results; noted here for the record
rather than patched, since patching isotope fractions risks subtly breaking the parser's
assumptions about the table's internal normalization.

### Elemental densities (`atomic_density`, atoms/cm^3)

Checked (September 2026) against NIST's X-Ray Mass Attenuation Coefficients reference
densities (Hubbell & Seltzer, NISTIR 5632,
<https://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html>), Z=1-92. Converting
`atomic_density` back to g/cm^3 via each element's mass and Avogadro's number and comparing
against that table: **median relative difference 0.3%** across the 75 elements with a
comparable (non-gas, non-zero) reference density.

Four elements differed by more than 5%. Three are allotrope-convention mismatches, not
errors — carbon (atom4.dat's 2.267 g/cm^3 is graphite's crystallographic density; NIST's
reference table uses a lower bulk value), phosphorus (atom4.dat uses white phosphorus,
1.82 g/cm^3; NIST uses red, ~2.2), and selenium (atom4.dat uses the stable grey allotrope,
4.79 g/cm^3; NIST's is lower). **Calcium had no such explanation** — atom4.dat carried
1.34 g/cm^3 against NIST's 1.55 g/cm^3 and multiple other literature sources in the same
1.54-1.59 g/cm^3 range — and has been corrected to match NIST's reference value
(`atomic_density` 2.0144e22 -> 2.3289e22 at/cm^3).

## `density.tab` — compound densities for thin-film thickness units

Bulk/theoretical densities, "compliments of Prof. Chris Palmstrom" per the file's own header.
Cross-checked against literature for the compound-semiconductor and oxide/nitride entries.

**One correction made**: the `GaP` entry was `0.2942e23 at/cm^3`, identical to the `InSb` row
below it — a copy-paste error in the original table. Corrected to `0.494e23 at/cm^3`, derived
from GaP's literature density of 4.14 g/cm^3 (cross-checked against the Ioffe Institute's
semiconductor database, ChemicalBook, and Wikipedia, all in agreement).

Everything else in the table (AlP, AlAs, AlSb, GaAs, GaSb, InP, InAs, InSb, SiO2, SiNx) matches
literature values within ~1%.

**Physics caveat, not a table defect**: SiO2 and SiNx here are theoretical/fully-dense values.
Real deposited thin films commonly run below them — sputtered SiO2 is often ~80-98% of this
table's value, PECVD Si3N4 as low as ~73% (hydrogen incorporation), LPCVD Si3N4 ~79-92%. Treat
density as a fit-checked parameter for real samples, not an assumed constant.

Al2O3 and TiO2 are not present in this table (not wrong, just absent) — add entries if you need
those materials' named-shorthand thickness units.

## `.adt` / R33 non-Rutherford cross-section files — not bundled

IBANDL (IAEA) third-party evaluations, A. F. Gurbich, Nucl. Instr. Meth. B136-138 (1998) 60.
Not vendored here — re-download from IBANDL if you need non-Rutherford cross sections. (Note:
this pyRUMP release doesn't wire tabulated-resonance cross sections into any user-facing
command yet regardless — see README's known limitations.)
