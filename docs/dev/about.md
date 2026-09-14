## Milestones

| | Milestone | Status |
|---|---|---|
| M0 | Reference oracle (pty driver + cffi library) | done |
| M1 | Elements, isotopes, compound densities | done |
| M2 | Stopping powers: ZBL85, Konac, Mylar, priority chain | done |
| M3 | STOP_SQRT polynomial refit, Bragg summation, session cache | done |
| M4 | Kinematics, geometry, cross-sections | done |
| M5 | Slab march → bricks | done |
| M6 | Brick → channel fill | done |
| M7 | Straggling (closed-form erf) | done |
| M8 | Detector convolution | done |
| M9 | Depth profiles (13 EQUATION forms) | done |
| M10 | Absorber, pile-up, fuzz, multiple scattering | done |
| M11 | File I/O (`.RBS` binary, ASCII) | done |
| M12 | Fitting (PERT) | done |
| M13 | CLI, plotting, `.lcm` subset | done |
| M14 | Interactive shell: buffers, SIM and PERT levels, macros | done |

**Known limitation**: non-Rutherford (tabulated-resonance) cross sections —
RUMP's `.adt`/R33 nuclear cross-section tables — have a complete, tested
reader (`pyrump.io.adt`) but aren't wired into simulation yet; pyRUMP
currently computes pure Rutherford + L'Ecuyer-screened scattering only. A
later milestone.


## Licensing and provenance

pyRUMP is MIT licensed, and an **independent reimplementation**: it is not
affiliated with, endorsed by, or derived from the RUMP source distribution.

RUMP and Genplot were trademarks of Computer Graphic Service, Ltd. (CGS). CGS
ceased operating as a business in June 2012 and `genplot.com` no longer
resolves; the authors stated at the time that GENPLOT and RUMP remain free to
download and use, which removes the trademark concern but **not** copyright in
the original source, which remains with its authors — hence the C tree is
still not redistributed here.

Four data tables are bundled with pyRUMP (`src/pyrump/data/`), independent of
CGS and checked against current CIAAW/NIST values and literature:

- `pscoef.dat` — the ZBL/TRIM `SCOEF` stopping-coefficient table
- `newstop.kal` — Konac/Kalbitzer stopping-power fits
- `atom4.dat` — elements and isotopes
- `density.tab` — compound densities

See `src/pyrump/data/SOURCES.md` for full provenance, verification notes, and
the one correction made (a data-entry error in the GaP density); citations are
in [References](#references).

Non-Rutherford cross-section tables (`*.adt`) are IBANDL evaluations and are
**not** bundled — obtain them separately from IBANDL if you need that data.


## References

- L. R. Doolittle, *Algorithms for the rapid simulation of Rutherford
  backscattering spectra*, Nucl. Instr. Meth. **B9** (1985) 344–351.
- L. R. Doolittle, *A new approach to Rutherford backscattering analysis*,
  Nucl. Instr. Meth. **B15** (1986) 227–231.
- J. F. Ziegler, J. P. Biersack, U. Littmark, *The Stopping and Range of Ions in
  Solids*, Pergamon (1985) — source of `pscoef.dat`, the ZBL/TRIM SCOEF table.
- G. Konac, S. Kalbitzer, Ch. Klatt, D. Niemann, R. Stoll, Nucl. Instr. Meth.
  **B136–138** (1998) 159–165 — source of `newstop.kal`.
- W.-K. Chu, J. W. Mayer, M.-A. Nicolet, *Backscattering Spectrometry*,
  Academic Press (1978) — the [ε] stopping cross-section factor and kinematics.
- S. Baker, R. D. Cousins, Nucl. Instr. Meth. **221** (1984) 437 — the Poisson
  fitting objective.
- J. L'Ecuyer et al., Nucl. Instr. Meth. **160** (1979) 337 — screening
  correction.
- J. F. Ziegler, Nucl. Instr. Meth. **B136–138** (1998) 141 — screening and
  cross-section formulae.
- V. Quillet, F. Abel, M. Schott, Nucl. Instr. Meth. **B83** (1993) 47 —
  screening and cross-section formulae.
- A. F. Gurbich, Nucl. Instr. Meth. **B136–138** (1998) 60 — non-Rutherford
  cross-section evaluations, as distributed via IBANDL.
</content>
