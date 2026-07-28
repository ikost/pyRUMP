"""Fitted stopping-power tables -- what the simulation actually evaluates.

Port of ``RbsStpfind`` + ``RbsGenStopp`` (stopping.c:270-568).

**This is the single most important indirection in RUMP.** The simulation never
calls Ziegler or Konac. At the start of a run it samples the chosen model at 201
points, fits a degree-5 polynomial per target element, and from then on evaluates
only that polynomial. Get the fit wrong by 0.3% and the whole depth scale is wrong
by 0.3%, with no other symptom.

Two consequences that are easy to miss:

* **The fit window depends on the beam energy** (``emin = 0.04*E``,
  ``emax = 1.15*E`` for STOP_SQRT). Change the beam and every coefficient changes.
* **The independent variable is sqrt(E)**, sampled uniformly in sqrt(E) -- not in
  E. So the low-energy end is oversampled relative to a linear grid.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum

import numpy as np

from .polyfit import fit_polynomial
from .registry import StoppingRegistry

#: Points sampled before fitting (``NPLOT``, stopping.c).
NPLOT = 201

#: Coefficients per fit: degree-5 polynomial (``NDEG``).
NDEG = 6


class StoppingType(IntEnum):
    """``STOPPING_TYPE`` -- which variable the polynomial is in."""

    LINEAR = 0
    SQRT = 1
    """The shipped default (stopping.c:73)."""


def fit_window(
    e_beam_MeV: float, kind: StoppingType = StoppingType.SQRT
) -> tuple[float, float, float]:
    """Return ``(emin, emax, cutoff)`` in MeV for a beam energy (stopping.c:308-321).

    The C's own comment here is "Need some experience in setting these limits".
    """
    if kind == StoppingType.LINEAR:
        return 0.08 * e_beam_MeV, 1.15 * e_beam_MeV, 0.03 * e_beam_MeV
    return 0.04 * e_beam_MeV, 1.15 * e_beam_MeV, 0.03 * e_beam_MeV


@dataclass(slots=True)
class StoppingTable:
    """Per-element fitted stopping polynomials for one projectile and beam energy."""

    z_beam: int
    m_beam: float
    e_beam_MeV: float
    emin: float
    emax: float
    cutoff: float
    kind: StoppingType
    coefficients: dict[int, np.ndarray]
    """``{Z_target: ascending-order coefficients in sqrt(E[keV])}``."""

    @classmethod
    def build(
        cls,
        registry: StoppingRegistry,
        z_beam: int,
        m_beam: float,
        e_beam_MeV: float,
        targets: list[int],
        kind: StoppingType = StoppingType.SQRT,
    ) -> "StoppingTable":
        emin, emax, cutoff = fit_window(e_beam_MeV, kind)

        # Grid is uniform in the *fit* variable, i.e. in sqrt(E) for STOP_SQRT
        # (stopping.c:465-471).
        fxmin, fxmax = emin * 1000.0, emax * 1000.0
        if kind == StoppingType.SQRT:
            fxmin, fxmax = np.sqrt(fxmin), np.sqrt(fxmax)

        scaled = fxmin + np.arange(NPLOT) * (fxmax - fxmin) / (NPLOT - 1.0)
        energies = scaled**2 if kind == StoppingType.SQRT else scaled.copy()

        coefficients: dict[int, np.ndarray] = {}
        for z_target in targets:
            # RbsGenStopp fits electronic + nuclear COMBINED (stopping.c:485).
            power = registry(z_beam, m_beam, z_target, energies).values
            # sigma = sqrt(|S|), so the weight 1/sigma^2 = 1/S is relative
            # (stopping.c:519-521).
            sigma = np.where(power != 0.0, np.sqrt(np.abs(power)), 1.0)
            coefficients[z_target] = fit_polynomial(scaled, power, sigma, order=NDEG - 1)

        return cls(
            z_beam=z_beam,
            m_beam=m_beam,
            e_beam_MeV=e_beam_MeV,
            emin=emin,
            emax=emax,
            cutoff=cutoff,
            kind=kind,
            coefficients=coefficients,
        )

    # ------------------------------------------------------------- evaluation
    def transform(self, energy_keV) -> np.ndarray:
        """``S_XFORM``: map energy to the polynomial's independent variable."""
        energy = np.asarray(energy_keV, dtype=np.float64)
        return np.sqrt(energy) if self.kind == StoppingType.SQRT else energy

    def __call__(self, z_target: int, energy_keV) -> np.ndarray:
        """Stopping power in eV/(1e15 atoms/cm^2).

        No range check: RUMP happily extrapolates beyond ``emax``, and the port
        must behave the same. Use :meth:`in_window` if you need to know.
        """
        x = self.transform(energy_keV)
        return np.polynomial.polynomial.polyval(x, self.coefficients[z_target])

    def derivative(
        self, z_target: int, energy_keV, order: int = 1, *, faithful: bool = True
    ) -> np.ndarray:
        """dS/dE or d2S/dE2 with respect to **E**, per ``DS_POWER``/``DDS_POWER``.

        .. warning::
           **RUMP's second-derivative macro is wrong, and pyRUMP reproduces it.**

           ``SQRT_DDS_POWER`` (stopping.h:47-49) reads ``p[2]`` where the
           mathematics requires ``p[1]``::

               (3.75 p5 e^4 + 2 p4 e^3 + 0.75 p3 e^2 - 0.25 p2) / e^3
                                                            ^^ should be p1

           The header's own comment two lines above documents the correct form,
           ``DDS = 1/(4e^4) * SUM i(i-2) a_i e^i``, whose ``i=1`` term is
           ``-p1 e``. Verified against the C: the oracle reproduces the macro
           exactly, while a numerical second derivative reproduces the corrected
           form exactly. The discrepancy is 35-50%.

           This feeds only the third-order term of the energy-loss expansion
           (creatr.c:1554), so its effect on depth scale is small but real.
           Pass ``faithful=False`` for the mathematically correct value.
        """
        coefficients = self.coefficients[z_target]
        x = self.transform(energy_keV)

        if self.kind == StoppingType.LINEAR:
            derived = coefficients
            for _ in range(order):
                derived = np.polynomial.polynomial.polyder(derived)
            return np.polynomial.polynomial.polyval(x, derived)

        # With e = sqrt(E):  dS/dE = S'(e) / (2e)   -- this macro is correct.
        first = np.polynomial.polynomial.polyder(coefficients)
        if order == 1:
            return np.polynomial.polynomial.polyval(x, first) / (2 * x)

        p = coefficients
        head = ((3.75 * p[5] * x + 2 * p[4]) * x + 0.75 * p[3]) * x * x
        tail = 0.25 * (p[2] if faithful else p[1])
        return (head - tail) / x**3

    def in_window(self, energy_keV) -> np.ndarray:
        energy = np.asarray(energy_keV, dtype=np.float64)
        return (energy >= self.emin * 1000.0) & (energy <= self.emax * 1000.0)

    def accepts(self, z_beam: int, m_beam: float, e_beam_MeV: float) -> bool:
        """Whether ``RbsStpfind`` would reuse this table for the given beam.

        Reuse conditions from stopping.c:274-279. Note the beam energy only has
        to *fit inside* the window -- it need not be the energy the table was
        built for.
        """
        if self.z_beam != z_beam:
            return False
        return 2 * self.emin <= e_beam_MeV <= self.emax

    def matches_mass(self, m_beam: float) -> bool:
        """Exact-isotope match, within the C's 0.2 amu tolerance."""
        return abs(self.m_beam - m_beam) <= 0.2

    def max_error(self, registry: StoppingRegistry, z_target: int) -> float:
        """Worst relative fit error, excluding the first and last 5 points.

        Mirrors the C's own quality metric (stopping.c:556-561), which it prints
        as "max error: %.2f%%".
        """
        fxmin, fxmax = self.emin * 1000.0, self.emax * 1000.0
        if self.kind == StoppingType.SQRT:
            fxmin, fxmax = np.sqrt(fxmin), np.sqrt(fxmax)
        scaled = fxmin + np.arange(NPLOT) * (fxmax - fxmin) / (NPLOT - 1.0)
        energies = scaled**2 if self.kind == StoppingType.SQRT else scaled

        reference = registry(self.z_beam, self.m_beam, z_target, energies).values
        fitted = np.polynomial.polynomial.polyval(scaled, self.coefficients[z_target])
        relative = np.abs((reference - fitted) / reference)
        return float(relative[5:-5].max())
