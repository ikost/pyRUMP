"""Inbound energy march and straggling accumulation.

Port of ``SimPrecal`` (creatr.c:1506-1568).

The beam energy at every slab interface is computed once, before any element is
considered, because the incoming path does not depend on what it eventually
scatters from. This is what makes the outward path the expensive part of the
simulation rather than the inward one.

Energy loss uses a **third-order Taylor expansion** of dE/da = -eps(E), not a
numerical ODE solve and not the surface approximation (creatr.c:1554)::

    de = secin*p0*(1 - secin*(0.5*p1 - 0.1666667*secin*(p1*p1 + p0*p2)))

Truncating after the first term would give the surface approximation; the extra
terms are what let RUMP use thick slabs and stay fast.

Two unit conventions are folded together and are easy to get wrong:

* internal energies are **keV**
* the stopping polynomials return **eV** through the slab (they are already
  multiplied by areal density)

The C hides the conversion inside the secant: ``secin = 1.0E-3 * samm->secin``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..stopping.bragg import evaluate_slab_loss
from ..stopping.table import NDEG, StoppingTable

#: Bohr straggling constant: 4*pi*(Z1 e^2)^2, keV^2 per (1e15 at/cm^2) per Z2
#: (creatr.c:1517). The literals are the C's own.
_STRAGGLE_CONSTANT = 12.56637e15
_ELECTRON_CHARGE_TERM = 1.4398e-10


@dataclass(slots=True)
class InboundPath:
    """Beam energy and accumulated straggling along the inward path."""

    energy: np.ndarray
    """``(n_slab+1,)`` beam energy in keV at each slab interface."""

    straggle: np.ndarray
    """``(n_slab+1,)`` cumulative Bohr variance in keV^2 at each interface."""

    reached: int
    """Number of slabs actually traversed before hitting the cutoff.

    Interfaces beyond this are **not** valid: the C leaves them stale from a
    previous iteration and relies on downstream cutoff checks (creatr.c:1544).
    """

    def valid(self) -> np.ndarray:
        return self.energy[: self.reached + 1]


def bohr_straggle_constant(z_beam: int, scale: float = 1.0) -> float:
    """``scon`` from creatr.c:1517.

    ``scale`` is RUMP's ``SIM STRAGGLE`` multiplier, which defaults to **0** --
    straggling is off unless explicitly enabled.
    """
    return scale * _STRAGGLE_CONSTANT * (z_beam * _ELECTRON_CHARGE_TERM) ** 2


def energy_loss_step(
    stopping: float, first_derivative: float, second_derivative: float, sec: float
) -> float:
    """Third-order energy loss through one slab, in keV (creatr.c:1554).

    ``sec`` must already carry the 1e-3 eV->keV factor.
    """
    return (
        sec
        * stopping
        * (
            1.0
            - sec
            * (
                0.5 * first_derivative
                # NOTE: literal 0.1666667, not 1/6 -- kept for fidelity.
                - 0.1666667
                * sec
                * (first_derivative * first_derivative + stopping * second_derivative)
            )
        )
    )


def march_inbound(
    table: StoppingTable,
    slab_coefficients: np.ndarray,
    slab_composition: np.ndarray,
    element_z: list[int],
    *,
    e0_keV: float,
    sec_in: float,
    cutoff_keV: float,
    straggle_scale: float = 0.0,
    z_beam: int = 2,
    e_scale: float = 1.0,
    first_slab: int = 0,
) -> InboundPath:
    """March the beam inward, recording energy and straggling per interface.

    This is inherently sequential -- each slab's loss depends on the energy left
    by the previous one -- so it stays a Python loop. It is O(n_slab) and run
    once per simulation, unlike the outward paths.

    ``first_slab`` is ``fsurf``: the march **starts there at the full beam
    energy** (creatr.c:1530). Absorber slabs sit between the sample and the
    detector, so the incoming beam never crosses them -- only the outgoing
    particle does. Marching through them on the way in would double-count their
    stopping and shift every edge.
    """
    n_slab = slab_coefficients.shape[0]
    energy = np.zeros(n_slab + 1, dtype=np.float64)
    straggle = np.zeros(n_slab + 1, dtype=np.float64)

    # The eV->keV conversion rides along on the secant, as in the C.
    sec = 1.0e-3 * sec_in
    scon = bohr_straggle_constant(z_beam, straggle_scale)

    z_weights = np.asarray(element_z, dtype=np.float64)
    running_z = 0.0

    energy[first_slab] = e0_keV
    current = e0_keV
    reached = first_slab

    for index in range(first_slab, n_slab):
        if current < cutoff_keV:
            break

        coefficients = slab_coefficients[index]
        x = table.transform(current * e_scale)
        powers = x ** np.arange(NDEG)

        p0 = float(coefficients @ powers)
        p1 = float(_derivative_at(table, coefficients, current * e_scale, 1) * e_scale)
        p2 = float(
            _derivative_at(table, coefficients, current * e_scale, 2) * e_scale * e_scale
        )

        current -= energy_loss_step(p0, p1, p2, sec)

        # Straggling accumulates across slabs: variance is proportional to the
        # running sum of N*Z2 over everything traversed so far.
        running_z += float(slab_composition[index] @ z_weights)
        straggle[index + 1] = scon * running_z

        energy[index + 1] = current
        reached = index + 1

    return InboundPath(energy=energy, straggle=straggle, reached=reached)


def _derivative_at(
    table: StoppingTable, coefficients: np.ndarray, energy_keV: float, order: int
) -> float:
    """Derivative of a slab-summed polynomial, using the table's conventions.

    Shares :meth:`StoppingTable.derivative`'s semantics, including RUMP's
    ``SQRT_DDS_POWER`` index bug for ``order=2``.
    """
    x = table.transform(energy_keV)
    first = np.polynomial.polynomial.polyder(coefficients)
    if order == 1:
        return float(np.polynomial.polynomial.polyval(x, first) / (2 * x))
    p = coefficients
    head = ((3.75 * p[5] * x + 2 * p[4]) * x + 0.75 * p[3]) * x * x
    return float((head - 0.25 * p[2]) / x**3)


def surface_energy_loss(
    table: StoppingTable,
    slab_coefficients: np.ndarray,
    energy_keV: float,
    sec_in: float,
) -> float:
    """First-order (surface-approximation) loss, for comparison and testing.

    This is what the expansion reduces to if truncated after the first
    derivative term (1985 paper, p. 345).
    """
    return 1.0e-3 * sec_in * float(
        evaluate_slab_loss(table, slab_coefficients, [energy_keV])[0]
    )
