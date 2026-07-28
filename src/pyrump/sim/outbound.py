"""Outbound path: flyout, the stopping cross-section factor, and the back layer.

Ports ``SimFlyout`` (creatr.c:1953), ``RbsEfact`` (creatr.c:2019) and
``SimBacklay`` (creatr.c:2062).

This is the expensive half of the simulation. The inbound march runs once, but a
particle scattered from slab *i* must be walked back out through slabs
*i, i-1, ..., 0*, for every slab and every isotope -- the O(N^2) term the 1985
paper's cost analysis is dominated by.

``SimFlyout`` returns two things:

* the energy that actually reaches the detector
* ``ratde``, the accumulated ratio of stopping powers before/after each slab.
  This is the spectrum-height compression factor: it accounts for the scattered
  beam's energy spread changing as it exits, and it is the same on both sides of
  an interface (1985 paper, p. 347).
"""

from __future__ import annotations

import numpy as np

from ..stopping.table import NDEG, StoppingTable


def _evaluate(table: StoppingTable, coefficients: np.ndarray, energy_keV: float) -> float:
    x = table.transform(energy_keV)
    return float(coefficients @ (x ** np.arange(NDEG)))


def _derivatives(
    table: StoppingTable, coefficients: np.ndarray, energy_keV: float
) -> tuple[float, float]:
    """First and second derivative, reproducing RUMP's ``SQRT_DDS_POWER`` bug."""
    x = table.transform(energy_keV)
    first = np.polynomial.polynomial.polyder(coefficients)
    d1 = float(np.polynomial.polynomial.polyval(x, first) / (2 * x))
    p = coefficients
    head = ((3.75 * p[5] * x + 2 * p[4]) * x + 0.75 * p[3]) * x * x
    d2 = float((head - 0.25 * p[2]) / x**3)
    return d1, d2


def flyout(
    table: StoppingTable,
    slab_coefficients_out: np.ndarray,
    from_slab: int,
    energy_keV: float,
    *,
    sec_out: float,
    cutoff_keV: float,
    e2_scale: float = 1.0,
    first_surface: int = 0,
) -> tuple[float, float]:
    """Walk a scattered particle out to the surface (creatr.c:1953).

    Returns ``(exit_energy, ratde)``.

    ``from_slab`` is the slab the particle scattered *out of*; the walk covers
    ``from_slab`` down to 0 inclusive. Passing -1 (scattering at the very
    surface) is a no-op with ``ratde = 1``.
    """
    if from_slab < 0:
        return energy_keV, 1.0

    ratde = 1.0
    energy = energy_keV
    sec = 1.0e-3 * sec_out

    for slab in range(from_slab, -1, -1):
        if energy < cutoff_keV:
            break
        # Absorber/dead layers sit in front of the sample and are not tilted
        # with it, so they are traversed at normal incidence (creatr.c:1971).
        step_sec = 1.0e-3 if slab < first_surface else sec

        coefficients = slab_coefficients_out[slab]
        p0 = _evaluate(table, coefficients, energy * e2_scale)
        d1, d2 = _derivatives(table, coefficients, energy * e2_scale)
        p1 = d1 * e2_scale
        p2 = d2 * e2_scale * e2_scale

        loss = (
            step_sec
            * p0
            * (
                1.0
                - step_sec
                * (0.5 * p1 - 0.1666667 * step_sec * (p1 * p1 + p0 * p2))
            )
        )
        energy -= loss
        if energy < cutoff_keV:
            break

        # Stopping power *after* the slab, for the dE'/dE compression ratio.
        p1_after = _evaluate(table, coefficients, energy * e2_scale)
        ratde = ratde * p0 / p1_after

    exit_energy = energy if energy >= cutoff_keV else cutoff_keV / 2.0
    return exit_energy, ratde


def epsilon_factor(
    table: StoppingTable,
    coefficients_in: np.ndarray,
    coefficients_out: np.ndarray,
    energy_keV: float,
    kinematic: float,
    *,
    sec_in: float,
    sec_out: float,
    e1_scale: float = 1.0,
    e2_scale: float = 1.0,
) -> float:
    """Surface stopping cross-section factor [eps], Chu et al. eq. 3.10.

    .. math::
        [\\varepsilon] = K\\,\\varepsilon(E)\\sec\\theta_{in}
                       + \\varepsilon(KE)\\sec\\theta_{out}

    Evaluated with slab-scaled coefficients, so the result is keV through the
    slab rather than a true cross-section (creatr.c:2019-2033).
    """
    eps_in = _evaluate(table, coefficients_in, energy_keV * e1_scale)
    eps_out = _evaluate(table, coefficients_out, energy_keV * e2_scale * kinematic)
    return (kinematic * eps_in * sec_in + eps_out * sec_out) * 1.0e-3


def back_layer_edge(
    energy_front: float, height_front: float, ratde: float, cutoff_keV: float
) -> tuple[float, float]:
    """Terminate a brick whose back edge falls below the cutoff (creatr.c:2062).

    The C calls this "a rather mediocre routine": it linearly extrapolates the
    last known front edge down to the cutoff, using a hard-coded 0.3 fudge on
    the height. It only affects the deep tail of the spectrum, but it produces
    the visible kink just above cutoff noted in the thesis (Fig. 4.2).
    """
    energy_back = cutoff_keV
    span = energy_front - energy_back
    if span <= 0:
        return energy_front, height_front
    # creatr.c:2069 -- the 0.3 is arbitrary and undocumented.
    height_back = height_front * (1.0 + 0.3 * span / max(energy_front, 1.0)) / ratde
    return energy_back, height_back
