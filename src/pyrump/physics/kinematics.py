"""Elastic-scattering kinematics.

Port of the kinematic factors in ``SimCideal`` (creatr.c:1659 backscatter,
creatr.c:1678 recoil). Non-relativistic throughout -- RUMP has no relativistic
correction anywhere.
"""

from __future__ import annotations

import math

import numpy as np


def kinematic_factor(m1: float, m2: float, scattering_angle_deg: float) -> float:
    """Backscattering kinematic factor K = E1/E0 (Chu et al. eq. 2.6).

    .. math::
        K = \\left[
            \\frac{\\sqrt{1 - (x\\sin\\phi)^2} + x\\cos\\phi}{1 + x}
        \\right]^2,\\quad x = m_1/m_2

    ``scattering_angle_deg`` is the **true** scattering angle. Raises if the
    projectile is heavier than the target and the angle exceeds the kinematic
    limit, where the square root turns imaginary.
    """
    x = m1 / m2
    phi = math.radians(scattering_angle_deg)
    # NOTE: the C computes sin/cos from the stored (supplement) phi and negates
    # the cosine; done here directly from the true angle, which is equivalent.
    inner = 1.0 - (x * math.sin(phi)) ** 2
    if inner < 0.0:
        limit = math.degrees(math.asin(1.0 / x))
        raise ValueError(
            f"scattering angle {scattering_angle_deg} deg exceeds the kinematic "
            f"limit {limit:.2f} deg for m1/m2 = {x:.3f}"
        )
    return ((math.sqrt(inner) + x * math.cos(phi)) / (1.0 + x)) ** 2


def recoil_factor(m1: float, m2: float, scattering_angle_deg: float) -> float:
    """Forward-recoil (ERD) kinematic factor (creatr.c:1678).

    .. math::
        K = \\frac{4 m_1 m_2 \\cos^2\\phi}{(m_1 + m_2)^2}

    The C writes this as ``4 cos^2 / (2 + x + 1/x)``, which is the same thing.
    """
    x = m1 / m2
    cos_phi = math.cos(math.radians(scattering_angle_deg))
    if cos_phi <= 0.0:
        raise ValueError(
            f"recoil requires a true scattering angle below 90 deg, got "
            f"{scattering_angle_deg}"
        )
    return 4.0 * cos_phi * cos_phi / (2.0 + x + 1.0 / x)


def edge_energy(e0_keV: float, m1: float, m2: float, scattering_angle_deg: float) -> float:
    """Surface-edge energy for a target nuclide: ``K * E0``.

    The sharpest available check on a simulated spectrum's energy calibration.
    """
    return e0_keV * kinematic_factor(m1, m2, scattering_angle_deg)


def kinematic_factors(
    m1: float, masses: np.ndarray, scattering_angle_deg: float
) -> np.ndarray:
    """Vectorised :func:`kinematic_factor` over target masses."""
    x = m1 / np.asarray(masses, dtype=np.float64)
    phi = math.radians(scattering_angle_deg)
    inner = 1.0 - (x * math.sin(phi)) ** 2
    if np.any(inner < 0.0):
        raise ValueError("one or more targets exceed the kinematic angle limit")
    return ((np.sqrt(inner) + x * math.cos(phi)) / (1.0 + x)) ** 2
