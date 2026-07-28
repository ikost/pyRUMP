"""Hard-coded Mylar stopping powers for the Z=93 pseudo-element.

Port of the fallback branch at stopping.c:490-500, described in the C itself as a
*"Kludge to support ancient Mylar hardcoded into DOS RUMP"*.

Mylar is outside Ziegler's Z=1..92 range, so stopping foils would otherwise have
no model at all. Only H, D and 4He are covered; anything else falls through.

The same three polynomials appear in three places in the distribution
(here, ``data/except.stp``, and ``MyXsect/Xsect.c``). These are the live ones --
``except.stp`` is never auto-loaded.
"""

from __future__ import annotations

import numpy as np

#: RUMP's Mylar pseudo-element.
MYLAR_Z = 93

#: Degree-5 polynomials in E[keV], giving eV/(1e15 atoms/cm^2) directly.
_H_IN_MYLAR = (13.37, -1.969e-2, +1.551e-5, -6.592e-9, +1.426e-12, -1.229e-16)
_D_IN_MYLAR = (10.37, -4.430e-4, -6.729e-6, +4.674e-9, -1.239e-12, +1.172e-16)
_HE_IN_MYLAR = (12.71, +5.538e-2, -5.557e-5, +2.440e-8, -5.175e-12, +4.274e-16)


def mylar_coefficients(z1: int, m1: float) -> tuple[float, ...] | None:
    """Select the polynomial for this projectile, or None if unsupported.

    Selection mirrors stopping.c:496: hydrogen with mass number 1 or 2, or
    helium-4. Note the C tests ``nint(m1) <= 2`` for hydrogen, so a mass number
    of 0 (natural abundance) also selects deuterium's curve -- preserved here.
    """
    mass_number = int(m1 + 0.5)
    if z1 == 1 and mass_number <= 2:
        return _H_IN_MYLAR if mass_number == 1 else _D_IN_MYLAR
    if z1 == 2 and mass_number == 4:
        return _HE_IN_MYLAR
    return None


def mylar_stopping(z1: int, m1: float, energy_keV) -> np.ndarray | None:
    """Combined stopping in Mylar, or None if the projectile is unsupported."""
    coefficients = mylar_coefficients(z1, m1)
    if coefficients is None:
        return None
    energy = np.atleast_1d(np.asarray(energy_keV, dtype=np.float64))
    # Horner, matching the C's nesting order exactly.
    result = np.full_like(energy, coefficients[5])
    for c in reversed(coefficients[:5]):
        result = result * energy + c
    return result
