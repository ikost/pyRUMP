"""Konac/Kalbitzer stopping powers.

Port of ``RbsCalcKalbitzerStop`` (stopping.c:790-825).

Two things about this model are easy to get wrong and both change results:

1. **It outranks Ziegler.** The priority chain tries Kalbitzer *before* ZBL
   (stopping.c:481 vs :483), so for H/D/3He/4He on C or Si -- the most common RBS
   cases by far -- RUMP is not using Ziegler at all.
2. **It already includes nuclear stopping.** The routine returns ``se + sn``, with
   ``sn`` taken from the same universal ZBL expression as :mod:`ziegler`. Callers
   receive a combined value, unlike ``zstop`` which separates the two.
"""

from __future__ import annotations

import numpy as np

from ..io.kalbitzer import KonacEntry
from ..model.element import Element

#: The literal in the C is `log(2.7182818 + E*beta)`, i.e. Euler's number appears
#: as an additive constant inside the log, not as a base (stopping.c:811).
_E_CONSTANT = 2.7182818


class KalbitzerStopping:
    """Konac fits for the specific ion/target pairs present in ``newstop.kal``."""

    def __init__(self, entries: list[KonacEntry], elements: list[Element]):
        self._entries = entries
        self._elements = elements

    def lookup(self, z1: int, m1: float, z2: int) -> KonacEntry | None:
        for entry in self._entries:
            if entry.matches(z1, m1, z2):
                return entry
        return None

    def __call__(self, z1: int, m1: float, z2: int, energy_keV) -> np.ndarray | None:
        """Combined electronic + nuclear stopping, or None if no entry matches.

        Result is in eV/(1e15 atoms/cm^2).
        """
        entry = self.lookup(z1, m1, z2)
        if entry is None:
            return None

        target = self._elements[z2 - 1].ziegler
        if target is None:
            return None
        m2 = target.mass_average

        energy = np.atleast_1d(np.asarray(energy_keV, dtype=np.float64))

        screen = z1**0.23 + z2**0.23
        epsilon0 = 32.53 * m2 / (z1 * z2 * (m1 + m2) * screen)
        sn0 = z1 * z2 * m1 * 8.462 / ((m1 + m2) * screen)

        # Electronic term. Note E is MeV/amu here, using the *table's* mass number.
        e_mev_amu = energy / 1000.0 / entry.m1
        a = entry.a
        denominator = (
            a[0]
            + a[1] * e_mev_amu**0.25
            + a[2] * e_mev_amu**0.5
            + a[3] * e_mev_amu**0.75
            + a[4] * e_mev_amu
            + a[5] * e_mev_amu ** (1 + entry.s)
        )
        se = (
            entry.scaling
            * e_mev_amu**entry.s
            * np.log(_E_CONSTANT + e_mev_amu * entry.beta)
            / denominator
        )

        # Universal ZBL nuclear stopping, same expression as ziegler.nuclear_stopping.
        epsilon = epsilon0 * energy
        helper = 0.01321 * epsilon**0.21226 + 0.19593 * np.sqrt(epsilon)
        sn = np.where(
            epsilon < 30.0,
            0.5 * np.log(1 + 1.1383 * epsilon) / (epsilon + helper),
            np.log(np.maximum(epsilon, 1.0 + 1e-30)) / (2 * epsilon),
        )
        return se + sn * sn0
