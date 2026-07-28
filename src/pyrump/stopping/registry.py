"""The stopping-power source priority chain.

Port of the selection logic in ``RbsGenStopp`` (stopping.c:474-515). RUMP tries
sources in a fixed order and takes the first that can answer:

1. a user DLL (``UserZStop``) -- e.g. a real SRIM installation; not supported here
2. **Konac/Kalbitzer** (``newstop.kal``) -- outranks Ziegler, and covers the most
   common RBS cases (H/D/He on C and Si)
3. **Ziegler ZBL85** (``pscoef.dat``) -- the general fallback for Z = 1..92
4. hard-coded **Mylar** polynomials for the Z=93 pseudo-element
5. failure

An ``except.stp`` expression override, if registered, is then applied *on top* of
whichever source won.

Everything downstream consumes **electronic + nuclear combined**, because that is
what ``RbsGenStopp`` fits (``spower[i] = se + sn``, stopping.c:485).
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import numpy as np

from ..model.element import Element
from .kalbitzer import KalbitzerStopping
from .mylar import MYLAR_Z, mylar_stopping
from .ziegler import ZieglerStopping


class StoppingSource(Enum):
    """Which model answered. RUMP prints a one-letter code for each."""

    KALBITZER = "k"
    ZIEGLER = "z"
    MYLAR = "m"
    USER = "u"


@dataclass(frozen=True, slots=True)
class StoppingResult:
    values: np.ndarray
    """Combined electronic + nuclear stopping, eV/(1e15 atoms/cm^2)."""

    source: StoppingSource


class StoppingRegistry:
    """Resolves (projectile, target, energy) to a stopping power."""

    def __init__(
        self,
        elements: list[Element],
        kalbitzer: KalbitzerStopping | None = None,
        ziegler: ZieglerStopping | None = None,
    ):
        self._elements = elements
        self._kalbitzer = kalbitzer
        self._ziegler = ziegler or ZieglerStopping(elements)

    def __call__(self, z1: int, m1: float, z2: int, energy_keV) -> StoppingResult:
        energy = np.atleast_1d(np.asarray(energy_keV, dtype=np.float64))

        if self._kalbitzer is not None:
            values = self._kalbitzer(z1, m1, z2, energy)
            if values is not None:
                return StoppingResult(values, StoppingSource.KALBITZER)

        if 1 <= z1 <= 92 and 1 <= z2 <= 92:
            electronic, nuclear = self._ziegler(z1, m1, z2, energy)
            return StoppingResult(electronic + nuclear, StoppingSource.ZIEGLER)

        if z2 == MYLAR_Z:
            values = mylar_stopping(z1, m1, energy)
            if values is not None:
                return StoppingResult(values, StoppingSource.MYLAR)

        raise ValueError(
            f"no stopping power available for Z1={z1} m1={m1} Z2={z2}; "
            "Ziegler covers Z=1..92 and Mylar (Z=93) only supports H, D and 4He"
        )

    def source_for(self, z1: int, m1: float, z2: int) -> StoppingSource:
        """Which source would answer, without evaluating it."""
        if self._kalbitzer is not None and self._kalbitzer.lookup(z1, m1, z2):
            return StoppingSource.KALBITZER
        if 1 <= z1 <= 92 and 1 <= z2 <= 92:
            return StoppingSource.ZIEGLER
        if z2 == MYLAR_Z and mylar_stopping(z1, m1, [1000.0]) is not None:
            return StoppingSource.MYLAR
        raise ValueError(f"no stopping source for Z1={z1} m1={m1} Z2={z2}")
