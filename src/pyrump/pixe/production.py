"""From ionisation to X-rays: line production cross sections.

A vacancy made in a subshell can move to a less tightly bound subshell of
the same shell before it decays (a Coster-Kronig transition), so the
vacancies a subshell ends up with are its own ionisations plus what flows
in from above. For the L shell, with the *direct* probabilities f_ij::

    n1 = sigma_L1
    n2 = sigma_L2 + f12 n1
    n3 = sigma_L3 + f13 n1 + f23 n2

(this is Maxwell et al. 1989's L3 expression, sigma3 + sigma2 f23 +
sigma1 (f13 + f12 f23), unfolded), and likewise for M1-M5. A line filling a
vacancy in subshell i is then produced with

    sigma_line = n_i * omega_i * F_line

where omega_i is the fluorescence yield and F_line the line's share of the
radiative decays (its radiative rate). Vacancies made in the K shell and
passed on to L by the K X-rays and Auger electrons are not followed: next to
direct L ionisation they are negligible for PIXE.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .atomic import AtomicData, XrayLine
from .xsect import IonisationTable

#: Subshells of each shell, most tightly bound first -- the order vacancies
#: flow in.
SHELLS = {"K": ("K",), "L": ("L1", "L2", "L3"), "M": ("M1", "M2", "M3", "M4", "M5")}

#: Lines produced this much more weakly than the strongest line of their
#: element and shell (K, L or M) are dropped: they would not show in any
#: spectrum. Per shell, because the shells' lines lie far apart in energy:
#: next to a 1.9 MeV He beam's L lines of Ni (0.85 keV, thousands of barns,
#: stopped by any filter) its Kβ (0.2 b) would be dropped, yet it is the
#: line a spectrum shows at 8.26 keV.
WEAK_LINE = 1e-4


@dataclass(frozen=True, slots=True)
class LineProduction:
    line: XrayLine
    sigma_barn: np.ndarray
    """Production cross section at each of the requested energies."""

    @property
    def family(self) -> str:
        """``"K"``, ``"L"`` or ``"M"``: which H value applies."""
        return self.line.shell[0]


def vacancies(
    z: int, family: str, energy_MeV: np.ndarray, ionisation: IonisationTable,
    atomic: AtomicData,
) -> dict[str, np.ndarray]:
    """Vacancies per subshell of one shell after Coster-Kronig transfer, in
    barn, at each energy. Subshells without cross sections give nothing."""
    available = set(ionisation.shells(z))
    subshells = SHELLS[family]
    result: dict[str, np.ndarray] = {}
    for j, shell in enumerate(subshells):
        if shell not in available:
            continue
        n = ionisation.sigma(z, shell, energy_MeV).copy()
        for i in range(j):
            upper = subshells[i]
            if upper in result:
                f = atomic.coster_kronig(z, f"{family}{i + 1}{j + 1}")
                n = n + f * result[upper]
        result[shell] = n
    return result


def line_production(
    z: int, energy_MeV, ionisation: IonisationTable, atomic: AtomicData
) -> list[LineProduction]:
    """Production cross section of every line of element ``z`` at each
    energy, the negligible ones (:data:`WEAK_LINE`) left out."""
    energy = np.atleast_1d(np.asarray(energy_MeV, dtype=np.float64))
    produced: list[LineProduction] = []
    for family in SHELLS:
        holes = vacancies(z, family, energy, ionisation, atomic)
        shell_lines = []
        for line in atomic.lines(z):
            if line.shell in holes:
                omega = atomic.fluorescence_yield(z, line.shell)
                shell_lines.append(LineProduction(line, holes[line.shell] * omega * line.rate))
        if shell_lines:
            strongest = max(float(p.sigma_barn.max()) for p in shell_lines)
            produced += [p for p in shell_lines if p.sigma_barn.max() >= WEAK_LINE * strongest]
    return produced
