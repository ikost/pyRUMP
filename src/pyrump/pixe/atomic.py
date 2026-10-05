"""X-ray atomic data for PIXE: lines, subshells, Coster-Kronig, attenuation.

Read from the tables in ``pyrump/data/pixe/``, which ``tools/pixe_atomic.py``
generates once from xraylib -- pyRUMP itself never imports xraylib. See that
script for what each table holds and why ``FLP13`` is left out.
"""

from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "pixe"


@dataclass(frozen=True, slots=True)
class XrayLine:
    z: int
    line: str
    """IUPAC name, e.g. ``KL3`` or ``M5N7``."""

    shell: str
    """The subshell whose vacancy the line fills: ``K``, ``L1`` ... ``M5``."""

    energy_keV: float
    rate: float
    """Share of the radiative decays of a ``shell`` vacancy that give this line."""


#: Lines grouped under their usual (Siegbahn) names, for markers and labels.
#: Each group is the lines a detector cannot separate.
LINE_GROUPS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("Kα", ("KL2", "KL3")),
    ("Kβ", ("KM2", "KM3")),
    ("Ll", ("L3M1",)),
    ("Lα", ("L3M4", "L3M5")),
    ("Lβ1", ("L2M4",)),
    ("Lβ2", ("L3N4", "L3N5")),
    ("Lγ1", ("L2N4",)),
    ("Mζ", ("M4N2", "M5N3")),
    ("Mα", ("M5N6", "M5N7")),
    ("Mβ", ("M4N6",)),
    ("Mγ", ("M3N5",)),
)

#: The groups marked by default -- the strongest line or two of each shell.
MAJOR_GROUPS = frozenset({"Kα", "Kβ", "Lα", "Lβ1", "Lβ2", "Lγ1", "Mα", "Mβ"})


@dataclass(frozen=True, slots=True)
class LineGroup:
    label: str
    energy_keV: float
    """Rate-weighted mean energy of the group's lines."""


def _rows(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        yield from csv.DictReader(line for line in handle if not line.startswith("#"))


class AtomicData:
    """The tables, indexed by atomic number."""

    def __init__(self, directory: Path = DATA_DIR):
        self._lines: dict[int, list[XrayLine]] = defaultdict(list)
        for row in _rows(directory / "xray_lines.csv"):
            z = int(row["Z"])
            self._lines[z].append(
                XrayLine(z, row["line"], row["shell"], float(row["energy_keV"]), float(row["rate"]))
            )

        self._shells: dict[tuple[int, str], tuple[float, float, float]] = {}
        for row in _rows(directory / "xray_shells.csv"):
            self._shells[int(row["Z"]), row["shell"]] = (
                float(row["edge_keV"]), float(row["fluor_yield"]), float(row["jump_factor"])
            )

        self._coster_kronig: dict[tuple[int, str], float] = {
            (int(row["Z"]), row["transition"]): float(row["probability"])
            for row in _rows(directory / "coster_kronig.csv")
        }

        energies: dict[int, list[float]] = defaultdict(list)
        mu: dict[int, list[float]] = defaultdict(list)
        for row in _rows(directory / "mass_attenuation.csv.gz"):
            z = int(row["Z"])
            energies[z].append(float(row["energy_keV"]))
            mu[z].append(float(row["mu_cm2_g"]))
        self._log_energy = {z: np.log(np.array(e)) for z, e in energies.items()}
        self._log_mu = {z: np.log(np.array(m)) for z, m in mu.items()}

    # -- lines --------------------------------------------------------------

    def lines(self, z: int) -> tuple[XrayLine, ...]:
        """Every tabulated line of element ``z``, K first."""
        return tuple(self._lines.get(z, ()))

    def line_groups(self, z: int, *, major_only: bool = True) -> list[LineGroup]:
        """The named line groups element ``z`` has, at their mean energies.

        The weighting uses the radiative rates only, which is exact within a
        group fed by one subshell (Kα, Lα, Mα ...) and close enough for a
        label otherwise.
        """
        by_name = {line.line: line for line in self.lines(z)}
        groups = []
        for label, names in LINE_GROUPS:
            if major_only and label not in MAJOR_GROUPS:
                continue
            members = [by_name[n] for n in names if n in by_name]
            weight = sum(line.rate for line in members)
            if weight > 0:
                energy = sum(line.rate * line.energy_keV for line in members) / weight
                groups.append(LineGroup(label, energy))
        return groups

    # -- subshells ----------------------------------------------------------

    def _shell(self, z: int, shell: str) -> tuple[float, float, float]:
        try:
            return self._shells[z, shell]
        except KeyError:
            raise KeyError(f"no {shell} subshell data for Z={z}") from None

    def edge(self, z: int, shell: str) -> float:
        """Absorption edge (binding energy) in keV."""
        return self._shell(z, shell)[0]

    def fluorescence_yield(self, z: int, shell: str) -> float:
        return self._shell(z, shell)[1]

    def jump_factor(self, z: int, shell: str) -> float:
        return self._shell(z, shell)[2]

    def coster_kronig(self, z: int, transition: str) -> float:
        """Direct Coster-Kronig probability, e.g. ``L13`` for f13; 0 where
        the table has none (the transition is energetically closed)."""
        return self._coster_kronig.get((z, transition), 0.0)

    # -- attenuation --------------------------------------------------------

    def mu(self, z: int, energy_keV) -> np.ndarray:
        """Mass attenuation coefficient in cm^2/g, log-log interpolated.

        Edges are sampled on both sides in the table, so interpolation never
        smooths one away. Energies outside the table are refused rather than
        extrapolated.
        """
        if z not in self._log_energy:
            raise KeyError(f"no attenuation data for Z={z}")
        log_e = np.log(np.asarray(energy_keV, dtype=np.float64))
        grid = self._log_energy[z]
        if np.any(log_e < grid[0]) or np.any(log_e > grid[-1]):
            raise ValueError(
                f"energy outside the attenuation table "
                f"({np.exp(grid[0]):g}-{np.exp(grid[-1]):g} keV)"
            )
        return np.exp(np.interp(log_e, grid, self._log_mu[z]))


@lru_cache(maxsize=1)
def atomic_data() -> AtomicData:
    """The shipped tables, read once per process."""
    return AtomicData()
