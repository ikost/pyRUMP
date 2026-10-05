"""Inner-shell ionisation cross sections: pyRUMP's ECPSSR tables.

``tools/ecpssr/make_tables.py`` computes them once (Brandt & Lapicki's
ECPSSR, pyRUMP's own implementation) for protons and 4He; pyRUMP ships the
result in ``data/pixe/ecpssr_*.csv.gz``. Values in between the tabulated
energies are interpolated linearly in ln(sigma) against ln(E); energies
outside the table are refused, never extrapolated.
"""

from __future__ import annotations

import csv
import gzip
from functools import lru_cache
from pathlib import Path

import numpy as np

from .atomic import DATA_DIR

#: The ions there are tables for: name -> (Z, mass number).
IONS = {"H1": (1, 1), "He4": (2, 4)}


def ion_name(z: int, mass: float) -> str:
    """The table for a beam of charge ``z`` and mass ``mass`` (amu)."""
    for name, (ion_z, mass_number) in IONS.items():
        if z == ion_z and round(mass) == mass_number:
            return name
    raise ValueError(
        f"no PIXE cross sections for a beam of Z={z}, mass {mass:g}: "
        "only protons and 4He are tabulated"
    )


class IonisationTable:
    """ECPSSR ionisation cross sections for one ion, by element and subshell."""

    def __init__(self, path: Path):
        curves: dict[tuple[int, str], tuple[list[float], list[float]]] = {}
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for row in csv.DictReader(line for line in handle if not line.startswith("#")):
                energies, sigmas = curves.setdefault((int(row["Z"]), row["shell"]), ([], []))
                energies.append(float(row["E_MeV"]))
                sigmas.append(float(row["sigma_barn"]))
        self._curves = {}
        for key, (energies, sigmas) in curves.items():
            e = np.array(energies)
            s = np.array(sigmas)
            order = np.argsort(e)
            e, s = e[order], s[order]
            # A few heavy-element K-shell values at the lowest energies are
            # exactly zero (the energy-loss factor closes the channel); they
            # stay zero rather than break the logarithm.
            positive = s > 0
            self._curves[key] = (np.log(e), np.log(np.where(positive, s, 1.0)), positive)
        self.e_min = min(np.exp(c[0][0]) for c in self._curves.values())
        self.e_max = max(np.exp(c[0][-1]) for c in self._curves.values())

    def shells(self, z: int) -> list[str]:
        """The subshells of ``z`` there are cross sections for."""
        return [shell for (zz, shell) in self._curves if zz == z]

    def sigma(self, z: int, shell: str, energy_MeV) -> np.ndarray:
        """Ionisation cross section in barn at the ion's total lab energy."""
        try:
            log_e, log_s, positive = self._curves[z, shell]
        except KeyError:
            raise KeyError(f"no {shell} ionisation cross sections for Z={z}") from None
        x = np.log(np.asarray(energy_MeV, dtype=np.float64))
        if np.any(x < log_e[0] - 1e-12) or np.any(x > log_e[-1] + 1e-12):
            raise ValueError(
                f"ion energy outside the cross-section table "
                f"({np.exp(log_e[0]):g}-{np.exp(log_e[-1]):g} MeV)"
            )
        sigma = np.exp(np.interp(x, log_e, log_s))
        # Zero below the last tabulated zero.
        if not positive.all():
            last_zero = log_e[~positive][-1]
            sigma = np.where(x <= last_zero, 0.0, sigma)
        return sigma


@lru_cache(maxsize=None)
def ionisation_table(ion: str) -> IonisationTable:
    """The shipped table for ``ion`` (``"H1"`` or ``"He4"``), read once."""
    if ion not in IONS:
        raise ValueError(f"unknown ion {ion!r}: expected one of {', '.join(IONS)}")
    return IonisationTable(DATA_DIR / f"ecpssr_{ion}.csv.gz")
