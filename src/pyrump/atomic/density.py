"""Compound density table and thickness units.

Reimplements ``SimLoadDensityTable`` (sim2.c:1415-1483) -- note the live parser is
in ``sim2.c``, not the dead ``tables.c``, which does not compile.

RUMP folds two different things into one lookup table, because both appear in the
same syntactic slot of a ``THICKNESS`` command (``thick 1000 A`` vs ``thick 151 ITO``):

* **thickness units** - five hard-coded entries (A, nm, um, /CM2, M/CM2)
* **compound densities** - name -> atomic density, read from ``density.tab``

Densities here are in **1e23 atoms/cm^3**, unlike :attr:`Element.atomic_density`
from ``atom4.dat`` which is absolute atoms/cm^3.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path

#: Fallback density when a compound is unknown: silicon, 0.4997e23 at/cm^3
#: (creatr.c). RUMP silently substitutes this rather than failing.
SILICON_DENSITY = 0.4997


class ThicknessKind(Enum):
    """What a thickness magnitude is measured in."""

    ANGSTROMS = "angstroms"
    ATOMIC = "atomic"
    """1e15 atoms/cm^2."""

    MOLECULAR = "molecular"
    """1e15 molecules/cm^2."""


@dataclass(frozen=True, slots=True)
class ThicknessUnit:
    name: str
    kind: ThicknessKind
    scale: float


#: The five permanent units prepended to every density table (sim2.c:1421-1427).
STATIC_UNITS: tuple[ThicknessUnit, ...] = (
    ThicknessUnit("A", ThicknessKind.ANGSTROMS, 1.0),
    ThicknessUnit("nm", ThicknessKind.ANGSTROMS, 10.0),
    ThicknessUnit("um", ThicknessKind.ANGSTROMS, 10000.0),
    ThicknessUnit("/CM2", ThicknessKind.ATOMIC, 1.0),
    ThicknessUnit("M/CM2", ThicknessKind.MOLECULAR, 1.0),
)


def parse_density_table(path: str | Path) -> dict[str, float]:
    """Parse ``density.tab`` into ``{compound_name: density_1e23_at_cm3}``.

    Blank lines and ``#`` comments are skipped; each remaining line is
    ``name  density``. Lookup is case-insensitive in RUMP, so keys are upper-cased.
    """
    out: dict[str, float] = {}
    for line in Path(path).read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) < 2:
            continue
        try:
            out[parts[0].upper()] = float(parts[1])
        except ValueError:
            continue
    return out


@dataclass(slots=True)
class DensityTable:
    """Combined unit + compound-density lookup, matching RUMP's single namespace."""

    compounds: dict[str, float]
    units: tuple[ThicknessUnit, ...] = STATIC_UNITS

    @classmethod
    def load(cls, path: str | Path) -> "DensityTable":
        return cls(compounds=parse_density_table(path))

    def unit(self, name: str) -> ThicknessUnit | None:
        for unit in self.units:
            if unit.name.upper() == name.upper():
                return unit
        return None

    def density(self, name: str, default: float = SILICON_DENSITY) -> float:
        """Look up a compound density in 1e23 at/cm^3, falling back to silicon."""
        return self.compounds.get(name.upper(), default)

    def __contains__(self, name: str) -> bool:
        return name.upper() in self.compounds or self.unit(name) is not None


def layer_atomic_density(
    composition: "np.ndarray", atomic_densities: "np.ndarray"
) -> float:
    """Average atomic density of a mixture, in 1e23 at/cm^3.

    Port of the ``IMPROVED`` branch of creatr.c:606-625, the post-1/97 default.
    Densities are combined as a composition-weighted average of **inverse**
    density -- "the idea of hard ball packing with weighted sum of cm^3/atom
    instead of atoms/cm^3" (creatr.c:589-592) -- then inverted:

    .. math::
        \\rho = \\left(\\frac{\\sum_i x_i / \\rho_i}{\\sum_i x_i}\\right)^{-1}

    Falls back to silicon (:data:`SILICON_DENSITY`) for an empty or degenerate
    composition, exactly as the C does.

    The pre-1/97 ``COMPATIBLE`` mode averaged densities directly; it is not
    implemented because nothing shipped since 1997 uses it.
    """
    import numpy as np

    composition = np.asarray(composition, dtype=np.float64)
    atomic_densities = np.asarray(atomic_densities, dtype=np.float64)

    total = composition.sum()
    usable = (composition != 0) & (atomic_densities > 0)
    if total <= 0 or not np.any(usable):
        return SILICON_DENSITY

    inverse = (composition[usable] / atomic_densities[usable]).sum() / total
    if inverse <= 0:
        return SILICON_DENSITY
    return float((1.0 / inverse) / 1e23)
