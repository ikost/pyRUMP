"""Element and isotope data model.

Mirrors the ``ATOMS``/``ISOTOPE`` structures of the legacy C (``rumpdata.h:60-91``),
which merge two independent sources into one record:

* ``atom4.dat``  - symbol, average mass, atomic density, isotope table
* ``pscoef.dat`` - Ziegler/ZBL stopping parameters (``z*`` fields)

They are kept as separate dataclasses here because the Ziegler block is only
meaningful for Z <= 92 and is loaded by a different parser.
"""

from __future__ import annotations

from dataclasses import dataclass, field

#: Maximum isotopes stored per element (``NISOT``, rumpdata.h:63).
NISOT = 6

#: RUMP models Mylar as a pseudo-element at Z=93 so that stopping foils can be
#: expressed as a normal layer composition. It has no Ziegler coefficients.
MYLAR_Z = 93


@dataclass(frozen=True, slots=True)
class Isotope:
    """One entry of an element's isotope table."""

    mass: float
    """Exact isotopic mass in amu."""

    fraction: float
    """Natural abundance as a fraction (not a percentage)."""

    @property
    def mass_number(self) -> int:
        return round(self.mass)


@dataclass(frozen=True, slots=True)
class ZieglerParameters:
    """Per-element ZBL stopping parameters from ``pscoef.dat``."""

    most_abundant_mass_number: int
    mass_most_abundant: float
    mass_average: float
    density_g_cm3: float
    atomic_density_e22: float
    """Target atomic density in 1e22 atoms/cm^3."""

    fermi_velocity: float
    """Fermi velocity of the solid, in units of the Bohr velocity v0."""

    lambda_screening: float
    """Lambda screening factor for heavy ions (``zlfctr``)."""

    proton_coefficients: tuple[float, ...]
    """The 8 Andersen-Ziegler proton stopping coefficients (``zpcoef``)."""


@dataclass(frozen=True, slots=True)
class Element:
    z: int
    symbol: str
    mass: float
    """Average (natural-abundance) atomic mass in amu."""

    atomic_density: float
    """Atomic density in atoms/cm^3.

    Note the units: ``atom4.dat`` stores this absolutely (e.g. 4.99e22 for Si),
    whereas the compound table ``density.tab`` uses 1e23 at/cm^3. See
    :mod:`pyrump.atomic.density`.
    """

    isotopes: tuple[Isotope, ...] = ()
    ziegler: ZieglerParameters | None = None

    def real_mass(self, mass_number: int = 0) -> float:
        """Resolve an integer mass number to an exact isotopic mass.

        Reproduces ``RbsGetRealMass`` (atomdo.c:172-197): ``mass_number == 0`` means
        natural abundance, an entry within 0.5 amu wins, and anything else falls back
        to the integer itself rather than raising.
        """
        if mass_number == 0:
            return self.mass
        for iso in self.isotopes:
            if abs(mass_number - iso.mass) < 0.5:
                return iso.mass
        return float(mass_number)

    @property
    def is_pseudo(self) -> bool:
        """True for RUMP's Mylar pseudo-element, which has no real Z."""
        return self.z == MYLAR_Z


@dataclass(slots=True)
class ElementRef:
    """An element optionally pinned to a specific isotope.

    RUMP writes these as ``Si``, ``28Si`` or ``Si+28``; ``mass_number == 0`` means
    natural abundance.
    """

    z: int
    mass_number: int = 0
    scale: float = 1.0
    """Per-element stopping-power fudge factor (``atom[].scale``)."""

    tags: dict[str, str] = field(default_factory=dict)
