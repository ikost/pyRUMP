"""The periodic table: element lookup, isotope resolution and reference parsing."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from ..io.atom4 import parse_atom4
from ..io.scoef import parse_pscoef
from ..model.element import Element, ElementRef

#: ``Si+28``, ``28Si`` and ``Si29`` all pin an isotope (atomdo.c:113-152).
_REF_SUFFIX = re.compile(r"^([A-Za-z]+)\+(\d+)$")
_REF_PREFIX = re.compile(r"^(\d+)([A-Za-z]+)$")
_REF_TRAILING = re.compile(r"^([A-Za-z]+)(\d*)$")


@dataclass(slots=True)
class PeriodicTable:
    """Element data merged from ``atom4.dat`` and ``pscoef.dat``."""

    elements: list[Element]
    _by_symbol: dict[str, Element] = field(default_factory=dict, repr=False)

    def __post_init__(self) -> None:
        self._by_symbol = {e.symbol.upper(): e for e in self.elements}

    @classmethod
    def load(cls, atom4_path: str | Path, pscoef_path: str | Path | None = None):
        """Load the table, optionally merging in Ziegler stopping parameters."""
        elements = parse_atom4(atom4_path)
        if pscoef_path is not None:
            ziegler = parse_pscoef(pscoef_path)
            elements = [
                Element(
                    z=e.z,
                    symbol=e.symbol,
                    mass=e.mass,
                    atomic_density=e.atomic_density,
                    isotopes=e.isotopes,
                    ziegler=ziegler.get(e.z),
                )
                for e in elements
            ]
        return cls(elements=elements)

    def __len__(self) -> int:
        return len(self.elements)

    def __iter__(self):
        return iter(self.elements)

    def by_z(self, z: int) -> Element:
        if not 1 <= z <= len(self.elements):
            raise KeyError(f"no element with Z={z}")
        return self.elements[z - 1]

    def by_symbol(self, symbol: str) -> Element:
        try:
            return self._by_symbol[symbol.upper()]
        except KeyError:
            raise KeyError(f"unknown element symbol {symbol!r}") from None

    def parse_ref(self, token: str) -> ElementRef:
        """Parse RUMP's element/isotope syntax into an :class:`ElementRef`.

        Accepts ``Si`` (natural abundance), ``Si+28``, ``28Si`` and ``Si29``.
        Reproduces ``RbsIdent`` (atomdo.c:113-152) but raises instead of printing.
        """
        token = token.strip()
        if not token:
            raise ValueError("empty element reference")

        for pattern, sym_group, iso_group in (
            (_REF_SUFFIX, 1, 2),
            (_REF_PREFIX, 2, 1),
            (_REF_TRAILING, 1, 2),
        ):
            match = pattern.match(token)
            if match:
                symbol = match.group(sym_group)
                raw = match.group(iso_group)
                mass_number = int(raw) if raw else 0
                return ElementRef(z=self.by_symbol(symbol).z, mass_number=mass_number)

        raise ValueError(f"illegal element/isotope specification: {token!r}")

    def real_mass(self, z: int, mass_number: int = 0) -> float:
        """Resolve (Z, mass number) to an exact mass; 0 means natural abundance."""
        return self.by_z(z).real_mass(mass_number)

    def mass_of(self, ref: ElementRef) -> float:
        return self.real_mass(ref.z, ref.mass_number)
