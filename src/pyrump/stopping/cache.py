"""Session cache of fitted stopping tables.

Port of the lookup half of ``RbsStpfind`` (stopping.c:270-296).

This is stateful behaviour that changes results, so it is not an optimisation --
it is part of the model. RUMP keeps every stopping table it has ever built for
the session and **reuses one whenever the new beam energy merely fits inside the
existing window** (``2*emin <= E <= emax``). Simulating at 3 MeV and then at
2 MeV therefore does *not* refit: the second run uses the 3 MeV window, and its
coefficients differ from what a fresh 2 MeV fit would produce.

Two match qualities, in order of preference:

* **exact** - same Z and mass within 0.2 amu; ``e_scale`` is 1
* **Z only** - same element, different isotope; the table is reused with
  ``e_scale = table_mass / beam_mass``, the Amsel energy-scaling trick that lets
  one table serve 3He/4He or H/D
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .registry import StoppingRegistry
from .table import StoppingTable, StoppingType


@dataclass(slots=True)
class StoppingTableCache:
    """Holds the tables built so far, reproducing RUMP's reuse rules."""

    registry: StoppingRegistry
    targets: list[int]
    kind: StoppingType = StoppingType.SQRT
    tables: list[StoppingTable] = field(default_factory=list)

    def find(
        self, z_beam: int, m_beam: float, e_beam_MeV: float
    ) -> tuple[StoppingTable, float] | None:
        """Return ``(table, e_scale)`` for an existing usable table, else None."""
        z_match: StoppingTable | None = None
        zm_match: StoppingTable | None = None

        for table in self.tables:
            if not table.accepts(z_beam, m_beam, e_beam_MeV):
                continue
            # Among candidates the C prefers the one with the lowest emin,
            # i.e. the widest window at the bottom end.
            if z_match is None or table.emin < z_match.emin:
                z_match = table
            if table.matches_mass(m_beam):
                if zm_match is None or table.emin < zm_match.emin:
                    zm_match = table

        if zm_match is not None:
            return zm_match, 1.0
        if z_match is not None:
            return z_match, z_match.m_beam / m_beam
        return None

    def get(
        self, z_beam: int, m_beam: float, e_beam_MeV: float
    ) -> tuple[StoppingTable, float]:
        """Fetch a usable table, fitting a new one only if none can be reused."""
        found = self.find(z_beam, m_beam, e_beam_MeV)
        if found is not None:
            return found

        table = StoppingTable.build(
            self.registry, z_beam, m_beam, e_beam_MeV, self.targets, self.kind
        )
        self.tables.append(table)
        return table, 1.0

    def clear(self) -> None:
        """Drop all cached tables, forcing a refit on the next request.

        RUMP has no equivalent short of restarting, but tests and batch runs
        need it to get reproducible fits.
        """
        self.tables.clear()
