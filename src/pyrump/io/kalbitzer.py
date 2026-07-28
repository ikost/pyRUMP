"""Parser for ``newstop.kal``, RUMP's Konac/Kalbitzer stopping fits.

Reimplements ``RbsLoadKalbitzerData`` (stopping.c:743-787): each data line carries
exactly 12 whitespace-separated fields::

    z1  m1  z2  scaling  s  a0 a1 a2 a3 a4 a5  beta

``#`` introduces a comment. The shipped file holds 8 records covering H, D, 3He
and 4He on carbon and silicon only.

Reference given in the file itself: G. Konac, S. Kalbitzer, Ch. Klatt, D. Niemann,
R. Stoll, *Nucl. Instr. Meth.* B136-138 (1998) 159-165.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

_N_FIELDS = 12


@dataclass(frozen=True, slots=True)
class KonacEntry:
    """One (projectile, target) Konac fit."""

    z1: int
    m1: int
    """Projectile mass *number*; matching is on ``round(m1)`` (stopping.c:800)."""

    z2: int
    scaling: float
    s: float
    a: tuple[float, float, float, float, float, float]
    beta: float

    def matches(self, z1: int, m1: float, z2: int) -> bool:
        return self.z1 == z1 and self.z2 == z2 and int(m1 + 0.5) == self.m1


def parse_kalbitzer(path: str | Path) -> list[KonacEntry]:
    """Parse ``newstop.kal`` into a list of :class:`KonacEntry`."""
    entries: list[KonacEntry] = []
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        fields = line.split()
        if len(fields) < _N_FIELDS:
            continue  # the C requires exactly 12 and skips anything shorter
        values = [float(v) for v in fields[:_N_FIELDS]]
        entries.append(
            KonacEntry(
                z1=int(values[0]),
                m1=int(values[1]),
                z2=int(values[2]),
                scaling=values[3],
                s=values[4],
                a=tuple(values[5:11]),  # type: ignore[arg-type]
                beta=values[11],
            )
        )
    return entries
