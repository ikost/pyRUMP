"""The brick: the simulation's intermediate representation.

From the 1985 paper (p. 344):

    Each simulated spectrum is made up of the superimposed contributions from
    each isotope of each sublayer in the sample. Any such contribution will be
    referred to as a **brick**.

A brick is a trapezoid in energy space -- front and back edge energies, front and
back heights -- plus, in the paper's design, an independently computed area. The
shipped C discards that area (see README.md, "RUMP quirks and defects found
while porting"), but pyRUMP carries it
so the parabolic mode can use it.

Stored as a structured array rather than a list of objects: the fill stage
operates on all bricks at once, and this is what makes that vectorisable.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

BRICK_DTYPE = np.dtype(
    [
        ("e_front", np.float64),
        ("e_back", np.float64),
        ("h_front", np.float64),
        ("h_back", np.float64),
        ("area", np.float64),
        ("sig_front", np.float64),
        ("sig_back", np.float64),
    ]
)


@dataclass(slots=True)
class Bricks:
    """A collection of bricks, one per (isotope, slab)."""

    data: np.ndarray

    @classmethod
    def empty(cls, n: int) -> "Bricks":
        return cls(data=np.zeros(n, dtype=BRICK_DTYPE))

    @classmethod
    def from_list(cls, rows: list[tuple]) -> "Bricks":
        return cls(data=np.array(rows, dtype=BRICK_DTYPE))

    def __len__(self) -> int:
        return int(self.data.size)

    def __getitem__(self, item):
        result = self.data[item]
        return Bricks(data=np.atleast_1d(result)) if isinstance(item, slice) else result

    @property
    def e_front(self) -> np.ndarray:
        return self.data["e_front"]

    @property
    def e_back(self) -> np.ndarray:
        return self.data["e_back"]

    @property
    def h_front(self) -> np.ndarray:
        return self.data["h_front"]

    @property
    def h_back(self) -> np.ndarray:
        return self.data["h_back"]

    @property
    def area(self) -> np.ndarray:
        return self.data["area"]

    @property
    def sig_front(self) -> np.ndarray:
        return self.data["sig_front"]

    @property
    def sig_back(self) -> np.ndarray:
        return self.data["sig_back"]

    @property
    def has_straggling(self) -> bool:
        return bool(np.any(self.sig_front) or np.any(self.sig_back))

    def trapezoid_area(self) -> np.ndarray:
        """Area of each brick treated as a plain trapezoid.

        This is what the shipped C's fill stage effectively integrates, as
        distinct from :attr:`area`, which is the paper's exact value.
        """
        return 0.5 * (self.h_front + self.h_back) * (self.e_front - self.e_back)

    def total_yield(self) -> float:
        return float(self.trapezoid_area().sum())

    def to_array(self) -> np.ndarray:
        """Plain ``(n, 7)`` float view, for comparison with oracle captures."""
        return np.column_stack(
            [
                self.e_front,
                self.e_back,
                self.h_front,
                self.h_back,
                self.area,
                self.sig_front,
                self.sig_back,
            ]
        )
