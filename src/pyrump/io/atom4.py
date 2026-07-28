"""Parser for RUMP's ``atom4.dat`` master element/isotope table.

Reimplements ``RbsLoadAtomicData`` (atomio.c:97-213). The file layout is:

1. ``#`` comment block
2. a line with the element count (93 in the shipped file)
3. element symbols packed into **fixed two-character cells**, wrapped over
   several lines
4. ``count`` average masses, whitespace-delimited across lines
5. ``count`` atomic densities in atoms/cm^3
6. a line ``n_with_isotopes  isotopes_per_element``
7. ``n_with_isotopes`` blocks of ``NISOT`` ``(fraction, mass)`` pairs

Element 93 is the Mylar pseudo-element, and the last ten elements carry no
isotope data.
"""

from __future__ import annotations

from pathlib import Path

from ..model.element import NISOT, Element, Isotope


def _symbols(stream: "_Tokenizer", count: int) -> list[str]:
    """Read `count` symbols from fixed two-character cells.

    The C walks a character pointer two bytes at a time: a leading space means a
    one-letter symbol, otherwise the second character is lower-cased
    (atomio.c:139-158). Chunking each line into 2-char cells and title-casing is
    equivalent and far clearer.
    """
    out: list[str] = []
    while len(out) < count:
        line = stream.raw_line()
        if line is None:
            raise ValueError("atom4: ran out of lines reading element symbols")
        # Do not strip: the leading space of " H" is part of the first cell.
        line = line.rstrip("\n")
        for i in range(0, len(line) - 1, 2):
            if len(out) == count:
                break
            cell = line[i : i + 2].strip()
            if cell:
                out.append(cell.capitalize())
    return out


class _Tokenizer:
    """Whitespace tokenizer that can also yield whole raw lines.

    The C mixes ``fgets`` (for symbols) with ``fscanf`` (for numbers), so the
    parser needs both views of the same stream.
    """

    def __init__(self, text: str):
        self._lines = text.splitlines(keepends=True)
        self._i = 0
        self._pending: list[str] = []

    def raw_line(self) -> str | None:
        if self._pending:
            raise RuntimeError("cannot take a raw line mid-token")
        if self._i >= len(self._lines):
            return None
        line = self._lines[self._i]
        self._i += 1
        return line

    def token(self) -> str:
        while not self._pending:
            line = self.raw_line()
            if line is None:
                raise ValueError("atom4: unexpected end of file")
            self._pending = line.split()
        return self._pending.pop(0)

    def number(self) -> float:
        return float(self.token())

    def integer(self) -> int:
        return int(self.token())

    def skip_comments(self) -> None:
        while self._i < len(self._lines) and self._lines[self._i].lstrip().startswith("#"):
            self._i += 1


def parse_atom4(path: str | Path) -> list[Element]:
    """Parse ``atom4.dat`` into a list of :class:`Element`, ordered by Z."""
    stream = _Tokenizer(Path(path).read_text())
    stream.skip_comments()

    count = stream.integer()
    if not 0 < count <= 115:
        raise ValueError(f"atom4: implausible element count {count}")

    symbols = _symbols(stream, count)
    masses = [stream.number() for _ in range(count)]
    densities = [stream.number() for _ in range(count)]

    n_with_isotopes, per_element = stream.integer(), stream.integer()
    if per_element != NISOT:
        raise ValueError(f"atom4: isotope table mismatch, expected {NISOT}, got {per_element}")

    isotopes: list[tuple[Isotope, ...]] = [() for _ in range(count)]
    for i in range(n_with_isotopes):
        entries = []
        for _ in range(NISOT):
            # File order is (fraction, mass) even though the C struct declares
            # mass first -- see the fscanf argument order at atomio.c:174.
            fraction = stream.number()
            mass = stream.number()
            if fraction > 0.0:
                entries.append(Isotope(mass=mass, fraction=fraction))
        isotopes[i] = tuple(entries)

    return [
        Element(
            z=i + 1,
            symbol=symbols[i],
            mass=masses[i],
            atomic_density=densities[i],
            isotopes=isotopes[i],
        )
        for i in range(count)
    ]
