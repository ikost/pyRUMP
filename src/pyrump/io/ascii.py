"""ASCII spectrum formats.

Reimplements the plain-text readers and writers of ``rump/rdwr.c``.

Three dialects, distinguished by extension in RUMP but sniffed here:

* **one column** (``.dat``, ``.asc``, ``.ascii``) -- counts in channel order
* **two column** -- ``channel value`` pairs
* **tab-delimited** (``.txt``, ``.xls``) -- as exported by spreadsheets

Any non-numeric leading lines become the identifier string (rdwr.c:1150-1170).

Also handles RUMP's own ``wrascii`` output, which prefixes a keyword header
block terminated by the literal line ``Swallow``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

#: Terminates the header block in ``wrascii`` output.
_SWALLOW = "swallow"


@dataclass(slots=True)
class AsciiSpectrum:
    counts: np.ndarray
    identifier: str = ""
    channels: np.ndarray | None = None
    metadata: dict[str, str] = field(default_factory=dict)


def _numbers(line: str) -> list[float] | None:
    """Parse a line as floats, or None if it is not numeric."""
    parts = line.replace(",", " ").split()
    if not parts:
        return None
    try:
        return [float(p) for p in parts]
    except ValueError:
        return None


def read_ascii(path: str | Path) -> AsciiSpectrum:
    """Read a one- or two-column ASCII spectrum, or ``wrascii`` output.

    The column count is detected from the data rather than the extension, so
    the same function handles every dialect.
    """
    lines = Path(path).read_text(errors="replace").splitlines()

    # wrascii output: keyword header, then "Swallow", then one count per line.
    metadata: dict[str, str] = {}
    for index, line in enumerate(lines):
        if line.strip().lower() == _SWALLOW:
            for header in lines[:index]:
                parts = header.split(None, 1)
                if len(parts) == 2:
                    metadata[parts[0]] = parts[1].strip()
            lines = lines[index + 1 :]
            break

    identifier_parts: list[str] = []
    rows: list[list[float]] = []
    for line in lines:
        values = _numbers(line)
        if values is None:
            # Leading non-numeric lines are the identifier; later ones are junk.
            if not rows:
                identifier_parts.append(line.strip())
            continue
        rows.append(values)

    if not rows:
        raise ValueError(f"{path}: no numeric data found")

    width = max(len(r) for r in rows)
    if width == 1:
        counts = np.array([r[0] for r in rows], dtype=np.float64)
        channels = None
    else:
        channels = np.array([r[0] for r in rows], dtype=np.float64)
        counts = np.array([r[1] for r in rows], dtype=np.float64)

    identifier = metadata.get("Ident", " ".join(identifier_parts).strip()).strip("'\"")
    return AsciiSpectrum(
        counts=counts, identifier=identifier, channels=channels, metadata=metadata
    )


def write_ascii(
    path: str | Path,
    counts: np.ndarray,
    *,
    identifier: str = "",
    two_column: bool = False,
    first_channel: int = 0,
) -> None:
    """Write a one- or two-column ASCII spectrum."""
    counts = np.asarray(counts, dtype=np.float64)
    lines: list[str] = []
    if identifier:
        lines.append(identifier)
    if two_column:
        lines.extend(
            f"{first_channel + i} {value:.6f}" for i, value in enumerate(counts)
        )
    else:
        lines.extend(f"{value:.6f}" for value in counts)
    Path(path).write_text("\n".join(lines) + "\n")
