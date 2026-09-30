"""ASCII spectrum formats.

Reimplements the plain-text readers and writers of ``rump/rdwr.c``.

Three dialects, distinguished by extension in RUMP but sniffed here:

* **one column** (``.dat``, ``.asc``, ``.ascii``) -- counts in channel order
* **two column** -- ``channel value`` pairs
* **tab-delimited** (``.txt``, ``.xls``) -- as exported by spreadsheets

Any non-numeric leading lines become the identifier string (rdwr.c:1150-1170).

Also handles RUMP's own ``wrascii`` output, which prefixes a keyword header
block terminated by the literal line ``Swallow``.

:func:`write_columns` is pyRUMP's own, write-only addition: a ``#``-commented
header over named columns, for spreadsheets and plotting tools rather than
for RUMP to read back.
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
    header: str | None = None,
) -> None:
    """Write a one- or two-column ASCII spectrum.

    With ``header``, writes RUMP's own ``WRASCII`` dialect instead: that text
    verbatim, then the literal line ``Swallow``, then one count per line, a
    trailing blank line, and nothing else -- ``identifier``/``two_column``/
    ``first_channel`` are ignored, since the header already carries the
    identifier and RUMP's ``WRASCII`` has no two-column form.
    """
    counts = np.asarray(counts, dtype=np.float64)
    if header is not None:
        stripped = header.rstrip("\n")
        body = "\n".join(f"{value:.6f}" for value in counts)
        Path(path).write_text(f"{stripped}\nSwallow\n{body}\n\n")
        return

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


def write_columns(
    path: str | Path,
    columns: dict[str, np.ndarray],
    *,
    comments: list[str] | None = None,
    delimiter: str = "\t",
    formats: dict[str, str] | None = None,
) -> None:
    """Write named columns under ``#``-prefixed comment lines.

    One row per element of the columns, which must all be the same length;
    the column names form the first uncommented line. ``formats`` maps a
    column name to a format spec (default ``.6f``). Origin, Excel, gnuplot,
    ``numpy.loadtxt`` and ``pandas.read_csv(comment="#")`` all skip the
    comments.
    """
    formats = formats or {}
    arrays = [np.asarray(values) for values in columns.values()]
    lengths = {array.size for array in arrays}
    if len(lengths) > 1:
        raise ValueError(f"columns differ in length: {sorted(lengths)}")
    specs = [formats.get(name, ".6f") for name in columns]

    lines = [f"# {line}".rstrip() for line in comments or []]
    lines.append(delimiter.join(columns))
    for row in zip(*arrays):
        lines.append(delimiter.join(format(value, spec) for value, spec in zip(row, specs)))
    Path(path).write_text("\n".join(lines) + "\n")
