"""Parser for RUMP's ``pscoef.dat`` Ziegler/ZBL stopping coefficient table.

Reimplements ``zread1`` (ziegler.c:81-129). The file is the TRIM/SRIM ``SCOEF.DAT``
table split into two 92-row blocks (Z = 1..92, so no Mylar):

* block 1 - ``Z, most_abundant_mass_number, mass_MAI, mass_average, density_g_cc,
  atomic_density_e22, fermi_velocity, lambda_factor``
* block 2 - ``Z`` followed by the 8 Andersen-Ziegler proton stopping coefficients

Lines starting with ``#`` or ``/*`` are comments. The row index must equal Z.

.. note::
   This table originates from Ziegler, Biersack & Littmark, *The Stopping and Range
   of Ions in Solids* (Pergamon, 1985) and carries no licence notice in the RUMP
   distribution. See README.md, "Notices and citations", for provenance and the
   plan to regenerate it from the published tables before release.
"""

from __future__ import annotations

from pathlib import Path

from ..model.element import ZieglerParameters

#: Ziegler's tables cover Z = 1..92 only; RUMP's ``zcheck`` rejects anything else.
ZIEGLER_MAX_Z = 92

_N_PROTON_COEFFICIENTS = 8


def _data_lines(path: Path):
    for line in path.read_text().splitlines():
        stripped = line.lstrip()
        if not stripped or stripped.startswith("#") or stripped.startswith("/*"):
            continue
        yield line


def parse_pscoef(path: str | Path) -> dict[int, ZieglerParameters]:
    """Parse ``pscoef.dat`` into ``{Z: ZieglerParameters}`` for Z = 1..92."""
    lines = list(_data_lines(Path(path)))
    if len(lines) < 2 * ZIEGLER_MAX_Z:
        raise ValueError(
            f"pscoef: expected {2 * ZIEGLER_MAX_Z} data lines, found {len(lines)}"
        )

    header, coefficients = lines[:ZIEGLER_MAX_Z], lines[ZIEGLER_MAX_Z : 2 * ZIEGLER_MAX_Z]
    out: dict[int, ZieglerParameters] = {}

    for index, (head_line, coef_line) in enumerate(zip(header, coefficients), start=1):
        head = head_line.split()
        z = int(head[0])
        if z != index:
            raise ValueError(f"pscoef: header block out of order at row {index} (got Z={z})")

        coef = coef_line.split()
        if int(coef[0]) != index:
            raise ValueError(f"pscoef: coefficient block out of order at row {index}")
        proton = tuple(float(v) for v in coef[1 : 1 + _N_PROTON_COEFFICIENTS])
        if len(proton) != _N_PROTON_COEFFICIENTS:
            raise ValueError(f"pscoef: row {index} has {len(proton)} proton coefficients")

        out[z] = ZieglerParameters(
            most_abundant_mass_number=int(head[1]),
            mass_most_abundant=float(head[2]),
            mass_average=float(head[3]),
            density_g_cm3=float(head[4]),
            # Stored in the file as 1e22 at/cm^3; the C scales it up on load
            # (ziegler.c:102). We keep the file's units and scale at the point of use.
            atomic_density_e22=float(head[5]),
            fermi_velocity=float(head[6]),
            lambda_screening=float(head[7]),
            proton_coefficients=proton,
        )

    return out
