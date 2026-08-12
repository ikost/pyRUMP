"""Helpers shared by the batch CLI and the interactive shell.

Locating the legacy data tables and building the periodic table / stopping
registry / density table is identical for both front ends, and the registry is
expensive enough that ``docs/usage.md`` warns to build it once and reuse it.
"""

from __future__ import annotations

import os
from pathlib import Path


def data_dir(explicit: str | None = None) -> Path:
    """Locate the legacy data tables.

    Searched in order: an explicit path, ``$PYRUMP_DATA``,
    ``$PYRUMP_C_REFERENCE/rump/data``, then ``./C-code/rump/data``. A directory
    counts only if it actually holds ``atom4.dat``.
    """
    candidates = []
    if explicit:
        candidates.append(Path(explicit))
    if os.environ.get("PYRUMP_DATA"):
        candidates.append(Path(os.environ["PYRUMP_DATA"]))
    if os.environ.get("PYRUMP_C_REFERENCE"):
        candidates.append(Path(os.environ["PYRUMP_C_REFERENCE"]) / "rump" / "data")
    candidates.append(Path.cwd() / "C-code" / "rump" / "data")
    for path in candidates:
        if (path / "atom4.dat").is_file():
            return path
    raise SystemExit(
        "Could not find the data tables (atom4.dat, pscoef.dat, newstop.kal).\n"
        "Pass --data DIR, or set PYRUMP_DATA."
    )


def load_tables(data: Path):
    """Build the periodic table, stopping registry and compound densities."""
    from pyrump.atomic.density import DensityTable
    from pyrump.atomic.tables import PeriodicTable
    from pyrump.io.kalbitzer import parse_kalbitzer
    from pyrump.stopping.kalbitzer import KalbitzerStopping
    from pyrump.stopping.registry import StoppingRegistry
    from pyrump.stopping.ziegler import ZieglerStopping

    table = PeriodicTable.load(data / "atom4.dat", data / "pscoef.dat")
    registry = StoppingRegistry(
        table.elements,
        kalbitzer=KalbitzerStopping(parse_kalbitzer(data / "newstop.kal"), table.elements),
        ziegler=ZieglerStopping(table.elements),
    )
    densities = DensityTable.load(data / "density.tab")
    return table, registry, densities


def read_spectrum(path: str | Path):
    """Read a spectrum by extension: ``.rbs`` and friends binary, else ASCII."""
    from pyrump.io.ascii import read_ascii
    from pyrump.io.rbs import read_rbs

    path = Path(path)
    if path.suffix.lower() in (".rbs", ".rump", ".frs", ".fres", ".pixe"):
        return read_rbs(path)
    return read_ascii(path)


def resolve_beam(table, spec: str) -> tuple[int, float]:
    """Resolve a beam specification such as ``He``, ``4He`` or ``Si+28``.

    Returns ``(z, mass)``. Trailing charge-state plus signs (``4He++``) are the
    caller's business -- RUMP carries them separately.
    """
    element = table.by_symbol(spec) if spec.isalpha() else None
    if element is not None:
        return element.z, element.mass
    reference = table.parse_ref(spec)
    return reference.z, table.real_mass(reference.z, getattr(reference, "mass_number", 0))
