"""Shared data-table lookup for tests.

Prefers the legacy ``C-code/`` tree when present (so local development
exercises the same tables the oracle tests compare against), falls back to
the tables bundled with the package otherwise, and returns ``None`` so
callers can skip cleanly if neither is available.

``tests/unit/test_io_adt.py`` and ``tests/unit/test_io_rbs.py`` need files
that aren't bundled (``.adt`` cross-section tables, and the ``Fixed/``
sample spectra) and keep their own ``C-code``-only lookups rather than using
this helper.
"""

from __future__ import annotations

from pathlib import Path

from pyrump.cli._common import data_dir as _resolve_data_dir


def data_dir() -> Path | None:
    """Locate the data tables, or ``None`` if unavailable (tests should skip)."""
    try:
        return _resolve_data_dir()
    except SystemExit:
        return None
