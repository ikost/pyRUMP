"""The shipped PIXE atomic-data tables and their loader.

Three kinds of check: physical invariants every table must satisfy; that the
loader returns exactly what the tables say; and -- only where xraylib is
installed -- that the tables are what ``tools/pixe_atomic.py`` generates and
that interpolation reproduces xraylib itself.
"""

from __future__ import annotations

import csv
import gzip
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

from pyrump.pixe.atomic import DATA_DIR, MAJOR_GROUPS, AtomicData, atomic_data

ROOT = Path(__file__).resolve().parents[2]
SHELLS = ["K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5"]
needs_xraylib = pytest.mark.skipif(
    importlib.util.find_spec("xraylib") is None, reason="xraylib not installed"
)


@pytest.fixture(scope="module")
def data() -> AtomicData:
    return atomic_data()


def _table(name: str) -> list[dict[str, str]]:
    opener = gzip.open if name.endswith(".gz") else open
    with opener(DATA_DIR / name, "rt", encoding="utf-8") as handle:
        return list(csv.DictReader(line for line in handle if not line.startswith("#")))


# -- physical invariants ------------------------------------------------------


def test_every_line_lies_below_its_own_edge(data):
    """A line fills a vacancy in its initial subshell, so its photon carries
    less than that subshell's binding energy."""
    for z in range(6, 93):
        for line in data.lines(z):
            assert line.energy_keV < data.edge(z, line.shell), (z, line.line)


def test_edges_are_ordered_and_yields_are_probabilities(data):
    for z in range(30, 93):  # every K-M5 subshell is bound from Zn on
        edges = [data.edge(z, shell) for shell in SHELLS]
        assert edges == sorted(edges, reverse=True), z
        for shell in SHELLS:
            assert 0.0 <= data.fluorescence_yield(z, shell) <= 1.0


def test_radiative_rates_of_a_subshell_add_up_to_at_most_one(data):
    for z in range(6, 93):
        totals: dict[str, float] = {}
        for line in data.lines(z):
            totals[line.shell] = totals.get(line.shell, 0.0) + line.rate
        assert all(total <= 1.0 + 1e-3 for total in totals.values()), (z, totals)


def test_coster_kronig_probabilities(data):
    """Direct f_ij only: each in [0, 1], and an L1 vacancy cannot move on
    with more than certainty."""
    for row in _table("coster_kronig.csv"):
        assert 0.0 <= float(row["probability"]) <= 1.0
    for z in range(30, 93):
        assert data.coster_kronig(z, "L12") + data.coster_kronig(z, "L13") <= 1.0
    assert data.coster_kronig(78, "M45") > 0
    assert data.coster_kronig(6, "L13") == 0.0  # not tabulated: closed


def test_attenuation_jumps_up_at_an_edge(data):
    edge = data.edge(14, "K")
    below, above = data.mu(14, [edge * 0.99, edge * 1.01])
    assert above > 5 * below


# -- the loader returns the tables ------------------------------------------


def test_loader_matches_the_tables(data):
    lines = {(int(r["Z"]), r["line"]): r for r in _table("xray_lines.csv")}
    pt_alpha = lines[78, "M5N7"]
    loaded = next(line for line in data.lines(78) if line.line == "M5N7")
    assert (loaded.shell, loaded.energy_keV, loaded.rate) == (
        pt_alpha["shell"], float(pt_alpha["energy_keV"]), float(pt_alpha["rate"])
    )
    shells = {(int(r["Z"]), r["shell"]): r for r in _table("xray_shells.csv")}
    row = shells[78, "L3"]
    assert data.edge(78, "L3") == float(row["edge_keV"])
    assert data.fluorescence_yield(78, "L3") == float(row["fluor_yield"])
    points = [r for r in _table("mass_attenuation.csv.gz") if r["Z"] == "29"]
    energies = np.array([float(r["energy_keV"]) for r in points])
    np.testing.assert_allclose(
        data.mu(29, energies), [float(r["mu_cm2_g"]) for r in points], rtol=1e-12
    )


def test_attenuation_outside_the_table_is_refused(data):
    with pytest.raises(ValueError, match="outside"):
        data.mu(14, 0.05)
    with pytest.raises(KeyError):
        data.mu(99, 5.0)


def test_line_groups(data):
    groups = {g.label: g.energy_keV for g in data.line_groups(29, major_only=False)}
    lines = {line.line: line.energy_keV for line in data.lines(29)}
    assert lines["KL2"] <= groups["Kα"] <= lines["KL3"]
    assert set(groups) >= {"Kα", "Kβ", "Lα"}
    major = {g.label for g in data.line_groups(78)}
    assert major <= MAJOR_GROUPS and {"Mα", "Lα"} <= major
    assert "Mζ" in {g.label for g in data.line_groups(78, major_only=False)}


# -- against xraylib, where installed ---------------------------------------


@needs_xraylib
def test_shipped_tables_are_what_the_generator_writes(tmp_path):
    """Regenerating must reproduce the shipped tables row for row (the
    header's date aside), so they can't drift from tools/pixe_atomic.py."""
    sys.path.insert(0, str(ROOT / "tools"))
    try:
        import pixe_atomic
    finally:
        sys.path.pop(0)
    pixe_atomic.write(tmp_path)
    for name in ("xray_lines.csv", "xray_shells.csv", "coster_kronig.csv",
                 "mass_attenuation.csv.gz"):
        opener = gzip.open if name.endswith(".gz") else open

        def rows(path):
            with opener(path, "rt", encoding="utf-8") as handle:
                return [line for line in handle if not line.startswith("#")]

        assert rows(tmp_path / name) == rows(DATA_DIR / name), name


@needs_xraylib
def test_interpolated_attenuation_reproduces_xraylib(data):
    import xraylib

    rng = np.random.default_rng(3)
    for z in (1, 6, 14, 25, 44, 48, 77, 78, 92):
        energies = np.exp(rng.uniform(np.log(0.11), np.log(99.0), 500))
        reference = np.array([xraylib.CS_Total(z, float(e)) for e in energies])
        np.testing.assert_allclose(data.mu(z, energies), reference, rtol=2e-3)
