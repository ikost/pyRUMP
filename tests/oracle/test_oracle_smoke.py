"""M0 acceptance: the legacy C binary is drivable and deterministic.

These tests do not exercise pyRUMP at all. They establish that the reference oracle
works, which is the precondition for every numerical comparison in later milestones.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
import driver  # noqa: E402

pytestmark = [
    pytest.mark.oracle,
    pytest.mark.skipif(not driver.available(), reason="legacy C tree not available"),
]

BARE_SI = "reset\nlayer 1\nthick 5000 A\ncomposition Si 1 /"


@pytest.fixture(scope="module")
def fixture_2a() -> Path:
    root = driver.reference_root()
    assert root is not None
    return root / "rump" / "data" / "Fixed" / "2A.rbs"


def test_binary_starts_and_loads_tables():
    oracle = driver.RumpOracle()
    result = oracle.run(["version"])
    assert "RBS Analysis and Simulation Package" in result.transcript
    # All four startup tables must load, or every downstream number is wrong.
    for table in ("atom4.dat", "pscoef.dat", "newstop.kal", "density.tab"):
        assert table in result.transcript, f"{table} did not load"


def test_simulation_produces_a_spectrum(tmp_path, fixture_2a):
    oracle = driver.RumpOracle()
    out = oracle.simulate_to_ascii(BARE_SI, tmp_path / "si.dat", data_file=fixture_2a)
    meta, counts = driver.read_wrascii(out)

    assert len(counts) == 2048, "should inherit the data buffer's channel count"
    assert meta["Spectrum"] == "RBS"

    nonzero = [i for i, c in enumerate(counts) if c > 0]
    assert nonzero, "simulation is empty - SIM 'recal' probably did not fire"
    # A 5000 A Si film gives a plateau, not a delta: expect a wide band of signal.
    assert len(nonzero) > 100
    assert max(counts) > 100


def test_simulation_is_deterministic(tmp_path, fixture_2a):
    """Byte-identical across runs, so frozen reference spectra are meaningful."""
    oracle = driver.RumpOracle()
    a = oracle.simulate_to_ascii(BARE_SI, tmp_path / "a.dat", data_file=fixture_2a)
    b = oracle.simulate_to_ascii(BARE_SI, tmp_path / "b.dat", data_file=fixture_2a)
    _, counts_a = driver.read_wrascii(a)
    _, counts_b = driver.read_wrascii(b)
    assert counts_a == counts_b


def test_thickness_scales_the_plateau(tmp_path, fixture_2a):
    """Sanity check that the engine responds to sample changes as physics requires."""
    oracle = driver.RumpOracle()
    thin = oracle.simulate_to_ascii(
        "reset\nlayer 1\nthick 2000 A\ncomposition Si 1 /",
        tmp_path / "thin.dat",
        data_file=fixture_2a,
    )
    thick = oracle.simulate_to_ascii(
        "reset\nlayer 1\nthick 8000 A\ncomposition Si 1 /",
        tmp_path / "thick.dat",
        data_file=fixture_2a,
    )
    _, cthin = driver.read_wrascii(thin)
    _, cthick = driver.read_wrascii(thick)

    # Same surface edge height, but a wider plateau and more total yield.
    assert sum(cthick) > sum(cthin) * 2
    width_thin = sum(1 for c in cthin if c > 1)
    width_thick = sum(1 for c in cthick if c > 1)
    assert width_thick > width_thin
