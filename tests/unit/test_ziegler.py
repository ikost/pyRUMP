"""M2 acceptance: ZBL85 stopping must match the C oracle.

Acceptance criterion from the plan: agreement to ~1e-6 relative for Z1 in {1,2},
all Z2 in 1..92, over E in logspace(10 keV, 10 MeV).

Tolerances are set by single precision, not by choice. The oracle is the `float`
build -- the `double` build is unusable because RUMP's table readers use
scanf "%f" against REAL fields (see tests/oracle/oracle.py). So the tables
themselves are float32, and a float64 reimplementation reading the same text
cannot agree better than ~1e-7 relative on inputs, which the ZBL exponentials
amplify.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pytest

from pyrump.atomic.tables import PeriodicTable
from pyrump.stopping.ziegler import MAX_Z, StoppingUnits, ZieglerStopping

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "oracle"))
import oracle as ora  # noqa: E402


def _data_dir() -> Path | None:
    env = os.environ.get("PYRUMP_C_REFERENCE")
    roots = [Path(env)] if env else []
    roots.append(Path(__file__).resolve().parents[2] / "C-code")
    for root in roots:
        if (root / "rump" / "data" / "atom4.dat").is_file():
            return root / "rump" / "data"
    return None


DATA = _data_dir()
pytestmark = pytest.mark.skipif(
    DATA is None or not ora.available(),
    reason="legacy tables or oracle library unavailable",
)

#: Single-precision tables are the floor; see module docstring.
RTOL = 2e-6


@pytest.fixture(scope="module")
def stopping() -> ZieglerStopping:
    assert DATA is not None
    table = PeriodicTable.load(DATA / "atom4.dat", DATA / "pscoef.dat")
    return ZieglerStopping(table.elements)


@pytest.fixture(scope="module")
def oracle() -> ora.Oracle:
    return ora.Oracle.load()


def _compare(stopping, oracle, z1, m1, z2, energies, rtol=RTOL):
    se, sn = stopping(z1, m1, z2, energies)
    for i, energy in enumerate(energies):
        ref_se, ref_sn = oracle.zstop(z1, m1, z2, float(energy))
        assert se[i] == pytest.approx(ref_se, rel=rtol), (
            f"electronic: Z1={z1} Z2={z2} E={energy}"
        )
        assert sn[i] == pytest.approx(ref_sn, rel=rtol), (
            f"nuclear: Z1={z1} Z2={z2} E={energy}"
        )


@pytest.mark.parametrize("z2", [1, 6, 8, 14, 26, 47, 79, 92])
def test_proton_stopping_matches_oracle(stopping, oracle, z2):
    energies = np.logspace(np.log10(10.0), np.log10(10_000.0), 40)
    _compare(stopping, oracle, 1, 1.00797, z2, energies)


@pytest.mark.parametrize("z2", [1, 6, 8, 14, 26, 47, 79, 92])
def test_helium_stopping_matches_oracle(stopping, oracle, z2):
    energies = np.logspace(np.log10(10.0), np.log10(10_000.0), 40)
    _compare(stopping, oracle, 2, 4.0026, z2, energies)


def test_full_z2_sweep_protons(stopping, oracle):
    """The M2 criterion proper: every target element, no exceptions."""
    energies = np.logspace(np.log10(10.0), np.log10(10_000.0), 25)
    for z2 in range(1, MAX_Z + 1):
        _compare(stopping, oracle, 1, 1.00797, z2, energies)


def test_full_z2_sweep_helium(stopping, oracle):
    energies = np.logspace(np.log10(10.0), np.log10(10_000.0), 25)
    for z2 in range(1, MAX_Z + 1):
        _compare(stopping, oracle, 2, 4.0026, z2, energies)


@pytest.mark.parametrize("z1, m1", [(3, 6.939), (8, 15.999), (29, 63.54), (79, 196.97)])
def test_heavy_ion_stopping_matches_oracle(stopping, oracle, z1, m1):
    """Exercises the Brandt-Kitagawa effective-charge branch."""
    energies = np.logspace(np.log10(100.0), np.log10(20_000.0), 30)
    _compare(stopping, oracle, z1, m1, 14, energies)


def test_heavy_ion_low_energy_branch(stopping, oracle):
    """Below yrmin histop switches to velocity-proportional stopping."""
    energies = np.logspace(np.log10(1.0), np.log10(200.0), 25)
    _compare(stopping, oracle, 29, 63.54, 14, energies)


def test_carbon_special_case(stopping, oracle):
    """Z2 == 6 takes a dedicated low-energy correction in histop."""
    energies = np.logspace(np.log10(1.0), np.log10(500.0), 25)
    _compare(stopping, oracle, 29, 63.54, 6, energies)


def test_out_of_range_z_returns_rump_dummies(stopping):
    """RUMP returns finite dummies rather than failing (ziegler.c:199-203)."""
    se, sn = stopping(2, 4.0026, 93, [2000.0])
    assert se[0] == 40.0 and sn[0] == 5.0


def test_reduced_energy_limit(stopping):
    with pytest.raises(ValueError, match="out of range"):
        stopping(1, 1.00797, 14, [200_000.0])


def test_known_literature_values(stopping):
    """Independent sanity check that does not go through the oracle."""
    se, _ = stopping(2, 4.0026, 14, [2000.0])  # 2 MeV He in Si
    assert se[0] == pytest.approx(48.9, rel=0.02)
    se, _ = stopping(1, 1.00797, 14, [2000.0])  # 2 MeV H in Si
    assert se[0] == pytest.approx(5.33, rel=0.02)


def test_units_conversion(stopping, oracle):
    for units in (StoppingUnits.EV_1E15_ATOMS, StoppingUnits.EV_PER_ANGSTROM):
        se, sn = stopping(2, 4.0026, 14, [2000.0], units=units)
        ref_se, ref_sn = oracle.zstop(2, 4.0026, 14, 2000.0, units=int(units))
        assert se[0] == pytest.approx(ref_se, rel=RTOL)
        assert sn[0] == pytest.approx(ref_sn, rel=RTOL)
