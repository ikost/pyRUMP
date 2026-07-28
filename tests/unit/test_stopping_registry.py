"""M2 acceptance: Kalbitzer, Mylar and the source priority chain."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from pyrump.atomic.tables import PeriodicTable
from pyrump.io.kalbitzer import parse_kalbitzer
from pyrump.stopping.kalbitzer import KalbitzerStopping
from pyrump.stopping.mylar import MYLAR_Z, mylar_stopping
from pyrump.stopping.registry import StoppingRegistry, StoppingSource
from pyrump.stopping.ziegler import ZieglerStopping


def _data_dir() -> Path | None:
    env = os.environ.get("PYRUMP_C_REFERENCE")
    roots = [Path(env)] if env else []
    roots.append(Path(__file__).resolve().parents[2] / "C-code")
    for root in roots:
        if (root / "rump" / "data" / "newstop.kal").is_file():
            return root / "rump" / "data"
    return None


DATA = _data_dir()
pytestmark = pytest.mark.skipif(DATA is None, reason="legacy data tables unavailable")


@pytest.fixture(scope="module")
def elements():
    assert DATA is not None
    return PeriodicTable.load(DATA / "atom4.dat", DATA / "pscoef.dat").elements


@pytest.fixture(scope="module")
def registry(elements) -> StoppingRegistry:
    assert DATA is not None
    entries = parse_kalbitzer(DATA / "newstop.kal")
    return StoppingRegistry(
        elements,
        kalbitzer=KalbitzerStopping(entries, elements),
        ziegler=ZieglerStopping(elements),
    )


def test_kalbitzer_table_parses():
    assert DATA is not None
    entries = parse_kalbitzer(DATA / "newstop.kal")
    assert len(entries) == 8, "shipped file holds 8 Konac records"
    assert all(len(e.a) == 6 for e in entries)
    # H, D, 3He, 4He on carbon and silicon.
    assert {e.z2 for e in entries} == {6, 14}
    assert {e.z1 for e in entries} == {1, 2}


@pytest.mark.parametrize(
    "z1, m1, z2, expected",
    [
        # Kalbitzer wins for the common RBS cases -- Ziegler is NOT used here.
        (2, 4.0026, 14, StoppingSource.KALBITZER),
        (2, 4.0026, 6, StoppingSource.KALBITZER),
        (1, 1.00797, 14, StoppingSource.KALBITZER),
        # No Konac entry -> Ziegler.
        (2, 4.0026, 79, StoppingSource.ZIEGLER),
        (1, 1.00797, 8, StoppingSource.ZIEGLER),
        (29, 63.54, 14, StoppingSource.ZIEGLER),
        # Mylar pseudo-element is outside Ziegler's range.
        (2, 4.0026, MYLAR_Z, StoppingSource.MYLAR),
        (1, 1.00797, MYLAR_Z, StoppingSource.MYLAR),
    ],
)
def test_priority_chain(registry, z1, m1, z2, expected):
    assert registry.source_for(z1, m1, z2) is expected
    result = registry(z1, m1, z2, [2000.0])
    assert result.source is expected
    assert result.values[0] > 0


def test_kalbitzer_differs_from_ziegler(registry, elements):
    """If these agreed, priority order would not matter -- it does."""
    ziegler = ZieglerStopping(elements)
    energies = np.array([500.0, 1000.0, 2000.0])
    chain = registry(2, 4.0026, 14, energies)
    se, sn = ziegler(2, 4.0026, 14, energies)

    assert chain.source is StoppingSource.KALBITZER
    relative = np.abs(chain.values - (se + sn)) / (se + sn)
    assert np.max(relative) > 0.01, "Konac and ZBL should differ by more than rounding"
    # ...but not wildly: both describe the same physics.
    assert np.max(relative) < 0.25


def test_kalbitzer_is_physically_sensible(registry):
    """He in Si peaks near 500 keV and falls at higher energy."""
    energies = np.array([250.0, 500.0, 1000.0, 2000.0, 3000.0])
    values = registry(2, 4.0026, 14, energies).values
    assert values.argmax() <= 1
    assert values[-1] < values[0]
    assert values[2] == pytest.approx(60.0, rel=0.15)


def test_mylar_polynomials():
    for z1, m1 in ((1, 1.0), (1, 2.0), (2, 4.0)):
        values = mylar_stopping(z1, m1, [500.0, 1000.0, 2000.0])
        assert values is not None
        assert np.all(values > 0)

    # Unsupported projectiles fall through rather than returning nonsense.
    assert mylar_stopping(2, 3.0, [1000.0]) is None
    assert mylar_stopping(6, 12.0, [1000.0]) is None


def test_mylar_matches_hardcoded_polynomial():
    """Spot-check against a direct evaluation of the stopping.c literals."""
    coefficients = [12.71, +5.538e-2, -5.557e-5, +2.440e-8, -5.175e-12, +4.274e-16]
    energy = 1500.0
    expected = sum(c * energy**i for i, c in enumerate(coefficients))
    assert mylar_stopping(2, 4.0, [energy])[0] == pytest.approx(expected, rel=1e-12)


def test_unsupported_combination_raises(registry):
    with pytest.raises(ValueError, match="no stopping power available"):
        registry(6, 12.0, MYLAR_Z, [1000.0])
