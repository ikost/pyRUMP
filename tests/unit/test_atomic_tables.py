"""M1 acceptance: element, isotope and density tables."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from pyrump.atomic.density import STATIC_UNITS, DensityTable
from pyrump.atomic.tables import PeriodicTable
from pyrump.io.scoef import ZIEGLER_MAX_Z, parse_pscoef
from pyrump.model.element import MYLAR_Z


def _reference_data() -> Path | None:
    env = os.environ.get("PYRUMP_C_REFERENCE")
    roots = [Path(env)] if env else []
    roots.append(Path(__file__).resolve().parents[2] / "C-code")
    for root in roots:
        data = root / "rump" / "data"
        if (data / "atom4.dat").is_file():
            return data
    return None


DATA = _reference_data()
pytestmark = pytest.mark.skipif(DATA is None, reason="legacy data tables not available")


@pytest.fixture(scope="module")
def table() -> PeriodicTable:
    assert DATA is not None
    return PeriodicTable.load(DATA / "atom4.dat", DATA / "pscoef.dat")


def test_element_count_and_bounds(table):
    assert len(table) == 93, "93 elements including the Mylar pseudo-element"
    assert table.by_z(1).symbol == "H"
    assert table.by_z(92).symbol == "U"
    assert table.by_z(MYLAR_Z).is_pseudo


@pytest.mark.parametrize(
    "z, symbol, mass, density",
    [
        (1, "H", 1.00797, 4.2716e22),
        (6, "C", 12.011, 1.1364e23),
        (14, "Si", 28.086, 4.9777e22),
        (79, "Au", 196.97, 5.9049e22),
    ],
)
def test_known_elements(table, z, symbol, mass, density):
    element = table.by_z(z)
    assert element.symbol == symbol
    assert element.mass == pytest.approx(mass, rel=1e-5)
    assert element.atomic_density == pytest.approx(density, rel=1e-4)


def test_symbols_are_unique_and_titlecased(table):
    symbols = [e.symbol for e in table]
    assert len(set(symbols)) == len(symbols)
    assert all(s == s.capitalize() for s in symbols)


def test_isotope_abundances_sum_to_one(table):
    """The M1 acceptance criterion: abundances sum to 1 within 1e-4."""
    checked = 0
    for element in table:
        if not element.isotopes:
            continue
        total = sum(i.fraction for i in element.isotopes)
        assert total == pytest.approx(1.0, abs=1e-4), f"Z={element.z} sums to {total}"
        checked += 1
    assert checked > 70, "expected isotope data for most elements"


def test_real_mass_resolution(table):
    silicon = table.by_z(14)
    assert silicon.real_mass(0) == pytest.approx(28.086, rel=1e-5)  # natural
    assert silicon.real_mass(28) == pytest.approx(27.977, rel=1e-4)
    assert silicon.real_mass(29) == pytest.approx(28.976, rel=1e-4)
    # Unknown isotopes fall back to the integer rather than raising (atomdo.c:194).
    assert silicon.real_mass(99) == 99.0


@pytest.mark.parametrize(
    "token, z, mass_number",
    [
        ("Si", 14, 0),
        ("si", 14, 0),
        ("Si+28", 14, 28),
        ("28Si", 14, 28),
        ("Si29", 14, 29),
        ("4He", 2, 4),
        ("Au", 79, 0),
    ],
)
def test_parse_element_reference(table, token, z, mass_number):
    ref = table.parse_ref(token)
    assert (ref.z, ref.mass_number) == (z, mass_number)


def test_parse_element_reference_rejects_junk(table):
    for bad in ("", "Zz", "+28", "123"):
        with pytest.raises((ValueError, KeyError)):
            table.parse_ref(bad)


def test_ziegler_block(table):
    assert DATA is not None
    ziegler = parse_pscoef(DATA / "pscoef.dat")
    assert len(ziegler) == ZIEGLER_MAX_Z
    assert set(ziegler) == set(range(1, ZIEGLER_MAX_Z + 1))
    assert all(len(z.proton_coefficients) == 8 for z in ziegler.values())

    # Mylar is outside Ziegler's range and must carry no coefficients.
    assert table.by_z(MYLAR_Z).ziegler is None

    silicon = ziegler[14]
    assert silicon.most_abundant_mass_number == 28
    assert silicon.density_g_cm3 == pytest.approx(2.321, rel=1e-2)
    assert silicon.atomic_density_e22 == pytest.approx(4.977, rel=1e-2)


def test_density_table():
    assert DATA is not None
    densities = DensityTable.load(DATA / "density.tab")
    assert densities.density("SiO2") == pytest.approx(0.660)
    assert densities.density("sio2") == pytest.approx(0.660), "lookup is case-insensitive"
    assert densities.density("GaAs") == pytest.approx(0.4428)
    # Unknown compounds silently fall back to silicon.
    assert densities.density("Unobtainium") == pytest.approx(0.4997)


def test_thickness_units():
    assert DATA is not None
    densities = DensityTable.load(DATA / "density.tab")
    assert len(STATIC_UNITS) == 5
    assert densities.unit("A").scale == 1.0
    assert densities.unit("nm").scale == 10.0
    assert densities.unit("um").scale == 10000.0
    assert densities.unit("bogus") is None
