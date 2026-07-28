"""M11b acceptance: non-Rutherford cross-section tables.

Every shipped ``.adt`` file must parse, across all three dialects, with the
right nuclide identities, angle and units.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from pyrump.io.adt import (
    AdtError,
    CrossSectionTable,
    Dialect,
    SigmaMode,
    parse_reaction,
    read_adt,
    read_adt_single,
)


def _data_dir() -> Path | None:
    env = os.environ.get("PYRUMP_C_REFERENCE")
    roots = [Path(env)] if env else []
    roots.append(Path(__file__).resolve().parents[2] / "C-code")
    for root in roots:
        if (root / "rump" / "data" / "oxy.adt").is_file():
            return root / "rump" / "data"
    return None


DATA = _data_dir()
pytestmark = pytest.mark.skipif(DATA is None, reason="legacy data tables unavailable")

#: (file, dialect, z1, m1, z2, m2, scattering angle, mode)
EXPECTED = [
    ("boron.adt", Dialect.DSIR_33A, 2, 4, 5, 11, 170.5, SigmaMode.BARNS),
    ("car.adt", Dialect.EARLY_RUMP, 2, 4, 6, 12, 170.0, SigmaMode.BARNS),
    ("car_pp.adt", Dialect.DSIR_33A, 1, 1, 6, 12, 160.0, SigmaMode.RELATIVE),
    ("carbon.adt", Dialect.DSIR_33A, 2, 4, 6, 12, 170.5, SigmaMode.BARNS),
    ("fluorine.adt", Dialect.DSIR_33A, 2, 4, 9, 19, 170.0, SigmaMode.BARNS),
    ("nitrogen.adt", Dialect.DSIR_33A, 2, 4, 7, 14, 165.0, SigmaMode.BARNS),
    ("oxy.adt", Dialect.EARLY_RUMP, 2, 4, 8, 16, 170.0, SigmaMode.BARNS),
    ("oxy2.adt", Dialect.DSIR_33A, 2, 4, 8, 16, 170.0, SigmaMode.RELATIVE),
    ("oxy_pp.adt", Dialect.DSIR_33A, 1, 1, 8, 16, 160.0, SigmaMode.RELATIVE),
    ("oxygen.adt", Dialect.DSIR_33A, 2, 4, 8, 16, 170.0, SigmaMode.BARNS),
    ("si.adt", Dialect.EARLY_RUMP, 2, 4, 14, 28, 170.0, SigmaMode.BARNS),
    ("silicon.adt", Dialect.DSIR_33A, 2, 4, 14, 28, 165.1, SigmaMode.BARNS),
]


@pytest.mark.parametrize(
    "name, dialect, z1, m1, z2, m2, angle, mode",
    EXPECTED,
    ids=[e[0] for e in EXPECTED],
)
def test_shipped_tables_parse(name, dialect, z1, m1, z2, m2, angle, mode):
    assert DATA is not None
    table = read_adt_single(DATA / name)

    assert table.dialect is dialect
    assert (table.z1, table.m1, table.z2, table.m2) == (z1, m1, z2, m2)
    assert table.scattering_angle == pytest.approx(angle, abs=0.05)
    assert table.mode is mode
    assert table.npt > 10


@pytest.mark.parametrize("name", [e[0] for e in EXPECTED])
def test_tables_are_strictly_sorted(name):
    """Interpolation slopes assume strict ordering; RUMP drops violators."""
    assert DATA is not None
    table = read_adt_single(DATA / name)
    assert np.all(np.diff(table.energy_keV) > 0)


@pytest.mark.parametrize("name", [e[0] for e in EXPECTED])
def test_cross_sections_are_positive(name):
    assert DATA is not None
    table = read_adt_single(DATA / name)
    assert np.all(table.sigma > 0)


def test_nitrogen_is_not_a_neutron():
    """Regression: 'N' and 'n' differ only in case.

    A case-insensitive element lookup with neutron at index 0 resolves nitrogen
    to Z=0 -- and the file parses cleanly, so nothing looks wrong.
    """
    assert DATA is not None
    table = read_adt_single(DATA / "nitrogen.adt")
    assert table.z2 == 7


def test_millibarns_are_converted():
    """mb/sr files must come out in barns/sr."""
    assert DATA is not None
    table = read_adt_single(DATA / "oxygen.adt")
    assert table.mode is SigmaMode.BARNS
    # Alpha-oxygen resonances are of order 0.01-1 b/sr, not 10-1000.
    assert table.sigma.max() < 10.0


def test_relative_mode_is_flagged_not_scaled():
    """'rtr' tables are ratios to Rutherford and must not be treated as barns."""
    assert DATA is not None
    table = read_adt_single(DATA / "oxy2.adt")
    assert table.mode is SigmaMode.RELATIVE
    # A ratio hovers around 1, so this is a cheap sanity check on the reading.
    assert 0.1 < np.median(table.sigma) < 20.0


def test_interpolation_and_coverage():
    assert DATA is not None
    table = read_adt_single(DATA / "oxygen.adt")

    # Exact at the knots -- the table is piecewise linear through them.
    assert table(table.energy_keV[5])[0] == pytest.approx(table.sigma[5])
    # Midpoint lies between neighbours.
    midpoint = 0.5 * (table.energy_keV[5] + table.energy_keV[6])
    low, high = sorted((table.sigma[5], table.sigma[6]))
    assert low <= table(midpoint)[0] <= high

    inside = table.covers([table.energy_keV[0], table.energy_keV[-1]])
    assert bool(inside[0]) and bool(inside[1])
    assert not bool(table.covers([table.energy_keV[0] - 1.0])[0])
    assert not bool(table.covers([table.energy_keV[-1] + 1.0])[0])


def test_slopes_match_finite_differences():
    assert DATA is not None
    table = read_adt_single(DATA / "carbon.adt")
    slopes = table.slopes()
    expected = np.diff(table.sigma) / np.diff(table.energy_keV)
    assert np.allclose(slopes[:-1], expected)
    assert slopes[-1] == 0.0  # flagged end-of-data, as RUMP stores it


# --------------------------------------------------------------- reactions


@pytest.mark.parametrize(
    "text, expected",
    [
        ("11B(a,a)11B", (2, 4, 5, 11)),
        ("12C(p,p)12C", (1, 1, 6, 12)),
        ("12C(p,p0)12C", (1, 1, 6, 12)),  # R33 state index
        ("16O(a,a)16O", (2, 4, 8, 16)),
        ("28Si(d,d)28Si", (1, 2, 14, 28)),
    ],
)
def test_parse_reaction(text, expected):
    assert parse_reaction(text) == expected


def test_excited_states_are_rejected():
    """p1 leaves the nucleus excited -- inelastic, so unusable."""
    with pytest.raises(AdtError, match="excited"):
        parse_reaction("12C(p,p1)12C")


def test_non_elastic_channels_are_rejected():
    with pytest.raises(AdtError, match="not elastic"):
        parse_reaction("14N(d,p)15N")
    for text in ("10B(n,a)7Li", "12C(p,g)13N"):
        with pytest.raises(AdtError):
            parse_reaction(text)


def test_unparseable_reaction():
    with pytest.raises(AdtError, match="cannot parse"):
        parse_reaction("nonsense")


# ------------------------------------------------------------- edge cases


def test_r33_format_exemplar_is_rejected_as_rump_rejects_it(tmp_path):
    """The shipped R33.Format uses 'Units: mb', which RUMP does not accept.

    Its accepted set is exactly {b/sr, mb/sr, rtr, rr, relative}
    (reswork.c:325-334), so RUMP cannot read its own bundled R33 example --
    it is a format exemplar from SigmaCalc, not a loadable table. Matching that
    refusal is fidelity, not a gap.
    """
    assert DATA is not None
    with pytest.raises(AdtError, match="invalid UNITS"):
        read_adt(DATA / "R33.Format")


def test_non_zero_q_is_rejected(tmp_path):
    path = tmp_path / "q.adt"
    path.write_text(
        "Version: R33\nReaction: 14N(a,a)14N\nQvalue: 2.5\nTheta: 170\n"
        "Units: mb/sr\nData:\n1000 0 1.0 0\n2000 0 2.0 0\n"
    )
    with pytest.raises(AdtError, match="non-zero Q"):
        read_adt(path)


def test_comma_separated_qvalue_list_is_accepted(tmp_path):
    """Several shipped files write 'Qvalue: 0.00, 0.00, ...'; atof reads one."""
    path = tmp_path / "q.adt"
    path.write_text(
        "Version: DSIR 33a\nReaction: 12C(p,p)12C\n"
        "Qvalue:      0.00,      0.00,      0.00\nTheta: 160\n"
        "Units: mb/sr\nData:\n1000 1.0\n2000 2.0\n"
    )
    table = read_adt_single(path)
    assert table.npt == 2


def test_multiple_angle_blocks(tmp_path):
    """EndData: recycles into header mode for another angle."""
    path = tmp_path / "multi.adt"
    path.write_text(
        "Version: DSIR 33a\nReaction: 16O(a,a)16O\nUnits: mb/sr\n"
        "Theta: 170\nData:\n1000 1.0\n2000 2.0\nEndData:\n"
        "Theta: 150\nData:\n1000 3.0\n2000 4.0\nEndData:\n"
    )
    tables = read_adt(path)
    assert len(tables) == 2
    assert tables[0].scattering_angle == pytest.approx(170.0)
    assert tables[1].scattering_angle == pytest.approx(150.0)
    with pytest.raises(AdtError, match="angle blocks"):
        read_adt_single(path)


def test_energy_and_sigma_factors(tmp_path):
    path = tmp_path / "scaled.adt"
    path.write_text(
        "Version: DSIR 33a\nReaction: 16O(a,a)16O\nUnits: b/sr\n"
        "Theta: 170\nEnfactors: 1000.0 5.0\nSigfactors: 2.0\n"
        "Data:\n1.0 1.0\n2.0 3.0\n"
    )
    table = read_adt_single(path)
    assert table.energy_keV == pytest.approx([1005.0, 2005.0])
    assert table.sigma == pytest.approx([2.0, 6.0])


def test_out_of_order_points_are_dropped(tmp_path):
    """RUMP deletes them "to maintain strict data sort" (reswork.c:428)."""
    path = tmp_path / "unsorted.adt"
    path.write_text(
        "Version: DSIR 33a\nReaction: 16O(a,a)16O\nUnits: b/sr\nTheta: 170\n"
        "Data:\n1000 1.0\n2000 2.0\n1500 9.0\n3000 3.0\n"
    )
    table = read_adt_single(path)
    assert table.energy_keV == pytest.approx([1000.0, 2000.0, 3000.0])
    assert 9.0 not in table.sigma


def test_early_rump_header_is_positional(tmp_path):
    path = tmp_path / "old.adt"
    path.write_text("# comment\n2 4.0 8 16.0 10.0 3\n1000 0.5\n2000 0.6\n3000 0.7\n")
    table = read_adt_single(path)
    assert table.dialect is Dialect.EARLY_RUMP
    assert (table.z1, table.m1, table.z2, table.m2) == (2, 4, 8, 16)
    assert table.phi == pytest.approx(10.0)
    assert table.mode is SigmaMode.BARNS  # always barns in this dialect
    assert table.npt == 3


def test_malformed_early_rump_header(tmp_path):
    path = tmp_path / "bad.adt"
    path.write_text("2 4.0 8\n1000 0.5\n")
    with pytest.raises(AdtError, match="Z1 M1 Z2 M2 Phi NPT"):
        read_adt(path)


def test_empty_file(tmp_path):
    path = tmp_path / "empty.adt"
    path.write_text("# only a comment\n")
    with pytest.raises(AdtError, match="no cross-section data"):
        read_adt(path)


def test_table_is_usable_as_a_callable():
    assert DATA is not None
    table: CrossSectionTable = read_adt_single(DATA / "si.adt")
    energies = np.linspace(table.energy_keV[0], table.energy_keV[-1], 50)
    values = table(energies)
    assert values.shape == energies.shape
    assert np.all(values > 0)
