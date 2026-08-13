"""M15 acceptance: element-edge physics shared by ELEMENT/MATRIX/INFO/WHATISIT."""

from __future__ import annotations

import pytest

from pyrump.analysis.elements import (
    cosines,
    cross_section_barns,
    kappa,
    locate_candidates,
    matrix_result,
    resolve_element,
)
from pyrump.model.detector import Measurement
from pyrump.model.geometry import Geometry
from pyrump.model.spectrum import Calibration, Spectrum
from pyrump.physics.xsec.rutherford import setup_scatter
from pyrump.shell.session import Buffer
from pyrump.sim.engine import Beam

from conftest import data_dir

DATA = data_dir()
needs_data = pytest.mark.skipif(DATA is None, reason="legacy data tables unavailable")


def test_kappa_returns_zero_rather_than_raising_past_the_kinematic_limit():
    """anlytc.c's RbsKappa returns 0.0 for a forbidden geometry, not an error."""
    # A light projectile can never backscatter off a much lighter target at a
    # near-180-degree angle -- kinematic_factor() itself raises for this case.
    with pytest.raises(ValueError):
        from pyrump.physics.kinematics import kinematic_factor
        kinematic_factor(m1=20.0, m2=1.0, scattering_angle_deg=170.0)
    assert kappa(20.0, 1.0, 170.0) == 0.0


def test_kappa_matches_kinematic_factor_when_allowed():
    from pyrump.physics.kinematics import kinematic_factor

    assert kappa(4.0026, 28.086, 170.0) == pytest.approx(
        kinematic_factor(4.0026, 28.086, 170.0)
    )


def test_cosines_from_geometry():
    geometry = Geometry(theta=0.0, phi=10.0)  # scattering angle 170 deg
    cosin, cosout = cosines(geometry)
    assert cosin == pytest.approx(1.0)  # cos(theta=0)
    # Cornell default: cosout = -cosin*cos_phi (RbsSetexp, anlytc.c:1182).
    assert cosout == pytest.approx(-cosin * geometry.cos_phi)


def test_cross_section_barns_is_unscreened():
    """anlytc.c's local RbsSigma has no L'Ecuyer screening term."""
    unscreened = cross_section_barns(2, 4.0026, 14, 28.086, 170.0, 500.0)
    screened = float(setup_scatter(2, 4.0026, 14, 28.086, 170.0, screening=True)(500.0)[0])
    plain = float(setup_scatter(2, 4.0026, 14, 28.086, 170.0, screening=False)(500.0)[0])
    assert unscreened == pytest.approx(plain)
    assert unscreened != pytest.approx(screened)


def _table_registry():
    from pyrump.cli._common import load_tables

    return load_tables(DATA)[:2]


@needs_data
def test_resolve_element_parses_isotope_syntax():
    table, _ = _table_registry()
    element, mass, mass_number = resolve_element(table, "28Si")
    assert element.symbol == "Si"
    assert mass_number == 28
    assert mass == pytest.approx(27.9769, abs=1e-3)


@needs_data
def test_matrix_result_height_for_a_known_case():
    """Cross-checked by hand: sigma * 1e-24 * 6.25e12 / cosin / epsilon."""
    table, registry = _table_registry()
    buffer = Buffer(
        spectrum=Spectrum.zeros(Calibration(kevch=5.0, npt=1024)),
        beam=Beam(e0_MeV=2.0, z=2, mass=4.0026),
        geometry=Geometry(theta=0.0, phi=10.0),
        measurement=Measurement(),
    )
    result = matrix_result(buffer, table, registry, "Au")
    assert result.k > 0.0
    assert result.height is not None
    assert result.height == pytest.approx(233.6, rel=1e-3)


@needs_data
def test_matrix_result_reports_forbidden_kinematics():
    """A heavy beam can't backscatter off a much lighter target near 180 deg."""
    table, registry = _table_registry()
    buffer = Buffer(
        spectrum=Spectrum.zeros(Calibration(kevch=5.0, npt=1024)),
        beam=Beam(e0_MeV=2.0, z=79, mass=196.97),
        geometry=Geometry(theta=0.0, phi=10.0),  # scattering angle 170 deg
        measurement=Measurement(),
    )
    result = matrix_result(buffer, table, registry, "Si")
    assert result.k == 0.0
    assert result.energy_keV is None
    assert result.height is None


@needs_data
def test_locate_candidates_finds_the_element_and_its_z_neighbors():
    table, registry = _table_registry()
    buffer = Buffer(
        spectrum=Spectrum.zeros(Calibration(kevch=5.0, npt=1024)),
        beam=Beam(e0_MeV=2.0, z=2, mass=4.0026),
        geometry=Geometry(theta=0.0, phi=10.0),
        measurement=Measurement(),
    )
    gold_edge = matrix_result(buffer, table, registry, "Au").energy_keV
    candidates = locate_candidates(buffer, table, gold_edge)
    symbols = [c.symbol for c in candidates]
    assert "Au" in symbols
    assert len(candidates) >= 3
    zs = [c.z for c in candidates]
    assert zs == sorted(zs)
