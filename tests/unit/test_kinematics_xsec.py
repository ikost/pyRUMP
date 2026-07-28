"""M4 acceptance: geometry, kinematics and cross-sections.

These are closed forms with no fitting or table lookup, so the tolerance is tight
(1e-10). Anything looser would be hiding a real disagreement.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

from pyrump.model.geometry import Geometry, GeometryKind
from pyrump.physics.kinematics import (
    edge_energy,
    kinematic_factor,
    kinematic_factors,
    recoil_factor,
)
from pyrump.physics.xsec.rutherford import (
    E2_OVER_2_SQUARED,
    E2_OVER_4_SQUARED,
    CrossSectionKind,
    setup_recoil,
    setup_scatter,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "oracle"))
import oracle as ora  # noqa: E402

RTOL = 1e-10

oracle_only = pytest.mark.skipif(
    not ora.available() or ora.data_dir() is None, reason="oracle unavailable"
)


@pytest.fixture(scope="module")
def oracle() -> ora.Oracle:
    return ora.Oracle.load()


# ------------------------------------------------------------------- geometry


def test_phi_is_the_supplement_of_the_scattering_angle():
    """A detector at 170 deg is entered as phi=10 -- RUMP's biggest gotcha."""
    geometry = Geometry(theta=0.0, phi=10.0)
    assert geometry.scattering_angle == 170.0
    assert geometry.cos_phi == pytest.approx(math.cos(math.radians(170.0)), rel=RTOL)


def test_cornell_geometry():
    geometry = Geometry(theta=0.0, phi=10.0, kind=GeometryKind.CORNELL)
    assert geometry.sec_in == pytest.approx(1.0, rel=RTOL)
    # sec_out = -sec_in/cos(true angle); cos(170 deg) is negative, so this is positive.
    expected = -1.0 / math.cos(math.radians(170.0))
    assert geometry.sec_out == pytest.approx(expected, rel=RTOL)
    geometry.validate()


def test_ibm_geometry():
    geometry = Geometry(theta=7.0, phi=10.0, kind=GeometryKind.IBM)
    assert geometry.sec_in == pytest.approx(1.0 / math.cos(math.radians(7.0)), rel=RTOL)
    assert geometry.sec_out == pytest.approx(
        1.0 / math.cos(math.radians(17.0)), rel=RTOL
    )


def test_general_geometry_uses_psi():
    geometry = Geometry(theta=0.0, phi=20.0, psi=20.0, kind=GeometryKind.GENERAL)
    assert geometry.sec_out == pytest.approx(1.0 / math.cos(math.radians(20.0)), rel=RTOL)


def test_geometry_rejects_grazing_exit():
    """Non-positive secants mean the path runs into the surface (creatr.c:423)."""
    with pytest.raises(ValueError, match="bad scattering geometry"):
        Geometry(theta=95.0, phi=10.0).validate()


def test_tilt_increases_path_length():
    flat = Geometry(theta=0.0, phi=10.0, kind=GeometryKind.IBM)
    tilted = Geometry(theta=60.0, phi=10.0, kind=GeometryKind.IBM)
    assert tilted.sec_in > flat.sec_in
    assert tilted.sec_in == pytest.approx(2.0, rel=1e-9)  # 1/cos(60)


# ----------------------------------------------------------------- kinematics


@pytest.mark.parametrize(
    "m1, m2, angle, expected",
    [
        # Evaluated independently from Chu eq. 2.6 with the masses in atom4.dat,
        # NOT copied from a published table -- tabulated K values use modern
        # atomic weights and differ in the 4th digit.
        (4.0026, 28.086, 170.0, 0.565737),  # He on Si
        (4.0026, 196.97, 170.0, 0.922491),  # He on Au
        (4.0026, 12.011, 170.0, 0.2526473),  # He on C
        (1.00797, 28.086, 170.0, 0.867165),  # H on Si
    ],
)
def test_kinematic_factor_against_closed_form(m1, m2, angle, expected):
    assert kinematic_factor(m1, m2, angle) == pytest.approx(expected, rel=1e-6)


def test_kinematic_factor_is_within_a_percent_of_published_tables():
    """Looser check against values as usually quoted in the RBS literature."""
    assert kinematic_factor(4.0026, 28.086, 170.0) == pytest.approx(0.566, rel=5e-3)
    assert kinematic_factor(4.0026, 196.97, 170.0) == pytest.approx(0.923, rel=5e-3)


def test_kinematic_factor_limits():
    """K -> 1 forward, and is minimal at exact backscatter."""
    assert kinematic_factor(4.0026, 28.086, 0.0) == pytest.approx(1.0, rel=RTOL)
    back = kinematic_factor(4.0026, 28.086, 180.0)
    near = kinematic_factor(4.0026, 28.086, 175.0)
    assert back < near


def test_kinematic_factor_equal_masses_at_90_degrees():
    """Equal masses transfer everything at 90 deg: K = 0."""
    assert kinematic_factor(4.0, 4.0, 90.0) == pytest.approx(0.0, abs=1e-12)


def test_kinematic_factor_rejects_impossible_angle():
    """A heavy projectile on a light target cannot scatter backwards."""
    with pytest.raises(ValueError, match="kinematic limit"):
        kinematic_factor(28.086, 4.0026, 170.0)


def test_recoil_factor():
    """4 m1 m2 cos^2 / (m1+m2)^2, peaking at zero degrees."""
    m1, m2 = 4.0026, 1.00797
    expected = 4 * m1 * m2 * math.cos(math.radians(30.0)) ** 2 / (m1 + m2) ** 2
    assert recoil_factor(m1, m2, 30.0) == pytest.approx(expected, rel=RTOL)
    assert recoil_factor(m1, m2, 0.0) > recoil_factor(m1, m2, 45.0)


def test_recoil_requires_forward_angle():
    with pytest.raises(ValueError, match="below 90 deg"):
        recoil_factor(4.0026, 1.00797, 120.0)


def test_edge_energy():
    """2 MeV He on Au at 170 deg puts the surface edge near 1846 keV."""
    assert edge_energy(2000.0, 4.0026, 196.97, 170.0) == pytest.approx(1846.0, rel=1e-3)


def test_vectorised_kinematic_factors_match_scalar():
    masses = np.array([12.011, 15.999, 28.086, 196.97])
    got = kinematic_factors(4.0026, masses, 170.0)
    for value, mass in zip(got, masses):
        assert value == pytest.approx(kinematic_factor(4.0026, mass, 170.0), rel=RTOL)


# -------------------------------------------------------------- cross sections


def test_rutherford_scales_as_inverse_energy_squared():
    xsec = setup_scatter(2, 4.0026, 79, 196.97, 170.0, screening=False)
    values = xsec([1000.0, 2000.0])
    assert values[0] / values[1] == pytest.approx(4.0, rel=RTOL)


def test_screening_reduces_the_cross_section():
    screened = setup_scatter(2, 4.0026, 79, 196.97, 170.0)
    bare = setup_scatter(2, 4.0026, 79, 196.97, 170.0, screening=False)
    assert screened.kind is CrossSectionKind.RUTHERFORD_SCREENED
    assert screened(1000.0)[0] < bare(1000.0)[0]
    # L'Ecuyer: sigma *= (1 - 0.049 Z1 Z2^{4/3} / E)
    factor = 1.0 - (0.049 * 2 * 79**1.3333) / 1000.0
    assert screened(1000.0)[0] == pytest.approx(bare(1000.0)[0] * factor, rel=RTOL)


def test_cross_section_constants_are_the_documented_values():
    assert E2_OVER_4_SQUARED == 1295.9358
    assert E2_OVER_2_SQUARED == 5183.7432


def test_recoil_selects_analytic_forms_in_their_windows():
    # 4He on 1H below 40 deg -> Ziegler's measured form.
    assert setup_recoil(2, 4.0026, 1, 1.00797, 30.0).kind is CrossSectionKind.ZIEGLER_H
    # ...but Rutherford outside it.
    assert setup_recoil(2, 4.0026, 1, 1.00797, 50.0).kind is CrossSectionKind.RUTHERFORD

    # 4He on 2H between 10 and 32 deg -> Quillet.
    assert setup_recoil(2, 4.0026, 1, 2.014, 20.0).kind is CrossSectionKind.QUILLET_D
    assert setup_recoil(2, 4.0026, 1, 2.014, 35.0).kind is CrossSectionKind.RUTHERFORD

    # Other targets are always plain Rutherford.
    assert setup_recoil(2, 4.0026, 6, 12.011, 20.0).kind is CrossSectionKind.RUTHERFORD


def test_quillet_is_constant_above_its_validity_limit():
    """Above 2.7 MeV the cross-section is frozen, not extrapolated (sigma.c:361)."""
    xsec = setup_recoil(2, 4.0026, 1, 2.014, 20.0)
    assert xsec(2700.0)[0] == pytest.approx(xsec(4000.0)[0], rel=RTOL)


def test_ziegler_h_ratio_freezes_but_energy_dependence_continues():
    """Above 4 MeV only the *ratio* is held; the result still falls as 1/E^2."""
    xsec = setup_recoil(2, 4.0026, 1, 1.00797, 30.0)
    high, higher = xsec(5000.0)[0], xsec(10000.0)[0]
    assert high / higher == pytest.approx(4.0, rel=1e-9)


# ------------------------------------------------------- comparison with the C


@oracle_only
@pytest.mark.parametrize(
    "z1, m1, z2, m2, angle",
    [
        (2, 4.0026, 14, 28.086, 170.0),
        (2, 4.0026, 79, 196.97, 170.0),
        (2, 4.0026, 8, 15.999, 160.0),
        (1, 1.00797, 14, 28.086, 150.0),
        (2, 4.0026, 6, 12.011, 175.0),
    ],
)
def test_scatter_matches_oracle(oracle, z1, m1, z2, m2, angle):
    energies = np.linspace(500.0, 3000.0, 40)
    ours = setup_scatter(z1, m1, z2, m2, angle)(energies)
    theirs = oracle.cross_section(z1, m1, z2, m2, angle, energies)
    assert np.allclose(ours, theirs, rtol=RTOL)


@oracle_only
@pytest.mark.parametrize(
    "z1, m1, z2, m2, angle",
    [
        (2, 4.0026, 1, 1.00797, 30.0),  # Ziegler H
        (2, 4.0026, 1, 2.014, 20.0),  # Quillet D
        (2, 4.0026, 1, 1.00797, 50.0),  # falls back to Rutherford
        (2, 4.0026, 6, 12.011, 25.0),  # plain Rutherford recoil
    ],
)
def test_recoil_matches_oracle(oracle, z1, m1, z2, m2, angle):
    energies = np.linspace(500.0, 2500.0, 40)
    ours = setup_recoil(z1, m1, z2, m2, angle)(energies)
    theirs = oracle.cross_section(z1, m1, z2, m2, angle, energies, recoil=True)
    assert np.allclose(ours, theirs, rtol=RTOL)


@oracle_only
def test_setup_constants_match_oracle(oracle):
    ours = setup_scatter(2, 4.0026, 79, 196.97, 170.0)
    theirs = oracle.sigma_constants(2, 4.0026, 79, 196.97, 170.0)
    assert ours.csigma == pytest.approx(theirs["csigma"], rel=RTOL)
    assert ours.csig_0 == pytest.approx(theirs["csig_0"], abs=1e-12)
    assert ours.csig_f == pytest.approx(theirs["csig_f"], rel=RTOL)
