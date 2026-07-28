"""M5 acceptance: pyRUMP's bricks must match the C's, brick for brick.

The oracle captures ``creatr.c``'s output by redirecting ``SimFillSpectrum``
(see tests/oracle/test_brick_capture.py), so this is a direct comparison of the
whole forward model up to the fill stage: slab discretization, inbound march,
outbound flyout, the [eps] factor, and the yield calculation.

Tolerances are set by float32 in the C, as everywhere else -- the stopping
coefficients alone are stored single-precision.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pytest

from pyrump.atomic.tables import PeriodicTable
from pyrump.io.kalbitzer import parse_kalbitzer
from pyrump.model.geometry import Geometry, GeometryKind
from pyrump.sim.engine import Beam, UniformSample, simulate_bricks
from pyrump.stopping.kalbitzer import KalbitzerStopping
from pyrump.stopping.registry import StoppingRegistry
from pyrump.stopping.ziegler import ZieglerStopping

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
    DATA is None or not ora.available(), reason="legacy tables or oracle unavailable"
)

ENERGY_RTOL = 1e-5
HEIGHT_RTOL = 5e-5

# Oracle brick columns: z mass efront eback hfront hback qqq sigf sigb
O_EFRONT, O_EBACK, O_HFRONT, O_HBACK, O_SIGF, O_SIGB = 2, 3, 4, 5, 7, 8


@pytest.fixture(scope="module")
def table() -> PeriodicTable:
    assert DATA is not None
    return PeriodicTable.load(DATA / "atom4.dat", DATA / "pscoef.dat")


@pytest.fixture(scope="module")
def registry(table) -> StoppingRegistry:
    assert DATA is not None
    return StoppingRegistry(
        table.elements,
        kalbitzer=KalbitzerStopping(parse_kalbitzer(DATA / "newstop.kal"), table.elements),
        ziegler=ZieglerStopping(table.elements),
    )


@pytest.fixture
def oracle() -> ora.Oracle:
    handle = ora.Oracle.load()
    handle.reset_stopping_tables()
    return handle


def _compare(oracle, registry, table, sample, beam, geometry):
    """Run both engines on the same problem and return (mine, theirs)."""
    mine = simulate_bricks(sample, beam, geometry, registry, table)

    oracle.set_beam(
        e0_MeV=beam.e0_MeV,
        zbeam=beam.z,
        mbeam=beam.mass,
        phi=geometry.phi,
        theta=geometry.theta,
        psi=geometry.psi,
        geom=int(geometry.kind),
        kevch=5.0,
        npt=1024,
    )
    oracle.set_sample(
        sample.thicknesses,
        sample.element_z,
        sample.compositions,
        sublayers=sample.sublayers,
        maxpth=sample.maxpth,
        straggle=sample.straggle,
    )
    theirs = oracle.simulate_bricks()
    return mine, theirs


def _assert_bricks_match(mine, theirs, *, straggling: bool = False):
    assert len(mine) == len(theirs), (
        f"brick count differs: {len(mine)} vs {len(theirs)}"
    )
    assert np.allclose(mine.e_front, theirs[:, O_EFRONT], rtol=ENERGY_RTOL)
    assert np.allclose(mine.e_back, theirs[:, O_EBACK], rtol=ENERGY_RTOL)
    assert np.allclose(mine.h_front, theirs[:, O_HFRONT], rtol=HEIGHT_RTOL)
    assert np.allclose(mine.h_back, theirs[:, O_HBACK], rtol=HEIGHT_RTOL)
    if straggling:
        assert np.allclose(mine.sig_front, theirs[:, O_SIGF], rtol=1e-4)
        assert np.allclose(mine.sig_back, theirs[:, O_SIGB], rtol=1e-4)


def test_bare_silicon(oracle, registry, table):
    sample = UniformSample([1000.0], [14], [[1.0]])
    mine, theirs = _compare(
        oracle, registry, table, sample, Beam(), Geometry(theta=0.0, phi=10.0)
    )
    _assert_bricks_match(mine, theirs)


@pytest.mark.parametrize("thickness", [200.0, 1000.0, 5000.0])
def test_thickness_sweep(oracle, registry, table, thickness):
    sample = UniformSample([thickness], [14], [[1.0]])
    mine, theirs = _compare(
        oracle, registry, table, sample, Beam(), Geometry(theta=0.0, phi=10.0)
    )
    _assert_bricks_match(mine, theirs)


@pytest.mark.parametrize("e0", [1.5, 2.0, 3.0])
def test_beam_energy_sweep(oracle, registry, table, e0):
    sample = UniformSample([1000.0], [14], [[1.0]])
    mine, theirs = _compare(
        oracle, registry, table, sample, Beam(e0_MeV=e0), Geometry(theta=0.0, phi=10.0)
    )
    _assert_bricks_match(mine, theirs)


@pytest.mark.parametrize("z_target", [6, 8, 22, 47, 79])
def test_element_sweep(oracle, registry, table, z_target):
    """Different Z exercises different cross-sections and isotope counts."""
    sample = UniformSample([1000.0], [z_target], [[1.0]])
    mine, theirs = _compare(
        oracle, registry, table, sample, Beam(), Geometry(theta=0.0, phi=10.0)
    )
    _assert_bricks_match(mine, theirs)


@pytest.mark.parametrize("phi", [5.0, 10.0, 20.0, 30.0])
def test_detector_angle_sweep(oracle, registry, table, phi):
    sample = UniformSample([1000.0], [14], [[1.0]])
    mine, theirs = _compare(
        oracle, registry, table, sample, Beam(), Geometry(theta=0.0, phi=phi)
    )
    _assert_bricks_match(mine, theirs)


@pytest.mark.parametrize("theta", [0.0, 30.0, 60.0])
def test_tilt_sweep_cornell(oracle, registry, table, theta):
    """Tilt changes both path secants and therefore the slab count."""
    sample = UniformSample([1000.0], [14], [[1.0]])
    mine, theirs = _compare(
        oracle,
        registry,
        table,
        sample,
        Beam(),
        Geometry(theta=theta, phi=10.0, kind=GeometryKind.CORNELL),
    )
    _assert_bricks_match(mine, theirs)


def test_ibm_geometry(oracle, registry, table):
    sample = UniformSample([1000.0], [14], [[1.0]])
    mine, theirs = _compare(
        oracle,
        registry,
        table,
        sample,
        Beam(),
        Geometry(theta=45.0, phi=10.0, kind=GeometryKind.IBM),
    )
    _assert_bricks_match(mine, theirs)


def test_compound_layer(oracle, registry, table):
    """SiO2: Bragg mixing plus two elements' isotopes."""
    sample = UniformSample([1000.0], [14, 8], [[1.0, 2.0]])
    mine, theirs = _compare(
        oracle, registry, table, sample, Beam(), Geometry(theta=0.0, phi=10.0)
    )
    _assert_bricks_match(mine, theirs)


def test_multilayer(oracle, registry, table):
    """A buried gold marker under silicon -- the classic depth-scale check."""
    sample = UniformSample(
        [500.0, 200.0, 1000.0],
        [14, 79],
        [[1.0, 0.0], [0.0, 1.0], [1.0, 0.0]],
        sublayers=[3, 1, 4],
    )
    mine, theirs = _compare(
        oracle, registry, table, sample, Beam(), Geometry(theta=0.0, phi=10.0)
    )
    _assert_bricks_match(mine, theirs)


def test_proton_beam(oracle, registry, table):
    sample = UniformSample([2000.0], [14], [[1.0]])
    mine, theirs = _compare(
        oracle,
        registry,
        table,
        sample,
        Beam(e0_MeV=2.0, z=1, mass=1.00797),
        Geometry(theta=0.0, phi=10.0),
    )
    _assert_bricks_match(mine, theirs)


def test_straggling_enabled(oracle, registry, table):
    sample = UniformSample([1000.0], [14], [[1.0]], straggle=1.0)
    mine, theirs = _compare(
        oracle, registry, table, sample, Beam(), Geometry(theta=0.0, phi=10.0)
    )
    _assert_bricks_match(mine, theirs, straggling=True)
    assert mine.has_straggling


def test_explicit_sublayers(oracle, registry, table):
    sample = UniformSample([1000.0], [14], [[1.0]], sublayers=[17])
    mine, theirs = _compare(
        oracle, registry, table, sample, Beam(), Geometry(theta=0.0, phi=10.0)
    )
    _assert_bricks_match(mine, theirs)
    assert len(mine) == 17 * 3  # three silicon isotopes


def test_bricks_tile_without_gaps(oracle, registry, table):
    """Within an isotope block, each back edge is the next front edge."""
    sample = UniformSample([1000.0], [14], [[1.0]])
    mine = simulate_bricks(
        sample, Beam(), Geometry(theta=0.0, phi=10.0), registry, table
    )
    block = mine[:6]
    assert np.allclose(block.e_back[:-1], block.e_front[1:], rtol=1e-12)
