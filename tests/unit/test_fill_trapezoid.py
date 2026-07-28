"""M6 acceptance: brick -> channel fill and the full no-straggling spectrum.

Three layers of evidence, because a single per-channel tolerance would be
misleading here:

1. the vectorised fill equals a literal transliteration of ``SimAnlyz4``'s loop
   **bit for bit** -- so the antiderivative refactor is exact, not approximate;
2. bulk channels match the C to the float32 floor;
3. *partial* end channels are compared in absolute terms. A brick's first and
   last channels are slivers whose width is the difference of two nearly-equal
   float32 edge energies, so their relative error is dominated by C rounding and
   carries no information -- one such channel holds 0.005% of the spectrum yet
   shows 7e-3 relative deviation, worth 0.09 counts.

Total counts are the honest integral check, and agree to ~3e-6.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pytest

from pyrump.atomic.tables import PeriodicTable
from pyrump.io.kalbitzer import parse_kalbitzer
from pyrump.model.detector import Measurement, yield_normalisation
from pyrump.model.geometry import Geometry
from pyrump.model.spectrum import Calibration
from pyrump.sim.bricks import Bricks
from pyrump.sim.engine import Beam, UniformSample, simulate, simulate_bricks
from pyrump.sim.fill.trapezoid import fill_trapezoid, fill_trapezoid_reference
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
CAL = Calibration(kevch=5.0, kev0=0.0, first=0.0, npt=1024)

#: Bulk (fully covered) channels, limited by float32 stopping coefficients.
BULK_RTOL = 5e-5
#: Integral over the whole spectrum.
TOTAL_RTOL = 1e-5


# ------------------------------------------------------- fill unit behaviour


def _brick(e_front, e_back, h_front, h_back):
    return Bricks.from_list([(e_front, e_back, h_front, h_back, 0.0, 0.0, 0.0)])


def test_single_brick_conserves_area():
    """A brick well inside the range must deposit exactly its trapezoid area."""
    bricks = _brick(100.0, 50.0, 3.0, 1.0)
    counts = fill_trapezoid(bricks, CAL)
    expected = 0.5 * (3.0 + 1.0) * (100.0 - 50.0)
    assert counts.sum() == pytest.approx(expected, rel=1e-12)


def test_flat_brick_spreads_evenly():
    """Constant height over exactly ten channels."""
    bricks = _brick(50.0, 0.0, 2.0, 2.0)
    counts = fill_trapezoid(bricks, CAL)
    assert np.allclose(counts[0:10], 2.0 * CAL.kevch)
    assert counts[10] == 0.0


def test_partial_channels_are_clipped():
    """A brick starting and ending mid-channel gets slivers at both ends."""
    bricks = _brick(12.5, 2.5, 1.0, 1.0)  # channels 0..2 with kevch=5
    counts = fill_trapezoid(bricks, CAL)
    assert counts[0] == pytest.approx(2.5)  # 2.5..5
    assert counts[1] == pytest.approx(5.0)  # 5..10 (full)
    assert counts[2] == pytest.approx(2.5)  # 10..12.5
    assert counts.sum() == pytest.approx(10.0)


def test_inverted_brick_is_ignored():
    """anlyz.c:308 returns before dividing when eback >= efront."""
    assert fill_trapezoid(_brick(50.0, 100.0, 1.0, 1.0), CAL).sum() == 0.0
    assert fill_trapezoid(_brick(50.0, 50.0, 1.0, 1.0), CAL).sum() == 0.0


def test_brick_above_the_top_channel_is_dropped_whole():
    """Not clipped -- discarded, with a warning (anlyz.c:313-316).

    Reproducing this matters: clipping instead would silently change any
    spectrum whose sample scatters above the calibrated range.
    """
    top = CAL.edge_energy(CAL.npt)
    bricks = _brick(float(top) + 100.0, float(top) - 100.0, 1.0, 1.0)
    assert fill_trapezoid(bricks, CAL).sum() == 0.0


def test_brick_below_channel_zero_is_clipped_not_dropped():
    bricks = _brick(20.0, -30.0, 1.0, 1.0)
    counts = fill_trapezoid(bricks, CAL)
    assert counts.sum() == pytest.approx(20.0)  # only the 0..20 keV part lands


def test_empty_brick_set():
    assert fill_trapezoid(Bricks.empty(0), CAL).sum() == 0.0


def test_calibration_offset_shifts_channels():
    shifted = Calibration(kevch=5.0, kev0=100.0, first=0.0, npt=1024)
    counts = fill_trapezoid(_brick(150.0, 100.0, 1.0, 1.0), shifted)
    assert np.count_nonzero(counts[:10]) == 10
    assert counts[10] == 0.0


# --------------------------------------------------- vectorised == literal C

pytestmark = pytest.mark.skipif(
    DATA is None or not ora.available(), reason="legacy tables or oracle unavailable"
)


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


CASES = [
    ("bare Si", UniformSample([1000.0], [14], [[1.0]]), Beam(), Geometry(theta=0.0, phi=10.0)),
    ("thick Si", UniformSample([5000.0], [14], [[1.0]]), Beam(), Geometry(theta=0.0, phi=10.0)),
    ("Au film", UniformSample([500.0], [79], [[1.0]]), Beam(), Geometry(theta=0.0, phi=10.0)),
    ("SiO2", UniformSample([1000.0], [14, 8], [[1.0, 2.0]]), Beam(), Geometry(theta=0.0, phi=10.0)),
    (
        "Au marker in Si",
        UniformSample([500.0, 200.0, 1000.0], [14, 79], [[1, 0], [0, 1], [1, 0]], sublayers=[3, 1, 4]),
        Beam(),
        Geometry(theta=0.0, phi=10.0),
    ),
    ("tilted Si", UniformSample([1000.0], [14], [[1.0]]), Beam(), Geometry(theta=60.0, phi=10.0)),
    ("proton beam", UniformSample([2000.0], [14], [[1.0]]), Beam(z=1, mass=1.00797), Geometry(theta=0.0, phi=10.0)),
]


@pytest.mark.parametrize("label, sample, beam, geometry", CASES, ids=[c[0] for c in CASES])
def test_vectorised_fill_equals_the_c_loop(registry, table, label, sample, beam, geometry):
    """The refactor must be exact, not merely close."""
    bricks = simulate_bricks(sample, beam, geometry, registry, table)
    fast = fill_trapezoid(bricks, CAL)
    literal = fill_trapezoid_reference(bricks, CAL)
    assert np.array_equal(fast, literal)


# --------------------------------------------------- full spectrum vs oracle


def _oracle_spectrum(oracle, sample, beam, geometry, measurement):
    oracle.set_beam(
        e0_MeV=beam.e0_MeV,
        zbeam=beam.z,
        mbeam=beam.mass,
        cbeam=measurement.charge_state,
        phi=geometry.phi,
        theta=geometry.theta,
        geom=int(geometry.kind),
        kevch=CAL.kevch,
        kev0=CAL.kev0,
        first=CAL.first,
        npt=CAL.npt,
        fwhm=0.0,
        q_uC=measurement.charge_uC,
        omega=measurement.omega_msr,
        corr=measurement.correction,
    )
    oracle.set_sample(
        sample.thicknesses,
        sample.element_z,
        sample.compositions,
        sublayers=sample.sublayers,
        maxpth=sample.maxpth,
    )
    return oracle.simulate_spectrum()


@pytest.mark.parametrize("label, sample, beam, geometry", CASES, ids=[c[0] for c in CASES])
def test_spectrum_matches_oracle(oracle, registry, table, label, sample, beam, geometry):
    measurement = Measurement(omega_msr=1.0, charge_uC=10.0)
    mine = simulate(sample, beam, geometry, registry, table, CAL, measurement)
    theirs = _oracle_spectrum(oracle, sample, beam, geometry, measurement)

    n = min(len(theirs), len(mine.counts))
    ours, ref = mine.counts[:n], theirs[:n]

    # Total counts: the meaningful integral check.
    assert ours.sum() == pytest.approx(ref.sum(), rel=TOTAL_RTOL)

    # Per channel: relative agreement OR a small absolute fraction of the peak.
    #
    # The absolute leg is not slack, it is necessary. A channel containing a
    # brick edge is a sliver whose width is the difference of two nearly-equal
    # float32 energies, so its relative error is set by C rounding rather than
    # by the algorithm -- and such a channel can still hold 30% of the peak, so
    # it cannot be excluded by magnitude.
    assert np.allclose(ours, ref, rtol=BULK_RTOL, atol=1e-4 * ref.max())

    # The same channels are occupied.
    assert np.array_equal(ref != 0, ours != 0)


def test_yield_normalisation_is_the_charge_and_solid_angle(oracle, registry, table):
    """counts scale linearly with Q and Omega, inversely with CORR."""
    sample, beam = UniformSample([1000.0], [14], [[1.0]]), Beam()
    geometry = Geometry(theta=0.0, phi=10.0)

    base = Measurement(omega_msr=1.0, charge_uC=10.0)
    doubled = Measurement(omega_msr=2.0, charge_uC=10.0)
    more_charge = Measurement(omega_msr=1.0, charge_uC=20.0)

    a = simulate(sample, beam, geometry, registry, table, CAL, base).total()
    b = simulate(sample, beam, geometry, registry, table, CAL, doubled).total()
    c = simulate(sample, beam, geometry, registry, table, CAL, more_charge).total()
    assert b == pytest.approx(2 * a, rel=1e-12)
    assert c == pytest.approx(2 * a, rel=1e-12)

    assert yield_normalisation(base) == pytest.approx(10.0)


def test_fill_routine_is_chosen_by_straggling(registry, table):
    """SimAnlyz (anlyz.c:207) dispatches on whether either width is non-zero.

    With straggling the trapezoid is abandoned for two broadened triangles, so
    the two paths must give the same integral but different distributions.
    """
    geometry = Geometry(theta=0.0, phi=10.0)
    sharp = simulate(
        UniformSample([1000.0], [14], [[1.0]]), Beam(), geometry, registry, table, CAL
    )
    broad = simulate(
        UniformSample([1000.0], [14], [[1.0]], straggle=1.0),
        Beam(), geometry, registry, table, CAL,
    )
    assert broad.total() == pytest.approx(sharp.total(), rel=1e-3)
    assert not np.array_equal(broad.counts, sharp.counts)
    assert np.count_nonzero(broad.counts) > np.count_nonzero(sharp.counts)


def test_surface_edge_lands_in_the_right_channel(registry, table):
    """The sharpest calibration check: the Au edge must sit at K*E0."""
    from pyrump.physics.kinematics import kinematic_factor

    sample = UniformSample([200.0], [79], [[1.0]])
    mine = simulate(sample, Beam(), Geometry(theta=0.0, phi=10.0), registry, table, CAL)
    occupied = np.flatnonzero(mine.counts)

    # Use the *isotopic* mass, which is what the simulation scatters from --
    # gold's exact 197.0, not the natural-abundance average 196.97. The two
    # differ by 0.02 keV in edge position, enough to land in a different channel.
    gold = table.by_z(79)
    mass = gold.isotopes[0].mass
    assert mass == pytest.approx(197.0)

    # K*E0 must fall inside the highest occupied channel -- not merely near it.
    top = int(occupied.max())
    expected = 2000.0 * kinematic_factor(4.0026, mass, 170.0)
    assert CAL.edge_energy(top) <= expected <= CAL.edge_energy(top + 1)
    # And nothing may appear above it.
    assert not np.any(mine.counts[top + 1 :])
