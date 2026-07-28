"""M10 acceptance: absorber, fuzz, pile-up and multiple scattering.

These are sample- and detector-level effects layered on top of the forward
model. Stage order follows ``SimCreateDetails`` (creatr.c:307-345): fill,
detector convolution, yield normalisation, pile-up, multiple-scattering tail.

The subtle one is the absorber. It sits between the sample and the *detector*,
so the beam reaches the sample undegraded — ``SimPrecal`` seeds the inbound
march at ``fsurf`` rather than at slab 0 (creatr.c:1533). Marching through the
absorber on the way in instead costs ~20% of the total yield.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.special import ndtri

from pyrump.atomic.tables import PeriodicTable
from pyrump.io.kalbitzer import parse_kalbitzer
from pyrump.model.detector import Measurement
from pyrump.model.geometry import Geometry
from pyrump.model.spectrum import Calibration
from pyrump.sim.absorber import AbsorberSpec, first_sample_slab, sample_depth
from pyrump.sim.engine import Beam, UniformSample, simulate
from pyrump.sim.fuzz import fuzz_steps, iteration_count, replica_thicknesses
from pyrump.sim.multiscatter import add_multiple_scattering
from pyrump.sim.pileup import new_pileup, old_pileup
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


# ------------------------------------------------------------------- fuzz


def test_fuzz_weights_sum_to_one():
    """Fuzzing redistributes yield; it must not create or destroy any."""
    for steps in (2, 3, 4, 5, 8, 9):
        weights = [s.weight for s in fuzz_steps(100.0, steps)]
        assert sum(weights) == pytest.approx(1.0, rel=1e-9), f"steps={steps}"


def test_fuzz_is_symmetric_about_nominal():
    """The weighted mean thickness offset must be zero."""
    for steps in (2, 3, 4, 5, 8):
        replicas = fuzz_steps(100.0, steps)
        mean = sum(s.weight * s.delta for s in replicas)
        assert mean == pytest.approx(0.0, abs=1e-9), f"steps={steps}"


def test_fuzz_spread_grows_with_amount():
    def spread(amount):
        replicas = fuzz_steps(amount, 8)
        return np.sqrt(sum(s.weight * s.delta**2 for s in replicas))

    assert spread(200.0) == pytest.approx(2 * spread(100.0), rel=1e-9)


def test_fuzz_disabled_gives_one_unperturbed_replica():
    for steps in (0, 1):
        replicas = fuzz_steps(100.0, steps)
        assert len(replicas) == 1
        assert replicas[0].delta == 0.0
        assert replicas[0].weight == 1.0
    # A zero amount also collapses, even with many steps requested.
    assert len(fuzz_steps(0.0, 9)) == 1


def test_fuzz_iteration_count_multiplies():
    """Every fuzzed layer multiplies the cost -- the most expensive option."""
    assert iteration_count([5]) == 5
    assert iteration_count([5, 5]) == 25
    assert iteration_count([3, 0, 4]) == 12
    assert iteration_count([0, 0]) == 1


def test_fuzz_replicas_enumerate_the_mixed_radix_counter():
    thicknesses = np.array([100.0, 200.0])
    replicas = replica_thicknesses(thicknesses, [10.0, 20.0], [3, 2])
    assert len(replicas) == 6
    assert sum(a for _, a in replicas) == pytest.approx(1.0, rel=1e-9)
    # Layer 1 is untouched in the first three replicas' inner index and so on;
    # what matters is that the weighted mean thickness is the nominal one.
    mean = sum(a * t for t, a in replicas)
    assert np.allclose(mean, thicknesses, atol=1e-9)


# --------------------------------------------------------------- absorber


def test_first_sample_slab():
    layer_index = np.array([0, 0, 1, 1, 1])
    assert first_sample_slab(layer_index, 0) == 0
    assert first_sample_slab(layer_index, 1) == 2
    # More absorber layers than exist: nothing is sample.
    assert first_sample_slab(layer_index, 5) == 5


def test_absorber_slabs_do_not_advance_depth():
    """Depth is measured from the sample surface, not the absorber's."""
    areal = np.array([50.0, 50.0, 100.0, 100.0])
    depth = sample_depth(areal, first_slab=2)
    assert depth[0] == 0.0 and depth[1] == 0.0
    assert depth[2] == 100.0
    assert depth[3] == 200.0


def test_fres_only_absorber_is_rejected():
    with pytest.raises(NotImplementedError, match="forward-recoil"):
        AbsorberSpec(layers=1, fres_only=True)


# ----------------------------------------------------------------- pile-up


def test_pileup_requires_tau_and_current():
    counts = np.linspace(0.0, 100.0, 200)
    assert np.array_equal(new_pileup(counts, tau_us=0.0, current_nA=10.0, charge_uC=10.0), counts)
    assert np.array_equal(new_pileup(counts, tau_us=5.0, current_nA=0.0, charge_uC=10.0), counts)


def test_pileup_moves_counts_upward():
    """Coincident events are recorded at their combined energy."""
    counts = np.zeros(400)
    counts[100:150] = 1000.0
    out = new_pileup(counts, tau_us=5.0, current_nA=500.0, charge_uC=10.0)
    # Something appears above the original band...
    assert out[200:].sum() > 0
    # ...and the real counts are reduced to pay for it.
    assert out[100:150].sum() < counts[100:150].sum()


def test_pileup_scales_with_rate():
    counts = np.zeros(400)
    counts[100:150] = 1000.0
    low = new_pileup(counts, tau_us=5.0, current_nA=100.0, charge_uC=10.0)
    high = new_pileup(counts, tau_us=5.0, current_nA=1000.0, charge_uC=10.0)
    assert high[200:].sum() > low[200:].sum()


def test_old_pileup_is_available_for_legacy_results():
    counts = np.zeros(200)
    counts[50] = 1000.0
    out = old_pileup(counts, current_nA=100.0)
    assert out[100] > 0, "counts should appear at twice the energy"


# ---------------------------------------------------------- multiple scatter


def test_multiple_scattering_off_by_default():
    counts = np.linspace(1.0, 100.0, 50)
    assert np.array_equal(
        add_multiple_scattering(counts, strength=0.0, charge_uC=10.0, omega_msr=1.0),
        counts,
    )


def test_multiple_scattering_adds_a_low_energy_tail():
    counts = np.zeros(200)
    counts[150] = 1e6
    out = add_multiple_scattering(counts, strength=100.0, charge_uC=10.0, omega_msr=1.0)
    # Everything below the peak gains; nothing above it does.
    assert np.all(out[:150] >= counts[:150])
    assert out[:150].sum() > 0
    assert np.allclose(out[151:], counts[151:])


def test_multiple_scattering_grows_with_strength():
    counts = np.zeros(200)
    counts[150] = 1e6
    weak = add_multiple_scattering(counts, strength=10.0, charge_uC=10.0, omega_msr=1.0)
    strong = add_multiple_scattering(counts, strength=100.0, charge_uC=10.0, omega_msr=1.0)
    assert strong[:150].sum() > weak[:150].sum()


# ------------------------------------------------------------ against the C

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
    handle.set_pileup(None)
    return handle


def _compare(oracle, registry, table, sample, measurement, *, theta=0.0, pileup=None):
    oracle.set_pileup(pileup)
    geometry = Geometry(theta=theta, phi=10.0)
    mine = simulate(sample, Beam(), geometry, registry, table, CAL, measurement)

    oracle.set_beam(
        e0_MeV=2.0, phi=10.0, theta=theta, kevch=CAL.kevch, kev0=CAL.kev0,
        first=CAL.first, npt=CAL.npt, fwhm=measurement.fwhm_keV,
        q_uC=measurement.charge_uC, omega=measurement.omega_msr, corr=1.0, cbeam=1,
        current_nA=measurement.current_nA, tau=measurement.tau_us,
    )
    oracle.set_sample(
        sample.thicknesses, sample.element_z, sample.compositions,
        straggle=sample.straggle, multiple=sample.multiple,
    )
    oracle.set_absorber(sample.absorber_layers)
    theirs = oracle.simulate_spectrum()

    n = min(len(theirs), len(mine.counts))
    return mine.counts[:n], theirs[:n]


ABSORBER_CASES = [
    ("one layer", UniformSample([200.0, 1000.0], [14], [[1.0], [1.0]], absorber_layers=1), 0.0),
    ("tilted 45", UniformSample([200.0, 1000.0], [14], [[1.0], [1.0]], absorber_layers=1), 45.0),
    ("tilted 60", UniformSample([200.0, 1000.0], [14], [[1.0], [1.0]], absorber_layers=1), 60.0),
    ("thick", UniformSample([800.0, 1000.0], [14], [[1.0], [1.0]], absorber_layers=1), 0.0),
    ("two layers", UniformSample([200.0, 300.0, 1000.0], [14], [[1.0]] * 3, absorber_layers=2), 0.0),
]


@pytest.mark.parametrize(
    "label, sample, theta", ABSORBER_CASES, ids=[c[0] for c in ABSORBER_CASES]
)
def test_absorber_matches_oracle(oracle, registry, table, label, sample, theta):
    measurement = Measurement(omega_msr=1.0, charge_uC=10.0, tau_us=0.0)
    ours, ref = _compare(oracle, registry, table, sample, measurement, theta=theta)
    assert ours.sum() == pytest.approx(ref.sum(), rel=1e-4)
    assert np.max(np.abs(ours - ref)) < 1e-3 * ref.max()


def test_absorber_attenuates_without_adding_yield(oracle, registry, table):
    """The absorber must lower the spectrum, not shift or add to it."""
    measurement = Measurement(omega_msr=1.0, charge_uC=10.0, tau_us=0.0)
    geometry = Geometry(theta=0.0, phi=10.0)

    bare = simulate(
        UniformSample([1000.0], [14], [[1.0]]),
        Beam(), geometry, registry, table, CAL, measurement,
    )
    absorbed = simulate(
        UniformSample([200.0, 1000.0], [14], [[1.0], [1.0]], absorber_layers=1),
        Beam(), geometry, registry, table, CAL, measurement,
    )
    # Same number of scattering slabs, so the same total yield...
    assert absorbed.total() == pytest.approx(bare.total(), rel=0.02)
    # ...but shifted down in energy by the exit-path loss.
    centroid = lambda s: float((np.arange(s.size) * s).sum() / s.sum())  # noqa: E731
    assert centroid(absorbed.counts) < centroid(bare.counts)


def test_pileup_matches_oracle(oracle, registry, table):
    sample = UniformSample([1000.0], [14], [[1.0]])
    measurement = Measurement(
        omega_msr=1.0, charge_uC=10.0, tau_us=5.0, current_nA=100.0
    )
    ours, ref = _compare(oracle, registry, table, sample, measurement, pileup="new")
    assert ours.sum() == pytest.approx(ref.sum(), rel=1e-4)
    assert np.max(np.abs(ours - ref)) < 1e-3 * ref.max()


def test_pileup_actually_changes_the_spectrum(oracle, registry, table):
    """Guard against the comparison passing because nothing happened."""
    sample = UniformSample([1000.0], [14], [[1.0]])
    _, without = _compare(
        oracle, registry, table, sample,
        Measurement(omega_msr=1.0, charge_uC=10.0, tau_us=0.0),
    )
    _, with_pileup = _compare(
        oracle, registry, table, sample,
        Measurement(omega_msr=1.0, charge_uC=10.0, tau_us=5.0, current_nA=100.0),
        pileup="new",
    )
    assert not np.allclose(without, with_pileup)


@pytest.mark.parametrize("strength", [1.0, 50.0, 200.0])
def test_multiple_scattering_matches_oracle(oracle, registry, table, strength):
    sample = UniformSample([1000.0], [14], [[1.0]], multiple=strength)
    measurement = Measurement(omega_msr=1.0, charge_uC=10.0, tau_us=0.0)
    ours, ref = _compare(oracle, registry, table, sample, measurement)
    assert ours.sum() == pytest.approx(ref.sum(), rel=1e-4)
    assert np.max(np.abs(ours - ref)) < 1e-3 * ref.max()


def test_multiple_scattering_actually_changes_the_spectrum(oracle, registry, table):
    measurement = Measurement(omega_msr=1.0, charge_uC=10.0, tau_us=0.0)
    _, plain = _compare(
        oracle, registry, table, UniformSample([1000.0], [14], [[1.0]]), measurement
    )
    _, tailed = _compare(
        oracle, registry, table,
        UniformSample([1000.0], [14], [[1.0]], multiple=200.0), measurement,
    )
    assert tailed.sum() > plain.sum() * 1.05


# ------------------------------------------------------------------- NDTRI


def test_rump_ndtri_matches_scipy(oracle):
    """RUMP's inverse normal CDF, used by fuzz, agrees with scipy.

    See docs/rump-quirks.md: gvcalc.c's NDTRI reads one element past the end of
    its Taylor coefficient array. In RUMP's own build the adjacent object makes
    that harmless; the oracle pads the array so the probe evaluates what the
    code intends rather than whatever follows it in memory.
    """
    probabilities = np.array([1e-6, 0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999])
    theirs = np.array([float(oracle.ndtri(p)[0]) for p in probabilities])
    assert np.allclose(theirs, ndtri(probabilities), atol=1e-8)
