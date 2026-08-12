"""M10 acceptance: absorber layers, fuzz, pile-up and multiple scattering.

Two things had to be right before any of these matched:

1. The inbound march **starts at** ``fsurf`` at full beam energy. Absorber
   layers sit between sample and detector, so the incoming beam never crosses
   them — only the outgoing particle does.
2. RUMP's ``NDTRI`` must be built unoptimised. ``poly_e`` reads one element past
   its coefficient array; at ``-O0`` that is harmless, at ``-O2`` the compiler
   exploits the UB and NDTRI returns ±0.15 for every mid-range argument. The
   shipped RUMP passes no ``-O`` flag, so it runs the benign version.
"""

from __future__ import annotations

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


from conftest import data_dir

DATA = data_dir()
CAL = Calibration(kevch=5.0, kev0=0.0, first=0.0, npt=1024)


# --------------------------------------------------------------- absorber


def test_first_sample_slab():
    layer_index = np.array([0, 0, 1, 1, 1, 2, 2])
    assert first_sample_slab(layer_index, 0) == 0
    assert first_sample_slab(layer_index, 1) == 2
    assert first_sample_slab(layer_index, 2) == 5
    # More absorber layers than exist: nothing is sample.
    assert first_sample_slab(layer_index, 9) == layer_index.size


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


# ------------------------------------------------------------------- fuzz


@pytest.mark.parametrize("steps", [2, 3, 4, 5, 7, 8])
def test_fuzz_weights_sum_to_one(steps):
    """Fuzzing redistributes yield; it must not create or destroy it."""
    total = sum(s.weight for s in fuzz_steps(100.0, steps))
    assert total == pytest.approx(1.0, rel=1e-12)


@pytest.mark.parametrize("steps", [2, 3, 4, 5, 7])
def test_fuzz_offsets_are_symmetric(steps):
    deltas = sorted(s.delta for s in fuzz_steps(100.0, steps))
    assert deltas == pytest.approx([-d for d in reversed(deltas)], abs=1e-9)


def test_fuzz_is_a_no_op_without_steps():
    for steps in (0, 1):
        result = fuzz_steps(100.0, steps)
        assert len(result) == 1
        assert result[0].delta == 0.0 and result[0].weight == 1.0
    assert fuzz_steps(0.0, 5)[0].delta == 0.0


def test_fuzz_spread_scales_with_amount():
    small = max(abs(s.delta) for s in fuzz_steps(100.0, 5))
    large = max(abs(s.delta) for s in fuzz_steps(300.0, 5))
    assert large == pytest.approx(3.0 * small, rel=1e-12)


def test_fuzz_uses_one_over_sqrt_two():
    """The C's .7071067 means the parameter is not the roughness sigma."""
    steps = fuzz_steps(100.0, 3)
    extreme = max(s.delta for s in steps)
    # k=0 of 3 steps: y = 0.5 * size, size = 1/5.
    assert extreme == pytest.approx(100.0 * abs(ndtri(0.1)) * 0.7071067, rel=1e-6)


def test_iteration_count_multiplies():
    """Every fuzzed layer multiplies the cost -- this is the expensive option."""
    assert iteration_count([0, 0]) == 1
    assert iteration_count([5, 0]) == 5
    assert iteration_count([5, 3]) == 15
    assert iteration_count([5, 3, 2]) == 30


def test_replica_enumeration():
    thicknesses = np.array([1000.0, 500.0])
    replicas = replica_thicknesses(thicknesses, [100.0, 50.0], [3, 2])
    assert len(replicas) == 6
    assert sum(a for _, a in replicas) == pytest.approx(1.0, rel=1e-12)

    # Each layer's offsets are exactly its own fuzz_steps deltas -- the mixed
    # radix enumeration pairs every step of layer 0 with every step of layer 1.
    for column, (amount, steps) in enumerate([(100.0, 3), (50.0, 2)]):
        nominal = thicknesses[column]
        expected = {round(nominal + s.delta, 9) for s in fuzz_steps(amount, steps)}
        assert {round(t[column], 9) for t, _ in replicas} == expected

    # Weighted mean thickness is preserved in both layers, so fuzzing perturbs
    # without biasing.
    for column, nominal in enumerate(thicknesses):
        mean = sum(t[column] * a for t, a in replicas)
        assert mean == pytest.approx(nominal, abs=1e-9)


# --------------------------------------------------------------- pile-up


def test_pileup_needs_tau_and_current():
    counts = np.linspace(0, 100, 50)
    assert np.array_equal(new_pileup(counts, tau_us=0, current_nA=10, charge_uC=10), counts)
    assert np.array_equal(new_pileup(counts, tau_us=5, current_nA=0, charge_uC=10), counts)


def test_pileup_moves_counts_upward():
    counts = np.zeros(200)
    counts[40:60] = 1000.0
    out = new_pileup(counts, tau_us=5.0, current_nA=50.0, charge_uC=10.0)
    # Pile-up sums two events, so intensity appears near twice the energy.
    assert out[100:130].sum() > 0
    assert out.size >= counts.size


def test_pileup_conserves_counts_approximately():
    """Every pile-up event removes two real ones, so the total drops slightly."""
    counts = np.zeros(200)
    counts[40:60] = 1000.0
    out = new_pileup(counts, tau_us=5.0, current_nA=50.0, charge_uC=10.0)
    assert out.sum() < counts.sum()
    assert out.sum() == pytest.approx(counts.sum(), rel=0.05)


# --------------------------------------------------- multiple scattering


def test_multiple_scattering_disabled_by_default():
    counts = np.linspace(1, 100, 50)
    assert np.array_equal(
        add_multiple_scattering(counts, strength=0.0, charge_uC=10, omega_msr=1), counts
    )


def test_multiple_scattering_adds_a_low_energy_tail():
    counts = np.zeros(200)
    counts[150:160] = 1000.0
    out = add_multiple_scattering(counts, strength=100.0, charge_uC=10, omega_msr=1)
    assert np.all(out >= counts - 1e-12)
    # Everything below the peak gains; nothing above it does.
    assert out[:150].sum() > 0
    assert np.allclose(out[160:], counts[160:])


# ------------------------------------------------------- against the C

pytestmark = pytest.mark.skipif(
    DATA is None or not ora.available() or ora.data_dir() is None, reason="legacy tables or oracle unavailable"
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


def _compare(oracle, registry, table, sample, measurement, theta=0.0, *, rtol=1e-4):
    mine = simulate(
        sample, Beam(), Geometry(theta=theta, phi=10.0), registry, table, CAL, measurement
    )
    oracle.set_beam(
        e0_MeV=2.0, phi=10.0, theta=theta, kevch=CAL.kevch, kev0=CAL.kev0,
        first=CAL.first, npt=CAL.npt, fwhm=measurement.fwhm_keV,
        q_uC=measurement.charge_uC, omega=measurement.omega_msr, corr=1.0,
        current_nA=measurement.current_nA, tau=measurement.tau_us,
    )
    oracle.set_sample(
        sample.thicknesses, sample.element_z, sample.compositions,
        straggle=sample.straggle, multiple=sample.multiple,
    )
    oracle.set_absorber(sample.absorber_layers)
    for i in range(len(sample.thicknesses)):
        oracle.set_layer_fuzz(i, 0.0, 0)
    if sample.fuzz_steps:
        amounts = sample.fuzz_amounts or [0.0] * len(sample.fuzz_steps)
        for i, steps in enumerate(sample.fuzz_steps):
            if steps > 1:
                oracle.set_layer_fuzz(i, amounts[i], steps)

    theirs = oracle.simulate_spectrum()
    n = min(len(theirs), len(mine.counts))
    ours, ref = mine.counts[:n], theirs[:n]
    assert ours.sum() == pytest.approx(ref.sum(), rel=rtol)
    assert np.max(np.abs(ours - ref)) < rtol * ref.max()
    return ours, ref


def test_ndtri_matches_the_c(oracle):
    """Guards the -O0 build requirement -- at -O2 this fails loudly."""
    p = np.array([0.05, 0.1, 0.2, 0.25, 0.3, 0.5, 0.7, 0.75, 0.9])
    assert np.allclose(oracle.ndtri(p), ndtri(p), atol=1e-8)


def test_absorber_matches_oracle(oracle, registry, table):
    sample = UniformSample(
        [200.0, 1000.0], [79, 14], [[1, 0], [0, 1]], absorber_layers=1
    )
    _compare(oracle, registry, table, sample, Measurement(omega_msr=1.0, charge_uC=10.0))


def test_absorber_with_tilt_matches_oracle(oracle, registry, table):
    """The absorber is not tilted with the sample -- exit is normal incidence."""
    sample = UniformSample(
        [200.0, 1000.0], [79, 14], [[1, 0], [0, 1]], absorber_layers=1
    )
    _compare(
        oracle, registry, table, sample,
        Measurement(omega_msr=1.0, charge_uC=10.0), theta=45.0,
    )


def test_absorber_attenuates(registry, table):
    """Physical check: an absorber shifts the sample's edge down in energy."""
    geometry = Geometry(theta=0.0, phi=10.0)
    measurement = Measurement(omega_msr=1.0, charge_uC=10.0)
    bare = simulate(
        UniformSample([1000.0], [14], [[1.0]]), Beam(), geometry, registry, table,
        CAL, measurement,
    )
    behind = simulate(
        UniformSample([200.0, 1000.0], [79, 14], [[1, 0], [0, 1]], absorber_layers=1),
        Beam(), geometry, registry, table, CAL, measurement,
    )
    assert np.flatnonzero(behind.counts).max() < np.flatnonzero(bare.counts).max()


@pytest.mark.parametrize("steps", [2, 3, 4, 5, 7])
def test_fuzz_matches_oracle(oracle, registry, table, steps):
    sample = UniformSample(
        [1000.0], [14], [[1.0]], fuzz_amounts=[150.0], fuzz_steps=[steps]
    )
    _compare(oracle, registry, table, sample, Measurement(omega_msr=1.0, charge_uC=10.0))


def test_fuzz_with_straggling_matches_oracle(oracle, registry, table):
    sample = UniformSample(
        [1000.0], [14], [[1.0]], straggle=1.0, fuzz_amounts=[150.0], fuzz_steps=[3]
    )
    _compare(oracle, registry, table, sample, Measurement(omega_msr=1.0, charge_uC=10.0))


@pytest.mark.parametrize("strength", [1.0, 5.0])
def test_multiple_scattering_matches_oracle(oracle, registry, table, strength):
    sample = UniformSample([1000.0], [14], [[1.0]], multiple=strength)
    _compare(oracle, registry, table, sample, Measurement(omega_msr=1.0, charge_uC=10.0))


def test_pileup_matches_oracle(oracle, registry, table):
    sample = UniformSample([1000.0], [14], [[1.0]])
    measurement = Measurement(
        omega_msr=1.0, charge_uC=10.0, current_nA=20.0, tau_us=5.0
    )
    # Pile-up is a full self-convolution in float32; a little looser than the
    # rest of the pipeline, and the thesis rates the model itself at "30%".
    _compare(oracle, registry, table, sample, measurement, rtol=1e-3)


def test_fuzz_broadens_an_edge(registry, table):
    """Roughness should soften the back edge without moving the total."""
    geometry = Geometry(theta=0.0, phi=10.0)
    measurement = Measurement(omega_msr=1.0, charge_uC=10.0)
    sharp = simulate(
        UniformSample([1000.0], [14], [[1.0]]), Beam(), geometry, registry, table,
        CAL, measurement,
    )
    rough = simulate(
        UniformSample([1000.0], [14], [[1.0]], fuzz_amounts=[300.0], fuzz_steps=[5]),
        Beam(), geometry, registry, table, CAL, measurement,
    )
    assert rough.total() == pytest.approx(sharp.total(), rel=0.02)
    assert np.count_nonzero(rough.counts) > np.count_nonzero(sharp.counts)


# --------------------------------------------------------------- guard tests
#
# Merged from a parallel M10 test file. These assert that an effect actually
# changes the spectrum -- without them a comparison can pass because both
# implementations did nothing.


def test_old_pileup_is_available_for_legacy_results():
    counts = np.zeros(200)
    counts[50] = 1000.0
    out = old_pileup(counts, current_nA=100.0)
    assert out[100] > 0, "counts should appear at twice the energy"


def test_pileup_scales_with_rate():
    counts = np.zeros(400)
    counts[100:150] = 1000.0
    low = new_pileup(counts, tau_us=5.0, current_nA=100.0, charge_uC=10.0)
    high = new_pileup(counts, tau_us=5.0, current_nA=1000.0, charge_uC=10.0)
    assert high[200:].sum() > low[200:].sum()


def test_multiple_scattering_grows_with_strength():
    counts = np.zeros(200)
    counts[150] = 1e6
    weak = add_multiple_scattering(counts, strength=10.0, charge_uC=10.0, omega_msr=1.0)
    strong = add_multiple_scattering(counts, strength=100.0, charge_uC=10.0, omega_msr=1.0)
    assert strong[:150].sum() > weak[:150].sum()


def test_pileup_actually_changes_the_spectrum(oracle, registry, table):
    """Guard against the comparison passing because nothing happened."""
    sample = UniformSample([1000.0], [14], [[1.0]])
    oracle.set_pileup(None)
    _, without = _compare(
        oracle, registry, table, sample,
        Measurement(omega_msr=1.0, charge_uC=10.0, tau_us=0.0),
    )
    oracle.set_pileup("new")
    _, with_pileup = _compare(
        oracle, registry, table, sample,
        Measurement(omega_msr=1.0, charge_uC=10.0, tau_us=5.0, current_nA=100.0),
    )
    oracle.set_pileup(None)
    assert not np.allclose(without, with_pileup)


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
