"""M5 infrastructure: brick-level capture from the unmodified simulation engine.

``creatr.c`` calls ``SimFillSpectrum`` once per brick, and that is a *function
pointer*. Pointing it at a recorder gives exact ground truth for the slab march
without patching a single line of the C.

These tests establish that the capture is trustworthy before pyRUMP's own slab
march is compared against it.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

from pyrump.atomic.tables import PeriodicTable
from pyrump.physics.kinematics import kinematic_factor

sys.path.insert(0, str(Path(__file__).parent))
import oracle as ora  # noqa: E402

pytestmark = [
    pytest.mark.oracle,
    pytest.mark.skipif(
        not ora.available() or ora.data_dir() is None, reason="oracle unavailable"
    ),
]

# Column layout of the captured brick array.
Z, MASS, EFRONT, EBACK, HFRONT, HBACK, QQQ, SIGF, SIGB = range(9)

BEAM_MASS = 4.0026
SCATTERING_ANGLE = 170.0  # entered as phi=10


@pytest.fixture
def oracle() -> ora.Oracle:
    handle = ora.Oracle.load()
    handle.reset_stopping_tables()
    handle.set_beam(e0_MeV=2.0, phi=10.0, theta=0.0, kevch=5.0, npt=1024)
    return handle


@pytest.fixture(scope="module")
def silicon_isotopes():
    data = ora.data_dir()
    assert data is not None
    element = PeriodicTable.load(data / "atom4.dat").by_z(14)
    # The engine emits one block per isotope, heaviest first.
    return sorted(element.isotopes, key=lambda i: -i.mass)


def test_capture_produces_one_block_per_isotope(oracle, silicon_isotopes):
    oracle.set_sample([1000.0], [14], [[1.0]])
    bricks = oracle.simulate_bricks()

    # 1000e15 at/cm^2 at the default maxpth=200 gives int(1 + 1000/200) = 6 slabs.
    assert len(bricks) == 6 * len(silicon_isotopes)
    # The z/mass columns identify the *scattered particle*, not the target --
    # they exist for the stopper-foil lookup (anlyz.c:180).
    assert np.all(bricks[:, Z] == 2)
    assert np.allclose(bricks[:, MASS], BEAM_MASS)


def test_surface_edges_match_kinematics(oracle, silicon_isotopes):
    """Independent confirmation of M4 against the running engine.

    The front energy of each block's first brick must be exactly K*E0.
    """
    oracle.set_sample([1000.0], [14], [[1.0]])
    bricks = oracle.simulate_bricks()
    per_block = len(bricks) // len(silicon_isotopes)

    for index, isotope in enumerate(silicon_isotopes):
        block = bricks[index * per_block : (index + 1) * per_block]
        expected = 2000.0 * kinematic_factor(BEAM_MASS, isotope.mass, SCATTERING_ANGLE)
        assert block[0, EFRONT] == pytest.approx(expected, abs=0.01)


def test_brick_heights_track_isotopic_abundance(oracle, silicon_isotopes):
    oracle.set_sample([1000.0], [14], [[1.0]])
    bricks = oracle.simulate_bricks()
    per_block = len(bricks) // len(silicon_isotopes)

    heights = [bricks[i * per_block, HFRONT] for i in range(len(silicon_isotopes))]
    abundances = [i.fraction for i in silicon_isotopes]
    ratios = [h / heights[0] for h in heights]
    expected = [a / abundances[0] for a in abundances]
    # Not exact: K and the stopping factor [eps] differ slightly per isotope.
    assert np.allclose(ratios, expected, rtol=0.01)


def test_bricks_descend_in_energy(oracle):
    """Deeper slabs scatter at lower energy; each brick spans front > back."""
    oracle.set_sample([1000.0], [14], [[1.0]])
    bricks = oracle.simulate_bricks()
    assert np.all(bricks[:, EFRONT] > bricks[:, EBACK])

    block = bricks[:6]
    assert np.all(np.diff(block[:, EFRONT]) < 0)
    # Bricks tile without gaps: one slab's back edge is the next slab's front.
    assert np.allclose(block[:-1, EBACK], block[1:, EFRONT], rtol=1e-6)


def test_yield_rises_with_depth(oracle):
    """Rutherford 1/E^2 means the plateau climbs as the beam slows."""
    oracle.set_sample([1000.0], [14], [[1.0]])
    block = oracle.simulate_bricks()[:6]
    assert np.all(np.diff(block[:, HFRONT]) > 0)


def test_slab_count_follows_maxpth(oracle, silicon_isotopes):
    n_isotopes = len(silicon_isotopes)
    for maxpth, expected_slabs in ((200.0, 6), (100.0, 11), (500.0, 3)):
        oracle.set_sample([1000.0], [14], [[1.0]], maxpth=maxpth)
        bricks = oracle.simulate_bricks()
        assert len(bricks) == expected_slabs * n_isotopes, f"maxpth={maxpth}"


def test_explicit_sublayers_override_maxpth(oracle, silicon_isotopes):
    oracle.set_sample([1000.0], [14], [[1.0]], sublayers=[4])
    bricks = oracle.simulate_bricks()
    assert len(bricks) == 4 * len(silicon_isotopes)


def test_straggling_is_zero_unless_enabled(oracle):
    """RUMP's SIM STRAGGLE defaults to 0, so sigf/sigb stay zero."""
    oracle.set_sample([1000.0], [14], [[1.0]])
    bricks = oracle.simulate_bricks()
    assert np.all(bricks[:, SIGF] == 0.0)
    assert np.all(bricks[:, SIGB] == 0.0)

    oracle.set_sample([1000.0], [14], [[1.0]], straggle=1.0)
    bricks = oracle.simulate_bricks()
    assert np.any(bricks[:, SIGB] > 0.0)


def test_straggling_grows_with_depth(oracle):
    oracle.set_sample([1000.0], [14], [[1.0]], straggle=1.0)
    block = oracle.simulate_bricks()[:6]
    assert np.all(np.diff(block[:, SIGB]) > 0)


def test_qqq_is_computed_even_though_the_fill_discards_it(oracle):
    """The 1985 paper's Rutherford integral is still calculated every run.

    ``SimAnlyz`` ignores the argument (the parabolic path is ``#if 0`` at
    anlyz.c:496), but ``SimPrecal`` fills ``layer[].qq`` regardless. Capturing it
    gives pyRUMP a reference for the opt-in parabolic mode.
    """
    oracle.set_sample([1000.0], [14], [[1.0]])
    bricks = oracle.simulate_bricks()
    assert np.all(bricks[:, QQQ] > 0.0)
    assert np.all(np.diff(bricks[:6, QQQ]) > 0)


def test_multilayer_sample(oracle):
    """Two layers of different composition produce distinguishable blocks."""
    oracle.set_sample(
        [500.0, 500.0],
        [14, 79],
        [[1.0, 0.0], [0.0, 1.0]],
        sublayers=[2, 2],
    )
    bricks = oracle.simulate_bricks()
    # Si has 3 isotopes, Au is monoisotopic: 3*2 + 1*2 = 8 bricks.
    assert len(bricks) == 8

    gold_edge = bricks[:, EFRONT].max()
    silicon_surface = 2000.0 * kinematic_factor(BEAM_MASS, 27.977, SCATTERING_ANGLE)
    gold_surface = 2000.0 * kinematic_factor(BEAM_MASS, 196.97, SCATTERING_ANGLE)

    # Gold sits far above silicon...
    assert gold_edge > silicon_surface + 500.0
    # ...but below K*E0, because it is buried under 500e15 at/cm^2 of silicon
    # and the beam arrives already slowed. That inbound loss is the whole point
    # of the depth scale.
    assert gold_edge < gold_surface
    assert gold_surface - gold_edge == pytest.approx(45.0, abs=15.0)


def test_simulation_is_deterministic(oracle):
    oracle.set_sample([1000.0], [14], [[1.0]])
    first = oracle.simulate_bricks()
    second = oracle.simulate_bricks()
    assert np.array_equal(first, second)


def test_full_spectrum_path_still_works(oracle):
    """With capture off, the real SimAnlyz runs and fills channels."""
    oracle.set_sample([1000.0], [14], [[1.0]])
    spectrum = oracle.simulate_spectrum()
    assert spectrum.size > 0
    assert np.any(spectrum > 0)
    # A 1000e15 at/cm^2 film gives a plateau, not a spike.
    assert np.count_nonzero(spectrum) > 20
