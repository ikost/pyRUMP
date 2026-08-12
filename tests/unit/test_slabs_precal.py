"""M5 (part 1) acceptance: slab discretization and the inbound energy march.

No oracle comparison yet -- reaching ``SimPrecal`` through the C requires the
whole ``SAMM``/``ALTBUF`` setup, so brick-level comparison waits for the
``ORACLE_DUMP`` hook. These tests pin the algorithm against physical invariants
and against the analytic limits the 1985 paper states, which catch a different
class of bug than an oracle diff does.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pyrump.atomic.tables import PeriodicTable
from pyrump.io.kalbitzer import parse_kalbitzer
from pyrump.sim.precal import (
    bohr_straggle_constant,
    energy_loss_step,
    march_inbound,
    surface_energy_loss,
)
from pyrump.sim.slabs import (
    DEFAULT_MAXPATH,
    build_uniform_grid,
    path_limited_maxpath,
    sublayer_count,
)
from pyrump.stopping.bragg import bragg_coefficients
from pyrump.stopping.kalbitzer import KalbitzerStopping
from pyrump.stopping.registry import StoppingRegistry
from pyrump.stopping.table import StoppingTable
from pyrump.stopping.ziegler import ZieglerStopping


from conftest import data_dir

DATA = data_dir()
pytestmark = pytest.mark.skipif(DATA is None, reason="legacy data tables unavailable")


@pytest.fixture(scope="module")
def registry() -> StoppingRegistry:
    assert DATA is not None
    table = PeriodicTable.load(DATA / "atom4.dat", DATA / "pscoef.dat")
    return StoppingRegistry(
        table.elements,
        kalbitzer=KalbitzerStopping(parse_kalbitzer(DATA / "newstop.kal"), table.elements),
        ziegler=ZieglerStopping(table.elements),
    )


@pytest.fixture(scope="module")
def silicon(registry) -> StoppingTable:
    return StoppingTable.build(registry, 2, 4.0026, 2.0, [14])


# --------------------------------------------------------------- slab grid


def test_sublayer_count_uses_truncation_not_rounding():
    """int(1 + t/maxpath): a layer exactly maxpath thick gets 2 slabs."""
    assert sublayer_count(200.0, 200.0) == 2
    assert sublayer_count(199.0, 200.0) == 1
    assert sublayer_count(1000.0, 200.0) == 6
    assert sublayer_count(0.0, 200.0) == 1


def test_sublayer_count_precedence():
    assert sublayer_count(1000.0, 200.0, explicit=3) == 3
    assert sublayer_count(1000.0, 200.0, sub_thickness=100.0) == 11
    # Explicit count beats explicit sub-thickness.
    assert sublayer_count(1000.0, 200.0, explicit=3, sub_thickness=100.0) == 3


def test_tilt_shortens_the_step():
    """Step size is path-length based, so tilting produces more, thinner slabs."""
    flat = path_limited_maxpath(DEFAULT_MAXPATH, 1.0, 1.0)
    tilted = path_limited_maxpath(DEFAULT_MAXPATH, 2.0, 3.0)
    assert flat == DEFAULT_MAXPATH
    assert tilted == pytest.approx(DEFAULT_MAXPATH / 3.0)

    grid_flat = build_uniform_grid([1000.0], [[1.0]], [14], sec_in=1.0, sec_out=1.0)
    grid_tilted = build_uniform_grid([1000.0], [[1.0]], [14], sec_in=2.0, sec_out=3.0)
    assert grid_tilted.n_slab > grid_flat.n_slab


def test_grid_conserves_thickness_and_composition():
    grid = build_uniform_grid([500.0, 1500.0], [[1.0, 0.0], [1.0, 2.0]], [14, 8])
    assert grid.total_thickness == pytest.approx(2000.0)
    # Layer 2 is SiO2: one third Si, two thirds O.
    layer2 = grid.composition[grid.layer_index == 1]
    fractions = layer2.sum(axis=0) / layer2.sum()
    assert fractions == pytest.approx([1 / 3, 2 / 3])
    # Every slab's per-element densities sum to its total.
    assert np.allclose(grid.composition.sum(axis=1), grid.areal_density)


def test_grid_skips_empty_layers():
    grid = build_uniform_grid([500.0, 0.0, 500.0], [[1.0]] * 3, [14])
    assert set(np.unique(grid.layer_index)) == {0, 2}


def test_grid_depth_is_monotonic():
    grid = build_uniform_grid([1000.0], [[1.0]], [14])
    assert np.all(np.diff(grid.depth) > 0)
    assert grid.depth[-1] == pytest.approx(1000.0)


def test_grid_rejects_shape_mismatch():
    with pytest.raises(ValueError, match="composition rows"):
        build_uniform_grid([100.0, 200.0], [[1.0]], [14])
    with pytest.raises(ValueError, match="composition columns"):
        build_uniform_grid([100.0], [[1.0, 2.0]], [14])


# ------------------------------------------------------------ energy march


def test_third_order_reduces_to_surface_approximation():
    """Zeroing the derivatives leaves the first-order term (1985 paper p.345)."""
    assert energy_loss_step(50.0, 0.0, 0.0, 1e-3) == pytest.approx(50.0 * 1e-3)


def test_third_order_correction_is_small_but_signed():
    """dS/dE < 0 above the stopping peak, so the correction raises the loss."""
    first_order = energy_loss_step(50.0, 0.0, 0.0, 1e-3)
    corrected = energy_loss_step(50.0, -0.02, 1e-4, 1e-3)
    assert corrected > first_order
    assert abs(corrected - first_order) / first_order < 0.01


def test_inbound_march_loses_energy_monotonically(registry, silicon):
    grid = build_uniform_grid([2000.0], [[1.0]], [14])
    coefficients = bragg_coefficients(silicon, grid.composition, grid.element_z)
    path = march_inbound(
        silicon,
        coefficients,
        grid.composition,
        grid.element_z,
        e0_keV=2000.0,
        sec_in=1.0,
        cutoff_keV=silicon.cutoff * 1000.0,
    )
    energies = path.valid()
    assert energies[0] == 2000.0
    assert np.all(np.diff(energies) < 0)
    assert path.reached == grid.n_slab


def test_inbound_march_matches_the_surface_approximation_for_a_thin_slab(
    registry, silicon
):
    """For a thin enough slab the higher-order terms must vanish."""
    grid = build_uniform_grid([1.0], [[1.0]], [14], explicit_sublayers=[1])
    coefficients = bragg_coefficients(silicon, grid.composition, grid.element_z)
    path = march_inbound(
        silicon,
        coefficients,
        grid.composition,
        grid.element_z,
        e0_keV=2000.0,
        sec_in=1.0,
        cutoff_keV=silicon.cutoff * 1000.0,
    )
    loss = path.energy[0] - path.energy[1]
    surface = surface_energy_loss(silicon, coefficients[0], 2000.0, 1.0)
    assert loss == pytest.approx(surface, rel=1e-4)


def test_inbound_march_converges_with_slab_count(registry, silicon):
    """Refining the grid must not move the exit energy much.

    This is the 1985 paper's own justification for maxpth=200: the third-order
    expansion is accurate on thick slabs, so few are needed.
    """
    exits = []
    for count in (1, 2, 5, 20, 100):
        grid = build_uniform_grid(
            [2000.0], [[1.0]], [14], explicit_sublayers=[count]
        )
        coefficients = bragg_coefficients(silicon, grid.composition, grid.element_z)
        path = march_inbound(
            silicon,
            coefficients,
            grid.composition,
            grid.element_z,
            e0_keV=2000.0,
            sec_in=1.0,
            cutoff_keV=silicon.cutoff * 1000.0,
        )
        exits.append(path.energy[path.reached])

    converged = exits[-1]
    # A single 2000e15 at/cm^2 slab is far beyond the recommended step and is
    # allowed to be poor; from maxpth-sized steps on it must be tight.
    assert abs(exits[3] - converged) / converged < 1e-5
    assert abs(exits[2] - converged) / converged < 1e-3


def test_tilt_increases_energy_loss(registry, silicon):
    grid = build_uniform_grid([1000.0], [[1.0]], [14])
    coefficients = bragg_coefficients(silicon, grid.composition, grid.element_z)
    losses = []
    for sec_in in (1.0, 2.0):
        path = march_inbound(
            silicon,
            coefficients,
            grid.composition,
            grid.element_z,
            e0_keV=2000.0,
            sec_in=sec_in,
            cutoff_keV=silicon.cutoff * 1000.0,
        )
        losses.append(path.energy[0] - path.energy[path.reached])
    # Twice the path, close to twice the loss (not exactly: S varies with E).
    assert 1.9 < losses[1] / losses[0] < 2.2


def test_march_stops_at_cutoff(registry, silicon):
    """Below cutoff the stopping fit is untrustworthy and the march bails."""
    grid = build_uniform_grid([50_000.0], [[1.0]], [14])
    coefficients = bragg_coefficients(silicon, grid.composition, grid.element_z)
    path = march_inbound(
        silicon,
        coefficients,
        grid.composition,
        grid.element_z,
        e0_keV=2000.0,
        sec_in=1.0,
        cutoff_keV=silicon.cutoff * 1000.0,
    )
    assert path.reached < grid.n_slab, "should have stopped early"
    assert path.energy[path.reached] < 2000.0


# -------------------------------------------------------------- straggling


def test_straggling_is_off_by_default():
    """RUMP's SIM STRAGGLE defaults to 0, disabling the model entirely."""
    assert bohr_straggle_constant(2, scale=0.0) == 0.0


def test_bohr_constant_scales_with_z_squared():
    assert bohr_straggle_constant(2) / bohr_straggle_constant(1) == pytest.approx(4.0)


def test_straggling_accumulates_with_depth(registry, silicon):
    grid = build_uniform_grid([2000.0], [[1.0]], [14])
    coefficients = bragg_coefficients(silicon, grid.composition, grid.element_z)
    path = march_inbound(
        silicon,
        coefficients,
        grid.composition,
        grid.element_z,
        e0_keV=2000.0,
        sec_in=1.0,
        cutoff_keV=silicon.cutoff * 1000.0,
        straggle_scale=1.0,
    )
    variance = path.straggle[: path.reached + 1]
    assert variance[0] == 0.0
    assert np.all(np.diff(variance) > 0), "variance must grow with depth"

    # Bohr variance is strictly linear in traversed areal density, so with equal
    # slabs the increments are identical. (Note it is linear in *depth*, not in
    # slab index fraction -- the grid need not divide evenly.)
    increments = np.diff(variance)
    assert np.allclose(increments, increments[0], rtol=1e-12)

    depth = np.concatenate([[0.0], grid.depth[: path.reached]])
    expected = variance[path.reached] * depth / depth[-1]
    assert np.allclose(variance, expected, rtol=1e-12)
