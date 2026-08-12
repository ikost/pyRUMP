"""M9 acceptance: depth-profile equations.

Every ``EQUATION`` form is compared brick-for-brick against the C. Three things
had to be right before any of them matched, and each was invisible with uniform
composition:

1. ``hfront`` is recomputed from *this* slab's areal density every iteration --
   only ``efront``/``ratde`` carry over (creatr.c:1804 sits outside the ``ok``
   guard). Reusing it shifts the whole spectrum by one slab.
2. Layer density comes from **inverse-density averaging** over the layer's own
   composition, not the 0.4997 silicon fallback. It sets the depth scale of
   every depth-dependent form.
3. A layer carrying an equation takes its sublayer count from the equation's
   recommended value, ignoring ``maxpth`` entirely.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

from pyrump.atomic.density import SILICON_DENSITY, layer_atomic_density
from pyrump.atomic.tables import PeriodicTable
from pyrump.io.kalbitzer import parse_kalbitzer
from pyrump.model.geometry import Geometry
from pyrump.profiles.equations import (
    INTEGRAL_FORMS,
    EquationType,
    Profile,
    mix_composition,
    recommended_sublayers,
    species_fraction,
)
from pyrump.sim.engine import Beam, UniformSample, simulate_bricks
from pyrump.sim.slabs import sublayer_count
from pyrump.stopping.kalbitzer import KalbitzerStopping
from pyrump.stopping.registry import StoppingRegistry
from pyrump.stopping.ziegler import ZieglerStopping

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "oracle"))
import oracle as ora  # noqa: E402


from conftest import data_dir

DATA = data_dir()
HEIGHT_RTOL = 1e-4

#: (oracle name, type, parameters) for every supported form.
FORMS = [
    ("constant", EquationType.CONSTANT, (0.1,)),
    ("linear", EquationType.LINEAR, (0.0, 0.2)),
    ("erfc", EquationType.ERFC, (0.3, 1e-14, 3600.0)),
    ("exponential", EquationType.EXPONENTIAL, (0.3, 1e-14, 1e-8)),
    ("semi-infinite", EquationType.SEMI_INFINITE, (0.4, 1e-14, 3600.0, 500.0)),
    ("thickfilm", EquationType.THICKFILM, (0.1, 0.3, 1e-14, 3600.0, 500.0)),
    ("timedependent", EquationType.TIMEDEPENDENT, (0.2, 1e-14, 3600.0, 1e-8)),
    ("gaussian", EquationType.GAUSSIAN, (100.0, 300.0, 200.0)),
    ("thinfilm", EquationType.THINFILM, (100.0, 1e-14, 3600.0)),
    ("buriedthinfilm", EquationType.BURIEDTHINFILM, (100.0, 1e-14, 3600.0, 400.0)),
    ("edgeworth", EquationType.EDGEWORTH, (100.0, 300.0, 150.0, 0.3, 0.2)),
]


# ------------------------------------------------------------ mixing rule


def test_mixing_normalises_both_compositions():
    """SiO2 written as (1,2) or (0.5,1) must give the same result."""
    a = mix_composition([0.5], matrix=np.array([1.0, 2.0]), species=np.array([0.0, 1.0]))
    b = mix_composition([0.5], matrix=np.array([0.5, 1.0]), species=np.array([0.0, 3.0]))
    assert np.allclose(a, b)


def test_mixing_endpoints():
    matrix = np.array([1.0, 0.0])
    species = np.array([0.0, 1.0])
    assert np.allclose(mix_composition([0.0], matrix, species)[0], [1.0, 0.0])
    assert np.allclose(mix_composition([1.0], matrix, species)[0], [0.0, 1.0])
    assert np.allclose(mix_composition([0.25], matrix, species)[0], [0.75, 0.25])


def test_mixing_clamps_above_one_but_allows_negative_then_clips():
    """The C caps the fraction at 1, lets negatives through, then clips density."""
    matrix = np.array([1.0, 0.0])
    species = np.array([0.0, 1.0])
    assert np.allclose(mix_composition([5.0], matrix, species)[0], [0.0, 1.0])
    # A negative fraction over-weights the matrix; the species side clips to 0.
    result = mix_composition([-1.0], matrix, species)[0]
    assert result[1] == 0.0
    assert result[0] == pytest.approx(2.0)


def test_mixing_rows_sum_to_one_for_valid_fractions():
    matrix = np.array([1.0, 2.0])
    species = np.array([3.0, 1.0])
    blended = mix_composition(np.linspace(0.0, 1.0, 11), matrix, species)
    assert np.allclose(blended.sum(axis=1), 1.0)


# --------------------------------------------------------- layer density


def test_inverse_density_average():
    """creatr.c averages cm^3/atom, not atoms/cm^3 -- 'hard ball packing'."""
    densities = np.array([4.9777e22, 5.9049e22])
    # Pure silicon.
    assert layer_atomic_density([1.0, 0.0], densities) == pytest.approx(0.49777, rel=1e-4)
    # A 50/50 mix sits between, but below the arithmetic mean.
    mixed = layer_atomic_density([1.0, 1.0], densities)
    arithmetic = 0.5 * (0.49777 + 0.59049)
    assert 0.49777 < mixed < 0.59049
    assert mixed < arithmetic


def test_density_falls_back_to_silicon():
    assert layer_atomic_density([0.0, 0.0], np.array([1e22, 1e22])) == SILICON_DENSITY


# ------------------------------------------------------- sublayer counts


def test_equation_overrides_maxpath():
    """Once an equation is attached, maxpth is ignored (creatr.c:700)."""
    # Without an equation a 1000e15 layer at maxpath 200 gives 6 sublayers.
    assert sublayer_count(1000.0, 200.0) == 6
    # With one, the equation's recommendation wins regardless of thickness.
    assert sublayer_count(1000.0, 200.0, equation_sublayers=20) == 20
    assert sublayer_count(50.0, 200.0, equation_sublayers=20) == 20
    # An explicit count still beats it.
    assert sublayer_count(1000.0, 200.0, explicit=7, equation_sublayers=20) == 7


@pytest.mark.parametrize(
    "equation, expected",
    [
        (EquationType.CONSTANT, 5),
        (EquationType.LINEAR, 10),
        (EquationType.ERFC, 10),
        (EquationType.GAUSSIAN, 20),
        (EquationType.THINFILM, 30),
        (EquationType.EDGEWORTH, 30),
    ],
)
def test_recommended_sublayers_match_eqlist(equation, expected):
    assert recommended_sublayers(equation) == expected


# ------------------------------------------------------ profile behaviour


def test_constant_is_flat():
    values = species_fraction(Profile(EquationType.CONSTANT, (0.3,)), 8, 1e-5)
    assert np.allclose(values, 0.3)


def test_linear_samples_sublayer_centres():
    """x = (i+0.5)/n, so a 0->1 ramp over 10 slabs starts at 0.05."""
    values = species_fraction(Profile(EquationType.LINEAR, (0.0, 1.0)), 10, 1e-5)
    assert values[0] == pytest.approx(0.05)
    assert values[-1] == pytest.approx(0.95)
    assert np.allclose(np.diff(values), 0.1)


def test_exponential_decays_with_depth():
    """f(x) = c0*exp(-x/(D/v)), sampled at sublayer centres.

    Note the surface value is *not* c0: the first centre already sits half a
    sublayer in, which here is half a decay length.
    """
    c0, diffusivity, velocity = 0.5, 1e-14, 1e-8
    thickness_cm, n = 2e-5, 20
    decay_length = diffusivity / velocity

    values = species_fraction(
        Profile(EquationType.EXPONENTIAL, (c0, diffusivity, velocity)), n, thickness_cm
    )
    centres = (np.arange(n) + 0.5) / n * thickness_cm
    assert np.allclose(values, c0 * np.exp(-centres / decay_length))
    assert np.all(np.diff(values) < 0)
    # log-linear in depth, which is the defining property.
    assert np.allclose(np.diff(np.log(values)), np.diff(np.log(values))[0])


def test_integral_forms_conserve_dose():
    """A Gaussian well inside the layer must deposit its full dose.

    Integral forms difference the cumulative normal at sublayer edges, so the
    total is exact however coarse the grid -- that is the point of them.
    """
    thickness = 1000.0  # 1e15 at/cm^2
    profile = Profile(EquationType.GAUSSIAN, (100.0, 1000.0, 400.0))
    for n in (20, 50, 200):
        values = species_fraction(
            profile, n, 2e-5, areal_thickness=thickness, species_density=0.59
        )
        deposited = (values * thickness / n).sum()
        assert deposited == pytest.approx(100.0, rel=0.02), f"n={n}"


def test_integral_forms_are_marked():
    assert EquationType.GAUSSIAN in INTEGRAL_FORMS
    assert EquationType.THINFILM in INTEGRAL_FORMS
    assert EquationType.LINEAR not in INTEGRAL_FORMS


def test_unsupported_forms_raise_rather_than_return_zero():
    for equation in (EquationType.SPLINE, EquationType.USEREQN):
        with pytest.raises(NotImplementedError, match="GENPLOT"):
            Profile(equation, (0.0,))


def test_none_gives_pure_matrix():
    values = species_fraction(Profile(EquationType.NONE), 5, 1e-5)
    assert np.allclose(values, 0.0)


# ----------------------------------------------------- against the C

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


@pytest.mark.parametrize("name, equation, params", FORMS, ids=[f[0] for f in FORMS])
def test_profile_matches_oracle(oracle, registry, table, name, equation, params):
    """Au diffusing into a Si matrix, one form at a time."""
    sample = UniformSample(
        [1000.0], [14, 79], [[1.0, 0.0]],
        profiles=[Profile(equation, params)], species=[[0.0, 1.0]],
    )
    mine = simulate_bricks(
        sample, Beam(), Geometry(theta=0.0, phi=10.0), registry, table
    )

    oracle.set_beam(e0_MeV=2.0, phi=10.0, theta=0.0, kevch=5.0, npt=1024)
    oracle.set_sample([1000.0], [14, 79], [[1.0, 0.0]])
    oracle.set_layer_equation(0, name, list(params), [0.0, 1.0])
    theirs = oracle.simulate_bricks()

    assert len(mine) == len(theirs), "sublayer count differs"
    assert np.allclose(mine.e_front, theirs[:, 2], rtol=1e-5)
    assert np.allclose(mine.e_back, theirs[:, 3], rtol=1e-5)
    assert np.allclose(mine.h_front, theirs[:, 4], rtol=HEIGHT_RTOL)
    assert np.allclose(mine.h_back, theirs[:, 5], rtol=HEIGHT_RTOL)


def test_sublayer_counts_match_the_c(oracle):
    for name, equation, _ in FORMS:
        assert oracle.equation_sublayers(name) == recommended_sublayers(equation), name


def test_profile_actually_varies_the_spectrum(registry, table):
    """A guard against the profile being silently ignored."""
    geometry = Geometry(theta=0.0, phi=10.0)
    flat = simulate_bricks(
        UniformSample(
            [1000.0], [14, 79], [[1.0, 0.0]],
            profiles=[Profile(EquationType.CONSTANT, (0.1,))], species=[[0.0, 1.0]],
        ),
        Beam(), geometry, registry, table,
    )
    graded = simulate_bricks(
        UniformSample(
            [1000.0], [14, 79], [[1.0, 0.0]],
            profiles=[Profile(EquationType.LINEAR, (0.0, 0.2))], species=[[0.0, 1.0]],
        ),
        Beam(), geometry, registry, table,
    )
    # Gold is the last block in each; a graded profile must not be flat.
    gold_flat = flat.h_front[-flat.__len__() // 4 :]
    gold_graded = graded.h_front[-graded.__len__() // 4 :]
    assert np.std(gold_graded) > np.std(gold_flat)
