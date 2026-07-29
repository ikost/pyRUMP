"""M12 acceptance: the fitting objective, windows and Levenberg-Marquardt.

Two things are tested very differently, and deliberately.

The **objective** defines what is being minimised, so it must match RUMP
exactly — it is compared to ``EvalChiPoisson`` channel by channel via the
oracle.

The **optimiser** only finds the minimum. scipy's LM and Bevington's CURFIT
reach the same place by different routes, so comparing iterates would be
meaningless. It is tested by recovering known parameters from synthetic data
instead.
"""

from __future__ import annotations

import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

from pyrump.atomic.tables import PeriodicTable
from pyrump.fit.lm import fit
from pyrump.fit.objective import chi_square, chi_square_exact, poisson_residuals
from pyrump.fit.parameters import (
    FitInputs,
    composition,
    pack,
    parameter,
    thickness,
    unpack,
)
from pyrump.fit.windows import MAX_ERROR_WINDOWS, Window, WindowSet
from pyrump.io.kalbitzer import parse_kalbitzer
from pyrump.model.detector import Measurement
from pyrump.model.geometry import Geometry
from pyrump.model.spectrum import Calibration
from pyrump.sim.engine import Beam, UniformSample, simulate
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
GEOMETRY = Geometry(theta=0.0, phi=10.0)


# ------------------------------------------------------------- the objective


def test_perfect_agreement_gives_zero():
    counts = np.linspace(1.0, 5000.0, 200)
    result = chi_square(counts, counts.copy())
    assert result.total == pytest.approx(0.0, abs=1e-20)
    assert np.allclose(result.residuals, 0.0)


def test_residual_sign_follows_theory_minus_data():
    """RUMP's convention: theory above data is a *positive* residual."""
    data = np.array([100.0, 100.0])
    theory = np.array([150.0, 60.0])
    residuals, _ = poisson_residuals(data, theory)
    assert residuals[0] > 0
    assert residuals[1] < 0


def test_zero_data_channel():
    """With no counts observed the log term vanishes: chi = sqrt(2t)."""
    residuals, _ = poisson_residuals(np.array([0.0]), np.array([8.0]))
    assert residuals[0] == pytest.approx(np.sqrt(16.0))


def test_non_positive_theory_is_invalid_and_contributes_nothing():
    data = np.array([10.0, 10.0, 10.0])
    theory = np.array([0.0, -1.0, 10.0])
    residuals, n_invalid = poisson_residuals(data, theory)
    assert n_invalid == 2
    assert residuals[0] == 0.0 and residuals[1] == 0.0


def test_series_branch_matches_the_exact_form():
    """The 0.9 < x < 1.1 substitution must agree to the C's stated accuracy."""
    data = np.full(2000, 4000.0)
    theory = data * np.linspace(0.9, 1.1, 2000)
    approximate = chi_square(data, theory).total
    exact = chi_square_exact(data, theory)
    assert approximate == pytest.approx(exact, rel=1e-4)


def test_series_is_continuous_at_the_branch_points():
    """No step at x = 0.9 or x = 1.1, where the implementation switches."""
    for edge in (0.9, 1.1):
        data = np.full(3, 1000.0)
        theory = data * np.array([edge - 1e-6, edge, edge + 1e-6])
        residuals, _ = poisson_residuals(data, theory)
        assert abs(residuals[1] - residuals[0]) < 1e-3
        assert abs(residuals[2] - residuals[1]) < 1e-3


def test_reduced_chi_square_divides_by_dof():
    data = np.full(100, 500.0)
    theory = np.full(100, 520.0)
    result = chi_square(data, theory, n_parameters=4)
    assert result.dof == 96
    assert result.reduced == pytest.approx(result.total / 96)


def test_shape_mismatch_rejected():
    with pytest.raises(ValueError, match="shape mismatch"):
        poisson_residuals(np.zeros(5), np.zeros(6))


@pytest.mark.oracle
@pytest.mark.skipif(not ora.available(), reason="oracle unavailable")
@pytest.mark.parametrize(
    "label, seed, scale",
    [
        ("typical", 1, 1.02),
        ("theory high", 2, 1.30),
        ("theory low", 3, 0.70),
        ("near unity", 4, 1.001),
        ("low counts", 5, 1.10),
    ],
)
def test_objective_matches_the_c(label, seed, scale):
    """Channel-by-channel against ``EvalChiPoisson``."""
    oracle = ora.Oracle.load()
    rng = np.random.default_rng(seed)
    expectation = np.linspace(3.0, 3000.0, 500)
    data = rng.poisson(expectation).astype(float)
    theory = expectation * scale

    mine = chi_square(data, theory, n_parameters=3)
    residuals, reduced, n_invalid = oracle.chi_poisson(data, theory, n_parameters=3)

    assert np.allclose(mine.residuals, residuals, atol=1e-4)
    assert mine.reduced == pytest.approx(reduced, rel=1e-5)
    assert mine.n_invalid == n_invalid


@pytest.mark.oracle
@pytest.mark.skipif(not ora.available(), reason="oracle unavailable")
def test_windowed_objective_matches_the_c():
    oracle = ora.Oracle.load()
    rng = np.random.default_rng(11)
    expectation = np.linspace(5.0, 2000.0, 400)
    data = rng.poisson(expectation).astype(float)
    theory = expectation * 1.05

    mask = np.zeros(400, dtype=bool)
    mask[50:150] = True
    mask[220:300] = True

    mine = chi_square(data, theory, valid=mask, n_parameters=2)
    residuals, reduced, _ = oracle.chi_poisson(data, theory, valid=mask, n_parameters=2)
    assert np.allclose(mine.residuals, residuals, atol=1e-4)
    assert mine.reduced == pytest.approx(reduced, rel=1e-5)
    assert mine.n_used == 180


# ----------------------------------------------------------------- windows


def test_window_mask():
    assert Window(2, 4).mask(8).tolist() == [False, False, True, True, True, False, False, False]


def test_inverted_window_rejected():
    with pytest.raises(ValueError, match="inverted"):
        Window(10, 5)


def test_windows_combine_by_union():
    windows = WindowSet(error=[Window(0, 2), Window(5, 6)])
    assert windows.mask(8).tolist() == [True] * 3 + [False, False] + [True, True] + [False]


def test_no_windows_uses_everything():
    assert WindowSet().mask(5).all()


def test_window_limit_enforced():
    with pytest.raises(ValueError, match="at most 10"):
        WindowSet(error=[Window(i, i) for i in range(MAX_ERROR_WINDOWS + 1)])


def test_normalisation_window_equalises_totals():
    data = np.full(100, 100.0)
    theory = np.full(100, 120.0)
    windows = WindowSet(normalisation=Window(10, 49))
    scale = windows.normalisation_factor(data, theory)
    assert scale == pytest.approx(1.2)
    assert (data[10:50] * scale).sum() == pytest.approx(theory[10:50].sum())


def test_normalisation_window_and_free_correction_are_rejected():
    """Degenerate: both scale the same thing (pert.c:1163)."""
    windows = WindowSet(normalisation=Window(0, 10))
    with pytest.raises(ValueError, match="degenerate"):
        windows.validate_against(varying_correction=True)
    windows.validate_against(varying_correction=False)  # fine on its own


# -------------------------------------------------------------- parameters


def _inputs(thick=1000.0, fwhm=15.0):
    return FitInputs(
        sample=UniformSample([thick], [14], [[1.0]]),
        beam=Beam(),
        geometry=GEOMETRY,
        calibration=CAL,
        measurement=Measurement(omega_msr=1.0, charge_uC=10.0, tau_us=0.0, fwhm_keV=fwhm),
    )


def test_parameter_round_trip():
    inputs = _inputs()
    params = [thickness(0), parameter("fwhm"), parameter("mev")]
    values = pack(params, inputs)
    assert values.tolist() == [1000.0, 15.0, 2.0]

    unpack(params, inputs, np.array([1500.0, 20.0, 3.0]))
    assert inputs.sample.thicknesses[0] == 1500.0
    assert inputs.measurement.fwhm_keV == 20.0
    assert inputs.beam.e0_MeV == 3.0


def test_setting_a_frozen_field_replaces_the_object():
    """Measurement and Calibration are frozen, so setters rebuild them."""
    inputs = _inputs()
    original = inputs.measurement
    unpack([parameter("fwhm")], inputs, np.array([99.0]))
    assert inputs.measurement is not original
    assert original.fwhm_keV == 15.0  # untouched


def test_composition_parameter():
    inputs = FitInputs(
        sample=UniformSample([1000.0], [14, 8], [[1.0, 2.0]]),
        beam=Beam(), geometry=GEOMETRY, calibration=CAL,
        measurement=Measurement(),
    )
    param = composition(0, 1)
    assert param.get(inputs) == 2.0
    param.set(inputs, 1.5)
    assert inputs.sample.compositions[0][1] == 1.5


def test_unknown_parameter_name_is_helpful():
    with pytest.raises(KeyError, match="unknown fit parameter"):
        parameter("wobble")


def test_parameter_names_follow_rump():
    for name in ("mev", "theta", "phi", "psi", "fwhm", "tau", "current",
                 "correction", "kev/ch", "kev(0)", "straggle", "multiple_scatter"):
        assert parameter(name).name == name


# ------------------------------------------------------------------- the fit

pytestmark_fit = pytest.mark.skipif(DATA is None, reason="legacy tables unavailable")


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
def runner(registry, table):
    def run(inputs: FitInputs) -> np.ndarray:
        return simulate(
            inputs.sample, inputs.beam, inputs.geometry, registry, table,
            inputs.calibration, inputs.measurement,
        ).counts
    return run


def _synthetic(runner, thick=1200.0, seed=3):
    truth = _inputs(thick=thick)
    counts = runner(truth)
    rng = np.random.default_rng(seed)
    return rng.poisson(np.maximum(counts, 0.0)).astype(float)


@pytestmark_fit
def test_recovers_a_known_thickness(runner):
    data = _synthetic(runner)
    inputs = _inputs(thick=900.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        result = fit(
            runner, data, inputs, [thickness(0)],
            windows=WindowSet(error=[Window(205, 240)]),
        )
    assert result.parameters["thickness[0]"] == pytest.approx(1200.0, rel=0.02)
    assert result.success


@pytestmark_fit
def test_recovers_two_parameters(runner):
    data = _synthetic(runner)
    inputs = _inputs(thick=900.0, fwhm=25.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        result = fit(
            runner, data, inputs, [thickness(0), parameter("fwhm")],
            windows=WindowSet(error=[Window(200, 260)]),
        )
    assert result.parameters["thickness[0]"] == pytest.approx(1200.0, rel=0.03)
    assert result.parameters["fwhm"] == pytest.approx(15.0, rel=0.15)
    assert result.correlation is not None
    assert result.correlation.shape == (2, 2)
    assert result.correlation[0, 0] == pytest.approx(1.0, abs=1e-9)


@pytestmark_fit
def test_reports_uncertainties(runner):
    data = _synthetic(runner)
    inputs = _inputs(thick=1100.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        result = fit(
            runner, data, inputs, [thickness(0)],
            windows=WindowSet(error=[Window(205, 240)]),
        )
    assert result.uncertainties["thickness[0]"] > 0
    assert result.dof > 0


@pytestmark_fit
def test_warns_when_the_window_is_mostly_empty(runner):
    """The manual's warning, surfaced.

    Channels where the model predicts zero contribute neither to chi-square nor
    to the gradient, so a window reaching well past the spectrum silently
    discards most of its own evidence.
    """
    data = _synthetic(runner)
    inputs = _inputs(thick=900.0)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = fit(
            runner, data, inputs, [thickness(0)],
            windows=WindowSet(error=[Window(200, 420)]),
        )
    assert result.n_invalid > 100
    assert any("non-positive predicted counts" in str(w.message) for w in caught)


@pytestmark_fit
def test_empty_window_is_rejected(runner):
    data = _synthetic(runner)
    inputs = _inputs()
    with pytest.raises(ValueError, match="select no channels"):
        fit(runner, data, inputs, [thickness(0)],
            windows=WindowSet(error=[Window(2000, 2010)]))


@pytestmark_fit
def test_normalisation_window_absorbs_a_charge_error(runner):
    """A 20% charge-integration error should not bias the thickness."""
    data = _synthetic(runner) * 1.2
    inputs = _inputs(thick=1000.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        result = fit(
            runner, data, inputs, [thickness(0)],
            windows=WindowSet(
                error=[Window(205, 240)], normalisation=Window(205, 240)
            ),
        )
    assert result.normalisation == pytest.approx(1 / 1.2, rel=0.1)
    assert result.parameters["thickness[0]"] == pytest.approx(1200.0, rel=0.05)
