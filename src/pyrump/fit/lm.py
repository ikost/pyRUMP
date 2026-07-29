"""Levenberg-Marquardt fitting of a simulation to a measured spectrum.

Corresponds to PERT's default ``THOMPSON`` method, which is
``CurveFit`` in ``genplot/curfit.c`` — a Bevington CURFIT-lineage LM. pyRUMP
uses ``scipy.optimize.least_squares`` instead of transliterating it.

That substitution is safe in a way the objective's was not. The objective
*defines* what is being minimised, so it must match RUMP exactly and is tested
against the C channel by channel. The optimiser only *finds* the minimum: two
correct LM implementations reach the same place by different routes, so
comparing iterates would be meaningless. What is compared instead is the
chi-square surface — evaluate the objective at the C's converged parameters and
check the value agrees.

RUMP's numerics, for reference (pert.c:96-98): ``EpsCrit = 1e-3``,
``MaxIterations = 10``, numerical derivatives at a 1% parameter step. The
iteration cap is low because each evaluation is a full simulation.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

from .objective import chi_square
from .parameters import FitInputs, Parameter, bounds, pack, unpack
from .windows import WindowSet

#: RUMP's defaults (pert.c:96-98).
DEFAULT_MAX_ITERATIONS = 10
DEFAULT_EPS = 1e-3
DEFAULT_DERIVATIVE_STEP = 0.01


@dataclass(slots=True)
class FitResult:
    """Outcome of a fit."""

    parameters: dict[str, float]
    values: np.ndarray
    reduced_chi_square: float
    chi_square: float
    dof: int
    n_evaluations: int
    success: bool
    message: str
    normalisation: float = 1.0
    """Data scale factor chosen by the normalisation window, if any."""

    covariance: np.ndarray | None = None
    uncertainties: dict[str, float] = field(default_factory=dict)
    correlation: np.ndarray | None = None

    n_invalid: int = 0
    """Windowed channels where the model predicted zero counts.

    These contribute **nothing** to the objective — Poisson likelihood is
    undefined at zero expectation — so they also carry no gradient. A large
    count means the fit is blind over much of its window, which is what the
    manual warns about:

        Make sure that the theory spectra never goes to zero within the error
        window since a zero will wreak havoc with the Poisson statistics.

    The classic symptom is a parameter that will not move: if the only channels
    that discriminate between two models are ones where the *current* model
    predicts zero, the objective discards exactly the evidence needed.
    """


def _covariance(jacobian: np.ndarray, reduced: float) -> np.ndarray | None:
    """Parameter covariance from the Jacobian at the solution.

    For a chi-square objective already scaled by its own errors, the covariance
    is ``(J^T J)^-1``. Returns None if the Jacobian is rank-deficient, which
    means at least one parameter is unconstrained by the data.
    """
    try:
        hessian = jacobian.T @ jacobian
        return np.linalg.inv(hessian)
    except np.linalg.LinAlgError:
        return None


def fit(
    simulate_fn,
    data: np.ndarray,
    inputs: FitInputs,
    parameters: list[Parameter],
    *,
    windows: WindowSet | None = None,
    max_iterations: int = DEFAULT_MAX_ITERATIONS,
    derivative_step: float = DEFAULT_DERIVATIVE_STEP,
    eps: float = DEFAULT_EPS,
) -> FitResult:
    """Fit ``parameters`` so the simulation matches ``data``.

    ``simulate_fn(inputs) -> np.ndarray`` runs the forward model and returns
    channel counts. ``inputs`` is mutated during the fit and left holding the
    best-fit values.

    Raises if a normalisation window is combined with a free ``correction``
    parameter, which is degenerate — RUMP rejects it too.
    """
    windows = windows or WindowSet()
    data = np.asarray(data, dtype=np.float64)

    varying_correction = any(p.name == "correction" for p in parameters)
    windows.validate_against(varying_correction)

    mask = windows.mask(data.size)
    if not np.any(mask):
        raise ValueError("the error windows select no channels")

    evaluations = 0
    last_normalisation = 1.0

    def residuals(values: np.ndarray) -> np.ndarray:
        nonlocal evaluations, last_normalisation
        evaluations += 1
        unpack(parameters, inputs, values)
        theory = np.asarray(simulate_fn(inputs), dtype=np.float64)

        n = min(theory.size, data.size)
        observed, expected = data[:n], theory[:n]
        window = mask[:n]

        # The normalisation window removes charge-integration error by scaling
        # the *data*, before the objective sees it.
        last_normalisation = windows.normalisation_factor(observed, expected)
        result = chi_square(
            observed * last_normalisation,
            expected,
            valid=window,
            n_parameters=len(parameters),
        )
        return result.residuals[window]

    start = pack(parameters, inputs)
    lower, upper = bounds(parameters)
    # scipy needs the start strictly inside the bounds.
    start = np.clip(start, lower + 1e-12, upper - 1e-12)

    outcome = least_squares(
        residuals,
        start,
        bounds=(lower, upper),
        # Relative step, matching RUMP's 1% numerical derivative.
        diff_step=derivative_step,
        xtol=eps,
        ftol=eps,
        max_nfev=max_iterations * (len(parameters) + 1) * 4,
    )

    unpack(parameters, inputs, outcome.x)
    final = np.asarray(simulate_fn(inputs), dtype=np.float64)
    n = min(final.size, data.size)
    scale = windows.normalisation_factor(data[:n], final[:n])
    summary = chi_square(
        data[:n] * scale, final[:n], valid=mask[:n], n_parameters=len(parameters)
    )

    covariance = _covariance(outcome.jac, summary.reduced)
    uncertainties: dict[str, float] = {}
    correlation = None
    if covariance is not None:
        sigma = np.sqrt(np.abs(np.diag(covariance)))
        uncertainties = {p.name: float(s) for p, s in zip(parameters, sigma)}
        outer = np.outer(sigma, sigma)
        with np.errstate(divide="ignore", invalid="ignore"):
            correlation = np.where(outer > 0, covariance / outer, 0.0)

    if summary.n_invalid:
        warnings.warn(
            f"{summary.n_invalid} of {summary.n_used} windowed channels have "
            "non-positive predicted counts; Poisson statistics are undefined "
            "there, so they contribute neither to chi-square nor to the "
            "gradient. Narrow the error window to where the model has counts.",
            RuntimeWarning,
            stacklevel=2,
        )

    return FitResult(
        parameters={p.name: float(v) for p, v in zip(parameters, outcome.x)},
        values=outcome.x,
        reduced_chi_square=summary.reduced,
        chi_square=summary.total,
        dof=summary.dof,
        n_evaluations=evaluations,
        success=bool(outcome.success),
        message=str(outcome.message),
        normalisation=scale,
        covariance=covariance,
        uncertainties=uncertainties,
        correlation=correlation,
        n_invalid=summary.n_invalid,
    )
