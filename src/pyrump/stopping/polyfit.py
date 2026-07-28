"""Weighted polynomial least squares, matching GENPLOT's ``FitPolynomial``.

Port of ``FitPolynomial`` (genplot/gptfit.c:1558). This is the fitter behind the
stopping-power refit, so its conventions matter:

* the weight is ``1/sigma**2``; ``RbsGenStopp`` passes ``sigma = sqrt(S)``, so the
  effective weight is ``1/S`` -- a *relative* error criterion, which is why the fit
  stays accurate across the two orders of magnitude that S spans
* ``x`` is rescaled to ``[-1, 1]`` and ``y`` to ``[0, 1]`` before solving, then the
  coefficients are transformed back

The C then solves normal equations; here we use a least-squares solve on the
design matrix, which is mathematically identical but better conditioned. The
*coefficients* may therefore differ in their last digits while the fitted
*function* agrees far more tightly -- so always compare functions, not
coefficients.
"""

from __future__ import annotations

import numpy as np


def fit_polynomial(
    x: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray | None = None,
    order: int = 5,
) -> np.ndarray:
    """Fit ``y(x)`` with a polynomial of the given order.

    Returns coefficients in **ascending** power order, in the original ``x``
    coordinate, so that ``sum(c[i] * x**i)`` reproduces the fit.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    n_terms = order + 1

    if x.size < n_terms:
        raise ValueError(f"need at least {n_terms} points to fit order {order}")
    xmin, xmax = float(x.min()), float(x.max())
    if xmin == xmax:
        raise ValueError("all x values are identical; no fit possible")

    # Mirror the C's conditioning transforms (gptfit.c:1594-1600).
    if order != 1:
        ax = 2.0 / (xmax - xmin)
        bx = -1.0 - ax * xmin
        ymin, ymax = float(y.min()), float(y.max())
    else:
        ax, bx = 1.0, 0.0
        ymin, ymax = 0.0, 1.0
    ynorm = 1.0 if ymax == ymin else ymax - ymin

    x_scaled = x * ax + bx
    y_scaled = (y - ymin) / ynorm

    weight = np.ones_like(y) if sigma is None else 1.0 / np.asarray(sigma, np.float64) ** 2
    # Least squares minimises the squared residual, so apply sqrt(weight).
    root_weight = np.sqrt(weight)

    design = np.vander(x_scaled, n_terms, increasing=True)
    coefficients, *_ = np.linalg.lstsq(
        design * root_weight[:, None], y_scaled * root_weight, rcond=None
    )

    return _denormalize(coefficients, ax, bx, ymin, ynorm)


def _denormalize(
    coefficients: np.ndarray, ax: float, bx: float, ymin: float, ynorm: float
) -> np.ndarray:
    """Undo the x and y rescaling, returning coefficients in original units.

    The fit is ``y_scaled = p(ax*x + bx)``, so
    ``y = ynorm * p(ax*x + bx) + ymin``.
    """
    # numpy.polynomial uses ascending order, matching our convention.
    from numpy.polynomial import polynomial as P

    inner = np.array([bx, ax])  # ax*x + bx
    result = np.zeros(1)
    for power, coefficient in enumerate(coefficients):
        if coefficient == 0.0:
            continue
        term = np.array([1.0])
        for _ in range(power):
            term = P.polymul(term, inner)
        result = P.polyadd(result, coefficient * term)

    result = result * ynorm
    result[0] += ymin

    # Pad back to the requested length; polyadd may have trimmed.
    out = np.zeros(len(coefficients))
    out[: len(result)] = result[: len(out)]
    return out
