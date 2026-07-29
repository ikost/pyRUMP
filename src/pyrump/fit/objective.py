r"""The fitting objective: Poisson maximum likelihood.

Port of ``EvalChiPoisson`` (genplot/curfit.c:557-614).

RBS spectra are counting measurements, so channel *i* holding :math:`n_i` counts
against an expectation :math:`t_i` is Poisson-distributed, not Gaussian. Using
ordinary least squares systematically mis-weights low-count channels — exactly
the ones that carry depth information in a tail. RUMP therefore uses the Poisson
maximum-likelihood ratio of Baker & Cousins, *Nucl. Instr. Meth.* **221** (1984)
437:

.. math::
    \chi^2 = 2\sum_i \left[ t_i - n_i + n_i \ln\frac{n_i}{t_i} \right]

which is asymptotically :math:`\chi^2`-distributed and reduces to the usual sum
of squares in the high-count limit.

Signed residuals
----------------

Levenberg-Marquardt wants a residual vector, not just a scalar, so the C splits
the sum into per-channel terms :math:`\chi_i` with
:math:`\chi^2 = \sum_i \chi_i^2` and a sign carrying the direction:

.. math::
    \chi_i = \operatorname{sign}(t_i - n_i)\,
             \sqrt{2 n_i \left(\tfrac{t_i}{n_i} - 1 - \ln\tfrac{t_i}{n_i}\right)}

Note the convention: **theory above data gives a positive residual**, the
opposite of the usual ``data - model``. Harmless as long as it is consistent,
but it matters if you compare residual vectors against another code.

Numerical stability
-------------------

Near :math:`t_i \approx n_i` the bracket is the difference of two nearly-equal
quantities and loses most of its significant digits. For
:math:`0.9 < t/n < 1.1` the C substitutes a series matched at the endpoints:

.. code-block:: c

    s = x - 1;
    chi = sqrt(data) * ((s*0.19547835 - 0.33469349)*s + 1.0) * s;

with a stated maximum error of 3.6e-6 (at x = 0.93 and 1.07). That is worth
keeping: it costs nothing and it removes a genuine cancellation problem, one the
leading coefficients confirm is the real expansion (the cubic term is -1/3).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

#: Bounds of the series-expansion branch (curfit.c:596-599).
_SERIES_LOW = 0.9
_SERIES_HIGH = 1.1

#: Series coefficients, verbatim (curfit.c:603).
_SERIES_A = 0.19547835
_SERIES_B = -0.33469349


@dataclass(slots=True)
class ChiSquare:
    """Outcome of one objective evaluation."""

    residuals: np.ndarray
    """Signed per-channel residuals; ``sum(residuals**2)`` is the total."""

    total: float
    """Unreduced chi-square."""

    reduced: float
    """``total / dof`` -- what RUMP reports as ``chisqr$``."""

    dof: int
    n_used: int
    """Channels that contributed, after windows and validity."""

    n_invalid: int
    """Channels with non-positive theory, where Poisson statistics do not apply."""


def poisson_residuals(data: np.ndarray, theory: np.ndarray) -> tuple[np.ndarray, int]:
    """Signed per-channel residuals, with the count of invalid channels.

    Reproduces ``EvalChiPoisson``'s three branches exactly:

    * ``theory <= 0`` — Poisson statistics are undefined, so the channel
      contributes nothing and is counted as invalid. The C prints
      "Poisson statistics invalid!" and returns the count.
    * ``data == 0`` — the log term vanishes and the residual is
      :math:`\\sqrt{2t}`.
    * otherwise — the exact form, or the series near :math:`t \\approx n`.
    """
    data = np.asarray(data, dtype=np.float64)
    theory = np.asarray(theory, dtype=np.float64)
    if data.shape != theory.shape:
        raise ValueError(f"shape mismatch: data {data.shape} vs theory {theory.shape}")

    residuals = np.zeros_like(data)
    invalid = theory <= 0.0
    n_invalid = int(np.count_nonzero(invalid))

    zero_data = (~invalid) & (data == 0.0)
    residuals[zero_data] = np.sqrt(2.0 * theory[zero_data])

    general = (~invalid) & (data != 0.0)
    if np.any(general):
        n = data[general]
        x = theory[general] / n
        root_n = np.sqrt(n)

        # Exact branch, evaluated where it is well conditioned.
        outer = (x > _SERIES_HIGH) | (x < _SERIES_LOW)
        safe = np.where(outer, x, 1.0)  # keep log() away from the series region
        magnitude = np.sqrt(np.maximum(2.0 * n * (safe - 1.0 - np.log(safe)), 0.0))
        exact = np.where(x > _SERIES_HIGH, magnitude, -magnitude)

        # Series branch near x == 1.
        s = x - 1.0
        series = root_n * ((s * _SERIES_A + _SERIES_B) * s + 1.0) * s

        residuals[general] = np.where(outer, exact, series)

    return residuals, n_invalid


def chi_square(
    data: np.ndarray,
    theory: np.ndarray,
    *,
    valid: np.ndarray | None = None,
    n_parameters: int = 0,
) -> ChiSquare:
    """Evaluate the Poisson objective over the selected channels.

    ``valid`` is the window mask; channels outside it are skipped entirely,
    which is how RUMP's error windows work (``nls->valid``, curfit.c:581).

    ``dof`` is ``n_used - n_parameters``, floored at 1 so the reduced value
    stays finite for an exactly-determined fit.
    """
    residuals, _ = poisson_residuals(data, theory)
    invalid = np.asarray(theory) <= 0.0

    if valid is not None:
        mask = np.asarray(valid, dtype=bool)
        if mask.shape != residuals.shape:
            raise ValueError("valid mask shape does not match the spectrum")
        residuals = np.where(mask, residuals, 0.0)
        n_used = int(np.count_nonzero(mask))
        # Count invalid channels only inside the window -- channels the fit
        # never looks at are not a problem worth reporting.
        invalid = invalid & mask
    else:
        n_used = int(residuals.size)
    n_invalid = int(np.count_nonzero(invalid))

    total = float(np.sum(residuals**2))
    dof = max(n_used - n_parameters, 1)
    return ChiSquare(
        residuals=residuals,
        total=total,
        reduced=total / dof,
        dof=dof,
        n_used=n_used,
        n_invalid=n_invalid,
    )


def chi_square_exact(data: np.ndarray, theory: np.ndarray) -> float:
    """The objective evaluated directly, with no series substitution.

    Reference implementation for testing: it says what the series branch is
    approximating, and is what the Baker-Cousins paper actually writes.
    """
    data = np.asarray(data, dtype=np.float64)
    theory = np.asarray(theory, dtype=np.float64)

    positive = theory > 0.0
    t = theory[positive]
    n = data[positive]

    # n*ln(n/t) -> 0 as n -> 0, which numpy would otherwise report as nan.
    log_term = np.zeros_like(n)
    counted = n > 0.0
    log_term[counted] = n[counted] * np.log(n[counted] / t[counted])
    return float(2.0 * np.sum(t - n + log_term))
