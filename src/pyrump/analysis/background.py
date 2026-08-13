"""BACKGROUND's weighted polynomial fit and subtraction.

Ports ``RbsBackground`` (anlytc.c:1011) using pyRUMP's own weighted
polynomial fitter (:func:`pyrump.stopping.polyfit.fit_polynomial`) rather than
reimplementing weighted least squares from scratch -- it already reproduces
``FitPolynomial``'s ``[-1,1]``/``[0,1]`` rescale-for-stability convention that
``RbsBackground`` itself needs (anlytc.c's own comment: "NEEDED TO NORMALIZE
THE X COORDINATES ... TO ELIMINATE INSTABILITY IN MATRIX INVERSION").
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..stopping.polyfit import fit_polynomial


@dataclass(frozen=True, slots=True)
class BackgroundFit:
    """The fitted background and the background-subtracted result."""

    channels: np.ndarray
    """Array indices ``i0..i3`` inclusive."""

    fit: np.ndarray
    """The fitted background curve, evaluated over :attr:`channels`."""

    stripped: np.ndarray
    """``counts[channels] - fit``."""

    order: int


def fit_background(
    counts: np.ndarray, i0: int, i1: int, i2: int, i3: int, order: int
) -> BackgroundFit:
    """``RbsBackground``'s fit (anlytc.c:1011-1080).

    Weighted least-squares polynomial fit (order 1-8) of the two flanking
    regions ``[i0,i1]`` and ``[i2,i3]`` -- weight ``1/sigma**2`` with
    ``sigma = sqrt(max(1,y))``, counting-statistics weighting -- evaluated
    across the whole ``[i0,i3]`` span, including the excluded "peak" region
    in between.
    """
    counts = np.asarray(counts, dtype=np.float64)
    if not (0 <= i0 < i1 < i2 < i3 < counts.size):
        raise ValueError(
            f"need 0 <= i0 < i1 < i2 < i3 < {counts.size}, got ({i0}, {i1}, {i2}, {i3})"
        )
    if not (1 <= order <= 8):
        raise ValueError("polynomial order must be between 1 and 8")

    flank_index = np.concatenate([np.arange(i0, i1 + 1), np.arange(i2, i3 + 1)])
    x = flank_index.astype(np.float64)
    y = counts[flank_index]
    sigma = np.sqrt(np.maximum(1.0, y))

    coefficients = fit_polynomial(x, y, sigma, order)  # ascending power order

    channels = np.arange(i0, i3 + 1, dtype=np.float64)
    fit = np.polyval(coefficients[::-1], channels)  # np.polyval wants descending
    stripped = counts[i0 : i3 + 1] - fit
    return BackgroundFit(channels=channels, fit=fit, stripped=stripped, order=order)
