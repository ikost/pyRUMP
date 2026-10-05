r"""The top-hat digital filter: background removal without a background model.

GUPIX's approach (Maxwell, Campbell & Teesdale, Nucl. Instr. Meth. B 43
(1989) 218; Maxwell & Campbell 2002), after Schamber (1977) and Statham.
Each channel is replaced by the mean of a central lobe of ``UW`` channels
minus the mean of two wings of ``LW`` channels on either side:

.. math::
    y'_i = \frac{1}{U}\sum_{\text{centre}} y_{i+j}
         - \frac{1}{2L}\sum_{\text{wings}} y_{i+j}

The kernel has zero area and is symmetric, so a constant or linearly sloping
continuum filters to exactly zero and any slowly varying one nearly so,
while peaks about one FWHM wide survive (with negative side lobes). Fitting
the filtered model to the filtered data therefore needs no model of the
continuum at all.

The lobe widths are given in units of the detector FWHM at each channel --
GUPIX's *variable* filter -- with Schamber and Statham's compromise of
``UW = 1 FWHM, LW = 0.5 FWHM`` as the default. A *two-region* filter uses
wider lobes above a split energy, where peaks are few, weak and isolated.

The variance of a filtered channel, for Poisson counts, is

.. math:: \sigma^2_i = \frac{1}{U^2}\sum_{\text{centre}} y_{i+j}
                     + \frac{1}{4L^2}\sum_{\text{wings}} y_{i+j}

with each channel's variance taken as its counts, but at least 1 -- the
usual floor, so that empty stretches still carry a weight.

Channels whose kernel would reach past either end of the spectrum are not
defined (NaN).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..model.spectrum import Calibration


@dataclass(frozen=True, slots=True)
class FilterSettings:
    """Lobe widths in units of the FWHM, and an optional wider upper region."""

    upper: float = 1.0
    """Central lobe width, in FWHM."""

    lower: float = 0.5
    """Width of each wing, in FWHM."""

    split_keV: float | None = None
    """Above this energy, ``split_upper``/``split_lower`` apply instead."""

    split_upper: float = 3.0
    split_lower: float = 0.5


def lobe_widths(
    calibration: Calibration, fwhm_keV: np.ndarray, settings: FilterSettings
) -> tuple[np.ndarray, np.ndarray]:
    """Per channel: half-width of the central lobe (it spans ``2u + 1``
    channels) and width of each wing, both in whole channels, at least 1."""
    energy = calibration.edge_energy(np.arange(calibration.npt)) + calibration.kevch / 2
    fwhm_channels = np.asarray(fwhm_keV, dtype=np.float64) / calibration.kevch
    upper = np.full(energy.size, settings.upper)
    lower = np.full(energy.size, settings.lower)
    if settings.split_keV is not None:
        above = energy > settings.split_keV
        upper[above] = settings.split_upper
        lower[above] = settings.split_lower
    half = np.maximum(np.rint((upper * fwhm_channels - 1) / 2), 0).astype(int)
    wing = np.maximum(np.rint(lower * fwhm_channels), 1).astype(int)
    return half, wing


def apply(counts, half: np.ndarray, wing: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Filter ``counts``; returns the filtered spectrum and its variance
    (from the counts themselves, as for measured data). Channels whose
    kernel runs off either end are NaN."""
    y = np.asarray(counts, dtype=np.float64)
    n = y.size
    cumulative = np.concatenate(([0.0], np.cumsum(y)))
    variance_source = np.concatenate(([0.0], np.cumsum(np.maximum(y, 1.0))))
    i = np.arange(n)
    left_outer = i - half - wing
    right_outer = i + half + wing
    valid = (left_outer >= 0) & (right_outer <= n - 1)
    filtered = np.full(n, np.nan)
    variance = np.full(n, np.nan)
    if not valid.any():
        return filtered, variance
    k = i[valid]
    h, w = half[valid], wing[valid]

    def window(sums, start, stop):  # sum over channels start..stop inclusive
        return sums[stop + 1] - sums[start]

    centre_width = 2 * h + 1
    centre = window(cumulative, k - h, k + h)
    wings = window(cumulative, k - h - w, k - h - 1) + window(cumulative, k + h + 1, k + h + w)
    filtered[valid] = centre / centre_width - wings / (2 * w)
    centre_var = window(variance_source, k - h, k + h)
    wings_var = (window(variance_source, k - h - w, k - h - 1)
                 + window(variance_source, k + h + 1, k + h + w))
    variance[valid] = centre_var / centre_width**2 + wings_var / (4.0 * w**2)
    return filtered, variance
