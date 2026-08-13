"""M15 acceptance: BACKGROUND's weighted polynomial fit and subtraction."""

from __future__ import annotations

import numpy as np
import pytest

from pyrump.analysis.background import fit_background


def test_fit_background_recovers_a_noiseless_linear_background():
    channels = np.arange(100, dtype=np.float64)
    background = 50.0 + 2.0 * channels
    counts = background.copy()
    counts[40:60] += 1000.0  # a "peak" excluded from the flanking fit regions

    fit = fit_background(counts, 10, 30, 70, 90, order=1)

    np.testing.assert_allclose(fit.fit, background[10:91], atol=1.0)
    stripped_at_peak = fit.stripped[40 - 10 : 60 - 10]
    np.testing.assert_allclose(stripped_at_peak, 1000.0, atol=5.0)
    # the flanking regions themselves should strip down close to zero
    stripped_flanks = np.concatenate([fit.stripped[: 30 - 10 + 1], fit.stripped[70 - 10 :]])
    np.testing.assert_allclose(stripped_flanks, 0.0, atol=5.0)


def test_fit_background_rejects_out_of_order_regions():
    counts = np.zeros(50)
    with pytest.raises(ValueError, match="i0 < i1"):
        fit_background(counts, 20, 10, 30, 40, order=2)


def test_fit_background_rejects_order_outside_1_to_8():
    counts = np.zeros(50)
    with pytest.raises(ValueError, match="between 1 and 8"):
        fit_background(counts, 0, 5, 10, 15, order=9)
    with pytest.raises(ValueError, match="between 1 and 8"):
        fit_background(counts, 0, 5, 10, 15, order=0)


def test_fit_background_channels_span_the_whole_region_including_the_peak():
    counts = np.full(50, 10.0)
    fit = fit_background(counts, 5, 10, 30, 35, order=1)
    np.testing.assert_array_equal(fit.channels, np.arange(5, 36, dtype=np.float64))
