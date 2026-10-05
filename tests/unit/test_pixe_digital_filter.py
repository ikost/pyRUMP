"""The top-hat digital filter: what it removes, what it keeps, its variance
and how its lobes follow the FWHM."""

from __future__ import annotations

import numpy as np
import pytest

from pyrump.model.spectrum import Calibration
from pyrump.pixe.digital_filter import FilterSettings, apply, lobe_widths


def _constant(n, half, wing):
    return np.full(n, half, dtype=int), np.full(n, wing, dtype=int)


def test_constant_and_linear_backgrounds_filter_to_zero():
    half, wing = _constant(200, 4, 3)
    i = np.arange(200, dtype=float)
    filtered, _ = apply(37.0 + 2.5 * i, half, wing)
    valid = np.isfinite(filtered)
    assert valid.sum() == 200 - 2 * (4 + 3)
    np.testing.assert_allclose(filtered[valid], 0.0, atol=1e-9)


def test_a_peak_survives_and_the_kernel_has_zero_area():
    half, wing = _constant(400, 4, 3)
    i = np.arange(400, dtype=float)
    peak = 1000 * np.exp(-0.5 * ((i - 200) / 4.0) ** 2)
    filtered, _ = apply(peak + 50.0, half, wing)
    assert np.nanargmax(filtered) == 200
    assert filtered[200] > 0 and np.nanmin(filtered) < 0  # negative side lobes
    assert np.nansum(filtered) == pytest.approx(0.0, abs=1e-6)  # zero-area kernel


def test_variance_of_a_flat_spectrum():
    """c counts per channel: c/(2h+1) from the centre, c/(2w) from the wings."""
    half, wing = _constant(100, 3, 2)
    _, variance = apply(np.full(100, 80.0), half, wing)
    expected = 80.0 / 7 + 80.0 / 4
    np.testing.assert_allclose(variance[np.isfinite(variance)], expected)
    _, empty = apply(np.zeros(100), half, wing)  # empty channels count as 1
    np.testing.assert_allclose(empty[np.isfinite(empty)], 1 / 7 + 1 / 4)


def test_lobes_follow_the_fwhm_and_the_split():
    calibration = Calibration(kevch=0.01, kev0=0.0, first=0.0, npt=1000)
    fwhm_keV = np.full(1000, 0.09)  # 9 channels
    half, wing = lobe_widths(calibration, fwhm_keV, FilterSettings())
    assert (2 * half[0] + 1, wing[0]) == (9, 4)  # UW = 1 FWHM, LW = 0.5 FWHM
    split = FilterSettings(split_keV=5.0, split_upper=3.0, split_lower=0.5)
    half, wing = lobe_widths(calibration, fwhm_keV, split)
    assert 2 * half[100] + 1 == 9 and 2 * half[900] + 1 == 27
