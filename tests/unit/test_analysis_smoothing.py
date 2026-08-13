"""M15 acceptance: SMOOTH's three algorithms (-sv, -conv, -fft)."""

from __future__ import annotations

import numpy as np
import pytest

from pyrump.analysis.smoothing import smooth_conv, smooth_fft, smooth_sv

# --------------------------------------------------------------------- -sv


def test_smooth_sv_matches_the_hand_computed_kernel():
    """A single spike, kernel [-3,12,17,12,-3]/35 applied around it."""
    counts = np.zeros(11)
    counts[5] = 10.0
    result = smooth_sv(counts, 2, 8)

    expected = {
        2: 0.0,
        3: -3 * 10 / 35,
        4: 12 * 10 / 35,
        5: 17 * 10 / 35,
        6: 12 * 10 / 35,
        7: -3 * 10 / 35,
        8: 0.0,
    }
    for index, value in expected.items():
        assert result[index] == pytest.approx(value)
    np.testing.assert_allclose(result[:2], counts[:2])
    np.testing.assert_allclose(result[9:], counts[9:])


def test_smooth_sv_reads_context_outside_the_selected_range():
    """Neighbor samples come from the whole buffer, not the selected window."""
    counts = np.zeros(11)
    counts[2] = 100.0  # just outside the selected [3, 9] range
    result = smooth_sv(counts, 3, 9)
    assert result[3] == pytest.approx(12 * 100 / 35)
    assert counts[2] == 100.0  # the input is not mutated


def test_smooth_sv_clamps_at_the_true_buffer_start():
    counts = np.zeros(11)
    counts[0] = 5.0
    result = smooth_sv(counts, 0, 6)
    # window at i=0 clamps offsets -2 and -1 to index 0: [5,5,5,0,0]
    expected = (-3 * 5 + 12 * 5 + 17 * 5 + 12 * 0 - 3 * 0) / 35
    assert result[0] == pytest.approx(expected)


def test_smooth_sv_needs_at_least_six_channels():
    counts = np.zeros(10)
    with pytest.raises(ValueError, match="at least 6"):
        smooth_sv(counts, 3, 6)  # high-low == 3


def test_smooth_sv_rejects_an_out_of_bounds_region():
    counts = np.zeros(10)
    with pytest.raises(ValueError, match="outside"):
        smooth_sv(counts, 0, 20)


# ------------------------------------------------------------------- -conv


def test_smooth_conv_reproduces_a_flat_plateau_after_the_first_pass():
    """Blurring a constant with an area-normalized kernel is a no-op.

    The *second* forced iteration (RbsSmooth_Conv always runs at least two,
    anlytc.c:881) actually makes things worse near the edges of the selected
    window: ``data`` outside the window stays permanently at 0, so that
    stale residual keeps leaking back in on every later pass -- a real
    property of the algorithm, not a bug, and exactly why the C's own
    "stop early if error stops improving" guard exists.
    """
    counts = np.full(200, 500.0)
    new_counts, rms_history = smooth_conv(counts, 50, 150, fwhm_keV=15.0, kevch=5.0)

    sigma = (15.0 / 2.0) / np.sqrt(np.log(2.0)) / 5.0
    iwidth = round(3 * sigma)
    # Well inside both the buffer's true edges and the selected window's own
    # boundary, where the second pass's edge leakage has fully decayed.
    lo, hi = 50 + 3 * iwidth, 150 - 3 * iwidth
    assert lo < hi
    np.testing.assert_allclose(new_counts[lo:hi], 500.0, atol=1e-6)
    assert rms_history[0] == pytest.approx(0.0, abs=1e-9)
    assert len(rms_history) == 2  # the second pass is worse, so it stops there


def test_smooth_conv_leaves_the_true_buffer_edges_untouched():
    """Channels within iwidth of the buffer's own start/end never get a value."""
    counts = np.full(200, 500.0)
    new_counts, _ = smooth_conv(counts, 0, 199, fwhm_keV=15.0, kevch=5.0)
    assert new_counts[0] == 0.0
    assert new_counts[-1] == 0.0


def test_smooth_conv_rejects_sigma_outside_1_to_50():
    counts = np.full(100, 10.0)
    with pytest.raises(ValueError, match="1-50"):
        smooth_conv(counts, 0, 99, fwhm_keV=0.01, kevch=5.0)


def test_smooth_conv_rejects_a_short_buffer():
    counts = np.full(10, 10.0)
    with pytest.raises(ValueError, match="at least 20"):
        smooth_conv(counts, 0, 9, fwhm_keV=15.0, kevch=5.0)


def test_smooth_conv_enforces_a_minimum_of_two_iterations():
    counts = np.full(200, 500.0)
    _, rms_history = smooth_conv(counts, 50, 150, fwhm_keV=15.0, kevch=5.0, n_iterations=1)
    assert len(rms_history) >= 2


# -------------------------------------------------------------------- -fft


def test_smooth_fft_preserves_a_pure_line():
    """A perfectly linear segment has no high-frequency content to filter."""
    n = 64
    y = 10.0 + 2.0 * np.arange(n)
    counts = np.zeros(100)
    counts[10 : 10 + n] = y
    result = smooth_fft(counts, 10, 10 + n - 1, width=5.0)
    np.testing.assert_allclose(result[10 : 10 + n], y, atol=1e-6)


def test_smooth_fft_needs_at_least_eight_channels():
    counts = np.zeros(20)
    with pytest.raises(ValueError, match="at least 8"):
        smooth_fft(counts, 0, 5, width=2.0)
