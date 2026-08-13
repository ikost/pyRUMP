"""SMOOTH's three algorithms: ``-sv``, ``-conv``, ``-fft``.

Ports ``RbsSmooth_SV``/``RbsSmooth_Conv``/``RbsSmooth_FFT`` (``anlytc.c``) and
``FFTSmooth`` (``genplot/fft_me.c``). ``FFT`` (the command) is ``SMOOTH -fft
-range`` with a forced range -- see ``smooth_fft`` below, used by both.

One deliberate departure from the C, not a faithful reproduction:
``RbsSmooth_Conv``'s inner ``conv()`` (anlytc.c:828) uses a fixed 25-slot ring
buffer to update its working array in place. That trick only reproduces its own
intended iterative update correctly for roughly ``2*iwidth+1 <= 25`` (i.e.
``sigma`` up to about 4 channels) -- well under the command's own allowed
``sigma`` range of 1-50. Nothing in the C comments this as deliberate; it reads
as an unexamined fixed-buffer limitation, not a documented behavior worth
bit-porting. :func:`smooth_conv` below implements the straightforward
(mathematically equivalent, for realistic FWHM/kevch ratios) sequential update
instead.
"""

from __future__ import annotations

import math

import numpy as np

#: RbsSmooth_SV's fixed 5-point kernel (Tracor manuals), anlytc.c:794.
SV_KERNEL = np.array([-3.0, 12.0, 17.0, 12.0, -3.0]) / 35.0


def smooth_sv(counts: np.ndarray, low: int, high: int) -> np.ndarray:
    """``RbsSmooth_SV`` (anlytc.c:776): fixed 5-point Savitzky-Golay smooth.

    Neighbor samples for channels near ``low``/``high`` are drawn from the
    *whole* buffer, not clamped to the selected sub-range -- only channels
    at the buffer's own true start/end are edge-extended (``anlytc.c``'s own
    ``max(0, ...)``/``min(..., npt-1)`` clamps).
    """
    counts = np.asarray(counts, dtype=np.float64)
    npt = counts.size
    if not (0 <= low <= high < npt):
        raise ValueError(f"region {low}-{high} outside 0-{npt - 1}")
    if high - low < 5:
        raise ValueError("SMOOTH -sv needs at least 6 channels (high-low >= 5)")

    index = np.arange(low, high + 1)
    accumulated = np.zeros(index.size, dtype=np.float64)
    for offset, weight in zip((-2, -1, 0, 1, 2), SV_KERNEL):
        accumulated += weight * counts[np.clip(index + offset, 0, npt - 1)]

    out = counts.copy()
    out[low : high + 1] = accumulated
    return out


def smooth_conv(
    counts: np.ndarray, low: int, high: int, fwhm_keV: float, kevch: float,
    n_iterations: int = 2,
) -> tuple[np.ndarray, list[float]]:
    """``RbsSmooth_Conv`` (anlytc.c:874): iterative Gaussian-residual smooth.

    Operates over the *whole* buffer (using the true buffer edges as the
    untouched boundary, exactly as the C's ``work``/``data`` arrays are sized
    to ``ibf->npt``, not to the selected region) and copies back only the
    requested ``[low, high]`` window -- so channels within ``iwidth`` of the
    buffer's own start/end are left at 0 if the selected window reaches them,
    a real (if surprising) consequence of the C's own array sizing, not
    something to special-case away.

    ``sigma`` uses RUMP's own convention, ``(FWHM/2)/sqrt(ln 2)/kevch`` --
    *not* the usual ``FWHM/(2*sqrt(2 ln 2))`` -- reproduced deliberately.

    Returns the modified full-length array and the per-iteration RMS history
    (for the caller to print, matching the C's live progress line).
    """
    full = np.asarray(counts, dtype=np.float64)
    npt = full.size
    if not (0 <= low <= high < npt):
        raise ValueError(f"region {low}-{high} outside 0-{npt - 1}")
    if npt < 20:
        raise ValueError("SMOOTH -conv needs at least 20 channels in the buffer")

    sigma = (fwhm_keV / 2.0) / math.sqrt(math.log(2.0)) / kevch
    if not (1.0 <= sigma <= 50.0):
        raise ValueError(
            f"SMOOTH -conv: characteristic width {sigma:.2f} channels is out of "
            "range 1-50 (check FWHM and keV/channel)"
        )
    iwidth = round(3.0 * sigma)
    if 2 * iwidth >= npt:
        raise ValueError("SMOOTH -conv: kernel too wide for this buffer")

    offsets = np.arange(-iwidth, iwidth + 1, dtype=np.float64)
    gauss = np.exp(-((offsets / sigma) ** 2))
    gauss /= gauss.sum()

    work = full.copy()
    data = np.zeros(npt, dtype=np.float64)
    valid = slice(iwidth, npt - iwidth)
    n_valid = npt - 2 * iwidth

    n_iterations = max(2, n_iterations)
    last_error = math.inf
    rms_history: list[float] = []
    for _ in range(n_iterations):
        residual = work - data
        # gauss is symmetric, so convolution and correlation coincide here --
        # this reproduces conv()'s sum_{j=-iwidth}^{iwidth} residual[i+j]*gauss[j+iwidth]
        # for every i in `valid`, since every tap i+j stays in [0, npt) there.
        blurred = np.convolve(residual, gauss, mode="same")
        data[valid] = data[valid] + blurred[valid]
        error = float(np.sqrt(np.sum((work[valid] - data[valid]) ** 2)) / n_valid)
        rms_history.append(error)
        if last_error <= error:
            break
        last_error = error

    out = full.copy()
    hi = min(high, npt - 1)
    out[low : hi + 1] = data[low : hi + 1]
    return out, rms_history


def _next_pow2(n: int) -> int:
    power = 1
    while power < n:
        power *= 2
    return power


def smooth_fft(counts: np.ndarray, low: int, high: int, width: float) -> np.ndarray:
    """``RbsSmooth_FFT``/``FFTSmooth`` (anlytc.c:952, genplot/fft_me.c:484).

    Self-contained to the selected ``[low, high]`` window (unlike ``-conv``):
    detrend -> zero-pad to the next power of 2 >= ``npt + 2*round(width)`` ->
    real FFT -> parabolic low-pass taper (cutoff set by ``width``, in
    channels) -> inverse FFT -> restore the linear trend.
    """
    full = np.asarray(counts, dtype=np.float64)
    npt_total = full.size
    if not (0 <= low <= high < npt_total) or high - low < 7:
        raise ValueError("SMOOTH -fft needs at least 8 channels in the selected range")

    y = full[low : high + 1].copy()
    n = y.size
    m = _next_pow2(n + 2 * round(width))

    a = y[0]
    b = (y[-1] - y[0]) / (n - 1)
    trend = a + b * np.arange(n)
    padded = np.zeros(m, dtype=np.float64)
    padded[:n] = y - trend

    spectrum = np.fft.rfft(padded)  # length m//2 + 1: bins 0..m//2
    k = np.arange(spectrum.size)
    factor = np.clip(1.0 - (width / m * k) ** 2, 0.0, None)
    factor[-1] = 0.0  # the Nyquist bin ("highest aliased point") is always zeroed
    spectrum = spectrum * factor

    smoothed = np.fft.irfft(spectrum, m)[:n] + trend

    out = full.copy()
    out[low : high + 1] = smoothed
    return out
