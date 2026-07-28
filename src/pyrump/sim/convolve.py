"""Detector-resolution convolution.

Port of ``ConvoluteDetector`` (creatr.c:1178-1298).

Applied **once to the finished spectrum**, not per slab. Since 7/9/94 the
detector width was deliberately removed from the per-brick straggling path
(anlyz.c:195-201), so straggling broadens each brick individually while detector
resolution smears the assembled result. The two are separate stages with
different mathematics, and a port that folds them together will not match.

The kernel is **channel-integrated**, not point-sampled: weight *j* is the
Gaussian's integral across channel *j*, computed from differences of the normal
CDF. That matters at coarse calibrations where a channel spans a sizeable
fraction of sigma.

.. warning::
   **The convolution is not count-conserving at either edge.** Contributions
   that would land outside the channel range are discarded: below channel 0 by
   the head loop (creatr.c:1233-1236) and above the last channel by the tail
   loop (creatr.c:1284-1291). This is defensible -- a real MCA cannot record
   counts in channels it does not have -- but it means the integral of a
   convolved spectrum is slightly less than the original whenever intensity sits
   within ~3 sigma of either end. Reproduced by default; ``mode="renormalize"``
   rescales the kernel per output channel to conserve counts instead.
"""

from __future__ import annotations

import numpy as np
from scipy.special import ndtr

from ..model.spectrum import Calibration, Spectrum

#: Maximum kernel half-width in channels (``LMAX``, creatr.c:1177).
LMAX = 100

#: FWHM -> sigma. The C's literal, not 2*sqrt(2*ln2) = 2.35482.
FWHM_TO_SIGMA = 2.355


def gaussian_kernel(fwhm_keV: float, kevch: float) -> np.ndarray:
    """Half-kernel of channel-integrated Gaussian weights (creatr.c:1197-1207).

    Element 0 is the central channel's weight; element *j* is the **one-sided**
    weight for channels *i ± j*, so the full kernel sums to
    ``gauss[0] + 2*sum(gauss[1:])``.

    Truncated at 3 sigma, and hard-clamped to :data:`LMAX` channels -- at fine
    calibrations the clamp bites and the C only warns.
    """
    if fwhm_keV <= 0.0:
        return np.array([1.0])

    sigma = fwhm_keV / FWHM_TO_SIGMA
    half_width = int(3.0 * sigma / kevch + 2)  # C truncation
    half_width = min(max(half_width, 1), LMAX)

    # The C works in units where the CDF argument is (j+0.5)*kevch/sigma; it
    # reaches that via scale = kevch/sigma*0.707107 and a 1.41421356 factor.
    edges = (np.arange(half_width) + 0.5) * kevch / sigma
    upper = ndtr(-edges)

    kernel = np.empty(half_width, dtype=np.float64)
    kernel[0] = 2.0 * (0.5 - upper[0])
    if half_width > 1:
        kernel[1:] = upper[:-1] - upper[1:]
    return kernel


def full_kernel(half: np.ndarray) -> np.ndarray:
    """Mirror a half-kernel into a symmetric one suitable for convolution."""
    if half.size == 1:
        return half.copy()
    return np.concatenate([half[:0:-1], half])


def convolve_detector(
    counts: np.ndarray,
    calibration: Calibration,
    fwhm_keV: float,
    *,
    mode: str = "rump",
) -> np.ndarray:
    """Apply detector resolution.

    ``mode="rump"`` reproduces the C exactly, dropping contributions that fall
    off either end. ``mode="renormalize"`` divides each output channel by the
    kernel weight actually available to it, conserving counts.
    """
    if fwhm_keV <= 0.0:
        return counts.copy()

    half = gaussian_kernel(fwhm_keV, calibration.kevch)
    kernel = full_kernel(half)

    if mode == "rump":
        # 'same' keeps the output aligned with the input and silently drops the
        # overhang at both ends -- precisely the C's behaviour.
        return np.convolve(counts, kernel, mode="same")

    if mode == "renormalize":
        # Conservation is a property of each *source* channel: scale channel j's
        # contribution by the kernel weight that actually lands inside the array
        # from j. Normalising the output instead would rescale received weight,
        # which does not conserve -- it inflates channels near the edges.
        #
        # This also compensates the 3-sigma truncation, so the total is exactly
        # preserved rather than merely edge-corrected.
        available = np.convolve(np.ones_like(counts), kernel, mode="same")
        return np.convolve(counts / available, kernel, mode="same")

    raise ValueError(f"unknown convolution mode {mode!r}")


def convolve_spectrum(
    spectrum: Spectrum, fwhm_keV: float, *, mode: str = "rump"
) -> Spectrum:
    """:func:`convolve_detector` lifted to a :class:`Spectrum`."""
    return Spectrum(
        counts=convolve_detector(
            spectrum.counts, spectrum.calibration, fwhm_keV, mode=mode
        ),
        calibration=spectrum.calibration,
    )
