r"""Brick -> channel fill with straggling.

Port of ``SimAnlyz3`` (anlyz.c:244) and ``SimStragf`` (anlyz.c:371).

With straggling on, RUMP gives up the trapezoid: each brick is split into two
triangles, and each triangle is convolved with a Gaussian analytically (1985
paper, p. 349, Fig. 5). ``SimStragf`` supplies the integral of that convolution.

**pyRUMP replaces ``SimStragf`` with a closed form.** The C evaluates it through a
hand-tuned rational approximation split across seven regimes -- the 1985 paper
explains why:

    Analytic expressions for the functions f and g ... are of limited utility.
    Direct evaluation is unnecessarily slow and often involves finding the small
    difference of large numbers. Single precision computation is inadequate ...
    I have therefore used IMSL routine IRATCU to fit the functions in terms of
    rational Chebyshev polynomials.

Those objections were about 1985 single-precision hardware. In float64 with
``scipy.special.erf`` the closed form is more accurate, vectorises, and avoids
transcribing a table of magic constants. The test suite measures the difference
against the C rather than assuming it away.

The function
------------

``SimStragf(x, s)`` is, in the C's own words, "the integral of the convolution
of a triangle with a Gaussian", where the triangle is 1 at :math:`x=0` falling
linearly to 0 at :math:`x=1`, offset so the result runs from :math:`-1/4` to
:math:`+1/4`:

.. math::
    R(x) = \tfrac12 \int_0^1 (1-u)\,\operatorname{erf}\!\big(\tfrac{x-u}{s}\big)\,du

Substituting :math:`w=(x-u)/s` and integrating by parts with

.. math::
    G_1(w) = \int\operatorname{erf} = w\operatorname{erf}(w) + \tfrac{e^{-w^2}}{\sqrt\pi}

    G_2(w) = \int w\operatorname{erf}(w)\,dw
           = \big(\tfrac{w^2}{2}-\tfrac14\big)\operatorname{erf}(w)
             + \tfrac{w e^{-w^2}}{2\sqrt\pi}

gives, with :math:`A=(x-1)/s` and :math:`B=x/s`:

.. math::
    R(x) = \tfrac{s^2}{2}\Big[G_2(B)-G_2(A) - A\big(G_1(B)-G_1(A)\big)\Big]

Note ``s`` follows the C's convention :math:`\exp(-t^2/s^2)`, i.e.
:math:`s = \sqrt2\,\sigma`.
"""

from __future__ import annotations

import numpy as np
from scipy.special import erf

from ...model.spectrum import Calibration
from ..bricks import Bricks

_INV_SQRT_PI = 1.0 / np.sqrt(np.pi)

#: sqrt(2) and 3*sqrt(2), as they appear in anlyz.c:256-258.
_SQRT2 = 1.414214
_THREE_SQRT2 = 4.2426402


def _g1(w: np.ndarray) -> np.ndarray:
    """Antiderivative of erf."""
    return w * erf(w) + _INV_SQRT_PI * np.exp(-np.minimum(w * w, 700.0))


def _g2(w: np.ndarray) -> np.ndarray:
    """Antiderivative of w*erf(w)."""
    return (0.5 * w * w - 0.25) * erf(w) + 0.5 * w * _INV_SQRT_PI * np.exp(
        -np.minimum(w * w, 700.0)
    )


def triangle_gaussian_integral(x, s: float) -> np.ndarray:
    r"""The convolution integral itself, in unscaled coordinates.

    ``x`` is measured in units of the triangle's base, so the triangle occupies
    :math:`[0, 1]`. Returns values in :math:`[-1/4, +1/4]`.
    """
    x = np.atleast_1d(np.asarray(x, dtype=np.float64))
    if s == 0.0:
        return np.where(x > 0.0, x - 0.5 * x * x - 0.25, -0.25)

    a = (x - 1.0) / s
    b = x / s
    return 0.5 * s * s * (_g2(b) - _g2(a) - a * (_g1(b) - _g1(a)))


def stragf(x, sig: float) -> np.ndarray:
    r"""Drop-in replacement for ``SimStragf`` (anlyz.c:371), same contract.

    .. important::
       ``SimStragf`` rescales its argument before doing anything else --
       ``newx = x * (1 + 3*sig)`` at anlyz.c:387. Its ``x`` is therefore in
       units of *the broadened width* :math:`|de| + 3\sqrt2\,\sigma`, not of the
       triangle base. Missing this makes the function look wrong by up to 0.23
       out of a total range of 0.5.

       The rescaling exactly undoes the ``fact`` normalisation applied by the
       caller (anlyz.c:257), so the composition is simply
       :math:`(E_j - E_{peak}) / |de|`.
    """
    x = np.atleast_1d(np.asarray(x, dtype=np.float64))
    return triangle_gaussian_integral(x * (1.0 + 3.0 * sig), sig)


def add_triangle(
    counts: np.ndarray,
    calibration: Calibration,
    height: float,
    energy: float,
    de: float,
    sigma: float,
) -> None:
    """Add one Gaussian-broadened triangle to the spectrum (``SimAnlyz3``).

    The triangle peaks at ``energy`` with the given ``height`` and reaches zero
    at ``energy + de``; ``de`` may be negative, which mirrors it. ``sigma`` is
    the true standard deviation in keV.
    """
    width = abs(de)
    if width == 0.0:
        return

    # anlyz.c:256-259. Note asig is in the exp(-t^2/s^2) convention.
    asig = _SQRT2 * sigma / width
    fact = 1.0 / (width + _THREE_SQRT2 * sigma)
    a = fact * calibration.kevch
    b = (
        (0.5 + calibration.first) * calibration.kevch + calibration.kev0 - energy
    ) * fact
    c = height * width

    # C truncation, not floor (anlyz.c:266-267).
    jmin = int(np.trunc(-(1.0 + b) / a + 1.0))
    jmax = int(np.trunc((1.0 - b) / a))

    if jmin < 0:
        jmin = 0
        running = float(stragf(b - a, asig)[0])
    else:
        running = -0.25
    jmax = min(jmax, calibration.npt - 1)

    if jmax >= jmin:
        channels = np.arange(jmin, jmax + 1)
        if de >= 0.0:
            values = stragf(a * channels + b, asig)
        else:
            # Mirrored: scan from x=1 to -1 (anlyz.c:288).
            values = -stragf(-a * channels - b, asig)
        previous = np.concatenate([[running], values[:-1]])
        np.add.at(counts, channels, c * (values - previous))
        running = float(values[-1])

    # anlyz.c:294-297: whatever is left goes into the channel past the end, so
    # the triangle's full area is conserved.
    if jmax >= 0 and jmax + 1 < calibration.npt:
        counts[jmax + 1] += c * (0.25 - running)


def fill_straggled(
    bricks: Bricks,
    calibration: Calibration,
    counts: np.ndarray | None = None,
) -> np.ndarray:
    """Fill bricks that carry straggling, as two triangles each.

    ``SimAnlyz`` (anlyz.c:207-215) splits the trapezoid at its two edges and
    convolves each with its own width -- explicitly "giving up the parabolic
    nature", and in fact giving up the trapezoid too.
    """
    if counts is None:
        counts = np.zeros(calibration.npt, dtype=np.float64)

    for index in range(len(bricks)):
        e_front = float(bricks.e_front[index])
        e_back = float(bricks.e_back[index])
        de = e_front - e_back
        add_triangle(
            counts,
            calibration,
            float(bricks.h_back[index]),
            e_back,
            de,
            float(np.sqrt(max(bricks.sig_back[index], 0.0))),
        )
        add_triangle(
            counts,
            calibration,
            float(bricks.h_front[index]),
            e_front,
            -de,
            float(np.sqrt(max(bricks.sig_front[index], 0.0))),
        )
    return counts
