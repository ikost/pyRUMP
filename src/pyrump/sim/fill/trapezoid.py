"""Brick -> channel fill without straggling.

Port of ``SimAnlyz4`` (anlyz.c:304-343).

The C walks each brick channel by channel, accumulating trapezoid areas and
clipping the two partial end channels. Since the brick top is a straight line,
that sum is exactly the integral of a linear function over each channel's
overlap with the brick, so pyRUMP evaluates the antiderivative at the clipped
boundaries and differences. Mathematically identical -- including the partial
end channels -- and it vectorises over every brick at once.

.. warning::
   **A brick reaching past the top channel is discarded entirely, not clipped**
   (anlyz.c:313-316). The C prints "Energy out of range on exit" and returns
   before adding anything. Reproduced, because clipping instead would silently
   change published spectra.
"""

from __future__ import annotations

import numpy as np

from ...model.spectrum import Calibration
from ..bricks import Bricks


def _c_trunc(values: np.ndarray) -> np.ndarray:
    """C's ``(int)`` cast: truncation toward zero, not floor.

    ``(int)(-2.7)`` is ``-2`` in C but ``-3`` under ``np.floor``. Every
    translated cast uses this to avoid off-by-one channel shifts.
    """
    return np.trunc(values).astype(np.int64)


def fill_trapezoid(
    bricks: Bricks,
    calibration: Calibration,
    counts: np.ndarray | None = None,
) -> np.ndarray:
    """Accumulate bricks into channels, returning the counts array.

    Values are areas in height x keV; RUMP's overall normalisation is applied
    later, in the equivalent of ``RbsNormK`` scaling.
    """
    if counts is None:
        counts = np.zeros(calibration.npt, dtype=np.float64)
    if len(bricks) == 0:
        return counts

    e_front = bricks.e_front
    e_back = bricks.e_back
    h_front = bricks.h_front
    h_back = bricks.h_back

    span = e_front - e_back
    # anlyz.c:308 bails on inverted or degenerate bricks before dividing.
    usable = span > 0
    if not np.any(usable):
        return counts

    slope = np.zeros_like(span)
    np.divide(h_front - h_back, span, out=slope, where=usable)

    # anlyz.c:311 and :330, including the C truncation semantics.
    k1 = _c_trunc(calibration.channel_of(e_front) + 1)
    k0 = _c_trunc(calibration.channel_of(e_back))

    # Bricks running off the top are dropped whole; ones entirely below the
    # first channel contribute nothing.
    usable &= (k1 < calibration.npt) & (k1 >= 0)
    k0 = np.maximum(k0, 0)
    usable &= k1 > k0
    if not np.any(usable):
        return counts

    index = np.flatnonzero(usable)
    widths = (k1 - k0)[index]

    # Flatten the ragged (brick, channel) pairs into one long index list.
    brick_id = np.repeat(index, widths)
    offsets = np.arange(widths.sum()) - np.repeat(
        np.cumsum(widths) - widths, widths
    )
    channel = np.repeat(k0[index], widths) + offsets

    # Clip each channel's span to the brick, exactly as the C's e0/e1 do.
    lower = np.maximum(calibration.edge_energy(channel), e_back[brick_id])
    upper = np.minimum(calibration.edge_energy(channel + 1), e_front[brick_id])

    # Integral of the linear brick top over [lower, upper]: the trapezoid rule
    # is exact here, which is why the C's loop and this agree identically.
    h_lower = h_back[brick_id] + slope[brick_id] * (lower - e_back[brick_id])
    h_upper = h_back[brick_id] + slope[brick_id] * (upper - e_back[brick_id])
    area = (upper - lower) * (h_lower + h_upper) * 0.5

    np.add.at(counts, channel, area)
    return counts


def fill_trapezoid_reference(
    bricks: Bricks, calibration: Calibration, counts: np.ndarray | None = None
) -> np.ndarray:
    """Literal transliteration of ``SimAnlyz4``'s loop.

    Kept as an executable statement of what the vectorised version must equal;
    the test suite asserts they agree bit for bit.
    """
    if counts is None:
        counts = np.zeros(calibration.npt, dtype=np.float64)

    for brick in range(len(bricks)):
        e_front = float(bricks.e_front[brick])
        e_back = float(bricks.e_back[brick])
        h_front = float(bricks.h_front[brick])
        h_back = float(bricks.h_back[brick])

        if e_back >= e_front:
            continue
        slope = (h_front - h_back) / (e_front - e_back)

        k1 = int(calibration.channel_of(e_front) + 1)
        if k1 >= calibration.npt or k1 < 0:
            continue

        zero_energy = calibration.edge_energy(0)
        if e_back < zero_energy:
            h_back = h_back + slope * (zero_energy - e_back)
            e_back = float(zero_energy)
        k0 = int(calibration.channel_of(e_back))

        e0, h0 = e_back, h_back
        e1 = float(calibration.edge_energy(k0 + 1))
        for channel in range(k0, k1):
            if e1 > e_front:
                e1 = e_front
            h1 = h_back + slope * (e1 - e_back)
            counts[channel] += (e1 - e0) * (h1 + h0) / 2.0
            h0, e0 = h1, e1
            e1 += calibration.kevch
    return counts
