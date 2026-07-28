"""Pulse pile-up.

Port of ``SimNewPileup`` (creatr.c:1331-1419), the Custer-thesis model, and the
legacy ``SimOldPileup``.

Pile-up happens when two events arrive within the detector's shaping time and
are recorded as one of their combined energy. The model convolves the spectrum
with itself, weighted by how long an event spends near each energy while its
pulse rises:

.. code-block:: text

    tails[i] = sum over j>i of counts[j]/j        # the 1/j is the dwell-time term
    pileup[i+j] += counts[i] * tails[j]

then scales by the probability of a coincidence, ``rate * tau``, and removes two
real counts for every pile-up count created.

Requires both ``tau`` (shaping time) and ``current``; RUMP silently does nothing
if either is non-positive.

.. note::
   The C's own comment flags a factor-of-two ambiguity: *"If the shaping is
   really exponential, then the factor should be multiplied by 2 to account for
   the total area under an exponential versus a triangle of the same timing."*
   The shipped code does **not** apply that factor, and neither does pyRUMP.
   The 1985 thesis rates this model's accuracy at "30% typical".
"""

from __future__ import annotations

import numpy as np


def new_pileup(
    counts: np.ndarray,
    *,
    tau_us: float,
    current_nA: float,
    charge_uC: float,
    max_channels: int | None = None,
) -> np.ndarray:
    """Apply the Custer pile-up model.

    Returns a new array which may be **longer** than the input: pile-up moves
    counts to higher energy, so the C doubles ``npt`` up to ``nptmax``
    (creatr.c:1345).
    """
    if tau_us <= 0 or current_nA <= 0:
        return counts.copy()

    counts = np.asarray(counts, dtype=np.float64)
    npt = counts.size
    limit = max_channels if max_channels is not None else 2 * npt
    nptmax = max(npt, min(2 * npt, limit))

    # tails[i] = sum_{j>i} counts[j]/j, accumulated downward. The 1/j accounts
    # for an event at energy 2E spending only half as long at each dE on the
    # way up (creatr.c:1351-1355).
    tails = np.zeros(npt, dtype=np.float64)
    index = np.arange(1, npt, dtype=np.float64)
    contribution = counts[1:] / index
    # tails[i] is the sum strictly above i, so reverse-cumsum then shift.
    upper = np.cumsum(contribution[::-1])[::-1]
    tails[1:] = np.concatenate([upper[1:], [0.0]])

    total = float(counts[1:].sum())
    if not np.any(tails):
        return counts.copy()

    # Highest channel whose tail is still zero -- above it there is nothing to
    # pile up with (creatr.c:1358).
    nonzero = np.flatnonzero(tails)
    valid = int(nonzero.min()) if nonzero.size else npt

    pileup = np.zeros(nptmax, dtype=np.float64)
    for i in range(valid, npt):
        if counts[i] == 0.0:
            continue
        span = min(npt, nptmax - i)
        if span <= 0:
            break
        pileup[i : i + span] += counts[i] * tails[:span]

    # rate * tau is the coincidence probability. The count total cancels between
    # normalising tails[] and the rate, so it never has to be formed
    # (creatr.c:1395-1400).
    seconds = (charge_uC / current_nA) * 1e9
    pileup *= (1.0 / seconds) * tau_us
    pileup_total = float(pileup.sum())

    # Each pile-up event consumes two real ones.
    scale = (total - 2.0 * pileup_total) / total if total else 1.0

    out = np.zeros(nptmax, dtype=np.float64)
    out[:npt] = counts * scale
    out += pileup
    return out


def old_pileup(
    counts: np.ndarray, *, current_nA: float, strength: float = 3e-9
) -> np.ndarray:
    """The legacy DOS-era model (thesis 4.11).

    "Integrate the yield with respect to energy, working down from high energy.
    At each point multiply the integral by a constant and add the result to the
    spectrum at twice the present energy."

    The author's own verdict was "difficult to justify... not quantitatively
    very accurate (30% typical)". Provided for reproducing old results only.
    """
    counts = np.asarray(counts, dtype=np.float64)
    out = counts.copy()
    running = 0.0
    factor = strength * current_nA
    for i in range(counts.size - 1, -1, -1):
        running += counts[i]
        target = 2 * i
        if target < out.size:
            out[target] += running * factor
    return out
