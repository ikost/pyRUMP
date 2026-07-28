"""Ad-hoc multiple-scattering tail.

Port of creatr.c:337-345.

.. warning::
   **This has no physical basis.** The C labels its own scale factor
   *"Ad-hoc scaling"*, and the 1985 paper lists multiple scattering among the
   effects the algorithms explicitly neglect. It is a one-parameter knob that
   adds a low-energy tail of roughly the right shape, nothing more.

   Provided for reproducing existing RUMP results. Treat any number it produces
   as qualitative.

The construction is a reverse cumulative sum -- each channel receives a share of
everything above it -- scaled by the spectrum's total intensity per unit charge
and solid angle:

.. code-block:: c

    sum = (total counts) / (q * omega) * 1.5E-9 * multiple;
    for (i = npt-1; i >= 0; i--) { tmp += counts[i]; counts[i] += tmp*sum; }
"""

from __future__ import annotations

import numpy as np

#: The C's ad-hoc constant (creatr.c:340).
_AD_HOC_SCALE = 1.5e-9


def add_multiple_scattering(
    counts: np.ndarray,
    *,
    strength: float,
    charge_uC: float,
    omega_msr: float,
) -> np.ndarray:
    """Add the empirical low-energy tail.

    ``strength`` is RUMP's ``SIM MULTIPLE_SCATTER`` value; 0 disables it, which
    is the default.
    """
    if strength == 0.0:
        return counts.copy()
    counts = np.asarray(counts, dtype=np.float64)
    if charge_uC == 0.0 or omega_msr == 0.0:
        return counts.copy()

    scale = counts.sum() / (charge_uC * omega_msr) * _AD_HOC_SCALE * strength
    # Reverse cumulative sum: channel i gains a share of everything at or above
    # it, including itself (the C adds counts[i] to tmp *before* using it).
    running = np.cumsum(counts[::-1])[::-1]
    return counts + running * scale
