"""Channel windows: which parts of a spectrum the fit actually sees.

Port of the window handling in ``pert.c``.

Two distinct kinds, easily confused:

**Error windows** select the channels the objective is summed over. RUMP allows
up to ten (``NUM_ERR_WINS``, pert.h:41), though the 1996 manual documents only
one. Everything outside every window is ignored entirely — not down-weighted.

**The normalisation window** is different in kind. Rather than selecting
channels, it removes a systematic: the total counts over that window are forced
to agree by scaling the *experimental* data's ``CORR`` factor before the
objective is evaluated. From the manual:

    To avoid errors from inaccuracies in the charge integration, a second
    window, the normalization window, is also defined. The total number of
    counts over the normalization window are made equal by varying the CORR
    factor of the experimental data before chi-square is computed.

Because that is itself a one-parameter fit, RUMP forbids combining it with
``CORR`` as a free parameter (pert.c:1163) — the two would be degenerate.

The manual also warns:

    Make sure that the theory spectra never goes to zero within the error
    window since a zero will wreak havoc with the Poisson statistics.

which is the ``theory <= 0`` branch of the objective.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

#: RUMP's limit on error windows (pert.h:41).
MAX_ERROR_WINDOWS = 10


@dataclass(frozen=True, slots=True)
class Window:
    """An inclusive channel range."""

    low: int
    high: int

    def __post_init__(self) -> None:
        if self.high < self.low:
            raise ValueError(f"window {self.low}..{self.high} is inverted")

    def mask(self, n_channels: int) -> np.ndarray:
        out = np.zeros(n_channels, dtype=bool)
        out[max(self.low, 0) : min(self.high + 1, n_channels)] = True
        return out


@dataclass(slots=True)
class WindowSet:
    """The error windows, plus an optional normalisation window."""

    error: list[Window] = field(default_factory=list)
    normalisation: Window | None = None

    def __post_init__(self) -> None:
        if len(self.error) > MAX_ERROR_WINDOWS:
            raise ValueError(
                f"RUMP allows at most {MAX_ERROR_WINDOWS} error windows, "
                f"got {len(self.error)}"
            )

    def mask(self, n_channels: int) -> np.ndarray:
        """Channels the objective is evaluated over.

        With no error windows the whole spectrum is used.
        """
        if not self.error:
            return np.ones(n_channels, dtype=bool)
        combined = np.zeros(n_channels, dtype=bool)
        for window in self.error:
            combined |= window.mask(n_channels)
        return combined

    def normalisation_factor(self, data: np.ndarray, theory: np.ndarray) -> float:
        """Scale for the *data* that equalises totals over the window.

        Returns 1.0 when no normalisation window is set. Also returns 1.0 rather
        than dividing by zero if the window contains no counts.
        """
        if self.normalisation is None:
            return 1.0
        mask = self.normalisation.mask(len(data))
        observed = float(np.sum(np.asarray(data)[mask]))
        expected = float(np.sum(np.asarray(theory)[mask]))
        if observed == 0.0:
            return 1.0
        return expected / observed

    def validate_against(self, varying_correction: bool) -> None:
        """Reject the combination RUMP forbids (pert.c:1163-1169)."""
        if varying_correction and self.normalisation is not None:
            raise ValueError(
                "a normalisation window cannot be used while CORRECTION is a free "
                "parameter -- they are degenerate, and RUMP rejects it too"
            )
