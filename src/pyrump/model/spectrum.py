"""Spectrum representation and channel calibration.

Mirrors the calibration fields of ``SPECTRUM`` (spectrum.h:35-65).

.. note::
   **RUMP's channel convention is inconsistent between its own comment and its
   own code, and the code wins.**

   The header comment at anlyz.c:150-153 says channel *N* covers
   ``(N-0.5)*kevch + kev0`` to ``(N+0.5)*kevch + kev0`` -- i.e. the channel is
   *centred* on ``N*kevch + kev0``.

   ``SimAnlyz4`` (anlyz.c:304) instead integrates channel *k* over
   ``[(k+first)*kevch + kev0, (k+1+first)*kevch + kev0)`` -- the nominal energy
   is the channel's *lower edge*, half a channel away from the comment.

   pyRUMP follows the code, since that is what produced every RUMP spectrum.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True, slots=True)
class Calibration:
    """Linear channel-to-energy mapping."""

    kevch: float = 5.0
    """keV per channel."""

    kev0: float = 0.0
    """Energy offset at channel zero."""

    first: float = 0.0
    """Channel number of the first stored data point."""

    npt: int = 1024
    """Number of channels."""

    def edge_energy(self, channel) -> np.ndarray:
        """Lower edge of a channel, in keV (anlyz.c:311)."""
        return (np.asarray(channel, dtype=np.float64) + self.first) * self.kevch + self.kev0

    def edges(self) -> np.ndarray:
        """All ``npt + 1`` channel boundaries."""
        return self.edge_energy(np.arange(self.npt + 1))

    def channel_of(self, energy) -> np.ndarray:
        """Fractional channel index for an energy (inverse of :meth:`edge_energy`)."""
        return (np.asarray(energy, dtype=np.float64) - self.kev0) / self.kevch - self.first


@dataclass(slots=True)
class Spectrum:
    """A channel spectrum plus the calibration needed to interpret it."""

    counts: np.ndarray
    calibration: Calibration = field(default_factory=Calibration)

    @classmethod
    def zeros(cls, calibration: Calibration) -> "Spectrum":
        return cls(
            counts=np.zeros(calibration.npt, dtype=np.float64), calibration=calibration
        )

    def __len__(self) -> int:
        return int(self.counts.size)

    @property
    def energies(self) -> np.ndarray:
        """Lower edge energy of each channel."""
        return self.calibration.edge_energy(np.arange(self.counts.size))

    def total(self) -> float:
        return float(self.counts.sum())
