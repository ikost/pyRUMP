r"""From line yields to a channel spectrum: the detector's response.

**Resolution.** A line of energy E is a Gaussian of variance

.. math:: \sigma^2(E) = \sigma_\text{noise}^2 + \varepsilon_\text{Si}\,F\,E

(electronic noise plus the statistics of charge creation, :math:`\varepsilon_\text{Si}`
= 3.64 eV per electron-hole pair, F the Fano factor). The detector is set up by
its FWHM at Mn Kα and F, which fix :math:`\sigma_\text{noise}`. Each line is
integrated over every channel's energy range, so the counts are exact however
narrow the peak.

**Si escape.** An X-ray above the Si K edge can eject a Si K photon that leaves
the crystal, and is then recorded 1.740 keV low. The escaping fraction, for
normal incidence (Reed & Ware, J. Phys. E 5 (1972) 582), is

.. math:: \eta = \tfrac12\,\omega_K\,(1 - 1/r)\,
          \big[1 - \tfrac{\mu_K}{\mu_E}\ln(1 + \tfrac{\mu_E}{\mu_K})\big]

with :math:`\omega_K` the Si K fluorescence yield, r the K jump ratio and
:math:`\mu_E`, :math:`\mu_K` silicon's attenuation at the line and at Si Kα.
That fraction moves from the line's peak to its escape peak.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.special import ndtr

from ..model.spectrum import Calibration
from .atomic import AtomicData, atomic_data
from .detector import PixeDetector
from .yields import LineYield

EPSILON_SI_EV = 3.64
SIGMA_PER_FWHM = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
#: Mn Kα1, the energy the FWHM is quoted at (xraylib's KL3).
MN_KALPHA_EV = 5898.8
SI_Z = 14


def resolution_sigma_keV(detector: PixeDetector, energy_keV) -> np.ndarray:
    """Gaussian sigma of a line, keV."""
    reference = (detector.fwhm_eV * SIGMA_PER_FWHM) ** 2
    fano = EPSILON_SI_EV * detector.fano
    noise = max(reference - fano * MN_KALPHA_EV, 0.0)
    energy_eV = np.asarray(energy_keV, dtype=np.float64) * 1000.0
    return np.sqrt(noise + fano * energy_eV) / 1000.0


def escape_fraction(energy_keV: float, detector: PixeDetector, atomic: AtomicData) -> float:
    """Fraction of a line recorded in its Si escape peak (0 for a non-Si
    crystal or a line below the Si K edge)."""
    if detector.crystal.material != "Si" or energy_keV <= atomic.edge(SI_Z, "K"):
        return 0.0
    omega = atomic.fluorescence_yield(SI_Z, "K")
    jump = atomic.jump_factor(SI_Z, "K")
    si_kalpha = _si_kalpha(atomic)
    ratio = float(atomic.mu(SI_Z, si_kalpha) / atomic.mu(SI_Z, energy_keV))
    return 0.5 * omega * (1 - 1 / jump) * (1 - ratio * np.log1p(1 / ratio))


def _si_kalpha(atomic: AtomicData) -> float:
    return next(g.energy_keV for g in atomic.line_groups(SI_Z) if g.label == "Kα")


@dataclass(slots=True)
class PixeSpectrum:
    """A simulated PIXE spectrum, with what went into it."""

    counts: np.ndarray
    by_element: dict[str, np.ndarray]
    lines: list[LineYield]
    calibration: Calibration
    detector: PixeDetector
    escape: bool = True

    def by_layer(self, layer: int) -> np.ndarray:
        """The part of the spectrum from one SIM layer (0 = the surface)."""
        lines = [replace(line, counts=line.by_layer.get(layer, 0.0)) for line in self.lines]
        return synthesize(lines, self.calibration, self.detector, escape=self.escape).counts


def synthesize(
    lines: list[LineYield],
    calibration: Calibration,
    detector: PixeDetector,
    *,
    escape: bool = True,
    atomic: AtomicData | None = None,
) -> PixeSpectrum:
    """Spread each line over the channels as the detector records it."""
    atomic = atomic or atomic_data()
    edges = calibration.edges()
    total = np.zeros(calibration.npt)
    by_element: dict[str, np.ndarray] = {}
    escape_shift = _si_kalpha(atomic)

    def add(energy: float, area: float, into: np.ndarray) -> None:
        if area <= 0:
            return
        sigma = float(resolution_sigma_keV(detector, energy))
        into += area * np.diff(ndtr((edges - energy) / sigma))

    for line in lines:
        element = by_element.setdefault(line.symbol, np.zeros(calibration.npt))
        eta = escape_fraction(line.energy_keV, detector, atomic) if escape else 0.0
        add(line.energy_keV, line.counts * (1 - eta), element)
        add(line.energy_keV - escape_shift, line.counts * eta, element)
    for spectrum in by_element.values():
        total += spectrum
    return PixeSpectrum(
        counts=total, by_element=by_element, lines=lines, calibration=calibration,
        detector=detector, escape=escape,
    )
