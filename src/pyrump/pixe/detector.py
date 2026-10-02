"""The X-ray detector, with pyRUMP's built-in defaults.

The defaults play the part :class:`~pyrump.model.detector.Measurement`'s do
for RBS: something sensible to simulate with until ``~/.pyrumprc`` or the
PIXE prompt says otherwise. Where a value could be read off a real
measurement, it was -- the MnPt example's SDD (an NEC RC43 endstation):

* energy calibration: the RC43 header of ``examples/MnPt.PIX``
  (10.097 eV/channel, offset -38.6 eV, channels numbered from 1, 2048 of
  them);
* resolution: 131 eV FWHM fitted to its 5.9 keV line, and a Fano factor of
  0.13 from the same fit together with Si Kα's 79 eV.

The geometry, window and crystal are typical SDD values, not measured ones.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..model.spectrum import Calibration

#: Channel calibration used when a spectrum's own header has none (see
#: :mod:`pyrump.io.oxford`), and for a PIXE simulation with no data.
DEFAULT_CALIBRATION = Calibration(kevch=0.01009699, kev0=-0.03864563, first=1.0, npt=2048)


@dataclass(frozen=True, slots=True)
class Absorber:
    """A uniform layer of one element between the sample and the crystal."""

    element: str
    """Element symbol."""

    thickness_um: float

    hole_percent: float = 0.0
    """Open area of a "funny filter", in percent; 0 for a solid foil."""


@dataclass(frozen=True, slots=True)
class PixeDetector:
    """Geometry and response of the X-ray detector."""

    angle_deg: float = 45.0
    """Detector axis to the sample normal, degrees (GUPIX's convention)."""

    solid_angle_msr: float = 1.0

    window: Absorber = Absorber("Be", 8.0)
    crystal: Absorber = Absorber("Si", 450.0)

    fwhm_eV: float = 131.0
    """Resolution at Mn Kα (5.899 keV)."""

    fano: float = 0.13

    filters: tuple[Absorber, ...] = ()
    """Absorbers in front of the window, in the order the X-rays meet them."""
