"""The X-ray detector, with pyRUMP's built-in defaults.

The defaults play the part :class:`~pyrump.model.detector.Measurement`'s do
for RBS: something sensible to simulate with until ``~/.pyrumprc`` or the
PIXE prompt says otherwise. They describe a real setup -- the NEC RC43
endstation the ``examples/MnPt`` pair was measured on:

* an Amptek silicon drift detector: 12.5 µm Be window, 25 mm^2 active
  area, a 500 µm crystal (Amptek's FAST SDD);
* looking at the sample through a tube from 7.125 in (181.0 mm) away, so
  Omega = 0.763 msr;
* at 45 deg to the normal of the untilted sample, in the plane the sample
  is tilted in: a THETA of -9 deg turns the sample towards it, 36 deg out;
* behind a Mylar filter against bremsstrahlung, specified as 125 µm but
  acting like 62 µm: the thick-target Si K yield gives 63.0 µm for a 258 nm
  SiO2-on-Si reference (2022) and 60.8 µm for ``examples/MnPt`` (2026)
  -- an effective value, which also takes up any error in the Si K
  cross section;
* resolution fitted to those spectra: Si Kα 78.5 eV, Mn Kα (as the Kα1/Kα2
  doublet) 122 +- 3 eV FWHM, hence a Fano factor of 0.104 and about 50 eV
  of electronic noise (the specification says 130 eV at 5.9 keV);
* energy calibration from the RC43 header of ``examples/MnPt.PIX``
  (10.097 eV/channel, offset -38.6 eV, channels numbered from 1).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from ..model.spectrum import Calibration

#: Channel calibration used when a spectrum's own header has none (see
#: :mod:`pyrump.io.oxford`), and for a PIXE simulation with no data.
DEFAULT_CALIBRATION = Calibration(kevch=0.01009699, kev0=-0.03864563, first=1.0, npt=2048)


@dataclass(frozen=True, slots=True)
class Compound:
    """An absorber material other than a pure element."""

    name: str
    atoms: dict[str, int]
    """Stoichiometry: element symbol -> atoms per formula unit."""

    density_g_cm3: float


#: Absorber materials known by name. Densities are NIST's (the X-ray mass
#: attenuation tables' material list).
COMPOUNDS = {
    "MYLAR": Compound("Mylar", {"C": 10, "H": 8, "O": 4}, 1.40),
    "KAPTON": Compound("Kapton", {"C": 22, "H": 10, "N": 2, "O": 5}, 1.42),
}


@dataclass(frozen=True, slots=True)
class Absorber:
    """A uniform layer between the sample and the crystal."""

    material: str
    """An element symbol, or a compound in :data:`COMPOUNDS` (``Mylar``)."""

    thickness_um: float

    hole_percent: float = 0.0
    """Open area of a "funny filter", in percent; 0 for a solid foil."""

    @property
    def compound(self) -> Compound | None:
        return COMPOUNDS.get(self.material.upper())


def disc_solid_angle_msr(area_mm2: float, distance_mm: float) -> float:
    """Solid angle of a round detector of ``area_mm2`` seen on axis from
    ``distance_mm`` away, in msr."""
    radius = math.sqrt(area_mm2 / math.pi)
    return 2 * math.pi * (1 - distance_mm / math.hypot(distance_mm, radius)) * 1000.0


@dataclass(frozen=True, slots=True)
class PixeDetector:
    """Geometry and response of the X-ray detector."""

    angle_deg: float = 45.0
    """Detector axis to the normal of the *untilted* sample, degrees."""

    tilt_sign: int = 1
    """How the sample's tilt (the RBS geometry's THETA) moves the exit
    angle: ``|angle + tilt_sign * THETA|``. 1 when a negative THETA turns
    the sample towards the detector, -1 when it turns it away, 0 when the
    tilt axis leaves the detector direction alone."""

    solid_angle_msr: float = field(default_factory=lambda: disc_solid_angle_msr(25.0, 180.975))

    window: Absorber = Absorber("Be", 12.5)
    crystal: Absorber = Absorber("Si", 500.0)

    fwhm_eV: float = 122.0
    """Resolution at Mn Kα (5.899 keV)."""

    fano: float = 0.104

    filters: tuple[Absorber, ...] = (Absorber("Mylar", 62.0),)
    """Absorbers in front of the window, in the order the X-rays meet them."""

    def exit_angle(self, theta_deg: float) -> float:
        """The X-rays' angle to the sample normal when the sample is tilted
        by ``theta_deg`` (the RBS THETA, signed)."""
        return abs(self.angle_deg + self.tilt_sign * theta_deg)
