"""Rutherford cross-sections, screening, and the built-in non-Rutherford forms.

Port of ``sigma.c``. Everything is **lab frame and non-relativistic**; RUMP has no
relativistic correction.

Structure follows the C: a setup step computes energy-independent constants, then
an evaluator applies the energy dependence. Cross-sections are barns/sr.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum

import numpy as np

#: ``(e^2/4)^2`` in keV^2-barn. Recomputed from CODATA in 2/2007 "on
#: recommendation of Nuno Barradas" (sigma.c:97).
E2_OVER_4_SQUARED = 1295.9358

#: ``(e^2/2)^2`` in keV^2-barn, for recoils (sigma.c:98).
E2_OVER_2_SQUARED = 5183.7432

#: Energy at which Quillet's D cross-section stops being valid (sigma.c:341).
QUILLET_MATCH_KEV = 2700.0

#: Above this the Ziegler H form is held constant in energy (sigma.c:303).
ZIEGLER_H_MAX_MEV = 4.0


class CrossSectionKind(Enum):
    """Which evaluator the setup selected -- mirrors ``sp->calc``."""

    RUTHERFORD = "rutherford"
    RUTHERFORD_SCREENED = "rutherford_screen"
    RUTHERFORD_OFFSET = "rutherford_offset"
    ZIEGLER_H = "ziegler_hscatt"
    QUILLET_D = "quillet_dscatt"


@dataclass(slots=True)
class CrossSection:
    """Configured cross-section, equivalent to RUMP's ``SP`` struct."""

    kind: CrossSectionKind
    csigma: float
    """Coefficient of the 1/E^2 term."""

    csig_0: float = 0.0
    """Constant term (manual-override path only)."""

    csig_f: float = 0.0
    """Screening rolloff coefficient, keV."""

    cos_phi: float = 1.0
    pf: tuple[float, ...] = field(default_factory=tuple)
    """Private angle-dependent parameters (Quillet)."""

    def __call__(self, energy_keV) -> np.ndarray:
        """Evaluate in barns/sr."""
        energy = np.atleast_1d(np.asarray(energy_keV, dtype=np.float64))

        if self.kind is CrossSectionKind.RUTHERFORD:
            return self.csigma / energy / energy
        if self.kind is CrossSectionKind.RUTHERFORD_SCREENED:
            return (self.csigma / energy / energy) * (1.0 - self.csig_f / energy)
        if self.kind is CrossSectionKind.RUTHERFORD_OFFSET:
            return self.csig_0 + self.csigma / energy / energy
        if self.kind is CrossSectionKind.ZIEGLER_H:
            return self._ziegler_h(energy)
        return self._quillet_d(energy)

    def _ziegler_h(self, energy: np.ndarray) -> np.ndarray:
        """1H(4He,1H)4He, Ziegler NIM B136-138 (1998) 141 refit (sigma.c:298).

        Valid for true scattering below 40 deg and below 4 MeV; above 4 MeV the
        *ratio* is frozen so the result falls off as pure 1/E^2.
        """
        c = (0.09720717, 1.359809, 0.1429851, 3.06073, 5.406188)
        mev = np.maximum(energy / 1000.0, 0.001)
        mev4 = np.minimum(mev, ZIEGLER_H_MAX_MEV)
        ratio = 1.0 + c[0] * mev4 ** c[1] + c[2] * mev4 ** c[3] * self.cos_phi ** c[4]
        return ratio / mev / mev * self.csigma

    def _quillet_d(self, energy: np.ndarray) -> np.ndarray:
        """2H(4He,2H)4He, Quillet NIM B83 (1993) 47 (sigma.c:337).

        Above 2.7 MeV the *cross-section itself* is held constant -- not scaled
        as Rutherford. The C offers a 1/E^2 alternative but it is compiled out
        (``USE_CONSTANT_CROSS_SECTION_EXTRAPOLOATION``, typo original).
        """
        a, b, c, d, k = self.pf
        clipped = np.minimum(energy, QUILLET_MATCH_KEV)
        sigma = 0.001 * (
            a / ((clipped - 2128.0) ** 2 + b) + c + clipped * (d + clipped * k)
        )
        return sigma * self.csigma


def setup_scatter(
    z1: int, m1: float, z2: int, m2: float, scattering_angle_deg: float,
    *, screening: bool = True,
) -> CrossSection:
    """RBS cross-section: Rutherford with L'Ecuyer screening (sigma.c:100).

    .. math::
        \\sigma = \\left(\\frac{Z_1 Z_2 e^2}{4E}\\right)^2
                  \\frac{4}{\\sin^4\\phi}
                  \\frac{\\left[\\sqrt{1-(x\\sin\\phi)^2}+\\cos\\phi\\right]^2}
                       {\\sqrt{1-(x\\sin\\phi)^2}}

    The screening correction is L'Ecuyer's ``(1 - 0.049 Z1 Z2^{4/3} / E)``
    (NIM 160 (1979) 337). Andersen screening is not implemented in RUMP.
    """
    x = m1 / m2
    phi = math.radians(scattering_angle_deg)
    sin_phi, cos_phi = math.sin(phi), math.cos(phi)

    sqirt = math.sqrt(1.0 - (x * sin_phi) ** 2)
    if sin_phi == 0.0:
        detail = (1.0 - x * x) ** 2
    else:
        detail = 4.0 / sin_phi**4 * (sqirt + cos_phi) ** 2 / sqirt

    return CrossSection(
        kind=CrossSectionKind.RUTHERFORD_SCREENED if screening else CrossSectionKind.RUTHERFORD,
        csigma=E2_OVER_4_SQUARED * detail * (z1 * z2) ** 2,
        # F&M 2.21; disabled by RUMP's NO_SIGMA_SCREEN switch.
        csig_f=(0.049 * z1 * z2**1.3333) if screening else 0.0,
        cos_phi=cos_phi,
    )


def setup_recoil(
    z1: int, m1: float, z2: int, m2: float, scattering_angle_deg: float,
    *, force_rutherford: bool = False,
) -> CrossSection:
    """ERD cross-section (sigma.c:144).

    .. math::
        \\sigma = \\left(\\frac{Z_1 Z_2 e^2 (m_1+m_2)}{2 m_2 E}\\right)^2
                  \\frac{1}{\\cos^3\\phi}

    Recoils get **no screening**. For 4He on H or D, RUMP substitutes measured
    analytic forms within their validity windows and silently falls back to
    Rutherford outside them.
    """
    phi = math.radians(scattering_angle_deg)
    cos_phi = math.cos(phi)
    if cos_phi <= 0.0:
        raise ValueError("recoil geometry needs a true scattering angle below 90 deg")

    csigma = E2_OVER_2_SQUARED * (z1 * z2 * (1.0 + m1 / m2)) ** 2 / cos_phi**3
    rutherford = CrossSection(
        kind=CrossSectionKind.RUTHERFORD, csigma=csigma, cos_phi=cos_phi
    )
    if force_rutherford or not (z1 == 2 and z2 == 1):
        return rutherford

    is_deuterium = round(m2) == 2
    if not is_deuterium:
        if scattering_angle_deg > 40.0:
            return rutherford  # outside Ziegler's angular range
        # csigma becomes a pure scaling factor here (default 1), because the
        # analytic form already carries the absolute magnitude (sigma.c:181).
        return CrossSection(
            kind=CrossSectionKind.ZIEGLER_H, csigma=1.0, cos_phi=cos_phi
        )

    if not (10.0 < scattering_angle_deg < 32.0):
        return rutherford  # outside Quillet's angular range

    p = scattering_angle_deg
    return CrossSection(
        kind=CrossSectionKind.QUILLET_D,
        csigma=1.0,
        cos_phi=cos_phi,
        pf=(
            -2.6e3 * p * p - 1.76e5 * p + 8.79e6,  # A
            0.18 * p * p - 10.0 * p + 1422.0,  # B
            -2.59 * p * p + 111.6 * p - 72.0,  # C
            3.09e-3 * p * p - 1.278e-1 * p + 5.83e-1,  # D
            -9.05e-7 * p * p + 3.645e-5 * p - 1.71e-4,  # K
        ),
    )
