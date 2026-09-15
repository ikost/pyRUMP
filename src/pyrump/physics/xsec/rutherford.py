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


class ScreeningModel(Enum):
    """Which screening correction to apply to the Rutherford cross-section.

    RUMP itself only ever had L'Ecuyer -- see ``setup_scatter``'s docstring.
    ANDERSEN is a pyRUMP-only addition with no C oracle to validate against.
    """

    NONE = "none"
    LECUYER = "lecuyer"
    ANDERSEN = "andersen"


class CrossSectionKind(Enum):
    """Which evaluator the setup selected -- mirrors ``sp->calc``."""

    RUTHERFORD = "rutherford"
    RUTHERFORD_SCREENED = "rutherford_screen"
    RUTHERFORD_ANDERSEN = "rutherford_andersen"
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
    """Screening rolloff coefficient, keV (L'Ecuyer)."""

    andersen_v1: float = 0.0
    """V1, the screening energy scale, keV (Andersen)."""

    andersen_k_cm: float = 1.0
    """E_CM / E_lab = m2/(m1+m2), energy-independent (Andersen)."""

    andersen_sin_half_theta_cm: float = 0.0
    """sin(theta_CM / 2), energy-independent (Andersen)."""

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
        if self.kind is CrossSectionKind.RUTHERFORD_ANDERSEN:
            return (self.csigma / energy / energy) * self._andersen_factor(energy)
        if self.kind is CrossSectionKind.RUTHERFORD_OFFSET:
            return self.csig_0 + self.csigma / energy / energy
        if self.kind is CrossSectionKind.ZIEGLER_H:
            return self._ziegler_h(energy)
        return self._quillet_d(energy)

    def _andersen_factor(self, energy: np.ndarray) -> np.ndarray:
        """Andersen screening correction, Eq. 4.14 of the SIMNRA manual.

        .. math::
            F_{Andersen} = \\frac{(1 + \\tfrac{1}{2}V_1/E_{CM})^2}
                {\\left\\{1 + V_1/E_{CM} +
                    \\left[V_1/(2 E_{CM}\\sin(\\theta_{CM}/2))\\right]^2\\right\\}^2}

        ``theta_CM`` and the ``E_CM/E_lab`` ratio are both energy-independent
        (non-relativistic two-body kinematics), so only ``V1/E_CM`` varies
        with the lab energy passed in here.
        """
        e_cm = energy * self.andersen_k_cm
        ratio = self.andersen_v1 / e_cm
        angular = self.andersen_v1 / (2.0 * e_cm * self.andersen_sin_half_theta_cm)
        return (1.0 + 0.5 * ratio) ** 2 / (1.0 + ratio + angular**2) ** 2

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
    *, screening: ScreeningModel = ScreeningModel.LECUYER,
) -> CrossSection:
    """RBS cross-section: Rutherford with a screening correction (sigma.c:100).

    .. math::
        \\sigma = \\left(\\frac{Z_1 Z_2 e^2}{4E}\\right)^2
                  \\frac{4}{\\sin^4\\phi}
                  \\frac{\\left[\\sqrt{1-(x\\sin\\phi)^2}+\\cos\\phi\\right]^2}
                       {\\sqrt{1-(x\\sin\\phi)^2}}

    ``screening=LECUYER`` (the default, and RUMP's own and only option) applies
    L'Ecuyer's ``(1 - 0.049 Z1 Z2^{4/3} / E)`` (NIM 160 (1979) 337) -- a first-order
    correction with no angular dependence.

    ``screening=ANDERSEN`` is a pyRUMP-only addition (RUMP never had it): Andersen,
    Besenbacher, Loftager & Moller's angle- *and* energy-dependent correction (Phys.
    Rev. A 21 (1980) 1891), transcribed from SIMNRA's manual Eq. 4.14. More accurate
    than L'Ecuyer at forward angles, where L'Ecuyer's angle-independence
    underestimates the needed correction -- but per SIMNRA's own manual, "may be
    inaccurate at energies below a few hundred keV and/or small scattering angles",
    which overlaps the low end of RUMP's own stopping-power fit window. No RUMP
    oracle exists to validate this against; treat it as literature-derived, not
    cross-checked bit-for-bit like the rest of this port.
    """
    x = m1 / m2
    phi = math.radians(scattering_angle_deg)
    sin_phi, cos_phi = math.sin(phi), math.cos(phi)

    sqirt = math.sqrt(1.0 - (x * sin_phi) ** 2)
    if sin_phi == 0.0:
        detail = (1.0 - x * x) ** 2
    else:
        detail = 4.0 / sin_phi**4 * (sqirt + cos_phi) ** 2 / sqirt

    csigma = E2_OVER_4_SQUARED * detail * (z1 * z2) ** 2

    if screening is ScreeningModel.ANDERSEN:
        # Non-relativistic lab -> CM angle transform, reusing `sqirt` (valid for
        # m1 < m2, i.e. x < 1, true of every RBS projectile/target pair).
        cos_theta_cm = cos_phi * sqirt - x * sin_phi * sin_phi
        sin_half_theta_cm = math.sqrt(max(0.0, 1.0 - cos_theta_cm) / 2.0)
        v1 = 0.04873 * z1 * z2 * math.sqrt(z1 ** (2.0 / 3.0) + z2 ** (2.0 / 3.0))
        return CrossSection(
            kind=CrossSectionKind.RUTHERFORD_ANDERSEN,
            csigma=csigma,
            andersen_v1=v1,
            andersen_k_cm=m2 / (m1 + m2),
            andersen_sin_half_theta_cm=sin_half_theta_cm,
            cos_phi=cos_phi,
        )

    lecuyer = screening is ScreeningModel.LECUYER
    return CrossSection(
        kind=CrossSectionKind.RUTHERFORD_SCREENED if lecuyer else CrossSectionKind.RUTHERFORD,
        csigma=csigma,
        # F&M 2.21; disabled by RUMP's NO_SIGMA_SCREEN switch.
        csig_f=(0.049 * z1 * z2**1.3333) if lecuyer else 0.0,
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
