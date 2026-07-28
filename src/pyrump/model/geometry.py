"""Scattering geometry: angles and the path-length secants.

Port of the geometry block of ``FillSimHeader`` (creatr.c:405-428).

**The angle convention is the single most confusing thing in RUMP's interface.**
``phi`` as stored in a spectrum and typed by the user is the *supplement* of the
scattering angle: a detector at 170 degrees is entered as ``phi = 10``. The
simulation converts once, at creatr.c:409-410, and everything downstream uses the
true angle. :class:`Geometry` stores the RUMP convention and exposes
:attr:`scattering_angle` for the physical one.

Three tilt conventions, differing only in how the exit path length is derived:

* ``CORNELL`` - tilt axis lies in the scattering plane
* ``IBM`` - tilt axis perpendicular to the scattering plane
* ``GENERAL`` - the exit angle ``psi`` is given explicitly
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import IntEnum


class GeometryKind(IntEnum):
    """``spectrum.h`` geometry codes; note GENERAL is -1, not 2."""

    GENERAL = -1
    CORNELL = 0
    IBM = 1


@dataclass(frozen=True, slots=True)
class Geometry:
    """Beam/detector geometry for one measurement."""

    theta: float
    """Incidence angle from the sample normal, in degrees (the tilt)."""

    phi: float
    """**Supplement** of the scattering angle, in degrees. 170 deg -> phi = 10."""

    psi: float = 0.0
    """Exit angle from the normal, in degrees. Used by GENERAL only."""

    kind: GeometryKind = GeometryKind.CORNELL

    @property
    def scattering_angle(self) -> float:
        """The true scattering angle in degrees (creatr.c:1629)."""
        return 180.0 - self.phi

    @property
    def sin_phi(self) -> float:
        """``sin`` of the stored phi (creatr.c:409)."""
        return math.sin(math.radians(self.phi))

    @property
    def cos_phi(self) -> float:
        """``cos`` of the *true* scattering angle.

        The C stores this negated (``samm->cosph = -COSD(phi)``, creatr.c:410)
        precisely because it wants cos of the true angle, and
        ``cos(180-x) = -cos(x)``.
        """
        return -math.cos(math.radians(self.phi))

    @property
    def sec_in(self) -> float:
        """Inbound path-length factor, ``1/cos(theta)``."""
        cos_theta = math.cos(math.radians(self.theta))
        if cos_theta == 0.0:
            return 0.0  # the C leaves secin at 0, then rejects it below
        return 1.0 / cos_theta

    @property
    def sec_out(self) -> float:
        """Outbound path-length factor, by tilt convention (creatr.c:413-421)."""
        if self.kind == GeometryKind.CORNELL:
            cos_phi = self.cos_phi
            return 0.0 if cos_phi == 0.0 else -self.sec_in / cos_phi
        if self.kind == GeometryKind.IBM:
            return 1.0 / math.cos(math.radians(self.theta + self.phi))
        return 1.0 / math.cos(math.radians(self.psi))

    def validate(self) -> None:
        """Reject geometries RUMP refuses (creatr.c:423-426).

        A non-positive secant means the beam or the exit path runs along or
        behind the surface, which the slab model cannot represent.
        """
        if self.sec_in <= 0 or self.sec_out <= 0:
            raise ValueError(
                f"bad scattering geometry: sec_in={self.sec_in:.4g}, "
                f"sec_out={self.sec_out:.4g} (theta={self.theta}, phi={self.phi}, "
                f"psi={self.psi}, {self.kind.name}) - check the angles"
            )
