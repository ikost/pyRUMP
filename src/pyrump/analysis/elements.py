"""Element-edge physics shared by ELEMENT, MATRIX, INFO, WHATISIT and WIDTH_THICK.

Ports ``anlytc.c``'s ``RbsSetexp``/``RbsKappa``/``RbsSigma``/``RbsStoper``/
``RbsEpsilon``, built on top of pyRUMP's already oracle-validated kinematics,
cross-section and stopping modules rather than reimplementing their raw C
formulas from scratch.

Two deliberate reuse decisions worth knowing:

* :func:`cross_section_barns` calls :func:`~pyrump.physics.xsec.rutherford.setup_scatter`
  with ``screening=False``. ``anlytc.c``'s own local ``RbsSigma`` (line 1277) is the
  plain Rutherford formula with no L'Ecuyer screening term, unlike the public
  ``sigma.c`` that ``setup_scatter`` normally ports with screening on.
* :func:`stopper` converts :class:`~pyrump.stopping.registry.StoppingRegistry`'s
  eV/(1e15 atoms/cm^2) convention to RUMP's raw eV*cm^2/atom via ``* 1e-15``,
  matching ``RbsStoper``'s own comment about where that factor is applied.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..atomic.tables import PeriodicTable
from ..model.geometry import Geometry
from ..physics.kinematics import kinematic_factor
from ..physics.xsec.rutherford import setup_scatter

#: RUMP's own literal constant (particles/uC), anlytc.c:274/321 -- kept as the
#: rounded value the C uses, not a more precise recomputation, so output
#: matches real RUMP's rather than merely being "more correct."
_PARTICLES_PER_UC_MATRIX = 6.25e12

#: 1 barn = 1e-24 cm^2. RbsSigma's raw C return value is absolute cm^2/sr;
#: cross_section_barns() reports barns (matching the "1E-24 cm2/sr" units
#: RUMP itself prints), so callers doing raw physics with it -- not just
#: printing it -- must convert back with this factor first.
_BARN_CM2 = 1e-24


def kappa(m1: float, m2: float, scattering_angle_deg: float) -> float:
    """``RbsKappa`` (anlytc.c:1240): the kinematic factor, or ``0.0`` if
    scattering cannot occur.

    :func:`~pyrump.physics.kinematics.kinematic_factor` raises past the
    kinematic limit; the C returns 0 instead, so translate that here rather
    than propagating the exception into every caller.
    """
    try:
        return kinematic_factor(m1, m2, scattering_angle_deg)
    except ValueError:
        return 0.0


def cosines(geometry: Geometry) -> tuple[float, float]:
    """``(cosin, cosout)`` from ``RbsSetexp`` (anlytc.c:1182).

    Raises :class:`ValueError` for a degenerate geometry (an undefined path
    length), matching ``RbsSetexp``'s own rejection of ``zbeam == 0``-style
    unusable setups -- callers should turn this into a ``CommandError``.
    """
    if geometry.sec_in == 0.0 or geometry.sec_out == 0.0:
        raise ValueError("geometry is degenerate: an incident or exit path is undefined")
    return 1.0 / geometry.sec_in, 1.0 / geometry.sec_out


def stopper(registry, z1: int, m1: float, z2: int, energy_MeV: float) -> float:
    """``RbsStoper`` (anlytc.c:1328): electronic+nuclear stopping, eV*cm^2/atom."""
    return float(registry(z1, m1, z2, energy_MeV * 1000.0).values[0]) * 1e-15


def epsilon(
    registry, z1: int, m1: float, z2: int, e0_MeV: float, k: float, cosin: float, cosout: float
) -> float:
    """``RbsEpsilon`` (anlytc.c:1362): combined inbound+outbound stopping factor.

    ``[e] = K*stopper(e0)/cosin + stopper(K*e0)/cosout``, eV*cm^2/atom.
    """
    return (
        k * stopper(registry, z1, m1, z2, e0_MeV) / cosin
        + stopper(registry, z1, m1, z2, e0_MeV * k) / cosout
    )


def cross_section_barns(
    z1: int, m1: float, z2: int, mass: float, scattering_angle_deg: float, energy_keV: float
) -> float:
    """``RbsSigma`` (anlytc.c's own, local, UNSCREENED) -- barns/sr, == 1e-24 cm^2/sr."""
    cross_section = setup_scatter(z1, m1, z2, mass, scattering_angle_deg, screening=False)
    return float(cross_section(energy_keV)[0])


def resolve_element(table: PeriodicTable, token: str) -> tuple[object, float, int]:
    """``(Element, real_mass, mass_number)`` for a RUMP element/isotope token.

    ``RbsQueryElement`` (anlytc.c:1123), minus the "no change, reuse the last
    element" default -- the shell has no re-prompt loop, so a missing/unknown
    token is always an error here.
    """
    ref = table.parse_ref(token)
    element = table.by_z(ref.z)
    mass = table.real_mass(ref.z, ref.mass_number)
    return element, mass, ref.mass_number


@dataclass(frozen=True, slots=True)
class MatrixResult:
    """MATRIX/INFO's shared report: an element's expected edge and height.

    All the scattering-dependent fields are ``None`` when ``k == 0.0``
    (kinematically forbidden) or, for ``cross_section``/``epsilon``/``height``,
    when the cross section itself is non-positive.
    """

    symbol: str
    z: int
    mass: float
    mass_number: int
    k: float
    energy_keV: float | None
    channel: float | None
    cross_section_barns: float | None
    epsilon_eVcm2: float | None
    height: float | None
    """counts/uC/keV/msr."""


def matrix_result(buffer, table: PeriodicTable, registry, token: str) -> MatrixResult:
    """MATRIX and INFO's shared computation (anlytc.c:258-275, 290-345)."""
    element, mass, mass_number = resolve_element(table, token)
    angle = buffer.geometry.scattering_angle
    z1, m1, e0 = buffer.beam.z, buffer.beam.mass, buffer.beam.e0_MeV

    k = kappa(m1, mass, angle)
    if k == 0.0:
        return MatrixResult(element.symbol, element.z, mass, mass_number, 0.0,
                             None, None, None, None, None)

    energy_keV = k * e0 * 1000.0
    channel = float(buffer.calibration.channel_of(energy_keV))

    sigma = cross_section_barns(z1, m1, element.z, mass, angle, e0 * 1000.0)
    if sigma <= 0.0:
        return MatrixResult(element.symbol, element.z, mass, mass_number, k,
                             energy_keV, channel, None, None, None)

    cosin, cosout = cosines(buffer.geometry)
    eps = epsilon(registry, z1, m1, element.z, e0, k, cosin, cosout)
    height = (sigma * _BARN_CM2) * _PARTICLES_PER_UC_MATRIX / cosin / eps
    return MatrixResult(element.symbol, element.z, mass, mass_number, k,
                         energy_keV, channel, sigma, eps, height)


@dataclass(frozen=True, slots=True)
class Candidate:
    """One WHATISIT match: an element and its predicted edge."""

    symbol: str
    z: int
    energy_keV: float
    channel: float


def locate_candidates(buffer, table: PeriodicTable, target_keV: float, n_neighbors: int = 2) -> list[Candidate]:
    """``RbsLocate`` (anlytc.c:1398): the element(s) whose surface edge is
    nearest ``target_keV``, plus ``n_neighbors`` on each side by Z.

    Elements with no natural isotopes (no defined mass) or a forbidden
    kinematic factor are skipped, matching the C's own scan.
    """
    angle = buffer.geometry.scattering_angle
    m1, e0 = buffer.beam.mass, buffer.beam.e0_MeV

    energies: dict[int, float] = {}
    for element in table:
        if not element.isotopes:
            continue
        k = kappa(m1, element.mass, angle)
        if k == 0.0:
            continue
        energies[element.z] = k * e0 * 1000.0

    if not energies:
        raise ValueError("no element gives a valid surface edge for this beam/geometry")

    best_z = min(energies, key=lambda z: abs(energies[z] - target_keV))
    zs = sorted(
        z for z in energies
        if best_z - n_neighbors <= z <= best_z + n_neighbors
    )
    return [
        Candidate(
            symbol=table.by_z(z).symbol, z=z, energy_keV=energies[z],
            channel=float(buffer.calibration.channel_of(energies[z])),
        )
        for z in zs
    ]
