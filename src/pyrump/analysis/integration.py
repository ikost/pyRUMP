"""INTEGRAL, THICKNESS and INTSET: RUMP's shared region-integration algorithm.

Ports ``RbsThickn`` (anlytc.c:1493), which all three commands call with a
different ``ourkey`` -- ``TH_INT`` (INTEGRAL, stop after gross/net), ``TH_THK``
(THICKNESS, continue into the surface/compensated thickness conversion), and
``TH_SET`` (INTSET, which only ever touches the two persistent mode flags and
never reaches this module at all -- see ``Session.integration_interp``/
``integration_qmode``).

Region arguments here are plain 0-based channel *indices* (matching
:meth:`~pyrump.model.spectrum.Calibration.channel_of`'s convention, already
used throughout the shell), not RUMP's own ``first``-relative "channel
number" -- the same simplification the existing ``DISPLAY``/``INTEGRAL``
commands already made.

Always reports normalized-yield units (``#/uC/msr``), matching the C's
non-``raw`` default -- pyRUMP's buffer model has no raw/normalized toggle to
key off of, and RUMP itself ships with ``raw=FALSE``.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..model.detector import norm_k
from .elements import cosines, cross_section_barns, kappa, resolve_element, stopper

#: RbsThickn's own literal constant (anlytc.c:1630), distinct from MATRIX/
#: INFO's 6.25e12 -- `gross` here has already had a `kevch` factor absorbed
#: via NormK, consuming three of the twelve powers of ten. Both constants are
#: correct for their own callers; kept separate deliberately.
_PARTICLES_PER_UC_THICKNESS = 6.25e9

#: Default target density (atoms/cm^3) when an element has none on record,
#: matching RbsThickn's own silent fallback (anlytc.c:1625).
_DEFAULT_DENSITY = 1.0e22

#: 1 barn = 1e-24 cm^2. cross_section_barns() reports barns; RbsThickn's raw
#: C formula needs absolute cm^2/sr, so convert before using it here.
_BARN_CM2 = 1e-24


def _integrate_boundary(y0: float, y1: float, dx: float) -> tuple[float, float]:
    """``RbsIntegrate``/``INTINT`` (anlytc.c:1739): the triangular-basis
    endpoint correction for interpolated-mode integration.

    Returns ``(interpolated_value_at_dx, gross_correction)`` -- the caller
    adds the correction onto its running ``gross`` sum.
    """
    dx1 = 1.0 - dx
    y = dx1 * y0 + dx * y1
    correction = 0.5 * (-dx1 * y0 + dx * y)
    return y, correction


def _quad(t: float) -> float:
    """RbsThickn's taper window, ``quad(x) = min(1, (1-x)**2)`` (anlytc.c:1510)."""
    return min(1.0, (1.0 - t) ** 2)


@dataclass(frozen=True, slots=True)
class ThicknessResult:
    gross_atoms_cm2: float
    gross_angstrom: float
    net_atoms_cm2: float
    net_angstrom: float
    density_g_cc: float
    compensated: ThicknessResult | None = None
    """The qmode>0 second pass, using the Chu et al. energy-loss-ratio method."""


@dataclass(frozen=True, slots=True)
class IntegrationResult:
    lo_channel: float
    hi_channel: float
    """Possibly rounded, in discrete mode, to the actual channels used."""

    ist: int
    iend: int
    gross: float
    net: float
    """``#/uC/msr``."""

    thickness: ThicknessResult | None = None


def _region_indices(npt: int, lo: float, hi: float, interp: bool):
    """Resolve real-valued region bounds to (ist, iend) plus the boundary
    fractions ``dxst``/``dxend`` interpolated mode needs, and the (possibly
    rounded) reported channel bounds -- ``RbsThickn``'s region-setup block
    (anlytc.c:1543-1571), in 0-based index space.
    """
    if interp:
        xst = max(lo, 0.0)
        xend = min(hi, npt - 1)
        ist = min(npt - 1, max(0, round(xst + 0.5)))
        iend = min(npt - 1, max(0, round(xend - 0.5)))
        dxst = ist - xst
        dxend = xend - iend
        return ist, iend, xst, xend, dxst, dxend
    ist = min(npt - 1, max(0, round(lo)))
    iend = min(npt - 1, max(0, round(hi)))
    return ist, iend, float(ist), float(iend), 0.0, 0.0


def integrate_region(
    buffer,
    lo: float,
    hi: float,
    *,
    interp: bool,
    registry=None,
    table=None,
    target_token: str | None = None,
    qmode: int = 0,
    alpha_override: float | None = None,
) -> IntegrationResult:
    """``RbsThickn``'s ``TH_INT``/``TH_THK`` body.

    ``target_token=None`` stops after gross/net (``INTEGRAL``). Given, it adds
    the surface-approximation thickness and, if ``qmode>0``, the compensated
    second pass (``THICKNESS``) -- both need ``registry``/``table``.
    """
    if hi < lo:  # OrderPair (anlytc.c:1550)
        lo, hi = hi, lo

    counts = buffer.spectrum.counts
    npt = counts.size
    ist, iend, xst, xend, dxst, dxend = _region_indices(npt, lo, hi, interp)

    gross = float(counts[ist : iend + 1].sum())

    if interp:
        y_before_ist = counts[ist - 1] if ist > 0 else counts[ist]
        y_after_iend = counts[iend + 1] if iend < npt - 1 else counts[iend]
        yst, correction_st = _integrate_boundary(counts[ist], y_before_ist, dxst)
        yend, correction_end = _integrate_boundary(counts[iend], y_after_iend, dxend)
        gross += correction_st + correction_end
        net = gross - 0.5 * (yst + yend) * (xend - xst)
    else:
        # The "+1" is a documented 1984 bugfix (anlytc.c:1586) -- reproduced.
        net = gross - (counts[ist] + counts[iend]) * (iend + 1 - ist) / 2.0

    kevch = buffer.calibration.kevch
    normalization = norm_k(
        buffer.measurement.omega_msr, buffer.measurement.charge_uC, kevch,
        buffer.measurement.charge_state, buffer.measurement.correction,
    ) * kevch
    gross *= normalization
    net *= normalization

    result = IntegrationResult(
        lo_channel=xst, hi_channel=xend, ist=ist, iend=iend, gross=gross, net=net,
    )
    if target_token is None:
        return result

    if iend == ist:
        raise ValueError("thickness aborted: no region specified")

    element, mass, _ = resolve_element(table, target_token)
    z1, m1, e0 = buffer.beam.z, buffer.beam.mass, buffer.beam.e0_MeV
    angle = buffer.geometry.scattering_angle

    sigma = cross_section_barns(z1, m1, element.z, mass, angle, e0 * 1000.0) * _BARN_CM2
    if sigma <= 0.0:
        raise ValueError(f"conventional backscattering does not occur for {element.symbol}")

    density = element.atomic_density if element.atomic_density > 0 else _DEFAULT_DENSITY
    density_g_cc = density / 6.022e23 * mass

    x = gross / sigma / _PARTICLES_PER_UC_THICKNESS * cosines(buffer.geometry)[0]
    x1 = x / density * 1e8
    x2 = x * net / gross if gross else 0.0
    x3 = x1 * net / gross if gross else 0.0
    thickness = ThicknessResult(
        gross_atoms_cm2=x, gross_angstrom=x1, net_atoms_cm2=x2, net_angstrom=x3,
        density_g_cc=density_g_cc,
    )

    if qmode == 0:
        return IntegrationResult(
            lo_channel=xst, hi_channel=xend, ist=ist, iend=iend,
            gross=gross, net=net, thickness=thickness,
        )

    # Compensated pass: Chu et al.'s energy-loss-ratio method (p.65).
    k = kappa(m1, mass, angle)
    cosin, cosout = cosines(buffer.geometry)
    if qmode == 2:
        if alpha_override is None:
            raise ValueError("INTSET QUERY mode needs an explicit alpha value")
        alpha = alpha_override
    else:
        alpha = stopper(registry, z1, m1, element.z, k * e0) \
            / stopper(registry, z1, m1, element.z, e0) * cosin / cosout

    # The "+1" here is anlytc.c's own (anlytc.c:1664); reproduced as given.
    surfi = float(buffer.calibration.channel_of(k * e0 * 1000.0)) + 1.0

    slope_per_channel = (kevch / 1000.0) / e0 / (alpha + k)
    slope = (counts[iend] - counts[ist]) / (iend - ist)
    intercept = counts[ist] - ist * slope

    comp_gross = 0.0
    comp_net = 0.0
    for i in range(ist, iend + 1):
        weight = _quad(slope_per_channel * (surfi - i))
        comp_gross += counts[i] * weight
        comp_net += (counts[i] - slope * i - intercept) * weight

    if interp:
        weight_ist = _quad(slope_per_channel * (surfi - ist))
        weight_ist_prev = _quad(slope_per_channel * (surfi - ist + 1))
        weight_iend = _quad(slope_per_channel * (surfi - iend))
        weight_iend_next = _quad(slope_per_channel * (surfi - iend - 1))
        y_before_ist = counts[ist - 1] if ist > 0 else counts[ist]
        y_after_iend = counts[iend + 1] if iend < npt - 1 else counts[iend]
        yst, correction_st = _integrate_boundary(
            counts[ist] * weight_ist, y_before_ist * weight_ist_prev, dxst
        )
        yend, correction_end = _integrate_boundary(
            counts[iend] * weight_iend, y_after_iend * weight_iend_next, dxend
        )
        comp_gross += correction_st + correction_end
        comp_net = comp_gross - 0.5 * (yst + yend) * (xend - xst)

    comp_gross *= normalization
    comp_net *= normalization

    xg = comp_gross / sigma / _PARTICLES_PER_UC_THICKNESS * cosin
    xg1 = xg / density * 1e8
    xg2 = xg * comp_net / comp_gross if comp_gross else 0.0
    xg3 = xg1 * comp_net / comp_gross if comp_gross else 0.0
    compensated = ThicknessResult(
        gross_atoms_cm2=xg, gross_angstrom=xg1, net_atoms_cm2=xg2, net_angstrom=xg3,
        density_g_cc=density_g_cc,
    )
    thickness = ThicknessResult(
        gross_atoms_cm2=x, gross_angstrom=x1, net_atoms_cm2=x2, net_angstrom=x3,
        density_g_cc=density_g_cc, compensated=compensated,
    )
    return IntegrationResult(
        lo_channel=xst, hi_channel=xend, ist=ist, iend=iend,
        gross=gross, net=net, thickness=thickness,
    )
