"""Detector and measurement parameters, and the absolute yield normalisation."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Measurement:
    """The quantities that set a spectrum's absolute scale."""

    omega_msr: float = 1.0
    """Detector solid angle, msr."""

    charge_uC: float = 10.0
    """Integrated beam charge."""

    correction: float = 1.0
    """RUMP's ``CORR`` factor -- the fudge that absorbs charge-integration error."""

    charge_state: int = 1
    """``cbeam``, the |charge| state of the beam."""

    current_nA: float = 0.0
    """Average beam current; only used by the pileup model."""

    fwhm_keV: float = 0.0
    """Detector resolution."""

    tau_us: float = 5.0
    """MCA shaping time constant; only used by the pileup model."""


def norm_k(
    omega_msr: float, charge_uC: float, kevch: float, charge_state: int, correction: float
) -> float:
    """``RbsNormK`` (bmanip.c:880).

    Falls back to 1.0 rather than dividing by zero, as the C does (with a
    "Couldn't normalize buffer" warning).
    """
    denominator = omega_msr * charge_uC * kevch
    if charge_state != 0 and denominator != 0.0 and correction != 0.0:
        return charge_state * correction / denominator
    return 1.0


def yield_normalisation(measurement: Measurement) -> float:
    """Factor converting the raw brick fill into real counts.

    ``SimCreateDetails`` builds the theory spectrum with ``q = cbeam``,
    ``omega = 1`` and ``corr = 1`` (creatr.c:287-289), then rescales by
    ``RbsNormK(ALTBUF) / RbsNormK(ibuf)`` (creatr.c:321) to put the data
    buffer's charge, solid angle and correction back in.

    Both ``RbsNormK`` calls share ``kevch``, so it cancels and the whole thing
    collapses to::

        omega * q / (cbeam * corr)

    Deriving it this way rather than reproducing the two-step dance avoids the
    intermediate's dependence on the theory buffer's transient state.
    """
    return (
        measurement.omega_msr
        * measurement.charge_uC
        / (measurement.charge_state * measurement.correction)
    )
