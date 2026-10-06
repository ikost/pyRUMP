"""Run the PIXE simulation for the shell: SIM sample in, PIXE spectrum out.

The sample is the one SIM edits and the RBS simulation uses -- all of it,
the substrate included with the thickness SIM gives it; the beam, tilt
and charge are the ACTIVE buffer's (one run, one charge integrator), the
live fraction its PIXE spectrum's; the detector, H values and escape
setting are the PIXE prompt's.
"""

from __future__ import annotations

from ..pixe.spectrum import PixeSpectrum, synthesize
from ..pixe.yields import Exposure, simulate_lines


def reference_buffer(session):
    """The buffer whose beam, tilt and charge the simulation uses: the
    ACTIVE data buffer, or the session defaults when there is none."""
    buffer = session.buffers.active_buffer
    if buffer is None or session.buffers.active == 0:
        return session.settings.experiment_defaults
    return buffer


def simulate(session) -> PixeSpectrum | None:
    """The PIXE simulation for the ACTIVE buffer, or ``None`` with no SIM
    sample. Raises ``ValueError`` when it can't be computed (an ion with no
    cross sections, an energy outside the tables ...)."""
    from ..script.lcm import to_sample

    if not session.script.layers:
        return None
    buffer = reference_buffer(session)
    sample = to_sample(session.script, session.table, session.densities)
    return simulate_with(
        session, sample, buffer.beam, buffer.geometry, buffer.measurement,
        session.pixe.h, getattr(buffer, "pixe", None),
    )


def simulate_with(session, sample, beam, geometry, measurement, h, data) -> PixeSpectrum:
    """The PIXE simulation for these inputs -- what :func:`simulate` runs
    for the ACTIVE buffer, and PERT's GO for each trial of a fit. ``data``,
    the measured PIXE spectrum (or ``None``), gives the live fraction and
    the channels; the detector and escape setting are the PIXE prompt's."""
    state = session.pixe
    exposure = Exposure(
        charge_uC=measurement.charge_uC,
        charge_state=measurement.charge_state,
        correction=measurement.correction,
        live_fraction=data.live_fraction if data is not None else 1.0,
        h=tuple(h),
    )
    lines = simulate_lines(
        sample, beam, geometry.theta, state.detector, exposure,
        session.registry, session.table, include_substrate=True,
        faithful=session.settings.faithful,
    )
    calibration = data.calibration if data is not None else state.calibration
    return synthesize(lines, calibration, state.detector, escape=state.escape)
