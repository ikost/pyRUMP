"""The fit stops when chi-square stops improving -- not when scipy judges a
step small next to the whole parameter vector.

A thickness in 1e15 at/cm^2 (~1000) varied with a composition (~1) made
every composition step look negligible, and the search stopped two steps in,
far from the minimum. RUMP stops on the relative chi-square change alone
(EpsCrit), which is what this checks: a W-Ta film started at twice the true
Ta must come back to it.
"""

from __future__ import annotations

import contextlib
import io

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from pyrump.model.detector import Measurement  # noqa: E402
from pyrump.model.geometry import Geometry, GeometryKind  # noqa: E402
from pyrump.model.spectrum import Calibration, Spectrum  # noqa: E402
from pyrump.shell.repl import execute_line  # noqa: E402
from pyrump.shell.session import Buffer, Session  # noqa: E402
from pyrump.sim.engine import Beam  # noqa: E402

from conftest import data_dir

DATA = data_dir()
needs_data = pytest.mark.skipif(DATA is None, reason="legacy data tables unavailable")

FILM = """Sim Reset
Layer 1
 Thick 1000 /cm2
 Composition W 1 Ta {ta} /
Next
 Thick 20000 /cm2
 Composition Si 1 /
"""


def run(session: Session, *lines: str) -> str:
    stack, out = ["rump"], io.StringIO()
    with contextlib.redirect_stdout(out):
        for line in lines:
            execute_line(session, line, stack)
    return out.getvalue()


@needs_data
def test_a_composition_varied_with_a_thickness_reaches_the_minimum(tmp_path):
    session = Session.create(str(DATA))
    session.buffers.load(Buffer(
        spectrum=Spectrum(counts=np.zeros(512), calibration=Calibration(kevch=4.0, npt=512)),
        beam=Beam(e0_MeV=2.0), geometry=Geometry(theta=0.0, phi=10.0, kind=GeometryKind.CORNELL),
        measurement=Measurement(omega_msr=2.0, charge_uC=10.0, fwhm_keV=15.0), name="wta",
    ), 1)
    session.buffers.active = 1
    (tmp_path / "truth.lcm").write_text(FILM.format(ta=1))
    (tmp_path / "guess.lcm").write_text(FILM.format(ta=2))
    run(session, f"sim get {tmp_path / 'truth.lcm'}")
    truth = session.simulation().spectrum.counts
    session.buffers[1].spectrum.counts = (
        np.random.default_rng(3).poisson(np.clip(truth, 0, None)).astype(float)
    )
    run(
        session, f"sim get {tmp_path / 'guess.lcm'}",
        "pert window 380 500", "pert composition 1 Ta", "pert thick 1", "pert go",
    )
    result = session.pert  # the fitted value is written back into the sample
    ta = session.script.layers[0].composition["Ta"]
    assert result is not None
    assert abs(ta - 1.0) < 0.15, ta  # it stopped at 1.15 before, from 2
