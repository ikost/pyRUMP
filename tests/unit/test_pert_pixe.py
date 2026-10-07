"""PERT fitting the PIXE spectrum together with the RBS one.

The acceptance tests are films whose elements RBS can hardly tell apart --
W-Ta (masses 184 and 181) and permalloy Ni-Fe -- simulated in RBS and PIXE
with a known ratio, given Poisson noise, and fitted from a wrong start, the
PIXE instrumental constant H included. With PIXE windows over the lines the
ratio must come back, and more tightly than from RBS alone.
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
from pyrump.pixe.data import PixeData  # noqa: E402
from pyrump.pixe.detector import DEFAULT_CALIBRATION  # noqa: E402
from pyrump.script.lcm import to_sample  # noqa: E402
from pyrump.shell import pixe_plotting, pixe_sim  # noqa: E402
from pyrump.shell.dispatch import CommandError  # noqa: E402
from pyrump.shell.repl import execute_line  # noqa: E402
from pyrump.shell.session import Buffer, Session  # noqa: E402
from pyrump.sim.engine import Beam  # noqa: E402

from conftest import data_dir

DATA = data_dir()
needs_data = pytest.mark.skipif(DATA is None, reason="legacy data tables unavailable")

FILM = """Sim Reset
Layer 1
 Thick 1000 /cm2
 Composition {first} 1 {second} {amount} /
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


def channel(energy_keV: float) -> int:
    """The PIXE channel number (as the .PIX file numbers them) of an energy."""
    c = DEFAULT_CALIBRATION
    return round((energy_keV - c.kev0) / c.kevch)


def film(tmp_path, first: str, second: str, *, truth=1.0, start=2.0, dose=1.0, seed=3):
    """A ``first 1 second <truth>`` film on Si, measured in RBS and PIXE,
    loaded with the sample at ``second <start>``. ``dose`` scales both
    measured spectra, as a charge reading that is off would."""
    session = Session.create(str(DATA))
    buffer = Buffer(
        spectrum=Spectrum(counts=np.zeros(512), calibration=Calibration(kevch=4.0, npt=512)),
        beam=Beam(e0_MeV=2.0), geometry=Geometry(theta=0.0, phi=10.0, kind=GeometryKind.CORNELL),
        measurement=Measurement(omega_msr=2.0, charge_uC=10.0, fwhm_keV=15.0), name="film",
    )
    buffer.pixe = PixeData(spectrum=Spectrum(
        counts=np.zeros(DEFAULT_CALIBRATION.npt), calibration=DEFAULT_CALIBRATION,
    ))
    session.buffers.load(buffer, 1)
    session.buffers.active = 1
    for name, amount in (("truth", truth), ("start", start)):
        (tmp_path / f"{name}.lcm").write_text(
            FILM.format(first=first, second=second, amount=amount)
        )
    run(session, f"sim get {tmp_path / 'truth.lcm'}")
    rbs = session.simulation().spectrum.counts
    sample = to_sample(session.script, session.table, session.densities)
    pixe = pixe_sim.simulate_with(
        session, sample, buffer.beam, buffer.geometry, buffer.measurement,
        (1.0, 1.0, 1.0), buffer.pixe,
    ).counts
    rng = np.random.default_rng(seed)
    buffer.spectrum.counts = rng.poisson(np.clip(rbs * dose, 0, None)).astype(float)
    buffer.pixe.spectrum.counts = rng.poisson(np.clip(pixe * dose, 0, None)).astype(float)
    run(session, f"sim get {tmp_path / 'start.lcm'}")
    return session


def fitted(output: str, name: str) -> tuple[float, float]:
    """``(value, uncertainty)`` of one parameter from GO's printout."""
    line = next(line for line in output.splitlines() if line.strip().startswith(name))
    parts = line.split()
    value = float(parts[len(name.split())])
    sigma = float(parts[parts.index("+/-") + 1])
    return value, sigma


# -- the acceptance tests -----------------------------------------------------


@needs_data
def test_w_ta_ratio_from_rbs_and_pixe(tmp_path):
    rbs_only = film(tmp_path, "W", "Ta")
    output = run(rbs_only, "pert window 380 500", "pert composition 1 Ta",
                 "pert thick 1", "pert go")
    _, sigma_rbs = fitted(output, "layer 1 composition Ta")

    both = film(tmp_path, "W", "Ta")
    both.pixe.h = (1.0, 1.25, 1.0)  # H_L off by 25%: the fit must find it
    output = run(
        both, "pert window 380 500", "pert composition 1 Ta", "pert thick 1",
        f"pert pixwin {channel(7.9)} {channel(10.0)}", "pert pixh l", "pert go",
    )
    ta, sigma = fitted(output, "layer 1 composition Ta")
    h, sigma_h = fitted(output, "PIXH L")
    assert abs(ta - 1.0) < 3 * sigma, (ta, sigma)
    assert abs(h - 1.0) < 3 * sigma_h, (h, sigma_h)
    assert sigma < sigma_rbs / 2, (sigma, sigma_rbs)
    assert "PIXE chi-square" in output and "RBS chi-square" in output
    assert both.pixe.h[1] == pytest.approx(h, rel=1e-4)  # written back


@needs_data
def test_permalloy_ni_fe_ratio_from_rbs_and_pixe(tmp_path):
    session = film(tmp_path, "Ni", "Fe", truth=0.25, start=0.6)
    session.pixe.h = (1.3, 1.0, 1.0)
    output = run(
        session, "pert window 300 420", "pert composition 1 Fe", "pert thick 1",
        f"pert pixwin {channel(6.2)} {channel(8.5)}", "pert pixh k", "pert go",
    )
    fe, sigma = fitted(output, "layer 1 composition Fe")
    assert abs(fe - 0.25) < 3 * sigma, (fe, sigma)


@needs_data
def test_a_normalisation_window_sets_the_pixe_dose_too(tmp_path):
    """Both spectra read 20% more ions than the charge says. The RBS
    normalisation window absorbs that; PIXE must get the same correction,
    or H would take the 20% instead."""
    session = film(tmp_path, "W", "Ta", dose=1.2)
    output = run(
        session, "pert window 380 500", "pert normalize 100 300",
        "pert composition 1 Ta", "pert thick 1",
        f"pert pixwin {channel(7.9)} {channel(10.0)}", "pert pixh l", "pert go",
    )
    h, sigma_h = fitted(output, "PIXH L")
    assert abs(h - 1.0) < 3 * sigma_h + 0.02, (h, sigma_h)


# -- the commands -------------------------------------------------------------


@needs_data
def test_pixe_windows_are_added_listed_and_cleared(tmp_path):
    session = film(tmp_path, "W", "Ta")
    run(session, "pert pixwin 786 900", "pert pixwin 920 994")
    assert [(w.low, w.high) for w in session.pert.pixe_windows] == [(786, 900), (920, 994)]
    assert "PIXE win    [1] 786-900  [2] 920-994" in run(session, "pert parms")
    run(session, "pert pixwin clear 1")
    assert [(w.low, w.high) for w in session.pert.pixe_windows] == [(920, 994)]
    run(session, "pert pixwin clear")
    assert session.pert.pixe_windows == []
    assert "(none -- RBS only)" in run(session, "pert pixwin")


@needs_data
def test_pixe_windows_and_h_round_trip_through_a_pert_file(tmp_path):
    session = film(tmp_path, "W", "Ta")
    run(session, "pert window 380 500", "pert pixwin 786 994", "pert pixh l 0.5 2",
        f"pert save {tmp_path / 'setup'}")
    text = (tmp_path / "setup.pert").read_text()
    assert "pixwin 786 994" in text and "pixh L 0.5 2" in text
    again = film(tmp_path, "W", "Ta")
    run(again, f"pert get {tmp_path / 'setup'}")
    assert [(w.low, w.high) for w in again.pert.pixe_windows] == [(786, 994)]
    assert [(v.name, v.bounds) for v in again.pert.varying] == [("PIXH L", (0.5, 2.0))]


@needs_data
def test_pixe_h_takes_k_l_or_m(tmp_path):
    session = film(tmp_path, "W", "Ta")
    with pytest.raises(CommandError, match="K, L or M"):
        run(session, "pert pixh x")


@needs_data
def test_go_needs_a_pixe_spectrum_for_pixe_windows(tmp_path):
    session = film(tmp_path, "W", "Ta")
    session.buffers.active_buffer.pixe = None
    with pytest.raises(CommandError, match="no PIXE spectrum"):
        run(session, "pert pixwin 786 994", "pert composition 1 Ta", "pert go")


@needs_data
def test_go_needs_pixe_windows_to_vary_h(tmp_path):
    session = film(tmp_path, "W", "Ta")
    with pytest.raises(CommandError, match="no PIXE windows"):
        run(session, "pert composition 1 Ta", "pert pixh l", "pert go")


@needs_data
def test_go_rejects_a_pixe_window_outside_the_spectrum(tmp_path):
    session = film(tmp_path, "W", "Ta")
    with pytest.raises(CommandError, match="outside the PIXE spectrum"):
        run(session, "pert pixwin 5000 6000", "pert composition 1 Ta", "pert go")


@needs_data
def test_pixe_compare_scores_the_pert_windows(tmp_path):
    session = film(tmp_path, "W", "Ta")
    run(session, "pert pixwin 786 900")
    mask = pixe_plotting.pert_windows(session, DEFAULT_CALIBRATION, DEFAULT_CALIBRATION.npt)
    first = round(DEFAULT_CALIBRATION.first)
    assert mask.sum() == 115 and mask[786 - first] and not mask[785 - first]


@needs_data
def test_pixe_in_pert_opens_the_pixe_prompt(tmp_path):
    """PIXWIN and PIXH need four letters, so PIX and PIXE reach the RUMP
    level's PIXE, as at every other prompt; the PERT setup stays."""
    session = film(tmp_path, "W", "Ta")
    stack = ["rump"]
    with contextlib.redirect_stdout(io.StringIO()):
        for line in ("pert", "pixw 786 994", "pixh l", "pixe"):
            execute_line(session, line, stack)
    assert stack == ["rump", "pixe"]
    assert [(w.low, w.high) for w in session.pert.pixe_windows] == [(786, 994)]


# -- PIXFIRST -----------------------------------------------------------------


@needs_data
def test_pixfirst_takes_the_ratio_from_pixe_when_rbs_is_off(tmp_path):
    """The RBS data two channels off the simulation's calibration -- a
    systematic misfit the counts can't average out. Fitted jointly, the RBS
    spectrum's far greater counts pull the Fe content with it; PIXFIRST
    takes it from the X-ray lines alone, and the thickness from RBS."""
    lines = ("pert window 300 420", "pert composition 1 Fe", "pert thick 1",
             f"pert pixwin {channel(6.2)} {channel(8.5)}", "pert pixh k")
    joint = film(tmp_path, "Ni", "Fe", truth=0.25, start=0.6)
    rbs = joint.buffers.active_buffer.spectrum
    rbs.counts = np.roll(rbs.counts, 2)
    fe_joint, _ = fitted(run(joint, *lines, "pert go"), "layer 1 composition Fe")

    first = film(tmp_path, "Ni", "Fe", truth=0.25, start=0.6)
    rbs = first.buffers.active_buffer.spectrum
    rbs.counts = np.roll(rbs.counts, 2)
    output = run(first, *lines, "pert pixfirst", "pert go")
    fe, sigma = fitted(output, "layer 1 composition Fe")
    assert abs(fe - 0.25) < 3 * sigma, (fe, sigma)
    assert abs(fe - 0.25) < abs(fe_joint - 0.25) / 3, (fe, fe_joint)
    assert "PIXE: layer 1 composition Fe, PIXH K" in output
    assert "RBS: layer 1 thickness" in output
    assert "PIXE: reduced chi-square" in output and "RBS: reduced chi-square" in output


@needs_data
def test_pixfirst_round_trips_and_is_listed(tmp_path):
    session = film(tmp_path, "W", "Ta")
    run(session, "pert window 380 500", "pert pixwin 786 994", "pert pixfirst",
        f"pert save {tmp_path / 'setup'}")
    assert "pixfirst" in (tmp_path / "setup.pert").read_text().split()
    again = film(tmp_path, "W", "Ta")
    run(again, f"pert get {tmp_path / 'setup'}")
    assert again.pert.pixe_first
    assert "PIXE fit    first" in run(again, "pert parms")
    run(again, "pert pixfirst off")
    assert not again.pert.pixe_first
    assert "PIXE fit    joint" in run(again, "pert parms")


@needs_data
def test_pixfirst_needs_pixe_windows(tmp_path):
    session = film(tmp_path, "W", "Ta")
    with pytest.raises(CommandError, match="no PIXE windows"):
        run(session, "pert composition 1 Ta", "pert pixfirst", "pert go")


@needs_data
def test_pixfirst_needs_something_for_pixe_to_fit(tmp_path):
    session = film(tmp_path, "W", "Ta")
    with pytest.raises(CommandError, match="nothing for the PIXE spectrum"):
        run(session, "pert pixwin 786 994", "pert thick 1", "pert pixfirst", "pert go")


@needs_data
def test_pixfirst_refuses_an_element_without_lines_in_the_windows(tmp_path):
    """Windows over the W and Ta L lines: Si's K lines (1.74 keV) lie far
    outside, so the PIXE spectrum can't fit Si -- as with O in an oxide."""
    session = film(tmp_path, "W", "Ta")
    with pytest.raises(CommandError, match="Si has no simulated counts"):
        run(session, "pert pixwin 786 994", "pert composition 2 Si", "pert pixfirst", "pert go")
