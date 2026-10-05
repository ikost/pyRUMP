"""The PIXE sub-processor: turning PIXE on and off, PIXE data on buffers,
pairing with RBS files, the detector settings and the PIXE window.

Driven through :func:`pyrump.shell.repl.execute_line`, like the other shell
tests, so the real dispatch and mode stack are exercised.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from pyrump.model.spectrum import Calibration, Spectrum  # noqa: E402
from pyrump.pixe.detector import DEFAULT_CALIBRATION, PixeDetector  # noqa: E402
from pyrump.shell import plotting  # noqa: E402
from pyrump.shell.commands.pixe import setup_lines  # noqa: E402
from pyrump.shell.dispatch import CommandError  # noqa: E402
from pyrump.shell.repl import execute_file, execute_line  # noqa: E402
from pyrump.shell.session import Buffer, Session  # noqa: E402

from conftest import data_dir

DATA = data_dir()
EXAMPLES = Path(__file__).resolve().parents[2] / "examples"


def _session() -> Session:
    if DATA is None:
        pytest.skip("legacy data tables unavailable")
    return Session.create(str(DATA))


@pytest.fixture
def session() -> Session:
    """A session with one synthetic RBS buffer, active, in slot 1."""
    built = _session()
    counts = np.full(64, 10.0)
    built.buffers.load(
        Buffer(spectrum=Spectrum(counts=counts, calibration=Calibration(npt=64)), name="rbs"), 1
    )
    built.buffers.active = 1
    yield built
    import matplotlib.pyplot as plt

    plt.close("all")


@pytest.fixture
def empty_session() -> Session:
    built = _session()
    yield built
    import matplotlib.pyplot as plt

    plt.close("all")


def run(session: Session, *lines: str, stack: list[str] | None = None) -> list[str]:
    stack = stack if stack is not None else ["rump"]
    for line in lines:
        execute_line(session, line, stack)
    return stack


# -- on and off ---------------------------------------------------------------


def test_entering_turns_pixe_on_and_return_keeps_it_on(session):
    stack = run(session, "pixe")
    assert stack == ["rump", "pixe"] and session.pixe.enabled
    run(session, "return", stack=stack)
    assert stack == ["rump"] and session.pixe.enabled


def test_disable_turns_pixe_off_and_leaves(session):
    stack = run(session, f"pixe get {EXAMPLES / 'MnPt.PIX'}", "pixe")
    assert session.pixe.figure is not None
    run(session, "disable", stack=stack)
    assert stack == ["rump"]
    assert not session.pixe.enabled and session.pixe.figure is None
    assert session.buffers[1].pixe is not None  # the data stays


def test_one_shots_do_not_turn_pixe_on(empty_session):
    stack = run(empty_session, "pixe angle 30", "pixe disable")
    assert stack == ["rump"] and not empty_session.pixe.enabled
    assert empty_session.pixe.detector.angle_deg == 30


def test_pyrumprc_block_ending_in_disable_stores_the_setup_only(empty_session, tmp_path):
    rc = tmp_path / "rc"
    rc.write_text("pixe\n angle 30\n solid 2.5\n window Be 12.5\n pair on\ndisable\nmev 2.5\n")
    stack = ["rump"]
    execute_file(empty_session, rc, stack)
    assert stack == ["rump"] and not empty_session.pixe.enabled
    detector = empty_session.pixe.detector
    assert (detector.angle_deg, detector.solid_angle_msr) == (30, 2.5)
    assert (detector.window.element, detector.window.thickness_um) == ("Be", 12.5)
    assert empty_session.pixe.pair


# -- defaults and SHOW --------------------------------------------------------


def test_built_in_defaults_without_a_pyrumprc(empty_session):
    assert empty_session.pixe.detector == PixeDetector()
    assert empty_session.pixe.calibration == DEFAULT_CALIBRATION
    assert not empty_session.pixe.pair


def test_show_prints_commands_that_rebuild_the_setup(empty_session, capsys):
    run(empty_session, "pixe angle 30", "pixe filter Al 25 10", "pixe calib 0.0101 -0.04",
        "pixe pair on")
    lines = setup_lines(empty_session)
    fresh = _session()
    run(fresh, *(f"pixe {line}" for line in lines))
    assert fresh.pixe.detector == empty_session.pixe.detector
    assert fresh.pixe.calibration == empty_session.pixe.calibration
    assert fresh.pixe.pair
    capsys.readouterr()
    run(empty_session, "pixe show")
    shown = capsys.readouterr().out
    assert all(line in shown for line in lines)


@pytest.mark.parametrize(
    "line, message",
    [
        ("pixe angle 95", "between -90 and 90"),
        ("pixe solid 0", "positive"),
        ("pixe window Xx 8", "unknown element"),
        ("pixe crystal Si -1", "positive"),
        ("pixe filter Al 25 100", "hole area"),
        ("pixe calib -0.01 0", "gain must be positive"),
        ("pixe pair maybe", "ON or OFF"),
    ],
)
def test_bad_settings_are_refused(empty_session, line, message):
    before = (empty_session.pixe.detector, empty_session.pixe.calibration)
    with pytest.raises(CommandError, match=message):
        run(empty_session, line)
    assert (empty_session.pixe.detector, empty_session.pixe.calibration) == before


# -- data on buffers ----------------------------------------------------------


def test_get_puts_the_spectrum_on_the_active_buffer(session):
    run(session, f"pixe get {EXAMPLES / 'MnPt'}")  # .PIX is tried for a bare name
    data = session.buffers[1].pixe
    assert data.path.name == "MnPt.PIX"
    assert data.spectrum.counts.size == 2048
    assert "+ PIXE" in session.buffers.listing()
    assert "PIXE       MnPt.PIX" in session.buffers[1].describe()


def test_get_with_no_data_buffer_makes_a_pixe_only_buffer(empty_session):
    run(empty_session, "mev 2.5", "charge 7", f"pixe get {EXAMPLES / 'MnPt.PIX'}")
    buffer = empty_session.buffers.active_buffer
    assert empty_session.buffers.active == 1
    assert buffer.pixe is not None and buffer.spectrum.total() == 0
    assert buffer.beam.e0_MeV == 2.5 and buffer.measurement.charge_uC == 7


def test_calib_sets_the_default_and_the_active_spectrum(session):
    run(session, f"pixe get {EXAMPLES / 'MnPt.PIX'}", "pixe calib 0.01 -0.05")
    calibration = session.buffers[1].pixe.calibration
    assert (calibration.kevch, calibration.kev0, calibration.first) == (0.01, -0.05, 1.0)
    assert (session.pixe.calibration.kevch, session.pixe.calibration.kev0) == (0.01, -0.05)


def _pair_folder(tmp_path: Path, rbs_name: str, pix_name: str | None, source: str) -> Path:
    shutil.copy(EXAMPLES / source, tmp_path / rbs_name)
    if pix_name:
        shutil.copy(EXAMPLES / "MnPt.PIX", tmp_path / pix_name)
    return tmp_path / rbs_name


def test_pair_on_brings_the_pix_along_with_xeq(empty_session, tmp_path):
    macro = _pair_folder(tmp_path, "S1.RBS", "s1.pix", "MnPt.RBS")  # case differs
    run(empty_session, "pixe pair on", f"xeq {macro}")
    assert empty_session.buffers.active_buffer.pixe.path.name == "s1.pix"


def test_pair_on_brings_the_pix_along_with_get(empty_session, tmp_path):
    spectrum = _pair_folder(tmp_path, "2A.rbs", "2A.PIX", "2A.rbs")
    run(empty_session, "pixe pair on", f"get {spectrum}")
    assert empty_session.buffers.active_buffer.pixe.path.name == "2A.PIX"


def test_pair_off_or_no_companion_leaves_the_buffer_alone(empty_session, tmp_path, capsys):
    macro = _pair_folder(tmp_path, "S1.RBS", "S1.PIX", "MnPt.RBS")
    run(empty_session, f"xeq {macro}")
    assert empty_session.buffers.active_buffer.pixe is None
    folder = tmp_path / "no-pix"
    folder.mkdir()
    lonely = _pair_folder(folder, "S2.RBS", None, "MnPt.RBS")
    run(empty_session, "pixe pair on", f"xeq {lonely}")
    assert empty_session.buffers.active_buffer.pixe is None
    assert "PAIR: no S2.PIX" in capsys.readouterr().out


# -- the PIXE window ----------------------------------------------------------


def test_the_pixe_window_follows_rbs_redraws(session):
    run(session, f"pixe get {EXAMPLES / 'MnPt.PIX'}")
    run(session, "plot 1")
    assert session.pixe.figure is None  # PIXE is off
    run(session, "pixe", "return", "plot 1")
    assert session.pixe.figure is not None
    assert session.figure is not session.pixe.figure


def test_region_and_figsave(session, tmp_path):
    stack = run(session, f"pixe get {EXAMPLES / 'MnPt.PIX'}", "pixe", "region 12 0.5")
    assert (session.pixe.emin, session.pixe.emax) == (0.5, 12)
    left, right = session.pixe.figure.axes[0].get_xlim()
    assert left >= 0.4 and right <= 12.1
    run(session, f"figsave {tmp_path / 'pixe'}", stack=stack)
    assert (tmp_path / "pixe.png").stat().st_size > 0
    with pytest.raises(CommandError, match="no PIXE channels"):
        run(session, "region 30 40", stack=stack)


def test_plot_needs_pixe_data(session):
    stack = run(session, "pixe")
    with pytest.raises(CommandError, match="no PIXE spectrum"):
        run(session, "plot", stack=stack)


def test_pump_services_both_windows(session):
    run(session, f"pixe get {EXAMPLES / 'MnPt.PIX'}", "pixe", "return", "plot 1")
    plotting._last_pump = 0.0
    plotting.pump(session)  # Agg: nothing to flush, but must not trip over two figures


# -- line markers -------------------------------------------------------------


def _labels(session) -> list[str]:
    return [text.get_text() for text in session.pixe.figure.axes[0].texts]


def test_markers_label_the_sim_samples_lines(session):
    run(session, f"sim get {EXAMPLES / 'MnPt.lcm'}", f"pixe get {EXAMPLES / 'MnPt.PIX'}",
        "pixe", "region 0.5 12")
    labels = " ".join(_labels(session))
    for expected in ("Mn Kα", "Mn Kβ", "Pt Lα", "Ru Lα", "Si Kα"):
        assert expected in labels
    assert "Pt Mα/Mβ" in labels  # too close to tell apart: one label
    assert "Ru Ll" not in labels  # a minor line: only with MARKERS ALL


def test_markers_all_and_off(session):
    stack = run(session, f"sim get {EXAMPLES / 'MnPt.lcm'}",
                f"pixe get {EXAMPLES / 'MnPt.PIX'}", "pixe", "region 0.5 12", "markers all")
    assert "Ru Ll" in " ".join(_labels(session))
    run(session, "markers off", stack=stack)
    assert not any(label.startswith(("Mn", "Pt", "Ru", "Si")) for label in _labels(session))
    with pytest.raises(CommandError, match="ON, ALL or OFF"):
        run(session, "markers some", stack=stack)


def test_markers_skip_rbs_absorber_layers(empty_session):
    from pyrump.shell.pixe_plotting import sample_elements

    run(empty_session, "sim layer 1", "sim thick 1000 A", "sim composition Al 1 /",
        "sim next", "sim thick 500 A", "sim composition Ti 1 /", "sim absorber 1")
    assert sample_elements(empty_session) == [(22, "Ti")]


def test_crowded_labels_spread_apart_but_keep_their_order():
    from pyrump.shell.pixe_plotting import _spread

    wanted = [2.0, 2.01, 2.02, 5.0, 9.99]
    positions = _spread(wanted, 0.1, 0.0, 10.0)
    assert positions == sorted(positions)
    assert all(b - a >= 0.1 - 1e-12 for a, b in zip(positions, positions[1:]))
    assert positions[3] == 5.0  # an uncrowded label stays where its line is
    assert np.mean(positions[:3]) == pytest.approx(np.mean(wanted[:3]))
    assert positions[-1] <= 10.0
