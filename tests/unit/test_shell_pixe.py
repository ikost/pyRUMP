"""The PIXE sub-processor: turning PIXE on and off, PIXE data on buffers,
pairing with RBS files, the detector settings and the PIXE window.

Driven through :func:`pyrump.shell.repl.execute_line`, like the other shell
tests, so the real dispatch and mode stack are exercised.
"""

from __future__ import annotations

import shutil
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from pyrump.model.geometry import Geometry, GeometryKind  # noqa: E402
from pyrump.model.spectrum import Calibration, Spectrum  # noqa: E402
from pyrump.pixe.detector import (  # noqa: E402
    DEFAULT_CALIBRATION, DEFAULT_GEOMETRY, PixeDetector,
)
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
    stack = run(session, f"pixe get {EXAMPLES / 'MnPt.PIX'}", "pixe", "plot")
    assert session.pixe.figure is not None
    run(session, "disable", stack=stack)
    assert stack == ["rump"]
    assert not session.pixe.enabled and session.pixe.figure is None
    assert session.buffers[1].pixe is not None  # the data stay


def test_one_shots_do_not_turn_pixe_on(empty_session):
    stack = run(empty_session, "pixe phi 30", "pixe disable")
    assert stack == ["rump"] and not empty_session.pixe.enabled
    assert empty_session.pixe.geometry.phi == 30


def test_pyrumprc_block_ending_in_disable_stores_the_setup_only(empty_session, tmp_path):
    rc = tmp_path / "rc"
    rc.write_text("pixe\n phi 30\n solid 2.5\n window Be 12.5\n pair on\ndisable\nmev 2.5\n")
    stack = ["rump"]
    execute_file(empty_session, rc, stack)
    assert stack == ["rump"] and not empty_session.pixe.enabled
    detector = empty_session.pixe.detector
    assert (empty_session.pixe.geometry.phi, detector.solid_angle_msr) == (30, 2.5)
    assert (detector.window.material, detector.window.thickness_um) == ("Be", 12.5)
    assert empty_session.pixe.pair


# -- defaults and SHOW --------------------------------------------------------


def test_built_in_defaults_without_a_pyrumprc(empty_session):
    assert empty_session.pixe.detector == PixeDetector()
    assert empty_session.pixe.calibration == DEFAULT_CALIBRATION
    assert empty_session.pixe.geometry == DEFAULT_GEOMETRY
    assert not empty_session.pixe.theta_rbs
    assert not empty_session.pixe.pair


def test_show_prints_commands_that_rebuild_the_setup(empty_session, capsys):
    run(empty_session, "pixe geometry general", "pixe phi 30", "pixe psi 20", "pixe theta rbs",
        "pixe filter 2 Al 25 hole 10%", "pixe calib 0.0101 -0.04", "pixe pair on")
    lines = setup_lines(empty_session)
    fresh = _session()
    run(fresh, *(f"pixe {line}" for line in lines))
    # SHOW prints the solid angle to 6 digits; everything else exactly.
    rebuilt, original = fresh.pixe.detector, empty_session.pixe.detector
    assert rebuilt.solid_angle_msr == pytest.approx(original.solid_angle_msr, rel=1e-6)
    assert replace(rebuilt, solid_angle_msr=original.solid_angle_msr) == original
    assert fresh.pixe.calibration == empty_session.pixe.calibration
    assert fresh.pixe.geometry == empty_session.pixe.geometry
    assert fresh.pixe.theta_rbs
    assert fresh.pixe.pair
    capsys.readouterr()
    run(empty_session, "pixe show")
    shown = capsys.readouterr().out
    assert all(line in shown for line in lines)


@pytest.mark.parametrize(
    "line, message",
    [
        ("pixe geometry tilted", "IBM, CORNELL or GENERAL"),
        ("pixe theta 95", "between -90 and 90"),
        ("pixe theta up", "degrees or RBS"),
        ("pixe phi 180", "from 0 up to 180"),
        ("pixe psi -5", "from 0 up to 90"),
        ("pixe solid 0", "positive"),
        ("pixe window Xx 8", "unknown element"),
        ("pixe crystal Si -1", "positive"),
        ("pixe filter 2 Al 25 100", "hole area"),
        ("pixe filter Al 25", "the filter's number first"),
        ("pixe filter 3 Al 25", "set filter 2 first"),
        ("pixe filter clear 4", "no filter 4"),
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


def _pixe_labels(session) -> list[str]:
    return [line.get_label() for line in session.pixe.figure.axes[0].lines
            if not line.get_label().startswith("_")]


def _two_buffers(session, tmp_path):
    """Buffer 1 with MnPt's RBS + PIXE pair, buffer 2 with a second copy."""
    for name in ("a", "b"):
        shutil.copy(EXAMPLES / "MnPt.RBS", tmp_path / f"{name}.RBS")
        shutil.copy(EXAMPLES / "MnPt.PIX", tmp_path / f"{name}.PIX")
    run(session, "pixe pair on", f"xeq {tmp_path / 'a.RBS'}", f"xeq {tmp_path / 'b.RBS'}")


def test_loading_and_entering_open_no_window(empty_session):
    run(empty_session, f"pixe get {EXAMPLES / 'MnPt.PIX'}", "pixe", "return",
        "pixe fwhm 125", "pixe pair on", f"xeq {EXAMPLES / 'MnPt.RBS'}")
    assert empty_session.pixe.figure is None


def test_plot_shows_the_same_buffer_in_both_windows(session):
    run(session, f"pixe get {EXAMPLES / 'MnPt.PIX'}", "plot 1")
    assert session.pixe.figure is None  # PIXE is off
    run(session, "pixe", "return", "plot 1")
    assert session.pixe.figure is not None and session.figure is not session.pixe.figure
    assert len(_pixe_labels(session)) == 1


def test_overlay_adds_the_second_buffers_pixe(empty_session, tmp_path):
    _two_buffers(empty_session, tmp_path)
    run(empty_session, "pixe", "return", "plot 1", "overlay 2")
    assert len(_pixe_labels(empty_session)) == 2


def test_a_buffer_without_pixe_is_noted(session):
    run(session, "pixe", "return", "plot 1")
    notes = [text.get_text() for text in session.pixe.figure.axes[0].texts]
    assert any("buffer 1: no PIXE spectrum" in note for note in notes)


def test_plot_0_and_splot_show_the_simulation(session):
    run(session, f"sim get {EXAMPLES / 'MnPt.lcm'}", f"pixe get {EXAMPLES / 'MnPt.PIX'}",
        "pixe", "return", "plot 0")
    assert _pixe_labels(session) == ["simulation"]
    run(session, "plot 1", "sim splot Mn", "sim splot 2")
    labels = _pixe_labels(session)
    assert "SIM(Mn)" in labels and "SIM(layer 2)" in labels


def test_compare_draws_data_simulation_and_residuals(session):
    run(session, f"sim get {EXAMPLES / 'MnPt.lcm'}", f"pixe get {EXAMPLES / 'MnPt.PIX'}",
        "pixe", "return", "compare")
    figure = session.pixe.figure
    assert len(figure.axes) == 2  # spectrum and residuals
    assert "simulation" in _pixe_labels(session)
    # On the residuals panel, as in the RBS COMPARE.
    assert any("chi-square" in text.get_text() for text in figure.axes[1].texts)
    run(session, "plot 1")  # back to one panel
    assert len(session.pixe.figure.axes) == 1


def test_compare_inside_the_pixe_prompt_leaves_rbs_alone(session):
    stack = run(session, f"sim get {EXAMPLES / 'MnPt.lcm'}", f"pixe get {EXAMPLES / 'MnPt.PIX'}",
                "pixe", "cmp")
    assert len(session.pixe.figure.axes) == 2
    assert session.figure is None
    run(session, "plot", stack=stack)
    assert len(session.pixe.figure.axes) == 1


def test_settings_redraw_an_open_window_with_a_fresh_simulation(session):
    stack = run(session, f"sim get {EXAMPLES / 'MnPt.lcm'}", f"pixe get {EXAMPLES / 'MnPt.PIX'}",
                "pixe", "plot 0")
    before = np.max(_sim_curve(session).get_ydata())
    run(session, "solid 1.5263", stack=stack)  # twice the default
    after = np.max(_sim_curve(session).get_ydata())
    assert after == pytest.approx(2 * before, rel=1e-3)
    run(session, "return", "charge 20", "pixe plot 0")  # an RBS parameter: the cache goes too
    assert np.max(_sim_curve(session).get_ydata()) == pytest.approx(4 * before, rel=1e-3)



def test_a_sim_edit_redraws_the_pixe_window_on_its_own(session):
    """LIVE: the PIXE prompt's own PLOT 0, with no RBS window open."""
    run(session, f"sim get {EXAMPLES / 'MnPt.lcm'}", f"pixe get {EXAMPLES / 'MnPt.PIX'}",
        "pixe", "plot 0")
    before = _sim_curve(session).get_ydata().copy()
    run(session, "sim", "layer 2", "thick 600 A")
    assert session.figure is None
    assert not np.array_equal(_sim_curve(session).get_ydata(), before)


def test_a_sim_edit_redraws_both_windows(session):
    """LIVE: the RBS window redraws, and the PIXE window follows it."""
    run(session, f"sim get {EXAMPLES / 'MnPt.lcm'}", f"pixe get {EXAMPLES / 'MnPt.PIX'}",
        "pixe", "return", "plot 0")
    before = _sim_curve(session).get_ydata().copy()
    run(session, "sim", "layer 2", "thick 600 A")
    assert session.traces[0].buffer is session.buffers.get(0)
    assert not np.array_equal(_sim_curve(session).get_ydata(), before)

def test_region_is_in_channels_with_energy_below_and_channels_on_top(session):
    stack = run(session, f"pixe get {EXAMPLES / 'MnPt.PIX'}", "pixe", "plot", "region 1200 50")
    assert (session.pixe.low, session.pixe.high) == (50, 1200)
    calibration = session.buffers[1].pixe.calibration
    left, right = session.pixe.figure.axes[0].get_xlim()
    assert left == pytest.approx(50 * calibration.kevch + calibration.kev0)
    assert right == pytest.approx(1201 * calibration.kevch + calibration.kev0)
    assert session.pixe.figure.axes[0].child_axes  # the channel axis on top
    with pytest.raises(CommandError, match="no PIXE channels"):
        run(session, "region 5000 6000", stack=stack)


def test_figsave_saves_both_windows_in_pair_mode(session, tmp_path):
    run(session, f"pixe get {EXAMPLES / 'MnPt.PIX'}", "pixe", "return", "plot 1")
    run(session, f"figsave {tmp_path / 'one'}")
    assert (tmp_path / "one_rbs.png").exists() and not (tmp_path / "one_pixe.png").exists()
    run(session, "pixe pair on", f"figsave {tmp_path / 'two'}")
    assert (tmp_path / "two_rbs.png").exists() and (tmp_path / "two_pixe.png").exists()
    stack = run(session, "pixe")
    run(session, f"figsave {tmp_path / 'three'}", stack=stack)
    assert (tmp_path / "three.png").exists() and not (tmp_path / "three_pixe.png").exists()
    assert not (tmp_path / "three_rbs.png").exists()


def test_plot_in_the_pixe_prompt_needs_something_to_show(empty_session):
    stack = run(empty_session, "pixe")
    with pytest.raises(CommandError, match="no active"):
        run(empty_session, "cmp", stack=stack)


def test_pump_services_both_windows(session):
    run(session, f"pixe get {EXAMPLES / 'MnPt.PIX'}", "pixe", "return", "plot 1")
    plotting._last_pump = 0.0
    plotting.pump(session)  # Agg: nothing to flush, but must not trip over two figures


# -- line markers -------------------------------------------------------------


def _labels(session) -> list[str]:
    return [text.get_text() for text in session.pixe.figure.axes[0].texts]


def test_markers_label_the_sim_samples_lines(session):
    run(session, f"sim get {EXAMPLES / 'MnPt.lcm'}", f"pixe get {EXAMPLES / 'MnPt.PIX'}",
        "pixe", "plot", "region 50 1190")
    labels = " ".join(_labels(session))
    for expected in ("Mn Kα", "Mn Kβ", "Pt Lα", "Ru Lα", "Si Kα"):
        assert expected in labels
    assert "Pt Mα/Mβ" in labels  # too close to tell apart: one label
    assert "Ru Ll" not in labels  # a minor line: only with MARKERS ALL


def test_markers_all_and_off(session):
    stack = run(session, f"sim get {EXAMPLES / 'MnPt.lcm'}",
                f"pixe get {EXAMPLES / 'MnPt.PIX'}", "pixe", "plot", "region 50 1190",
                "markers all")
    assert "Ru Ll" in " ".join(_labels(session))
    run(session, "markers off", stack=stack)
    assert not any(label.startswith(("Mn", "Pt", "Ru", "Si")) for label in _labels(session))
    with pytest.raises(CommandError, match="ON, ALL or OFF"):
        run(session, "markers some", stack=stack)


def test_element_lists_and_marks_the_x_ray_lines_with_pixe_on(session, capsys):
    """ELEMENT adds each element's X-ray lines -- energy, and the channel on
    the PIXE plot's channel axis -- and labels them on the PIXE plot even
    with MARKERS OFF; with PIXE off it is RUMP's ELEMENT."""
    run(session, "element Fe")
    assert "X-rays" not in capsys.readouterr().out
    stack = run(session, f"pixe get {EXAMPLES / 'MnPt.PIX'}", "pixe", "calib 0.01 0",
                "markers off", "plot", "region 50 1190")
    capsys.readouterr()
    assert not any(label.startswith("Fe") for label in _labels(session))
    run(session, "return", stack=stack)
    run(session, "element Fe Pt")
    out = capsys.readouterr().out
    assert "X-rays  Kα 6.400 keV ch 640.0   Kβ 7.058 keV ch 705.8" in out
    assert "Lα 9.435 keV ch 943.5" in out
    assert "Kα 66." not in out  # beyond the spectrum: not listed
    labels = " ".join(_labels(session))
    assert "Fe Kα" in labels and "Pt Lα" in labels


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


# -- the simulation in the shell ----------------------------------------------


def _sim_curve(session):
    ax = session.pixe.figure.axes[0]
    return next((line for line in ax.lines if line.get_label() == "simulation"), None)


def test_simulation_alone_without_pixe_data(empty_session):
    run(empty_session, f"sim get {EXAMPLES / 'MnPt.lcm'}", "mev 1.9", "beam 4He", "pixe",
        "plot 0")
    assert _sim_curve(empty_session) is not None


def test_h_scales_the_simulation(session):
    stack = run(session, f"sim get {EXAMPLES / 'MnPt.lcm'}",
                f"pixe get {EXAMPLES / 'MnPt.PIX'}", "pixe", "cmp")
    before = np.array(_sim_curve(session).get_ydata())
    run(session, "h 2 2 2", stack=stack)
    after = np.array(_sim_curve(session).get_ydata())
    np.testing.assert_allclose(after, 2 * before, rtol=1e-9)
    run(session, "h K 1", stack=stack)
    assert session.pixe.h == (1.0, 2.0, 2.0)
    with pytest.raises(CommandError, match="positive"):
        run(session, "h 0 1 1", stack=stack)


def test_lines_lists_film_and_substrate_lines(session, capsys):
    run(session, f"sim get {EXAMPLES / 'MnPt.lcm'}", f"pixe get {EXAMPLES / 'MnPt.PIX'}")
    capsys.readouterr()
    run(session, "pixe lines")
    out = capsys.readouterr().out
    assert "Mn" in out and "KL3" in out and "Pt" in out
    assert "M5N7" in out  # Pt Ma, now that there are M-shell cross sections
    assert "Si " in out  # the Si substrate's lines too


def test_simulation_errors_are_reported_not_raised(session, capsys):
    run(session, f"sim get {EXAMPLES / 'MnPt.lcm'}", f"pixe get {EXAMPLES / 'MnPt.PIX'}",
        "beam 7Li", "pixe", "return", "plot 1", "overlay 0")
    out = capsys.readouterr().out
    assert "only protons and 4He" in out
    assert _sim_curve(session) is None  # the data are still drawn
    assert len(_pixe_labels(session)) == 1
    with pytest.raises(CommandError, match="only protons and 4He"):
        run(session, "pixe lines")


# -- detector geometry and materials ------------------------------------------


def test_numbered_filters(empty_session, capsys):
    """Numbered from the sample outwards; setting n replaces it or adds it
    after the last, CLEAR n removes one and the rest move up."""
    run(empty_session, "pixe filter 1 mylar 50 hole 60%", "pixe filter 2 Al 10 45%",
        "pixe filter 3 KAPTON 7.5")

    def stack():
        return [(f.material, f.thickness_um, f.hole_percent)
                for f in empty_session.pixe.detector.filters]

    assert stack() == [("Mylar", 50.0, 60.0), ("Al", 10.0, 45.0), ("Kapton", 7.5, 0.0)]
    run(empty_session, "pixe filter 2 Be 25")
    assert stack()[1] == ("Be", 25.0, 0.0)
    run(empty_session, "pixe filter clear 2")
    assert stack() == [("Mylar", 50.0, 60.0), ("Kapton", 7.5, 0.0)]
    capsys.readouterr()
    run(empty_session, "pixe filter")
    out = capsys.readouterr().out
    assert "filter 1 Mylar 50 hole 60%" in out and "filter 2 Kapton 7.5" in out
    assert "transmission" in out
    run(empty_session, "pixe filter clear")
    assert stack() == []
    with pytest.raises(CommandError, match="Mylar, Kapton"):
        run(empty_session, "pixe filter 1 Teflon 10")
    with pytest.raises(CommandError, match="unknown element"):
        run(empty_session, "pixe crystal Mylar 500")  # a crystal is an element


def test_solid_from_area_and_distance(empty_session, capsys):
    run(empty_session, "pixe solid 25 7.125 in")
    assert empty_session.pixe.detector.solid_angle_msr == pytest.approx(0.7633, rel=1e-3)
    run(empty_session, "pixe solid 25 180.975")
    assert empty_session.pixe.detector.solid_angle_msr == pytest.approx(0.7633, rel=1e-3)
    run(empty_session, "pixe solid 1.5")
    assert empty_session.pixe.detector.solid_angle_msr == 1.5
    with pytest.raises(CommandError, match="MM or IN"):
        run(empty_session, "pixe solid 25 7 ft")


def test_pixe_tilt_is_its_own_until_theta_rbs(empty_session, capsys):
    # The RBS buffer's THETA leaves PIXE at normal incidence ...
    run(empty_session, "theta -9")
    capsys.readouterr()
    run(empty_session, "pixe show")
    assert "THETA 0: beam 0 deg, X-rays 45 deg" in capsys.readouterr().out
    # ... until THETA RBS links it: IBM, |-9 + 45| = 36.
    run(empty_session, "pixe theta rbs")
    assert "THETA -9 from the RBS buffer: beam 9 deg, X-rays 36 deg" in capsys.readouterr().out
    run(empty_session, "theta 9", "pixe show")
    assert "beam 9 deg, X-rays 54 deg" in capsys.readouterr().out
    # Linked, it follows any RBS geometry -- PERT THETA's trials too.
    trial = Geometry(theta=-4.0, phi=10.0, kind=GeometryKind.CORNELL)
    assert empty_session.pixe.geometry_for(trial) == replace(DEFAULT_GEOMETRY, theta=-4.0)
    # A number unlinks it.
    run(empty_session, "pixe theta -20")
    assert not empty_session.pixe.theta_rbs
    assert "THETA -20: beam 20 deg, X-rays 25 deg" in capsys.readouterr().out


def test_pair_on_makes_the_tilt_follow_rbs(empty_session, capsys):
    run(empty_session, "theta -9", "pixe theta 5", "pixe pair on")
    assert "THETA -9 from the RBS buffer (PAIR ON): beam 9 deg, X-rays 36 deg" in (
        capsys.readouterr().out)
    trial = Geometry(theta=-4.0, phi=10.0)
    assert empty_session.pixe.geometry_for(trial).theta == -4.0
    # A number would break the pairing: refused, nothing changes.
    with pytest.raises(CommandError, match="PAIR OFF first"):
        run(empty_session, "pixe theta 12")
    assert empty_session.pixe.geometry.theta == 5
    # PAIR OFF: the PIXE prompt's own THETA again.
    run(empty_session, "pixe pair off")
    assert "THETA 5: beam 5 deg, X-rays 50 deg" in capsys.readouterr().out
    assert empty_session.pixe.geometry_for(trial).theta == 5


def test_pixe_geometry_follows_rbs_rules(empty_session, capsys):
    run(empty_session, "pixe theta 25", "pixe phi 35", "pixe geometry cornell")
    assert "X-rays 42.06 deg" in capsys.readouterr().out
    run(empty_session, "pixe geometry general", "pixe psi 53.6")
    assert "X-rays 53.6 deg" in capsys.readouterr().out
    assert empty_session.pixe.geometry == Geometry(
        theta=25.0, phi=35.0, psi=53.6, kind=GeometryKind.GENERAL)


def test_pixe_angles_shadow_rbs_ones_at_the_pixe_prompt(empty_session):
    stack = run(empty_session, "pixe", "theta 12", "phi 30")
    assert stack == ["rump", "pixe"]
    assert (empty_session.pixe.geometry.theta, empty_session.pixe.geometry.phi) == (12, 30)
    rbs = empty_session.settings.experiment_defaults.geometry
    assert (rbs.theta, rbs.phi) != (12, 30)


def test_defaults_describe_the_rc43_setup():
    detector = PixeDetector()
    assert (detector.window.material, detector.window.thickness_um) == ("Be", 12.5)
    assert (detector.fwhm_eV, detector.fano) == (122.0, 0.104)
    assert detector.solid_angle_msr == pytest.approx(0.763, rel=1e-3)
    assert [(f.material, f.thickness_um) for f in detector.filters] == [("Mylar", 62.0)]


def test_lin_is_the_linear_axis_and_lines_must_be_typed_in_full(session, capsys):
    stack = run(session, f"sim get {EXAMPLES / 'MnPt.lcm'}", f"pixe get {EXAMPLES / 'MnPt.PIX'}",
                "pixe")
    for abbreviation in ("lin", "line", "linear"):
        run(session, "log", abbreviation, stack=stack)
        assert session.pixe.plot.yscale == "linear", abbreviation
    capsys.readouterr()
    run(session, "lines", stack=stack)
    assert "KL3" in capsys.readouterr().out


# -- exports --------------------------------------------------------------------


def _columns(path: Path) -> tuple[list[str], np.ndarray]:
    lines = [l for l in path.read_text().splitlines() if not l.startswith("#")]
    delimiter = "," if path.suffix == ".csv" else "\t"
    names = lines[0].split(delimiter)
    return names, np.array([[float(v) for v in l.split(delimiter)] for l in lines[1:]])


def test_rump_level_exports_add_the_pixe_file_in_pair_mode(session, tmp_path):
    run(session, f"sim get {EXAMPLES / 'MnPt.lcm'}", f"pixe get {EXAMPLES / 'MnPt.PIX'}")
    run(session, f"export {tmp_path / 'off'}", f"ec {tmp_path / 'offc'}")
    assert sorted(p.name for p in tmp_path.iterdir()) == ["off.txt", "offc.txt"]
    run(session, "pixe pair on", f"export {tmp_path / 'on'}", f"ec {tmp_path / 'onc.csv'}")
    assert (tmp_path / "on_pixe.txt").exists() and (tmp_path / "onc_pixe.csv").exists()
    names, rows = _columns(tmp_path / "on_pixe.txt")
    assert names == ["channel", "energy_keV", "counts", "error"]
    assert rows[0, 0] == 1  # numbered as in the .PIX file
    np.testing.assert_allclose(rows[:, 2], session.buffers[1].pixe.spectrum.counts)


def test_pixe_prompt_exports_only_the_pixe_file(session, tmp_path):
    stack = run(session, f"sim get {EXAMPLES / 'MnPt.lcm'}", f"pixe get {EXAMPLES / 'MnPt.PIX'}",
                "pixe pair on", "pixe")
    run(session, f"export {tmp_path / 'a'}", f"exportcmp {tmp_path / 'b'}", stack=stack)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["a_pixe.txt", "b_pixe.txt"]
    names, rows = _columns(tmp_path / "b_pixe.txt")
    assert names[:6] == ["channel", "energy_keV", "counts", "simulation", "diff", "residual"]
    assert {"sim_Mn", "sim_Pt", "sim_Ru", "sim_Si"} <= set(names)
    elements = rows[:, [names.index(n) for n in names if n.startswith("sim_")]].sum(axis=1)
    np.testing.assert_allclose(elements, rows[:, names.index("simulation")], rtol=1e-5, atol=1e-5)


def test_pixe_export_needs_pixe_data(session, tmp_path):
    stack = run(session, "pixe")
    with pytest.raises(CommandError, match="no PIXE spectrum"):
        run(session, f"export {tmp_path / 'x'}", stack=stack)


def test_show_reports_the_dose_shared_with_rbs(empty_session, capsys):
    """CORRECTION divides the dose for PIXE as for RBS -- SHOW says so."""
    run(empty_session, "charge 10", "correction 0.8", "pixe show")
    output = capsys.readouterr().out
    assert "dose, shared with RBS: CHARGE 10 uC / CORRECTION 0.8 = 12.5 uC" in output
    assert "7.802e+13 ions" in output  # 12.5 uC / 1.602e-13 uC, charge state 1


def test_help_explains_fwhm_and_fano(empty_session, capsys):
    run(empty_session, "pixe help fwhm", "pixe help fano")
    output = capsys.readouterr().out
    assert "full width at half maximum of a line at" in output
    assert "FWHM(E)^2 = noise^2 + 2.355^2 x 3.64 eV x F x E" in output
