"""The interactive shell: buffers, plot state, mode stack and XEQ.

Commands are driven through :func:`pyrump.shell.repl.execute_line`, the same
entry point the prompt and ``XEQ`` use, so these tests exercise the real
dispatch path rather than calling handlers directly.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from pyrump.io.ascii import write_ascii  # noqa: E402
from pyrump.model.spectrum import Calibration, Spectrum  # noqa: E402
from pyrump.shell.commands.rump import Quit  # noqa: E402
from pyrump.shell.dispatch import CommandError  # noqa: E402
from pyrump.shell.repl import execute_file, execute_line  # noqa: E402
from pyrump.shell.session import Buffer, BufferSet, PlotState, Session  # noqa: E402
from pyrump.shell import plotting  # noqa: E402


from conftest import data_dir

DATA = data_dir()
needs_data = pytest.mark.skipif(DATA is None, reason="legacy data tables unavailable")


def make_buffer(total: float = 100.0, channels: int = 64, name: str = "test") -> Buffer:
    counts = np.full(channels, total / channels, dtype=float)
    return Buffer(
        spectrum=Spectrum(counts=counts, calibration=Calibration(npt=channels)),
        name=name,
    )


@pytest.fixture
def session(tmp_path) -> Session:
    """A session with the real tables, and one synthetic buffer in slot 1."""
    if DATA is None:
        pytest.skip("legacy data tables unavailable")
    built = Session.create(str(DATA))
    built.buffers.load(make_buffer(), 1)
    built.buffers.active = 1
    return built


def run(session: Session, *lines: str) -> None:
    stack = ["rump"]
    for line in lines:
        execute_line(session, line, stack)


# -- BufferSet, no data tables needed --------------------------------------


def test_buffers_start_with_only_the_simulation_slot():
    buffers = BufferSet()
    assert len(buffers) == 1
    assert buffers.get(0) is None


def test_first_free_skips_occupied_slots():
    buffers = BufferSet()
    assert buffers.load(make_buffer(), None) == 1
    assert buffers.load(make_buffer(), None) == 2
    buffers.release(1)
    assert buffers.load(make_buffer(), None) == 1


def test_release_moves_active_to_a_surviving_buffer():
    buffers = BufferSet()
    buffers.load(make_buffer(), 1)
    buffers.load(make_buffer(), 2)
    buffers.active = 1
    buffers.release(1)
    assert buffers.active == 2


def test_buffer_zero_cannot_be_released():
    with pytest.raises(KeyError, match="simulation"):
        BufferSet().release(0)


def test_find_path_matches_a_resolved_path(tmp_path):
    buffers = BufferSet()
    target = tmp_path / "a.rbs"
    target.write_text("")
    buffer = make_buffer()
    buffer.path = target
    buffers.load(buffer, 1)
    assert buffers.find_path(Path(target)) == 1
    assert buffers.find_path(tmp_path / "other.rbs") is None


def test_plot_state_region_defaults_to_the_whole_spectrum():
    assert PlotState().region(512) == (0, 511)
    assert PlotState(low=10, high=20).region(512) == (10, 20)
    # A region past the end is clipped rather than erroring.
    assert PlotState(low=10, high=9999).region(512) == (10, 511)


def test_plot_state_rejects_an_empty_region():
    with pytest.raises(ValueError, match="empty plot region"):
        PlotState(low=100, high=100).region(512)


# -- dispatch through the REPL ---------------------------------------------


@needs_data
def test_unknown_command_is_reported(session):
    with pytest.raises(CommandError, match="unrecognized command: wiggle"):
        run(session, "wiggle 3")


@needs_data
def test_comments_and_blank_lines_are_ignored(session):
    run(session, "", "   ", "/* just a comment", "# hash", "! bang")


@needs_data
def test_quit_propagates(session):
    with pytest.raises(Quit):
        run(session, "quit")


@needs_data
def test_region_updates_the_plot_state(session):
    run(session, "region 100 400")
    assert (session.plot.low, session.plot.high) == (100, 400)


@needs_data
def test_region_rejects_an_inverted_range(session):
    with pytest.raises(CommandError, match="empty region"):
        run(session, "region 400 100")


@needs_data
def test_yield_scale_commands(session):
    run(session, "sqrt")
    assert session.plot.yscale == "sqrt"
    run(session, "log")
    assert session.plot.yscale == "log"
    run(session, "linear")
    assert session.plot.yscale == "linear"


@needs_data
def test_normalize_and_raw_toggle(session):
    run(session, "normalize")
    assert session.plot.normalized is True
    run(session, "raw")
    assert session.plot.normalized is False


@needs_data
def test_faithful_toggle(session):
    assert session.settings.faithful is True
    run(session, "faithful off")
    assert session.settings.faithful is False
    run(session, "faithful on")
    assert session.settings.faithful is True


@needs_data
def test_faithful_persists_through_a_pyrumprc_style_macro(session, tmp_path):
    rc = tmp_path / ".pyrumprc"
    rc.write_text("faithful off\n")
    execute_file(session, rc)
    assert session.settings.faithful is False

    fresh = Session.create(str(DATA))
    execute_file(fresh, rc)
    assert fresh.settings.faithful is False


_SIM_SAMPLE = (
    "Sim Reset\nLayer 1\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"
)


@needs_data
def test_mev_sets_the_session_default_with_no_buffer_active():
    """MEV before any GET can't touch a real buffer, so it sets the
    session's fallback experiment settings instead of erroring."""
    fresh = Session.create(str(DATA))
    assert fresh.buffers.active_buffer is None
    run(fresh, "mev 3.5")
    assert fresh.settings.experiment_defaults.beam.e0_MeV == 3.5


@needs_data
def test_sim_only_simulation_uses_the_session_default_energy(tmp_path):
    """A sample can be explored with PLOT 0 before any GET, using the
    default MEV set beforehand -- e.g. from ~/.pyrumprc."""
    fresh = Session.create(str(DATA))
    sample = tmp_path / "defaults.lcm"
    sample.write_text(_SIM_SAMPLE)
    run(fresh, "mev 3.5", f"sim get {sample}")

    buffer = fresh.simulation()
    assert buffer.beam.e0_MeV == 3.5


@needs_data
def test_defaults_are_dormant_once_a_real_buffer_is_active():
    """Once real data is loaded, MEV/etc. go back to targeting it, and the
    session defaults are left untouched -- never a silent override."""
    fresh = Session.create(str(DATA))
    run(fresh, "mev 3.5")  # sets the default, no buffer active yet

    fresh.buffers.load(make_buffer(), 1)
    fresh.buffers.active = 1
    run(fresh, "mev 7.0")

    assert fresh.buffers[1].beam.e0_MeV == 7.0
    assert fresh.settings.experiment_defaults.beam.e0_MeV == 3.5


@needs_data
def test_experiment_defaults_persist_through_a_pyrumprc_style_macro(tmp_path):
    rc = tmp_path / ".pyrumprc"
    rc.write_text("mev 3.5\ntheta 5\n")
    fresh = Session.create(str(DATA))
    execute_file(fresh, rc)
    assert fresh.settings.experiment_defaults.beam.e0_MeV == 3.5
    assert fresh.settings.experiment_defaults.geometry.theta == 5.0


@needs_data
def test_ascii_load_picks_up_the_session_defaults(tmp_path):
    """ASCII carries no metadata of its own, so it should pick up the
    session's defaults rather than the code's hardcoded ones -- but its
    real channel count must never be overwritten by the placeholder."""
    fresh = Session.create(str(DATA))
    run(fresh, "mev 3.5", "theta 5")

    counts = np.zeros(200)
    path = tmp_path / "bare.dat"
    write_ascii(path, counts)
    run(fresh, f"get {path}")

    buffer = fresh.buffers.require_active()
    assert buffer.beam.e0_MeV == 3.5
    assert buffer.geometry.theta == 5.0
    assert buffer.n_channels == 200  # the real file's count, not the default's

    # A further MEV on this now-active buffer must not alias back into the
    # session default (Beam is mutable, so this is the real risk).
    run(fresh, "mev 9.0")
    assert buffer.beam.e0_MeV == 9.0
    assert fresh.settings.experiment_defaults.beam.e0_MeV == 3.5


@needs_data
def test_rbs_load_is_unaffected_by_session_defaults(tmp_path):
    from pyrump.io.rbs import RbsSpectrum, write_rbs
    from pyrump.model.detector import Measurement
    from pyrump.model.geometry import Geometry, GeometryKind

    fresh = Session.create(str(DATA))
    run(fresh, "mev 3.5")

    path = tmp_path / "real.rbs"
    write_rbs(
        path,
        RbsSpectrum(
            counts=np.zeros(100), calibration=Calibration(npt=100),
            geometry=Geometry(theta=0.0, phi=10.0, kind=GeometryKind.CORNELL),
            measurement=Measurement(), e0_MeV=2.5, zbeam=2, mbeam=4.0026,
        ),
    )
    run(fresh, f"get {path}")

    assert fresh.buffers.require_active().beam.e0_MeV == 2.5


# -- SIM MAXPTH --------------------------------------------------------------


_MAXPTH_SAMPLE = "Sim Reset\nLayer 1\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"


@needs_data
def test_maxpth_with_no_argument_shows_the_current_value(session, tmp_path, capsys):
    """Every other layer/sample setter needs a value; MAXPTH is sample-wide
    state, so "what is it now" is always answerable -- it must not fall
    through to the generic per-layer setter path and raise a raw IndexError
    from the missing float(rest[0])."""
    sample = tmp_path / "s.lcm"
    sample.write_text(_MAXPTH_SAMPLE)
    run(session, f"sim get {sample}")

    capsys.readouterr()
    run(session, "sim maxpth")
    assert "maxpth 200" in capsys.readouterr().out


@needs_data
def test_maxpth_with_a_value_sets_it_and_echoes_it_back(session, tmp_path, capsys):
    sample = tmp_path / "s.lcm"
    sample.write_text(_MAXPTH_SAMPLE)
    run(session, f"sim get {sample}")

    capsys.readouterr()
    run(session, "sim maxpth 500")
    assert "maxpth 500" in capsys.readouterr().out
    assert session.script.maxpth == 500.0

    capsys.readouterr()
    run(session, "sim maxpth")
    assert "maxpth 500" in capsys.readouterr().out


@needs_data
def test_abbreviations_work_through_the_repl(session):
    run(session, "reg 100 400", "sq", "norm")
    assert (session.plot.low, session.plot.high) == (100, 400)
    assert session.plot.yscale == "sqrt"
    assert session.plot.normalized is True


@needs_data
def test_metadata_setters_rebuild_the_frozen_dataclasses(session):
    run(session, "theta 7", "phi 15", "fwhm 20", "charge 5")
    buffer = session.buffers[1]
    assert buffer.geometry.theta == 7.0
    assert buffer.geometry.phi == 15.0
    assert buffer.measurement.fwhm_keV == 20.0
    assert buffer.measurement.charge_uC == 5.0


@needs_data
def test_setting_a_parameter_marks_the_simulation_stale(session):
    session.dirty = False
    run(session, "fwhm 22")
    assert session.dirty is True


@needs_data
def test_integral_sums_a_channel_range(session, capsys):
    """RbsThickn/TH_INT: gross/net in normalized #/uC/msr, not a raw sum."""
    run(session, "integral 0 63")
    out = capsys.readouterr().out
    assert "Gross:       10.00" in out
    assert "Net:        0.00" in out


@needs_data
def test_integral_clamps_channels_outside_the_spectrum(session, capsys):
    """RBSINDEX clamps out-of-range channels rather than rejecting them --
    matches the C, and INTEGRAL 0 9999 lands on the same result as 0 63."""
    run(session, "integral 0 9999")
    out = capsys.readouterr().out
    assert "Gross:       10.00" in out
    assert "Net:        0.00" in out


@needs_data
def test_copy_and_move_buffers(session, capsys):
    run(session, "copy 1 3")
    assert session.buffers.get(3) is not None
    # A copy is independent of its source.
    session.buffers[3].spectrum.counts[0] = 999.0
    assert session.buffers[1].spectrum.counts[0] != 999.0
    run(session, "move 1 3")
    assert session.buffers[1].spectrum.counts[0] == 999.0


# -- the mode stack --------------------------------------------------------


@needs_data
def test_sim_and_pert_push_and_pop_modes(session):
    stack = ["rump"]
    execute_line(session, "sim", stack)
    assert stack == ["rump", "sim"]
    execute_line(session, "return", stack)
    assert stack == ["rump"]
    execute_line(session, "pert", stack)
    assert stack == ["rump", "pert"]
    execute_line(session, "return", stack)
    assert stack == ["rump"]


@needs_data
def test_an_unknown_sim_command_falls_through_and_auto_returns(session):
    """sim.htm: RUMP is returned to automatically on a command SIM lacks."""
    stack = ["rump"]
    execute_line(session, "sim", stack)
    execute_line(session, "region 10 20", stack)
    assert stack == ["rump"]
    assert (session.plot.low, session.plot.high) == (10, 20)


@needs_data
def test_sim_accepts_a_one_shot_command_from_the_rump_level(session):
    run(session, "sim thick 500 A", "sim composition Si 1 /")
    assert session.script.layers[0].thickness == 500.0
    assert session.script.layers[0].composition == {"Si": 1.0}


@needs_data
def test_sim_editing_marks_the_simulation_stale(session):
    run(session, "sim thick 500 A")
    assert session.dirty is True


# -- XEQ -------------------------------------------------------------------


@needs_data
def test_xeq_runs_a_macro(session, tmp_path):
    macro = tmp_path / "m.cmd"
    macro.write_text("/* set up the view\nregion 50 250\nsqrt\n")
    run(session, f"xeq {macro}")
    assert (session.plot.low, session.plot.high) == (50, 250)
    assert session.plot.yscale == "sqrt"


@needs_data
def test_a_macro_and_typed_lines_reach_the_same_state(session, tmp_path):
    macro = tmp_path / "m.cmd"
    macro.write_text("region 50 250\nsqrt\nnormalize\n")
    execute_file(session, macro)
    from_macro = (session.plot.low, session.plot.high, session.plot.yscale,
                  session.plot.normalized)

    fresh = Session.create(str(DATA))
    fresh.buffers.load(make_buffer(), 1)
    fresh.buffers.active = 1
    run(fresh, "region 50 250", "sqrt", "normalize")
    assert from_macro == (fresh.plot.low, fresh.plot.high, fresh.plot.yscale,
                          fresh.plot.normalized)


@needs_data
def test_xeq_reports_the_failing_line(session, tmp_path):
    macro = tmp_path / "m.cmd"
    macro.write_text("region 50 250\nwiggle\n")
    with pytest.raises(CommandError, match=r"m\.cmd:2: unrecognized command"):
        run(session, f"xeq {macro}")


@needs_data
def test_xeq_on_a_missing_file(session):
    with pytest.raises(CommandError, match="no such command file"):
        run(session, "xeq nowhere.cmd")


@needs_data
@pytest.mark.parametrize("suffix", [".rbs", ".RBS"])
def test_xeq_finds_a_bare_name_with_an_rbs_extension(session, tmp_path, suffix):
    """Real RBS acquisition software often writes its output as a plain
    EMPTY/SWALLOW command macro under a ``.rbs``/``.RBS`` extension, meant
    to be replayed with XEQ rather than read with GET -- so a bare XEQ
    argument should find one of those, the same way it already finds a
    bare ``.cmd``."""
    macro = tmp_path / f"m{suffix}"
    macro.write_text("region 50 250\nsqrt\n")
    run(session, f"xeq {tmp_path / 'm'}")
    assert (session.plot.low, session.plot.high) == (50, 250)
    assert session.plot.yscale == "sqrt"


@needs_data
def test_xeq_of_genuine_binary_rbs_data_fails_cleanly(session, tmp_path):
    """A real spectrum file belongs with GET, not XEQ -- trying to XEQ one
    should give a clear error naming the mistake, not a raw decode
    traceback."""
    from pyrump.io.rbs import RbsSpectrum, write_rbs
    from pyrump.model.detector import Measurement
    from pyrump.model.geometry import Geometry
    from pyrump.model.spectrum import Calibration as RbsCalibration

    path = tmp_path / "data.rbs"
    write_rbs(
        path,
        RbsSpectrum(
            counts=np.arange(64, dtype=float),
            calibration=RbsCalibration(npt=64),
            geometry=Geometry(theta=0.0, phi=10.0),
            measurement=Measurement(),
        ),
    )
    with pytest.raises(CommandError, match="not a text command file"):
        run(session, f"xeq {path}")


@needs_data
def test_a_self_calling_macro_is_stopped(session, tmp_path):
    macro = tmp_path / "loop.cmd"
    macro.write_text(f"xeq {tmp_path / 'loop.cmd'}\n")
    with pytest.raises(CommandError, match="nested more than"):
        run(session, f"xeq {macro}")


@needs_data
def test_record_writes_a_replayable_macro(session, tmp_path):
    log = tmp_path / "session.cmd"
    run(session, f"record {log}", "region 50 250", "sqrt", "record off")
    written = log.read_text().splitlines()
    assert "region 50 250" in written
    assert "sqrt" in written

    fresh = Session.create(str(DATA))
    fresh.buffers.load(make_buffer(), 1)
    fresh.buffers.active = 1
    execute_file(fresh, log)
    assert (fresh.plot.low, fresh.plot.high) == (50, 250)
    assert fresh.plot.yscale == "sqrt"


# -- Analysis (anlytc.c's command family) ------------------------------------


@needs_data
def test_element_reports_the_surface_edge(session, capsys):
    run(session, "element Si Au")
    out = capsys.readouterr().out
    assert "Si" in out and "Au" in out
    assert "Channel=" in out


@needs_data
def test_element_rejects_an_unknown_symbol(session):
    with pytest.raises(CommandError, match="unknown element"):
        run(session, "element Xx")


@needs_data
def test_matrix_reports_expected_energy_and_height(session, capsys):
    run(session, "matrix Au")
    out = capsys.readouterr().out
    assert "Au expected at" in out
    assert "height" in out


@needs_data
def test_matrix_rejects_forbidden_kinematics(session):
    session.buffers.active_buffer.beam.z = 79
    session.buffers.active_buffer.beam.mass = 196.97
    with pytest.raises(CommandError, match="cannot occur"):
        run(session, "matrix Si")


@needs_data
def test_matrix_marks_an_existing_plot(session):
    # A wide enough calibration to actually cover Au's predicted edge energy
    # (~1845 keV), which the default 5 keV/ch test buffer doesn't reach.
    wide = make_buffer(channels=64)
    wide.spectrum.calibration = Calibration(kevch=50.0, npt=64)
    session.buffers.load(wide, 2)
    run(session, "plot 2", "matrix Au")
    ax = session.figure.axes[0]
    assert any(line.get_marker() == "+" for line in ax.lines)
    assert any(text.get_text() == "Au" for text in ax.texts)


@needs_data
def test_matrix_does_not_mark_with_no_plot_up(session, capsys):
    run(session, "matrix Au")
    assert session.figure is None


@needs_data
def test_whatisit_finds_candidates_near_a_channel(session, capsys):
    run(session, "whatisit 20")
    out = capsys.readouterr().out
    assert "near channel 20" in out


@needs_data
def test_whatisit_warns_with_no_plot_to_mark(session, capsys):
    run(session, "whatisit 20")
    out = capsys.readouterr().out
    assert "Plot device not enabled" in out


@needs_data
def test_whatisit_marks_an_existing_plot(session, capsys):
    run(session, "plot 1", "whatisit 20")
    out = capsys.readouterr().out
    assert "Plot device not enabled" not in out
    ax = session.figure.axes[0]
    assert len(ax.texts) >= 1
    assert any(len(line.get_xdata()) == 2 for line in ax.lines)


@needs_data
def test_info_reports_a_full_element_summary(session, capsys):
    run(session, "info Si")
    out = capsys.readouterr().out
    assert "Density:" in out
    assert "Abundance:" in out


@needs_data
def test_integral_honors_intset_interp_mode(session, capsys):
    run(session, "intset interp")
    run(session, "integral 0 63")
    out = capsys.readouterr().out
    assert "Interpolated integration" in out


@needs_data
def test_intset_question_mark_prints_the_mode_table(session, capsys):
    run(session, "intset ?")
    out = capsys.readouterr().out
    assert "current mode" in out


@needs_data
def test_intset_rejects_an_unrecognized_mode(session):
    with pytest.raises(CommandError, match="unrecognized mode"):
        run(session, "intset bogus")


@needs_data
def test_thickness_reports_atoms_and_angstroms(session, capsys):
    run(session, "thickness 0 63 Si")
    out = capsys.readouterr().out
    assert "Atoms/cm**2" in out
    assert "Angstroms" in out


@needs_data
def test_thickness_query_mode_needs_an_alpha_argument(session):
    run(session, "intset query")
    with pytest.raises(CommandError, match="alpha"):
        run(session, "thickness 0 63 Si")


@needs_data
def test_background_creates_a_new_cropped_buffer(session, capsys):
    run(session, "background 0 10 50 63 2 -noplot")
    out = capsys.readouterr().out
    assert "result in buffer" in out
    assert session.buffers.get(2) is not None


@needs_data
def test_background_inplace_modifies_the_active_buffer(session, capsys):
    active_before = session.buffers.active
    run(session, "background 0 10 50 63 2 -inplace -noplot")
    out = capsys.readouterr().out
    assert "subtracted in place" in out
    assert session.buffers.active == active_before


@needs_data
def test_smooth_default_mode_preserves_a_flat_spectrum(session):
    before = session.buffers.active_buffer.spectrum.counts.copy()
    run(session, "smooth")
    after = session.buffers.active_buffer.spectrum.counts
    np.testing.assert_allclose(after, before)


@needs_data
def test_fft_matches_smooth_dash_fft_dash_range(session):
    run(session, "smooth -fft -range 10 50 4")
    via_smooth = session.buffers.active_buffer.spectrum.counts.copy()
    session.buffers.active_buffer.spectrum.counts = np.full(64, 100.0 / 64)
    run(session, "fft 10 50 4")
    via_fft = session.buffers.active_buffer.spectrum.counts
    np.testing.assert_allclose(via_smooth, via_fft)


@needs_data
def test_width_thick_reports_a_thickness(session, capsys):
    run(session, "width_thick 10 40 Si")
    out = capsys.readouterr().out
    assert "Areal density:" in out
    assert "Thickness:" in out


@needs_data
def test_calibrate_updates_the_active_buffers_calibration(session, capsys):
    run(session, "calibrate 100 Si 500 Au")
    out = capsys.readouterr().out
    assert "Conversion:" in out
    assert "keV(0)" in out


@needs_data
def test_calibrate_rejects_near_identical_channels(session):
    with pytest.raises(CommandError, match="DIFFERENT channel"):
        run(session, "calibrate 100 Si 101 Au")


@needs_data
def test_calibrate_rejects_the_same_element_twice(session):
    with pytest.raises(CommandError, match="DIFFERENT elements"):
        run(session, "calibrate 100 Si 500 Si")


@needs_data
def test_offset_changes_only_kev0(session, capsys):
    run(session, "offset 12")
    calibration = session.buffers[1].calibration
    assert calibration.kev0 == 12.0
    assert calibration.kevch == 5.0  # untouched -- CONVERSION's job, not OFFSET's
    assert "offset = 12" in capsys.readouterr().out


@needs_data
def test_offset_with_no_argument_reports_the_current_value(session, capsys):
    run(session, "offset")
    assert "offset = 0" in capsys.readouterr().out


@needs_data
def test_offset_chains_into_a_further_command(session):
    run(session, "offset 5 fwhm 20")
    buffer = session.buffers[1]
    assert buffer.calibration.kev0 == 5.0
    assert buffer.measurement.fwhm_keV == 20.0


@needs_data
def test_compare_respects_region(session, tmp_path):
    """COMPARE = ``PLOT NOW ... OV THEORY`` in the original (cmds.htm), so it
    should inherit REGION the same way PLOT/OVERLAY do."""
    sample = tmp_path / "region_compare.lcm"
    sample.write_text(
        "Sim Reset\nLayer 1\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"
    )
    run(session, f"sim get {sample}", "region 20 40", "compare")

    top = session.figure.axes[0]
    line = top.lines[0]
    assert len(line.get_xdata()) == 21  # channels 20..40 inclusive


@needs_data
def test_compare_legend_shows_buffer_names_not_generic_labels(session, tmp_path):
    """PLOT's legend shows the buffer's own name; COMPARE should match
    instead of hard-coding "data"/"simulation"."""
    sample = tmp_path / "labels_compare.lcm"
    sample.write_text(
        "Sim Reset\nLayer 1\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"
    )
    run(session, f"sim get {sample}", "compare")

    top = session.figure.axes[0]
    labels = top.get_legend_handles_labels()[1]
    assert labels == ["test", "SIM"]  # the active buffer's name, then SIM's


@needs_data
def test_structlabel_off_by_default_keeps_sim_label(session, tmp_path):
    sample = tmp_path / "structlabel_off.lcm"
    sample.write_text(
        "Sim Reset\nLayer 1\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"
    )
    run(session, f"sim get {sample}", "compare")

    labels = session.figure.axes[0].get_legend_handles_labels()[1]
    assert labels == ["test", "SIM"]


@needs_data
def test_structlabel_on_shows_the_sample_structure(session, tmp_path):
    sample = tmp_path / "structlabel_on.lcm"
    sample.write_text(
        "Sim Reset\nLayer 1\n Thick 30 A\n Composition Ru 1 /\n"
        "Next\n Thick 150 A\n Composition Mn 3 Pt 1 /\n"
        "Next\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"
    )
    run(session, f"sim get {sample}", "structlabel on", "compare")

    labels = session.figure.axes[0].get_legend_handles_labels()[1]
    assert labels == ["test", "Si [500/cm2] - Mn3Pt [150A] - Ru [30A]"]

    run(session, "sim splot")
    trace_labels = [t.label for t in session.traces]
    assert trace_labels == ["Si [500/cm2] - Mn3Pt [150A] - Ru [30A]"]

    run(session, "structlabel off", "compare")
    labels = session.figure.axes[0].get_legend_handles_labels()[1]
    assert labels == ["test", "SIM"]


_SPLOT_SAMPLE = (
    "Sim Reset\nLayer 1\n Thick 30 A\n Composition Ru 1 /\n"
    "Next\n Thick 150 A\n Composition Mn 3 Pt 1 /\n"
    "Next\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"
)


@needs_data
def test_splot_element_overlays_only_that_elements_contribution(session, tmp_path):
    # A wide enough calibration to actually catch these edges (~hundreds of
    # keV), which the default 5 keV/ch test buffer doesn't reach.
    session.buffers.get(1).spectrum.calibration = Calibration(kevch=50.0, npt=64)
    sample = tmp_path / "splot_element.lcm"
    sample.write_text(_SPLOT_SAMPLE)
    run(session, f"sim get {sample}", "plot 1", "sim splot Ru")
    trace = session.traces[-1]
    assert trace.label == "SIM(Ru)"
    full = session.simulation()
    assert trace.buffer.spectrum.total() < full.spectrum.total()


@needs_data
def test_splot_layer_overlays_only_that_layers_contribution(session, tmp_path):
    session.buffers.get(1).spectrum.calibration = Calibration(kevch=50.0, npt=64)
    sample = tmp_path / "splot_layer.lcm"
    sample.write_text(_SPLOT_SAMPLE)
    run(session, f"sim get {sample}", "plot 1", "sim splot 2")
    trace = session.traces[-1]
    assert trace.label == "SIM(layer 2)"
    full = session.simulation()
    assert trace.buffer.spectrum.total() < full.spectrum.total()


@needs_data
def test_splot_rejects_unknown_element(session, tmp_path):
    sample = tmp_path / "splot_unknown.lcm"
    sample.write_text(
        "Sim Reset\nLayer 1\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"
    )
    run(session, f"sim get {sample}", "plot 1")
    with pytest.raises(CommandError, match="does not exist in target"):
        run(session, "sim splot Au")


@needs_data
def test_splot_rejects_out_of_range_layer(session, tmp_path):
    sample = tmp_path / "splot_layer_range.lcm"
    sample.write_text(
        "Sim Reset\nLayer 1\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"
    )
    run(session, f"sim get {sample}", "plot 1")
    with pytest.raises(CommandError, match="doesn't exist"):
        run(session, "sim splot 5")


@needs_data
def test_compfrac_on_shows_atomic_fraction(session, tmp_path):
    sample = tmp_path / "compfrac_on.lcm"
    sample.write_text(
        "Sim Reset\nLayer 1\n Thick 150 A\n Composition Mn 3 Pt 1 /\n"
        "Next\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"
    )
    run(session, f"sim get {sample}", "structlabel on", "compfrac on", "compare")

    labels = session.figure.axes[0].get_legend_handles_labels()[1]
    assert labels == ["test", "Si [500/cm2] - Mn0.75Pt0.25 [150A]"]

    run(session, "compfrac off", "compare")
    labels = session.figure.axes[0].get_legend_handles_labels()[1]
    assert labels == ["test", "Si [500/cm2] - Mn3Pt [150A]"]


@needs_data
def test_plot_and_compare_reuse_the_same_figure(session, tmp_path):
    """PLOT (1 panel) -> COMPARE (2 panels) -> PLOT (1 panel again) must draw
    into the same window throughout, so a user's dragged window position
    survives switching between them."""
    sample = tmp_path / "reuse_figure.lcm"
    sample.write_text(
        "Sim Reset\nLayer 1\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"
    )

    run(session, "plot")
    figure = session.figure
    number = figure.number

    run(session, f"sim get {sample}", "compare")
    assert session.figure is figure
    assert session.figure.number == number
    assert len(session.figure.axes) == 2

    run(session, "plot")
    assert session.figure is figure
    assert session.figure.number == number
    assert len(session.figure.axes) == 1


@needs_data
def test_a_hand_resized_window_survives_switching_to_compare(session, tmp_path):
    """A window the user dragged bigger (or smaller) must stay that size
    across a panel-count change, not snap back to COMPARE's own default --
    only a genuinely new figure (e.g. after the user closes the window)
    should ever apply a default size."""
    sample = tmp_path / "resize_survives.lcm"
    sample.write_text(
        "Sim Reset\nLayer 1\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"
    )

    run(session, "plot")
    session.figure.set_size_inches(14, 10)

    run(session, f"sim get {sample}", "compare")
    assert session.figure.get_size_inches() == pytest.approx((14, 10))

    run(session, "plot")
    assert session.figure.get_size_inches() == pytest.approx((14, 10))


@needs_data
def test_profile_prints_the_verbatim_stub_message(session, capsys):
    before = session.buffers.active_buffer.spectrum.counts.copy()
    run(session, "profile")
    out = capsys.readouterr().out
    assert "not implemented" in out
    np.testing.assert_array_equal(session.buffers.active_buffer.spectrum.counts, before)


@needs_data
def test_cursor_prints_the_verbatim_stub_message(session, capsys):
    run(session, "cursor")
    out = capsys.readouterr().out
    assert "Cursor not enabled" in out


@needs_data
def test_cursor_reads_clicked_points_until_stopped(session, monkeypatch, capsys):
    run(session, "plot 1")
    points = iter([(20.0, 5.0), (40.0, 1.5), None])
    monkeypatch.setattr(plotting, "read_cursor_point", lambda figure: next(points))
    run(session, "cursor")
    out = capsys.readouterr().out
    assert "Cursor not enabled" not in out
    assert out.count("Channel:") == 2
    assert "Counts:" in out


@needs_data
def test_cursor_reports_yield_when_normalized(session, monkeypatch, capsys):
    run(session, "normalize", "plot 1")
    points = iter([(20.0, 5.0), None])
    monkeypatch.setattr(plotting, "read_cursor_point", lambda figure: next(points))
    run(session, "cursor")
    out = capsys.readouterr().out
    assert "Yield:" in out
    assert "/uC/keV/msr" in out


@pytest.mark.parametrize(
    "abbreviation, expected",
    [
        ("cur", "CURSOR"), ("el", "ELEMENT"), ("mat", "MATRIX"), ("what", "WHATISIT"),
        ("inf", "INFO"), ("int", "INTEGRAL"), ("thic", "THICKNESS"),
        ("back", "BACKGROUND"), ("smo", "SMOOTH"), ("wid", "WIDTH_THICK"),
        ("pro", "PROFILE"), ("intset", "INTSET"), ("cal", "CALIBRATE"),
        ("dis", "DISPLAY"), ("fft", "FFT"),
    ],
)
def test_analysis_command_abbreviations_resolve_correctly(abbreviation, expected):
    """Table order regression guard -- cheap insurance against a future
    reordering (like DISPLAY's move to its real cmlist position) silently
    breaking which command a short abbreviation lands on."""
    from pyrump.shell.commands.rump import TABLE

    matched = TABLE.match(abbreviation)
    assert matched is not None
    assert matched.name == expected


def test_compare_shows_its_cmp_synonym_at_every_level():
    """CMP is registered as its own hidden entry (kept out of the listing to
    avoid a redundant "synonym for COMPARE" row), so without note_synonym a
    bare HELP/`?` at any of the three levels would leave it looking
    undiscoverable."""
    from pyrump.shell.commands.pert import TABLE as PERT_TABLE
    from pyrump.shell.commands.rump import TABLE as RUMP_TABLE
    from pyrump.shell.commands.sim import TABLE as SIM_TABLE

    for table in (RUMP_TABLE, SIM_TABLE, PERT_TABLE):
        assert table.match("COMPARE").display == "COMPARE / CMP"
        assert table.match("COMP") is None or table.match("COMP").name != "COMPARE"
