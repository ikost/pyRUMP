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
from pyrump.physics.xsec.rutherford import ScreeningModel  # noqa: E402
from pyrump.shell.commands.rump import Quit  # noqa: E402
from pyrump.shell.dispatch import CommandError  # noqa: E402
from pyrump.shell.repl import execute_file, execute_line  # noqa: E402
from pyrump.shell.session import Buffer, BufferSet, PlotState, Session  # noqa: E402


from conftest import data_dir

DATA = data_dir()
needs_data = pytest.mark.skipif(DATA is None, reason="legacy data tables unavailable")


def make_buffer(total: float = 100.0, channels: int = 64, name: str = "test") -> Buffer:
    counts = np.full(channels, total / channels, dtype=float)
    return Buffer(
        spectrum=Spectrum(counts=counts, calibration=Calibration(npt=channels, kevch=5.0)),
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


def test_clear_also_blanks_buffer_zero():
    buffers = BufferSet()
    buffers.set(0, make_buffer())
    buffers.load(make_buffer(), 1)
    buffers.active = 1
    buffers.clear()
    assert buffers.get(0) is None
    assert buffers.active == 0


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
def test_expand_narrows_the_plotted_region(session):
    run(session, "plot 1", "expand 10 40")
    assert (session.plot.low, session.plot.high) == (10, 40)


@needs_data
def test_expand_swaps_limits_given_high_to_low(session):
    run(session, "plot 1", "expand 40 10")
    assert (session.plot.low, session.plot.high) == (10, 40)


@needs_data
def test_rejected_expand_leaves_the_region_alone(session):
    run(session, "plot 1", "region 10 40")
    with pytest.raises(CommandError, match="subset"):
        run(session, "expand 5 50")
    with pytest.raises(CommandError, match="empty region"):
        run(session, "expand 20 20")
    assert (session.plot.low, session.plot.high) == (10, 40)
    run(session, "replot")


@needs_data
def test_expand_with_nothing_plotted_stores_nothing(session):
    with pytest.raises(CommandError, match="nothing plotted"):
        run(session, "expand 10 40")
    assert (session.plot.low, session.plot.high) == (None, None)


@needs_data
def test_yield_scale_commands(session):
    run(session, "sqrt")
    assert session.plot.yscale == "sqrt"
    run(session, "log")
    assert session.plot.yscale == "log"
    run(session, "linear")
    assert session.plot.yscale == "linear"


@needs_data
def test_wrascii_writes_rumps_own_text_format(session, tmp_path):
    """bmanip.c:595-650's keyword header, not a bare counts dump -- verified
    byte-for-byte against the real binary in test_io_rbs's oracle tests."""
    out = tmp_path / "buffer.dat"
    run(session, f"wrascii {out}")

    text = out.read_text()
    lines = text.splitlines()
    assert lines[0] == "Empty File 'test'"
    assert lines[1] == "Spectrum    RBS"
    assert "Swallow" in lines
    assert lines[lines.index("Swallow") + 1] == "1.562500"  # 100/64 counts

    from pyrump.io.ascii import read_ascii

    result = read_ascii(out)
    assert result.metadata["Spectrum"] == "RBS"
    assert np.allclose(result.counts, session.buffers[1].spectrum.counts)


@needs_data
def test_wrascii_defaults_to_dat_extension(session, tmp_path):
    run(session, f"wrascii {tmp_path / 'buffer'}")
    assert (tmp_path / "buffer.dat").exists()


@needs_data
def test_wrascii_warns_about_an_rbs_extension(session, tmp_path, capsys):
    """.rbs/.RBS makes GET mistake this output for an RC43 macro and refuse
    it (see file-formats.md), so it's worth a warning rather than a silent
    footgun."""
    target = tmp_path / "buffer.rbs"
    run(session, f"wrascii {target}")
    assert target.exists()  # still written -- just warned about
    assert "GET will refuse this" in capsys.readouterr().out


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
def test_faithful_toggle_changes_the_simulated_spectrum(session, tmp_path):
    """Regression for the FAITHFUL command being a no-op end to end: the
    setting must reach sim/precal.py and sim/outbound.py, and toggling it
    must invalidate the cached buffer 0 (session.touch()), not just leave
    session.settings.faithful set with nothing rereading it."""
    # A wide enough calibration to actually catch these edges (~hundreds of
    # keV), which the default 5 keV/ch test buffer doesn't reach.
    session.buffers.get(1).spectrum.calibration = Calibration(kevch=50.0, npt=64)
    sample = tmp_path / "faithful_toggle.lcm"
    sample.write_text(_SPLOT_SAMPLE)
    run(session, f"sim get {sample}")

    run(session, "faithful on")
    on_first = session.simulation().spectrum.counts.copy()

    run(session, "faithful off")
    off = session.simulation().spectrum.counts.copy()
    assert not np.array_equal(on_first, off)

    run(session, "faithful on")
    on_again = session.simulation().spectrum.counts.copy()
    assert np.array_equal(on_first, on_again)


@needs_data
def test_faithful_persists_through_a_pyrumprc_style_macro(session, tmp_path):
    rc = tmp_path / ".pyrumprc"
    rc.write_text("faithful off\n")
    execute_file(session, rc)
    assert session.settings.faithful is False

    fresh = Session.create(str(DATA))
    execute_file(fresh, rc)
    assert fresh.settings.faithful is False


@needs_data
def test_screening_select(session, capsys):
    assert session.settings.screening is ScreeningModel.LECUYER
    run(session, "screening andersen")
    assert session.settings.screening is ScreeningModel.ANDERSEN
    assert "screening is ANDERSEN" in capsys.readouterr().out
    run(session, "screening none")
    assert session.settings.screening is ScreeningModel.NONE
    run(session, "screening lecuyer")
    assert session.settings.screening is ScreeningModel.LECUYER


@needs_data
def test_screening_change_invalidates_the_cached_simulation(session, tmp_path):
    """SCREENING had the same staleness bug FAITHFUL did: changing it never
    called session.touch(), so buffer 0 kept the cross-section computed under
    the old model until something unrelated marked it dirty."""
    session.buffers.get(1).spectrum.calibration = Calibration(kevch=50.0, npt=64)
    sample = tmp_path / "screening_touch.lcm"
    sample.write_text(_SPLOT_SAMPLE)
    run(session, f"sim get {sample}")

    lecuyer = session.simulation().spectrum.counts.copy()
    run(session, "screening andersen")
    andersen = session.simulation().spectrum.counts.copy()
    assert not np.array_equal(lecuyer, andersen)


@needs_data
def test_screening_with_no_argument_shows_the_current_value(session, capsys):
    capsys.readouterr()
    run(session, "screening")
    assert "screening is LECUYER" in capsys.readouterr().out


@needs_data
def test_screening_rejects_an_unknown_model(session):
    with pytest.raises(CommandError, match="unknown screening model"):
        run(session, "screening bogus")


@needs_data
def test_screening_persists_through_a_pyrumprc_style_macro(session, tmp_path):
    rc = tmp_path / ".pyrumprc"
    rc.write_text("screening andersen\n")
    execute_file(session, rc)
    assert session.settings.screening is ScreeningModel.ANDERSEN

    fresh = Session.create(str(DATA))
    execute_file(fresh, rc)
    assert fresh.settings.screening is ScreeningModel.ANDERSEN


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
def test_element_labels_its_energy_in_keV(session, capsys):
    run(session, "element Au")
    out = capsys.readouterr().out
    assert " keV  Channel=" in out


@needs_data
def test_element_does_not_mark_with_no_plot_up(session):
    run(session, "element Si Au")
    assert session.figure is None


@needs_data
def test_element_ticks_each_edge_on_an_existing_plot(session):
    # Wide enough to cover both Si (~1130 keV) and Au (~1845 keV) at 2 MeV.
    wide = make_buffer(channels=64)
    wide.spectrum.calibration = Calibration(kevch=50.0, npt=64)
    session.buffers.load(wide, 2)
    run(session, "plot 2", "element Si Au")
    ax = session.figure.axes[0]
    assert {text.get_text() for text in ax.texts} >= {"Si", "Au"}
    ticks = [line for line in ax.lines if len(line.get_xdata()) == 2]
    assert len(ticks) == 2


@needs_data
def test_element_labels_an_isotope_with_its_mass_number(session):
    wide = make_buffer(channels=64)
    wide.spectrum.calibration = Calibration(kevch=50.0, npt=64)
    session.buffers.load(wide, 2)
    run(session, "plot 2", "element 28Si")
    assert any(text.get_text() == "$^{28}$Si" for text in session.figure.axes[0].texts)


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
def test_matrix_expands_y_range_to_fit_a_tall_marker(session):
    # Au's predicted height dwarfs this buffer's flat, low-count trace, so
    # draw()'s frozen y range (set by PLOT, which disables autoscale) would
    # clip the marker unless mark_matrix() expands it to fit.
    wide = make_buffer(total=1.0, channels=64)
    wide.spectrum.calibration = Calibration(kevch=50.0, npt=64)
    session.buffers.load(wide, 2)
    run(session, "plot 2")
    _, before_top = session.figure.axes[0].get_ylim()
    run(session, "matrix Au")
    ax = session.figure.axes[0]
    marker_y = next(line.get_ydata()[0] for line in ax.lines if line.get_marker() == "+")
    assert marker_y > before_top
    assert ax.get_ylim()[1] >= marker_y


@needs_data
def test_matrix_respects_an_explicit_yhigh(session):
    wide = make_buffer(total=1.0, channels=64)
    wide.spectrum.calibration = Calibration(kevch=50.0, npt=64)
    session.buffers.load(wide, 2)
    run(session, "plot 2", "blowup 10")
    run(session, "matrix Au")
    assert session.figure.axes[0].get_ylim()[1] == pytest.approx(10.0)


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
def test_whatisit_marks_a_comparison_without_wiping_it(session, tmp_path):
    sample = tmp_path / "compare_whatisit.lcm"
    sample.write_text("Sim Reset\nLayer 1\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n")
    run(session, f"sim get {sample}", "compare", "whatisit 20")
    assert len(session.figure.axes) == 2
    top = session.figure.axes[0]
    assert top.get_legend_handles_labels()[1] == ["test", "SIM"]
    assert any(len(line.get_xdata()) == 2 for line in top.lines)


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
def test_slope_changes_only_kevch(session, capsys):
    run(session, "slope 4.5")
    calibration = session.buffers[1].calibration
    assert calibration.kevch == 4.5
    assert calibration.kev0 == 0.0  # untouched -- CONVERSION's job, not SLOPE's
    assert "slope = 4.5" in capsys.readouterr().out


@needs_data
def test_slope_with_no_argument_reports_the_current_value(session, capsys):
    run(session, "slope")
    assert "slope = 5" in capsys.readouterr().out


@needs_data
def test_slope_chains_into_a_further_command(session):
    run(session, "slope 4.5 fwhm 20")
    buffer = session.buffers[1]
    assert buffer.calibration.kevch == 4.5
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
def test_compare_shows_goodness_of_fit_over_region_with_no_pert_window(session, tmp_path):
    """With nothing selected in PERT, the chi-square text should fall back
    to the plot's own REGION -- 21 channels (20..40 inclusive) here."""
    sample = tmp_path / "gof_region.lcm"
    sample.write_text(
        "Sim Reset\nLayer 1\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"
    )
    run(session, f"sim get {sample}", "region 20 40", "compare")

    bottom = session.figure.axes[1]
    texts = [t.get_text() for t in bottom.texts]
    assert any("reduced chi-square" in t and "(21 dof)" in t for t in texts)


@needs_data
def test_compare_goodness_of_fit_uses_pert_error_window(session, tmp_path):
    """With a PERT error window set, the chi-square text should match what
    GO itself would report: dof = the window's channel count minus however
    many parameters are currently selected to vary."""
    sample = tmp_path / "gof_window.lcm"
    sample.write_text(
        "Sim Reset\nLayer 1\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"
    )
    run(session, f"sim get {sample}", "pert", "window 10 19", "thick 1", "return", "compare")

    bottom = session.figure.axes[1]
    texts = [t.get_text() for t in bottom.texts]
    # 10 channels (10..19 inclusive) minus 1 varying parameter.
    assert any("reduced chi-square" in t and "(9 dof)" in t for t in texts)


_COMPARE_SAMPLE = (
    "Sim Reset\nLayer 1\n Thick 30 A\n Composition Ru 1 /\n"
    "Next\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"
)


@needs_data
def test_splot_after_compare_adds_to_the_comparison(session, tmp_path):
    """COMPARE is ``PLOT NOW ... OV THEORY`` and SPLOT overlays (RbsSplot's
    ``RbsPlot(PLT_OV, ...)``), so SPLOT Ru after COMPARE must keep the data,
    the simulation and the residuals panel, adding its curve on top."""
    sample = tmp_path / "compare_splot.lcm"
    sample.write_text(_COMPARE_SAMPLE)
    run(session, f"sim get {sample}", "compare", "sim splot Ru")

    assert len(session.figure.axes) == 2
    top, bottom = session.figure.axes
    assert top.get_legend_handles_labels()[1] == ["test", "SIM", "SIM(Ru)"]
    assert any("reduced chi-square" in t.get_text() for t in bottom.texts)


@needs_data
def test_bare_splot_after_compare_replaces_the_simulation_in_place(session, tmp_path):
    sample = tmp_path / "compare_bare_splot.lcm"
    sample.write_text(_COMPARE_SAMPLE)
    run(session, f"sim get {sample}", "compare", "sim splot")

    assert len(session.figure.axes) == 2
    assert session.figure.axes[0].get_legend_handles_labels()[1] == ["test", "SIM"]


@needs_data
def test_overlay_after_compare_adds_to_the_comparison(session, tmp_path):
    sample = tmp_path / "compare_overlay.lcm"
    sample.write_text(_COMPARE_SAMPLE)
    session.buffers.load(make_buffer(name="other"), 2)
    run(session, f"sim get {sample}", "compare", "overlay 2")

    assert len(session.figure.axes) == 2
    assert session.figure.axes[0].get_legend_handles_labels()[1] == ["test", "SIM", "other"]


@needs_data
def test_region_after_compare_redraws_the_comparison(session, tmp_path):
    sample = tmp_path / "compare_region.lcm"
    sample.write_text(_COMPARE_SAMPLE)
    run(session, f"sim get {sample}", "compare", "region 20 40")

    assert len(session.figure.axes) == 2
    top, bottom = session.figure.axes
    assert len(top.lines[0].get_xdata()) == 21
    assert any("(21 dof)" in t.get_text() for t in bottom.texts)


@needs_data
def test_log_after_compare_rescales_the_top_panel(session, tmp_path):
    sample = tmp_path / "compare_log.lcm"
    sample.write_text(_COMPARE_SAMPLE)
    run(session, f"sim get {sample}", "compare", "log")

    assert len(session.figure.axes) == 2
    assert session.figure.axes[0].get_yscale() == "log"


@needs_data
def test_compare_honours_a_log_scale_set_beforehand(session, tmp_path):
    """The residuals are in sigma and go negative, so only the top panel
    goes log."""
    sample = tmp_path / "log_compare.lcm"
    sample.write_text(_COMPARE_SAMPLE)
    run(session, f"sim get {sample}", "log", "compare")

    top, bottom = session.figure.axes
    assert (top.get_yscale(), bottom.get_yscale()) == ("log", "linear")


@needs_data
@pytest.mark.parametrize("command", ["plot 1", "compare"])
def test_log_drops_a_zero_counts_floor_instead_of_warning(session, tmp_path, recwarn, command):
    sample = tmp_path / "log_floor.lcm"
    sample.write_text(_COMPARE_SAMPLE)
    run(session, f"sim get {sample}", "counts 0 5", "log", command)

    bottom, top = session.figure.axes[0].get_ylim()
    assert 0 < bottom < top == 5
    assert not [w for w in recwarn if "non-positive ylim" in str(w.message)]


@needs_data
def test_normalize_applies_to_compare_and_flags_the_goodness_of_fit(session, tmp_path):
    from pyrump.model.detector import yield_normalisation

    sample = tmp_path / "normalize_compare.lcm"
    sample.write_text(_COMPARE_SAMPLE)
    run(session, f"sim get {sample}", "compare")
    top, bottom = session.figure.axes
    raw = top.lines[0].get_ydata().copy()
    assert not any("normalized" in t.get_text() for t in bottom.texts)

    run(session, "normalize")
    top, bottom = session.figure.axes
    factor = yield_normalisation(session.buffers[1].measurement)
    np.testing.assert_allclose(top.lines[0].get_ydata(), raw / factor)
    assert top.get_ylabel() == "Yield (counts/msr/uC)"
    assert any(
        "reduced chi-square" in t.get_text() and "normalized yield" in t.get_text()
        for t in bottom.texts
    )


@needs_data
def test_plot_after_compare_returns_to_a_single_panel(session, tmp_path):
    sample = tmp_path / "compare_then_plot.lcm"
    sample.write_text(_COMPARE_SAMPLE)
    run(session, f"sim get {sample}", "compare", "plot 1", "sim splot Ru")

    assert len(session.figure.axes) == 1
    assert [t.label for t in session.traces] == ["test", "SIM(Ru)"]


_EC_SAMPLE = "Sim Reset\nLayer 1\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"


def _ec_sample(session, tmp_path, text: str = _EC_SAMPLE) -> Path:
    """Write the sample, and widen buffer 1 to 50 keV/ch so its 64 channels
    reach the Si edge (~1.1 MeV) -- the default 5 keV/ch stops at 320 keV,
    where the simulation is all zeros."""
    session.buffers.get(1).spectrum.calibration = Calibration(kevch=50.0, npt=64)
    sample = tmp_path / "ec.lcm"
    sample.write_text(text)
    return sample


def _read_export(path):
    """(comment lines without their "# ", column names, data rows) of an
    EXPORTCMP file, split on its own delimiter."""
    lines = path.read_text().splitlines()
    comments = [line[2:] for line in lines if line.startswith("#")]
    body = [line for line in lines if not line.startswith("#")]
    delimiter = "," if path.suffix == ".csv" else "\t"
    names = body[0].split(delimiter)
    rows = np.array([[float(v) for v in line.split(delimiter)] for line in body[1:]])
    return comments, names, rows


@needs_data
def test_exportcmp_writes_data_simulation_and_difference_as_columns(session, tmp_path):
    sample = _ec_sample(session, tmp_path)
    out = tmp_path / "fit.txt"
    run(session, f"sim get {sample}", f"ec {out}")

    _, names, rows = _read_export(out)
    assert names == ["channel", "energy_keV", "counts", "simulation", "diff", "residual"]
    data = session.buffers[1]
    theory = session.simulation().spectrum.counts
    np.testing.assert_array_equal(rows[:, 0], np.arange(data.n_channels))
    np.testing.assert_allclose(rows[:, 1], data.spectrum.energies, atol=1e-4)
    np.testing.assert_allclose(rows[:, 2], data.spectrum.counts, atol=1e-6)
    np.testing.assert_allclose(rows[:, 3], theory, atol=1e-6)
    np.testing.assert_allclose(rows[:, 4], data.spectrum.counts - theory, atol=1e-5)


@needs_data
def test_exportcmp_residual_is_the_poisson_residual_and_nan_without_a_model(session, tmp_path):
    from pyrump.fit.objective import poisson_residuals

    sample = _ec_sample(session, tmp_path)
    out = tmp_path / "fit.txt"
    run(session, f"sim get {sample}", f"ec {out}")

    _, _, rows = _read_export(out)
    residual = rows[:, 5]
    # Against the full-precision spectra, not the file's 6-decimal columns:
    # in the edge tail, where the simulation is tiny, that rounding alone
    # shifts the residual.
    counts = session.buffers[1].spectrum.counts
    theory = session.simulation().spectrum.counts
    modelled = theory > 0
    assert modelled.any() and (~modelled).any()
    expected, _ = poisson_residuals(counts[modelled], theory[modelled])
    np.testing.assert_allclose(residual[modelled], expected, atol=1e-5)
    assert np.isnan(residual[~modelled]).all()


@needs_data
def test_exportcmp_gof_matches_compare_over_region(session, tmp_path):
    sample = _ec_sample(session, tmp_path)
    out = tmp_path / "fit.txt"
    run(session, f"sim get {sample}", "region 20 40", "compare", f"ec {out}")

    shown = next(
        t.get_text() for t in session.figure.axes[1].texts if "reduced chi-square" in t.get_text()
    )
    comments, _, _ = _read_export(out)
    gof = next(line for line in comments if line.startswith("GOF"))
    assert gof == f"GOF            {shown} over REGION channels 20-40"


@needs_data
def test_exportcmp_gof_names_the_pert_error_window(session, tmp_path):
    sample = _ec_sample(session, tmp_path)
    out = tmp_path / "fit.txt"
    run(session, f"sim get {sample}", "pert", "window 10 19", "thick 1", "return", f"ec {out}")

    comments, _, _ = _read_export(out)
    gof = next(line for line in comments if line.startswith("GOF"))
    assert "(9 dof)" in gof
    assert gof.endswith("over PERT error windows 10-19 (1 varying parameter)")


@needs_data
def test_exportcmp_writes_raw_counts_under_normalize(session, tmp_path):
    sample = _ec_sample(session, tmp_path)
    raw, normalized = tmp_path / "raw.txt", tmp_path / "norm.txt"
    run(session, f"sim get {sample}", f"ec {raw}", "normalize", f"ec {normalized}")

    raw_comments, _, raw_rows = _read_export(raw)
    norm_comments, _, norm_rows = _read_export(normalized)
    np.testing.assert_array_equal(raw_rows, norm_rows)
    gof = next(line for line in norm_comments if line.startswith("GOF"))
    assert gof == next(line for line in raw_comments if line.startswith("GOF"))
    assert "normalized" not in gof


@needs_data
def test_exportcmp_csv_is_comma_separated_and_a_bare_name_gets_txt(session, tmp_path):
    sample = _ec_sample(session, tmp_path)
    run(session, f"sim get {sample}", f"ec {tmp_path / 'fit.csv'}", f"ec {tmp_path / 'fit'}")

    _, names, _ = _read_export(tmp_path / "fit.csv")
    assert names[0] == "channel" and len(names) == 6
    assert (tmp_path / "fit.txt").exists()
    assert "\t" in (tmp_path / "fit.txt").read_text().splitlines()[-1]


@needs_data
def test_exportcmp_header_describes_physics_and_sample(session, tmp_path):
    sample = _ec_sample(session, tmp_path, _SPLOT_SAMPLE)
    out = tmp_path / "fit.txt"
    run(session, f"sim get {sample}", "screening andersen", f"ec {out}")

    comments, _, _ = _read_export(out)
    rows = {line.split()[0]: line for line in comments}
    assert rows["Physics"].split() == ["Physics", "FAITHFUL", "on", "SCREENING", "ANDERSEN"]
    assert rows["Sample"].endswith("Si [500/cm2] - Mn3Pt [150A] - Ru [30A]")
    assert rows["Areal"].endswith("(1e15 at/cm2)") and "Si [500.0/cm2]" in rows["Areal"]
    assert rows["Normalisation"].startswith("Normalisation  divide counts by")


@needs_data
def test_exportcmp_refreshes_a_stale_simulation(session, tmp_path):
    sample = _ec_sample(session, tmp_path)
    out = tmp_path / "fit.txt"
    run(session, f"sim get {sample}", "compare")
    stale = session.buffers.get(0).spectrum.counts.copy()
    run(session, "sim thick 900 /cm2", f"ec {out}")

    _, _, rows = _read_export(out)
    fresh = session.simulation().spectrum.counts
    assert not np.allclose(stale, fresh)
    np.testing.assert_allclose(rows[:, 3], fresh, atol=1e-6)


@needs_data
def test_exportcmp_from_pert_stays_in_pert(session, tmp_path):
    sample = _ec_sample(session, tmp_path)
    stack = ["rump"]
    for line in (f"sim get {sample}", "pert", f"ec {tmp_path / 'fit.txt'}"):
        execute_line(session, line, stack)
    assert stack == ["rump", "pert"]
    assert (tmp_path / "fit.txt").exists()


@needs_data
def test_exportcmp_refuses_the_simulation_buffer_as_data(session, tmp_path):
    sample = _ec_sample(session, tmp_path)
    run(session, f"sim get {sample}", "compare")
    session.buffers.active = 0
    with pytest.raises(CommandError, match="data buffer"):
        run(session, f"ec {tmp_path / 'fit.txt'}")


@needs_data
def test_export_writes_the_active_buffer_with_poisson_errors(session, tmp_path):
    out = tmp_path / "spectrum.txt"
    run(session, f"export {out}")

    comments, names, rows = _read_export(out)
    assert names == ["channel", "energy_keV", "counts", "error"]
    data = session.buffers[1]
    np.testing.assert_array_equal(rows[:, 0], np.arange(data.n_channels))
    np.testing.assert_allclose(rows[:, 1], data.spectrum.energies, atol=1e-4)
    np.testing.assert_allclose(rows[:, 2], data.spectrum.counts, atol=1e-6)
    np.testing.assert_allclose(rows[:, 3], np.sqrt(data.spectrum.counts), atol=1e-6)
    keys = [line.split()[0] for line in comments]
    assert "buffer export" in comments[0]
    assert {"Data", "Beam", "Geometry", "Conversion", "Detector", "Dose"} <= set(keys)
    assert not {"Physics", "Sample", "GOF"} & set(keys)


@needs_data
def test_export_writes_raw_counts_under_normalize(session, tmp_path):
    raw, normalized = tmp_path / "raw.txt", tmp_path / "norm.txt"
    run(session, f"expo {raw}", "normalize", f"expo {normalized}")
    np.testing.assert_array_equal(_read_export(raw)[2], _read_export(normalized)[2])


@needs_data
def test_export_csv_is_comma_separated_and_a_bare_name_gets_txt(session, tmp_path):
    run(session, f"export {tmp_path / 'spectrum.csv'}", f"export {tmp_path / 'spectrum'}")

    _, names, _ = _read_export(tmp_path / "spectrum.csv")
    assert names == ["channel", "energy_keV", "counts", "error"]
    assert (tmp_path / "spectrum.txt").exists()


@needs_data
def test_export_of_buffer_0_writes_a_fresh_simulation_and_its_sample(session, tmp_path):
    sample = _ec_sample(session, tmp_path)
    out = tmp_path / "sim.txt"
    run(session, f"sim get {sample}", "compare")
    stale = session.buffers.get(0).spectrum.counts.copy()
    session.buffers.active = 0
    run(session, "sim thick 900 /cm2", f"export {out}")

    comments, _, rows = _read_export(out)
    fresh = session.simulation().spectrum.counts
    assert not np.allclose(stale, fresh)
    np.testing.assert_allclose(rows[:, 2], fresh, atol=1e-6)
    rows_by_key = {line.split()[0]: line for line in comments}
    assert rows_by_key["Data"].endswith("buffer 0)")
    assert rows_by_key["Sample"].endswith("Si [900/cm2]")
    assert "Physics" in rows_by_key and "Areal" in rows_by_key


@needs_data
def test_export_from_sim_and_pert_stays_at_that_level(session, tmp_path):
    sample = _ec_sample(session, tmp_path)
    stack = ["rump"]
    for line in ("sim", f"get {sample}", f"export {tmp_path / 'a.txt'}"):
        execute_line(session, line, stack)
    assert stack == ["rump", "sim"]
    for line in ("return", "pert", f"export {tmp_path / 'b.txt'}"):
        execute_line(session, line, stack)
    assert stack == ["rump", "pert"]
    assert (tmp_path / "a.txt").exists() and (tmp_path / "b.txt").exists()


def test_export_abbreviations_leave_expand_and_exportcmp_alone():
    from pyrump.shell.commands.rump import TABLE, cmd_expand, cmd_export, cmd_exportcmp

    assert TABLE.match("exp").handler is cmd_expand
    assert TABLE.match("expo").handler is cmd_export
    assert TABLE.match("export").handler is cmd_export
    assert TABLE.match("exportc").handler is cmd_exportcmp


def test_newall_leaves_no_active_buffer_for_a_following_pert_go(session, tmp_path):
    """NEWALL must blank buffer 0 too, so a stale simulation left over from
    before the reset can't masquerade as PERT GO's "observed" data."""
    sample = tmp_path / "reset.lcm"
    sample.write_text(
        "Sim Reset\nLayer 1\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"
    )
    run(session, f"sim get {sample}", "pert", "thick 1", "return", "compare", "newall")
    assert session.buffers.get(0) is None
    with pytest.raises(CommandError, match="no active buffer"):
        run(session, "pert go")


def test_figsave_with_nothing_plotted_is_rejected(session):
    with pytest.raises(CommandError, match="nothing plotted yet"):
        run(session, "figsave out.png")


@needs_data
def test_figsave_writes_an_image_file(session, tmp_path):
    sample = tmp_path / "figsave.lcm"
    sample.write_text(
        "Sim Reset\nLayer 1\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"
    )
    run(session, f"sim get {sample}", "compare", f"figsave {tmp_path / 'out.png'}")
    out = tmp_path / "out_rbs.png"
    assert out.exists()
    assert out.stat().st_size > 0
    assert not (tmp_path / "out.png").exists()


@needs_data
def test_figsave_writes_png_at_300_dpi(session, tmp_path):
    from matplotlib.image import imread
    sample = tmp_path / "figsave_dpi.lcm"
    sample.write_text(
        "Sim Reset\nLayer 1\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"
    )
    out = tmp_path / "dpi_rbs.png"
    run(session, f"sim get {sample}", "compare", f"figsave {tmp_path / 'dpi.png'}")
    width_in, height_in = session.figure.get_size_inches()
    height_px, width_px = imread(out).shape[:2]
    assert (width_px, height_px) == (round(width_in * 300), round(height_in * 300))


@needs_data
def test_figsave_defaults_to_png_with_no_extension(session, tmp_path):
    sample = tmp_path / "figsave_noext.lcm"
    sample.write_text(
        "Sim Reset\nLayer 1\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"
    )
    run(session, f"sim get {sample}", "compare", f"figsave {tmp_path / 'out'}")
    assert (tmp_path / "out_rbs.png").exists()


@needs_data
def test_figsave_does_not_double_an_rbs_suffix(session, tmp_path):
    sample = tmp_path / "figsave_rbs.lcm"
    sample.write_text(
        "Sim Reset\nLayer 1\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"
    )
    run(session, f"sim get {sample}", "compare", f"figsave {tmp_path / 'fit_rbs.pdf'}")
    assert (tmp_path / "fit_rbs.pdf").exists()
    assert not (tmp_path / "fit_rbs_rbs.pdf").exists()


@needs_data
def test_hcopy_is_a_synonym_for_figsave(session, tmp_path):
    sample = tmp_path / "hcopy.lcm"
    sample.write_text(
        "Sim Reset\nLayer 1\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"
    )
    run(session, f"sim get {sample}", "compare", f"hcopy {tmp_path / 'hc.png'}")
    assert (tmp_path / "hc_rbs.png").exists()


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
def test_plot_legend_shows_a_sanitized_id_not_the_raw_buffer_label(session, tmp_path):
    """Reproduces a WRASCII macro's FILENAME line stamping a full path into
    buffer.name (see cmd_filename, and GO's own header/footer) -- the legend
    must show the same clean ID, not that path."""
    sample = tmp_path / "labels_filename.lcm"
    sample.write_text(
        "Sim Reset\nLayer 1\n Thick 500 /cm2\n Composition Si 1 /\nMaxpth 200\n"
    )
    run(session, r"filename c:\RBS\data\2026\08\MA8410.RBS", f"sim get {sample}", "compare")

    labels = session.figure.axes[0].get_legend_handles_labels()[1]
    assert labels == ["MA8410", "SIM"]


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
    assert trace_labels == ["test", "Si [500/cm2] - Mn3Pt [150A] - Ru [30A]"]

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
def test_splot_element_keeps_its_own_label_under_structlabel(session, tmp_path):
    # STRUCTLABEL's structure-of-the-whole-sample legend text is right for a
    # plain SPLOT, but must not paper over a selective one's own name.
    session.buffers.get(1).spectrum.calibration = Calibration(kevch=50.0, npt=64)
    sample = tmp_path / "splot_element_structlabel.lcm"
    sample.write_text(_SPLOT_SAMPLE)
    run(session, f"sim get {sample}", "structlabel on", "plot 1", "sim splot Ru")
    trace = session.traces[-1]
    assert trace.label == "SIM(Ru)"


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
def test_splot_element_does_not_replace_an_existing_full_sim_overlay(session, tmp_path):
    session.buffers.get(1).spectrum.calibration = Calibration(kevch=50.0, npt=64)
    sample = tmp_path / "splot_coexist.lcm"
    sample.write_text(_SPLOT_SAMPLE)
    run(session, f"sim get {sample}", "plot 1", "ov 0", "sim splot Ru")
    labels = [t.label for t in session.traces]
    assert "SIM" in labels
    assert "SIM(Ru)" in labels


@needs_data
def test_splot_element_updates_in_place_on_repeat(session, tmp_path):
    session.buffers.get(1).spectrum.calibration = Calibration(kevch=50.0, npt=64)
    sample = tmp_path / "splot_repeat.lcm"
    sample.write_text(_SPLOT_SAMPLE)
    run(session, f"sim get {sample}", "plot 1", "sim splot Ru", "sim splot Ru")
    labels = [t.label for t in session.traces]
    assert labels.count("SIM(Ru)") == 1


@needs_data
def test_splot_element_and_layer_coexist(session, tmp_path):
    session.buffers.get(1).spectrum.calibration = Calibration(kevch=50.0, npt=64)
    sample = tmp_path / "splot_both.lcm"
    sample.write_text(_SPLOT_SAMPLE)
    run(session, f"sim get {sample}", "plot 1", "sim splot Ru", "sim splot 2")
    labels = [t.label for t in session.traces]
    assert "SIM(Ru)" in labels
    assert "SIM(layer 2)" in labels


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
    assert labels == ["test", "Si [500/cm2] - Mn0.750 Pt0.250 [150A]"]

    run(session, "compfrac off", "compare")
    labels = session.figure.axes[0].get_legend_handles_labels()[1]
    assert labels == ["test", "Si [500/cm2] - Mn3Pt [150A]"]


def test_compfrac_on_shows_atomic_fraction_in_show(session, tmp_path, capsys):
    sample = tmp_path / "compfrac_show.lcm"
    sample.write_text(
        "Sim Reset\nLayer 1\n Thick 150 A\n Composition Mn 3 Pt 1 /\nMaxpth 200\n"
    )
    run(session, f"sim get {sample}", "compfrac on")
    capsys.readouterr()
    run(session, "sim show")
    assert "Mn 0.750 Pt 0.250" in capsys.readouterr().out

    run(session, "compfrac off")
    capsys.readouterr()
    run(session, "sim show")
    assert "Mn 3 Pt 1" in capsys.readouterr().out


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
def test_cursor_reports_channel_energy_and_counts(session, capsys):
    run(session, "cursor 20")
    out = capsys.readouterr().out
    counts = session.buffers.active_buffer.spectrum.counts[20]
    assert "buffer 1 (test)" in out
    assert "Channel:     20" in out
    assert "Energy:    100.0 keV" in out  # (20 + first=0) * kevch=5.0 + kev0=0.0
    assert f"Counts: {counts:10.4f}" in out


@needs_data
def test_cursor_names_whichever_buffer_is_actually_active(session, capsys):
    """OVERLAY/SPLOT never change the active buffer, so CURSOR keeps reading
    whatever PLOT last pointed it at -- printing which one that is up front,
    so a stack of overlays never leaves that a guess."""
    session.buffers.load(make_buffer(name="second"), 2)
    run(session, "plot 2", "cursor 5")
    out = capsys.readouterr().out
    assert "buffer 2 (second)" in out


@needs_data
def test_cursor_snaps_a_fractional_channel_to_the_nearest_sample(session, capsys):
    run(session, "cursor 20.4")
    assert "Channel:     20" in capsys.readouterr().out
    run(session, "cursor 20.6")
    assert "Channel:     21" in capsys.readouterr().out


@needs_data
def test_cursor_rejects_a_channel_outside_the_buffer(session):
    with pytest.raises(CommandError, match="outside the buffer"):
        run(session, "cursor -1")
    with pytest.raises(CommandError, match="outside the buffer"):
        run(session, "cursor 9999")


@needs_data
def test_cursor_reports_yield_when_normalized(session, capsys):
    run(session, "normalize", "cursor 20")
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


# -- LIVE: the simulation on the plot follows the sample -------------------------


def _live_sample(session, tmp_path, text: str = _SPLOT_SAMPLE) -> None:
    # Wide enough to catch the edges (see the SPLOT tests above).
    session.buffers.get(1).spectrum.calibration = Calibration(kevch=50.0, npt=64)
    sample = tmp_path / "live.lcm"
    sample.write_text(text)
    run(session, f"sim get {sample}")


def _drawn(session, label: str) -> np.ndarray:
    """The y data of the curve labelled ``label`` on the (top) plot panel."""
    handles, labels = session.figure.axes[0].get_legend_handles_labels()
    return np.asarray(handles[labels.index(label)].get_ydata()).copy()


@needs_data
def test_live_redraws_the_simulation_after_a_sim_edit(session, tmp_path):
    _live_sample(session, tmp_path)
    run(session, "plot 0")
    before = _drawn(session, "SIM")

    run(session, "sim", "layer 1", "thick 60 A")

    assert session.traces[0].buffer is session.buffers.get(0)
    assert not np.array_equal(_drawn(session, "SIM"), before)


@needs_data
def test_live_keeps_the_comparison_and_updates_its_residuals(session, tmp_path):
    _live_sample(session, tmp_path, _COMPARE_SAMPLE)
    run(session, "compare")

    def chi_square_text():
        return [t.get_text() for t in session.figure.axes[1].texts if "chi-square" in t.get_text()]

    before = chi_square_text()

    run(session, "sim", "layer 1", "thick 60 A")

    assert len(session.figure.axes) == 2
    assert session.traces[1].buffer is session.buffers.get(0)
    assert chi_square_text() != before


@needs_data
def test_live_recomputes_once_for_a_whole_xeq(session, tmp_path, monkeypatch):
    import pyrump.sim.engine as engine

    _live_sample(session, tmp_path)
    run(session, "plot 0")
    calls = []
    real = engine.simulate
    monkeypatch.setattr(engine, "simulate", lambda *a, **k: calls.append(1) or real(*a, **k))
    macro = tmp_path / "edits.cmd"
    macro.write_text("sim\nlayer 1\nthick 40 A\nthick 50 A\nthick 60 A\nreturn\n")

    run(session, f"xeq {macro}")

    assert len(calls) == 1
    assert session.traces[0].buffer is session.buffers.get(0)


@needs_data
def test_live_recomputes_splot_curves_and_drops_a_vanished_element(session, tmp_path, capsys):
    _live_sample(session, tmp_path)
    run(session, "plot 1", "sim splot Mn", "sim splot 1")
    layer_before = _drawn(session, "SIM(layer 1)")

    run(session, "sim", "layer 1", "thick 60 A")
    assert not np.array_equal(_drawn(session, "SIM(layer 1)"), layer_before)

    capsys.readouterr()
    run(session, "sim", "layer 2", "composition Pt 1 /")
    labels = session.figure.axes[0].get_legend_handles_labels()[1]
    assert labels == ["test", "SIM(layer 1)"]
    assert "SPLOT Mn removed" in capsys.readouterr().out


@needs_data
def test_live_never_reopens_a_closed_window(session, tmp_path):
    import matplotlib.pyplot as plt

    _live_sample(session, tmp_path)
    run(session, "plot 0")
    plt.close(session.figure)
    open_before = plt.get_fignums()

    run(session, "sim", "layer 1", "thick 60 A")

    assert plt.get_fignums() == open_before


@needs_data
def test_live_off_leaves_the_plot_until_live_is_back_on(session, tmp_path, capsys):
    _live_sample(session, tmp_path)
    run(session, "plot 0", "live off")
    assert "live off" in capsys.readouterr().out
    shown = session.traces[0].buffer

    run(session, "sim", "layer 1", "thick 60 A")
    assert session.traces[0].buffer is shown

    run(session, "live")
    assert session.traces[0].buffer is session.buffers.get(0)
    assert session.traces[0].buffer is not shown


@needs_data
def test_live_works_at_the_sim_prompt_without_leaving_it(session, tmp_path):
    stack = ["rump"]
    for line in ("sim", "live off"):
        execute_line(session, line, stack)
    assert stack == ["rump", "sim"]
    assert session.live is False


@needs_data
def test_replot_brings_the_simulation_up_to_date_even_with_live_off(session, tmp_path):
    _live_sample(session, tmp_path)
    run(session, "live off", "plot 0")
    before = _drawn(session, "SIM")
    run(session, "sim", "layer 1", "thick 60 A")
    assert np.array_equal(_drawn(session, "SIM"), before)

    run(session, "replot")

    assert session.traces[0].buffer is session.buffers.get(0)
    assert not np.array_equal(_drawn(session, "SIM"), before)


@needs_data
def test_live_reports_a_broken_sample_once_and_keeps_the_plot(session, tmp_path, capsys):
    _live_sample(session, tmp_path)
    run(session, "plot 0")
    before = _drawn(session, "SIM")
    capsys.readouterr()

    run(session, "sim", "reset", "reset")

    assert capsys.readouterr().out.count("simulation: no sample described") == 1
    assert np.array_equal(_drawn(session, "SIM"), before)
