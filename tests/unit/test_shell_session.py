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

from pyrump.model.spectrum import Calibration, Spectrum  # noqa: E402
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
    run(session, "integral 0 63")
    assert "integral 100.0 counts" in capsys.readouterr().out


@needs_data
def test_integral_rejects_channels_outside_the_spectrum(session):
    with pytest.raises(CommandError, match="outside"):
        run(session, "integral 0 9999")


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
