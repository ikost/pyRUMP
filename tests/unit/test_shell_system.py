"""The general system command tier: navigation, listing, viewing, macros.

Ported from ``LexSystem`` (C-code/lexp/system.c:175), which RUMP's main loop
falls through to after its own tables miss (rump.c:283-302).

Everything here is filesystem-portable -- ``tmp_path`` and ``pathlib`` only, no
shell, no hard-coded separators -- because these commands are expected to behave
identically on Linux, macOS and Windows.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

from pyrump.shell.commands import system
from pyrump.shell.commands.rump import TABLE as RUMP_TABLE
from pyrump.shell.commands.system import TABLE as SYSTEM_TABLE
from pyrump.shell.dispatch import CommandError
from pyrump.shell.repl import execute_line
from pyrump.shell.session import Session


@pytest.fixture(autouse=True)
def _restore_cwd():
    """These commands really chdir, so put the process back afterwards."""
    origin = Path.cwd()
    yield
    os.chdir(origin)


@pytest.fixture
def session() -> Session:
    """A session with no data tables -- none of this tier needs them."""
    return Session(table=None, registry=None, densities=None, data=Path("."))


def run(session: Session, *lines: str, stack: list[str] | None = None) -> list[str]:
    stack = ["rump"] if stack is None else stack
    for line in lines:
        execute_line(session, line, stack)
    return stack


@pytest.fixture
def tree(tmp_path: Path) -> Path:
    """A small directory tree to navigate."""
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "a.rbs").write_bytes(b"\0" * 8)
    (tmp_path / "data" / "b.rbs").write_bytes(b"\0" * 16)
    (tmp_path / "data" / "notes.txt").write_text("line one\nline two\n")
    (tmp_path / "empty").mkdir()
    return tmp_path


# -- abbreviations ---------------------------------------------------------


@pytest.mark.parametrize(
    "token, expected",
    [
        ("ls", "LS"),
        ("dir", "DIRECTORY"),
        ("direct", "DIRECTORY"),
        ("sl", "SL"),
        ("ll", "LL"),
        ("cd", "CD"),
        ("chdir", "CHDIR"),
        ("pushd", "PUSHDIR"),
        ("pushdir", "PUSHDIR"),
        ("popd", "POPDIR"),
        ("whe", "WHERE"),
        ("pwd", "PWD"),
        ("ty", "TYPE"),
        ("cat", "CAT"),
        ("more", "MORE"),
        ("cls", "CLS"),
        ("xe", "XEQ"),          # minlen 2 in the C, not 3 (system.c:201)
        ("call", "CALL"),
        ("exe", "EXECUTE"),
        ("echo", "ECHO"),
        ("quiet", "QUIET"),
        ("script", "SCRIPT"),
        ("logf", "LOGFILE"),
        ("rec", "RECORD"),
    ],
)
def test_system_abbreviations(token, expected):
    assert SYSTEM_TABLE.match(token).name == expected


@pytest.mark.parametrize("token", ["c", "l", "p", "pu", "pop", "wh", "sc", "lo"])
def test_too_short_system_tokens_do_not_match(token):
    assert SYSTEM_TABLE.match(token) is None


def test_synonyms_stay_out_of_the_listing():
    listed = {command.name for command in SYSTEM_TABLE.visible()}
    assert {"DIRECTORY", "LS", "CD", "WHERE", "TYPE", "XEQ", "SCRIPT"} <= listed
    assert not listed & {"SL", "LL", "CHDIR", "PWD", "CAT", "MORE", "RECORD"}


def test_log_still_means_the_log_axis_not_the_logfile():
    """The collision the C avoided by requiring four characters for LOGFILE."""
    assert RUMP_TABLE.match("log").name == "LOG"
    assert RUMP_TABLE.match("logf") is None
    assert SYSTEM_TABLE.match("logf").name == "LOGFILE"


def test_display_beats_directory_because_rump_is_tried_first(session, tree):
    """'dis' is RUMP's DISPLAY; 'dir' falls through to the system tier."""
    assert RUMP_TABLE.match("dis").name == "DISPLAY"
    assert RUMP_TABLE.match("dir") is None


# -- pwd / cd --------------------------------------------------------------


def test_pwd_reports_the_working_directory(session, tree, capsys):
    os.chdir(tree)
    run(session, "pwd")
    assert str(Path.cwd()) in capsys.readouterr().out


def test_cd_changes_directory(session, tree):
    os.chdir(tree)
    run(session, "cd data")
    assert Path.cwd() == (tree / "data").resolve()


def test_cd_accepts_an_absolute_path(session, tree):
    run(session, f"cd {tree / 'data'}")
    assert Path.cwd() == (tree / "data").resolve()


def test_cd_with_no_argument_goes_home(session, tree):
    """system.c:447-453 -- bare CD falls back to HOME."""
    os.chdir(tree)
    run(session, "cd")
    assert Path.cwd() == Path.home().resolve()


def test_cd_expands_a_tilde(session, tree):
    os.chdir(tree)
    run(session, "cd ~")
    assert Path.cwd() == Path.home().resolve()


def test_cd_to_a_missing_directory_is_a_clean_error(session, tree):
    os.chdir(tree)
    with pytest.raises(CommandError, match="no such directory"):
        run(session, "cd nowhere")


def test_cd_onto_a_file_is_a_clean_error(session, tree):
    os.chdir(tree)
    with pytest.raises(CommandError, match="not a directory|no such directory"):
        run(session, "cd data/a.rbs")


# -- pushdir / popdir ------------------------------------------------------


def test_pushdir_and_popdir_round_trip(session, tree):
    os.chdir(tree)
    origin = Path.cwd()
    run(session, "pushdir data")
    assert Path.cwd() == (tree / "data").resolve()
    run(session, "popdir")
    assert Path.cwd() == origin


def test_pushdir_nests(session, tree):
    os.chdir(tree)
    origin = Path.cwd()
    run(session, "pushdir data", "pushdir ..", "popdir")
    assert Path.cwd() == (tree / "data").resolve()
    run(session, "popdir")
    assert Path.cwd() == origin


def test_popdir_on_an_empty_stack_errors(session, tree):
    os.chdir(tree)
    with pytest.raises(CommandError, match="nothing on the directory stack"):
        run(session, "popdir")


def test_a_failed_pushdir_does_not_push(session, tree):
    os.chdir(tree)
    with pytest.raises(CommandError):
        run(session, "pushdir nowhere")
    assert session.directory_stack == []


# -- ls --------------------------------------------------------------------


def test_ls_lists_the_working_directory(session, tree, capsys):
    os.chdir(tree / "data")
    run(session, "ls")
    out = capsys.readouterr().out
    assert "a.rbs" in out and "b.rbs" in out and "notes.txt" in out


def test_ls_expands_a_glob_itself(session, tree, capsys):
    """No OS shell is involved, so this works identically on Windows."""
    os.chdir(tree / "data")
    run(session, "ls *.rbs")
    out = capsys.readouterr().out
    assert "a.rbs" in out and "b.rbs" in out
    assert "notes.txt" not in out


def test_ls_takes_a_directory_argument(session, tree, capsys):
    os.chdir(tree)
    run(session, "ls data")
    assert "a.rbs" in capsys.readouterr().out


def test_ls_marks_directories(session, tree, capsys):
    os.chdir(tree)
    run(session, "ls")
    assert "data/" in capsys.readouterr().out


def test_ls_of_an_empty_directory(session, tree, capsys):
    os.chdir(tree)
    run(session, "ls empty")
    assert "empty" in capsys.readouterr().out


def test_ls_of_a_missing_name_errors(session, tree):
    os.chdir(tree)
    with pytest.raises(CommandError, match="no such file or directory"):
        run(session, "ls nowhere.rbs")


def test_ls_with_a_glob_matching_nothing_is_not_an_error(session, tree, capsys):
    os.chdir(tree / "data")
    run(session, "ls *.zzz")
    assert "empty" in capsys.readouterr().out


def test_ll_shows_sizes(session, tree, capsys):
    os.chdir(tree / "data")
    run(session, "ll *.rbs")
    out = capsys.readouterr().out
    assert "16" in out and "a.rbs" in out


# -- type ------------------------------------------------------------------


def test_type_prints_a_file_without_blocking(session, tree, capsys):
    """Under pytest there is no tty, so it must never wait for a keypress."""
    os.chdir(tree / "data")
    run(session, "type notes.txt")
    out = capsys.readouterr().out
    assert "line one" in out and "line two" in out
    assert "-- more --" not in out


def test_type_on_a_directory_errors(session, tree):
    os.chdir(tree)
    with pytest.raises(CommandError, match="is a directory"):
        run(session, "type data")


def test_type_on_a_missing_file_errors(session, tree):
    os.chdir(tree)
    with pytest.raises(CommandError, match="no such file"):
        run(session, "cat nowhere.txt")


# -- cls -------------------------------------------------------------------


def test_cls_uses_ansi_when_available(session, capsys, monkeypatch):
    monkeypatch.setattr(system, "enable_ansi", lambda: True)
    run(session, "cls")
    assert "\x1b[2J" in capsys.readouterr().out


def test_cls_falls_back_to_the_windows_console(session, monkeypatch):
    """The pre-Windows-10 path, exercised without a real console."""
    called = []
    monkeypatch.setattr(system, "enable_ansi", lambda: False)
    monkeypatch.setattr(system.os, "system", lambda cmd: called.append(cmd))
    run(session, "cls")
    assert called == ["cls"]


def test_enable_ansi_is_true_off_windows(monkeypatch):
    monkeypatch.setattr(system.os, "name", "posix")
    assert system.enable_ansi() is True


# -- macros and session logging -------------------------------------------


def test_script_records_a_replayable_macro(session, tmp_path, capsys):
    log = tmp_path / "session.cmd"
    run(session, f"script {log}", "pwd", "script off")
    assert "pwd" in log.read_text().splitlines()


def test_logfile_and_record_are_synonyms_for_script(session, tmp_path):
    log = tmp_path / "a.cmd"
    run(session, f"logfile {log}", "pwd", "script off")
    assert "pwd" in log.read_text()
    other = tmp_path / "b.cmd"
    run(session, f"record {other}", "pwd", "script off")
    assert "pwd" in other.read_text()


def test_echo_and_quiet(session, capsys):
    run(session, "echo")
    assert session.echo is True
    run(session, "quiet")
    assert session.echo is False


def test_xeq_runs_a_macro_at_minlen_two(session, tree, tmp_path):
    macro = tmp_path / "m.cmd"
    macro.write_text(f"cd {tree / 'data'}\n")
    run(session, f"xe {macro}")
    assert Path.cwd() == (tree / "data").resolve()


# -- interaction with the mode stack --------------------------------------


def test_system_commands_are_reachable_from_sim(session, tree):
    """SIM does not know PWD, so RUMP is returned to -- then LexSystem runs it."""
    os.chdir(tree)
    stack = ["rump", "sim"]
    run(session, "cd data", stack=stack)
    assert stack == ["rump"]
    assert Path.cwd() == (tree / "data").resolve()


def test_an_unknown_command_still_errors(session):
    with pytest.raises(CommandError, match="unrecognized command"):
        run(session, "wiggle")


# -- the regression this change is really about ----------------------------


def _write_spectrum(path: Path) -> None:
    """A minimal valid .rbs, so GET works without the atomic data tables."""
    import numpy as np

    from pyrump.io.rbs import RbsSpectrum, write_rbs
    from pyrump.model.detector import Measurement
    from pyrump.model.geometry import Geometry
    from pyrump.model.spectrum import Calibration

    write_rbs(
        path,
        RbsSpectrum(
            counts=np.arange(64, dtype=float),
            calibration=Calibration(npt=64),
            geometry=Geometry(theta=0.0, phi=10.0),
            measurement=Measurement(),
            identifier="regression fixture",
        ),
    )


def test_the_same_file_reached_by_different_relative_paths_is_one_buffer(
    session, tmp_path
):
    """Without absolute Buffer.path, CD makes GET load the same file twice.

    ``find_path`` resolves against the working directory, so a buffer holding
    the relative "data/a.rbs" stops matching once that directory *is* the
    working directory -- and the file is silently read into a second buffer.
    """
    (tmp_path / "data").mkdir()
    _write_spectrum(tmp_path / "data" / "a.rbs")

    os.chdir(tmp_path)
    run(session, "get data/a.rbs")
    assert session.buffers.active == 1

    run(session, "cd data", "get a.rbs")
    assert session.buffers.active == 1
    assert session.buffers.get(2) is None, "the file was loaded a second time"

    # And by yet another spelling of the same path.
    run(session, "get ./a.rbs")
    assert session.buffers.get(2) is None


def test_buffer_listings_stay_short_after_moving(session, tmp_path, capsys):
    (tmp_path / "data").mkdir()
    _write_spectrum(tmp_path / "data" / "a.rbs")
    os.chdir(tmp_path)
    run(session, "get data/a.rbs")

    capsys.readouterr()
    run(session, "buffers")
    assert "data" in capsys.readouterr().out

    run(session, "cd data")
    capsys.readouterr()
    run(session, "buffers")
    # Relative to the new working directory, it is just the file name.
    assert "a.rbs" in capsys.readouterr().out
