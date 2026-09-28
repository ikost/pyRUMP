"""Reading a command line at the prompt.

:func:`~pyrump.shell.repl._read_line` has two paths, and which one runs is
decided by ``sys.stdin.isatty()``. On a terminal it is plain :func:`input`, so
readline keeps history and completion and the GUI toolkit's own input hook
keeps the plot window alive. Off a terminal -- a pipe, which is how mintty,
Git Bash and IDE run panels present stdin on Windows -- that hook is never
reached, so the read moves to a thread and the wait is spent pumping the
figure instead of blocking on it.
"""

from __future__ import annotations

import threading

import pytest

from pyrump.shell import plotting, repl
from pyrump.shell.session import Session


@pytest.fixture
def session() -> Session:
    return Session(table=None, registry=None, densities=None, data=None)


@pytest.fixture(autouse=True)
def _fresh_reader():
    """The stdin reader thread is module state, cached across calls -- each
    test needs its own, or it reads the previous test's fake stdin."""
    repl._stdin_lines = None
    yield
    repl._stdin_lines = None


class FakeStdin:
    """A non-tty stdin handing out fixed lines, then EOF.

    ``gate``, when given, holds the first ``readline`` until the test releases
    it, which is how the pumping-while-waiting test forces the wait loop to
    time out at least once without depending on timing.
    """

    def __init__(self, lines, gate=None):
        self._lines = list(lines)
        self._gate = gate

    def isatty(self) -> bool:
        return False

    def readline(self) -> str:
        if self._gate is not None:
            self._gate.wait(5)
            self._gate = None
        return self._lines.pop(0) if self._lines else ""


def test_a_terminal_still_goes_through_input(session, monkeypatch):
    """On a tty nothing changes: readline's editing, history and completion
    are worth more than pumping, and the toolkit's input hook covers the
    window anyway."""
    monkeypatch.setattr(repl.sys, "stdin", type("T", (), {"isatty": lambda self: True})())
    monkeypatch.setattr("builtins.input", lambda prompt="": f"typed after {prompt!r}")

    assert repl._read_line(session, "Feed me! ") == "typed after 'Feed me! '"


def test_a_pipe_reads_lines_in_order(session, monkeypatch, capsys):
    monkeypatch.setattr(repl.sys, "stdin", FakeStdin(["PLOT\n", "QUIT\n"]))

    assert repl._read_line(session, "> ") == "PLOT"
    assert repl._read_line(session, "> ") == "QUIT"
    assert capsys.readouterr().out == "> > "


def test_a_pipe_strips_a_windows_line_ending(session, monkeypatch):
    """A macro or here-doc written on Windows arrives with CRLF when stdin
    escaped newline translation; the CR is not part of the command."""
    monkeypatch.setattr(repl.sys, "stdin", FakeStdin(["REPLOT\r\n"]))

    assert repl._read_line(session, "> ") == "REPLOT"


def test_end_of_input_raises_eof_every_time(session, monkeypatch):
    """The caller's loop treats EOFError as "stdin closed, leave" -- if the
    sentinel were consumed once, a second prompt would block forever."""
    monkeypatch.setattr(repl.sys, "stdin", FakeStdin([]))

    with pytest.raises(EOFError):
        repl._read_line(session, "> ")
    with pytest.raises(EOFError):
        repl._read_line(session, "> ")


def test_the_plot_window_is_pumped_while_the_prompt_waits(session, monkeypatch):
    """The whole point of the non-tty path: a session sitting at the prompt
    with a figure open must keep that figure's event loop running."""
    gate = threading.Event()
    pumps = []

    def fake_pump(_session):
        pumps.append(_session)
        gate.set()  # let the fake stdin deliver its line now

    monkeypatch.setattr(plotting, "pump", fake_pump)
    monkeypatch.setattr(repl.sys, "stdin", FakeStdin(["PLOT\n"], gate=gate))

    assert repl._read_line(session, "> ") == "PLOT"
    assert pumps and pumps[0] is session
