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
    """On a tty nothing is threaded: readline's editing, history and
    completion are worth more than pumping, and the toolkit's input hook
    covers the plot window anyway. Pinned on a non-Windows platform so it
    tests the routing, not read_prompt's own platform split."""
    monkeypatch.setattr(repl.sys, "platform", "linux")
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

# -- read_prompt: who draws the prompt ---------------------------------------


def test_the_prompt_is_written_by_us_on_windows(monkeypatch, capsys):
    """pyreadline3 must be handed an empty prompt, with the text already on
    screen -- that is the point: it never computes coordinates of its own, so
    it cannot get them wrong."""
    seen = []
    monkeypatch.setattr(repl.sys, "platform", "win32")
    monkeypatch.setattr("builtins.input", lambda prompt="": seen.append(prompt) or "PLOT")

    assert repl.read_prompt("Feed me! ") == "PLOT"
    assert seen == [""]
    assert capsys.readouterr().out == "Feed me! "


def test_gnu_readline_still_gets_the_prompt(monkeypatch, capsys):
    """GNU readline and editline place the prompt correctly and need it to
    lay the line out, so Linux and macOS keep the plain input() call."""
    seen = []
    monkeypatch.setattr(repl.sys, "platform", "linux")
    monkeypatch.setattr("builtins.input", lambda prompt="": seen.append(prompt) or "QUIT")

    assert repl.read_prompt("Yes Master? ") == "QUIT"
    assert seen == ["Yes Master? "]
    assert capsys.readouterr().out == ""


def test_a_terminal_draws_its_prompt_through_read_prompt(session, monkeypatch):
    """The tty branch of _read_line must go through read_prompt, or Windows
    loses the fix the moment it is reached from the REPL."""
    calls = []
    monkeypatch.setattr(repl.sys, "stdin", type("T", (), {"isatty": lambda self: True})())
    monkeypatch.setattr(repl, "read_prompt", lambda prompt: calls.append(prompt) or "PLOT")

    assert repl._read_line(session, "Next? ") == "PLOT"
    assert calls == ["Next? "]


# -- tab completion of paths with spaces ------------------------------------


def _complete(line: str) -> list[str]:
    """What completion offers at the end of ``line``, as readline would call
    it: ``begidx`` after the last space, since spaces are word delimiters."""
    begidx = max(line.rfind(" "), line.rfind("\t")) + 1
    return repl._argument_completions(line, begidx)


@pytest.fixture
def spaced(tmp_path, monkeypatch):
    (tmp_path / "My Data").mkdir()
    (tmp_path / "My Data" / "MnPt.lcm").write_text("")
    (tmp_path / "plain").mkdir()
    monkeypatch.chdir(tmp_path)
    return tmp_path


def test_a_directory_with_a_space_completes_inside_an_open_quote(spaced):
    assert _complete("cd My") == ['"My Data/']


def test_completion_carries_on_inside_an_open_quote(spaced):
    # readline only replaces what follows the last space: 'Da' -> 'Data/'
    assert _complete('cd "My Da') == ["Data/"]
    assert _complete('get "My Data/Mn') == ['Data/MnPt.lcm"']


def test_a_file_completes_with_its_quote_closed_and_tokenizes_whole(spaced):
    from pyrump.shell.dispatch import tokenize

    line = 'sim get "My Data/Mn'
    completed = line[: line.rfind(" ") + 1] + _complete(line)[0]
    assert tokenize(completed) == ["sim", "get", "My Data/MnPt.lcm"]


def test_completion_keeps_backslash_escaping_if_that_is_how_it_was_typed(spaced):
    assert _complete(r"cd My\ D") == ["Data/"]
    assert _complete(r"get My\ Data/M") == ["Data/MnPt.lcm"]
    assert _complete(r"cd M") == ['"My Data/']  # no escape typed yet: quote it


def test_a_name_without_spaces_completes_unquoted(spaced):
    assert _complete("cd pl") == ["plain/"]


def test_a_new_empty_argument_lists_the_directory(spaced):
    assert _complete("cd ") == ['"My Data/', "plain/"]


# -- the separator after a completed directory --------------------------------


@pytest.fixture
def windows_separators(monkeypatch):
    monkeypatch.setattr(repl.os, "sep", "\\")
    monkeypatch.setattr(repl.os, "altsep", "/")


@pytest.fixture
def years(tmp_path, monkeypatch):
    for year in ("2025", "2026"):
        (tmp_path / "data" / year).mkdir(parents=True)
    monkeypatch.chdir(tmp_path)


def test_windows_completes_a_directory_with_a_backslash(spaced, windows_separators):
    assert _complete("cd pl") == ["plain\\"]
    assert _complete("cd My") == ['"My Data\\']


def test_windows_respells_a_typed_forward_slash(years, windows_separators):
    assert _complete("cd data/202") == ["data\\2025\\", "data\\2026\\"]


def test_windows_cannot_respell_what_is_before_the_last_space(spaced, windows_separators):
    # readline replaces only from the last space on; '"My ' stays as typed
    assert _complete('get "My Data/Mn') == ['Data\\MnPt.lcm"']


def test_windows_escaped_spaces_complete_with_a_forward_slash(spaced, windows_separators):
    # "Data\" followed by a typed space would read back as an escaped space
    assert _complete(r"cd My\ D") == ["Data/"]


def test_posix_never_completes_with_a_backslash(years):
    assert _complete("cd data/202") == ["data/2025/", "data/2026/"]
