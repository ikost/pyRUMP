"""The read-eval-print loop, the mode stack, and macro execution.

RUMP nests three command levels -- RUMP, SIM and PERT -- each with its own
prompt. A command the sub-level does not recognise falls through to its parent
and the sub-level is left automatically (sim.htm: "RUMP will be returned to
automatically on a command it does not understand"). That fall-through is
reproduced here by :func:`execute_line` walking the mode stack outwards.

``XEQ`` feeds a file through exactly the same :func:`execute_line`, so a macro
and a typed session cannot drift apart.
"""

from __future__ import annotations

import os
import queue
import random
import sys
import threading
from pathlib import Path

from .dispatch import ArgReader, CommandError, strip_comment, tokenize
from .session import Session, XeqFrame

#: RUMP rotates through these (rump.c:128). tests/oracle/driver.py matches them
#: when driving the real binary, so they are worth keeping verbatim.
PROMPTS = (
    "Your wish? ", "Ready if you are! ", "Yes Master? ",
    "You called? ", "Yes dear? ", "At your service! ",
    "Beam me up Scottie! ", "Up periscope! ", "Feed me! ",
    "Here I am! ", "Next? ", "Hey, man, what next? ",
    "Next command: ", "Whoopee! ", "Go for it! ",
    "I want a cookie!! ", "Igen, uram! ",
)

#: How deep XEQ may nest before we assume a macro calls itself.
MAX_XEQ_DEPTH = 20

HISTORY = Path.home() / ".pyrump_history"
RC_FILE = Path.home() / ".pyrumprc"


def tables_for(mode: str):
    """The command table for a mode name."""
    from .commands import pert, pixe, rump, sim, system

    return {
        "rump": rump.TABLE,
        "sim": sim.TABLE,
        "pert": pert.TABLE,
        "pixe": pixe.TABLE,
        "system": system.TABLE,
    }[mode]


def prompt_for(session, stack: list[str], plain: bool) -> str:
    mode = stack[-1]
    if mode == "sim":
        return "SIM Command: "
    if mode == "pert":
        return "PERT Command: "
    if mode == "pixe":
        return "PIXE Command: "
    return "pyrump> " if plain else random.choice(PROMPTS)


def _invoke(session: Session, command, rest: list[str], stack: list[str]) -> None:
    """Run one resolved command, translating its exceptions for the caller."""
    from .commands.rump import EnterMode, Return

    args = ArgReader(rest, command=command.name.lower())
    try:
        command.handler(session, args)
    except EnterMode as entered:
        stack.append(entered.name)
    except Return:
        if len(stack) > 1:
            stack.pop()
    except CommandError:
        raise
    except KeyError as error:
        raise CommandError(str(error).strip("'")) from None
    except (ValueError, OSError) as error:
        raise CommandError(f"{command.name.lower()}: {error}") from None


def execute_line(session: Session, line: str, stack: list[str]) -> None:
    """Run one command line against the innermost mode that understands it.

    Falls outwards through the mode stack, popping levels as it goes -- RUMP's
    automatic return from SIM and PERT -- and finally to the general system
    commands, which run in place without changing the stack (see the note
    below).

    A line typed at the prompt -- not one inside a macro, so a whole XEQ
    counts as one -- ends by redrawing the simulation if it changed (LIVE,
    :func:`~pyrump.shell.plotting.follow_simulation`), even if it failed
    part-way.
    """
    try:
        _execute_line(session, line, stack)
    finally:
        if session.xeq_depth == 0:
            from . import plotting

            plotting.follow_simulation(session)


def _execute_line(session: Session, line: str, stack: list[str]) -> None:
    text = strip_comment(line)
    if not text:
        return
    tokens = tokenize(text)
    if not tokens:
        return

    session.write_log(text)
    if session.echo:
        print(f"> {text}")

    name, rest = tokens[0], tokens[1:]

    for depth in range(len(stack) - 1, -1, -1):
        command = tables_for(stack[depth]).match(name)
        if command is None:
            continue
        del stack[depth + 1 :]  # auto-return out of the levels we fell through
        _invoke(session, command, rest, stack)
        return

    # The system tier sits below every mode, and unlike a mode's own tables,
    # matching here does not leave SIM/PERT -- sim2.c:473 and pert2.c:198 both
    # just "continue" their loop after a general command runs, the same way
    # the top-level RUMP loop does at rump.c:283-291.
    command = tables_for("system").match(name)
    if command is not None:
        _invoke(session, command, rest, stack)
        return

    raise CommandError(f"unrecognized command: {name}")


#: Extensions tried, in order, for a bare (extension-less) XEQ/CALL/EXECUTE
#: argument: RUMP's own list (lexp.c:168, ``LexMacroExtList``) -- ``.xeq``
#: is what SNAPSHOT writes, being no file type Windows would run, as it does
#: ``.cmd`` -- then ``.rbs``/``.RBS``, because real RBS acquisition software
#: often writes its output as a plain ``EMPTY``/``SWALLOW`` command macro
#: under that extension (see the README's Buffers section) -- meant to be
#: replayed with XEQ, not read with GET, despite the misleading name.
MACRO_EXTENSIONS = (".mac", ".MAC", ".xeq", ".XEQ", ".cmd", ".CMD", ".rbs", ".RBS")


def execute_file(session: Session, path: Path, stack: list[str] | None = None) -> Path:
    """Run a macro file. Aborts at the first failing line, naming it.

    Returns the file actually run (a bare name gains an extension, and one
    named inside another macro may be found in that macro's folder).
    """
    found = session.locate(Path(path), MACRO_EXTENSIONS)
    if found is None:
        raise CommandError(f"no such command file: {path}")
    path = found

    depth = session.xeq_depth
    if depth >= MAX_XEQ_DEPTH:
        raise CommandError(f"XEQ nested more than {MAX_XEQ_DEPTH} deep: {path}")

    if stack is None:
        stack = ["rump"]

    try:
        text = path.read_text()
    except UnicodeDecodeError as error:
        raise CommandError(
            f"{path} is not a text command file ({error}) -- "
            "binary spectrum data belongs with GET, not XEQ"
        ) from None

    frame = XeqFrame(lines=text.splitlines(), directory=path.resolve().parent)
    session.xeq_depth = depth + 1
    session.xeq_stack.append(frame)
    try:
        while frame.index < len(frame.lines):
            number = frame.index + 1
            raw = frame.lines[frame.index]
            frame.index += 1
            try:
                execute_line(session, raw, stack)
            except CommandError as error:
                raise CommandError(f"{path}:{number}: {error}") from None
    finally:
        session.xeq_depth = depth
        session.xeq_stack.pop()
    return path


def _setup_readline(session, stack: list[str]) -> None:
    """History and tab-completion, when readline is available."""
    try:
        import readline
    except ImportError:  # pragma: no cover - Windows without pyreadline
        return

    try:
        readline.read_history_file(HISTORY)
    except (OSError, ValueError):
        pass
    readline.set_history_length(2000)

    def complete(text: str, state: int):
        buffer = readline.get_line_buffer()
        if buffer[: readline.get_begidx()].strip():
            # Past the command word: complete paths, which is what GET, CD, XEQ
            # and SIM GET all want.
            options = _argument_completions(
                buffer[: readline.get_endidx()], readline.get_begidx()
            )
        else:
            options = tables_for(stack[-1]).completions(text)
            if stack[-1] != "rump":
                options += tables_for("rump").completions(text)
            options += tables_for("system").completions(text)
        return options[state] if state < len(options) else None

    readline.set_completer(complete)
    readline.set_completer_delims(" \t\n")
    is_libedit = getattr(readline, "backend", None) == "editline" or "libedit" in (
        readline.__doc__ or ""
    )
    if is_libedit:
        # macOS ships Python's readline module backed by libedit, which uses
        # editline's bind syntax instead of GNU readline's "tab: complete".
        readline.parse_and_bind("bind ^I rl_complete")
    else:
        readline.parse_and_bind("tab: complete")


def _argument_start(line: str) -> tuple[int, str | None]:
    """Where the last argument of ``line`` starts, and its quote if still open.

    Walks the line the way :func:`~pyrump.shell.dispatch.tokenize` does, so
    ``"My Da`` (an open quote) and ``My\\ Da`` (an escaped space) are each one
    argument. A line ending in whitespace starts a new, empty one.
    """
    index, length = 0, len(line)
    start = length
    while index < length:
        while index < length and line[index].isspace():
            index += 1
        if index >= length:
            return length, None
        start = index
        if line[index] in ("'", '"'):
            end = line.find(line[index], index + 1)
            if end == -1:
                return start, line[index]
            index = end + 1
        else:
            while index < length and not line[index].isspace():
                index += 2 if line.startswith("\\ ", index) else 1
    return start, None


def _argument_completions(line: str, begidx: int) -> list[str]:
    """Path completions for the argument at the end of ``line`` (the buffer up
    to the cursor), as replacements for ``line[begidx:]``.

    Readline splits words at every space (``set_completer_delims``), so for a
    name with spaces ``begidx`` falls mid-argument; the whole argument is
    completed here and only the part after ``begidx`` handed back. A name
    with spaces completes as ``"My Data/`` -- the quote left open on a
    directory so completion can carry on into it, which ``tokenize`` accepts
    as running to end of line -- or as ``My\\ Data/`` if that is how the user
    started writing it.
    """
    start, quote = _argument_start(line)
    typed = line[start:]
    raw = typed[1:] if quote else typed.replace("\\ ", " ")
    escaped = quote is None and "\\ " in typed

    def spelt(path: str) -> str:
        if escaped:
            return path.replace(" ", "\\ ")
        opening = quote or ('"' if " " in path else "")
        closing = "" if not opening or path.endswith(("/", os.sep)) else opening
        return opening + path + closing

    already = line[start:begidx]
    return [
        full[len(already):]
        for full in map(spelt, _path_completions(raw))
        if full.startswith(already)
    ]


def _path_completions(text: str) -> list[str]:
    """Filesystem completions for a partially typed path."""
    try:
        expanded = Path(text).expanduser()
        if text.endswith(("/", os.sep)):
            base, prefix = expanded, ""
        else:
            base, prefix = expanded.parent, expanded.name
        if not base.is_dir():
            return []
        # Keep the directory part the user typed, so the completion substitutes
        # cleanly into the line.
        head = text[: len(text) - len(prefix)]
        return sorted(
            head + child.name + ("/" if child.is_dir() else "")
            for child in base.iterdir()
            if child.name.startswith(prefix)
        )
    except OSError:
        return []


def _save_history() -> None:
    try:
        import readline

        readline.write_history_file(HISTORY)
    except (ImportError, OSError):
        pass


#: How long :func:`_read_line` waits on stdin before pumping the plot window
#: again. ``plotting.pump`` self-throttles, so this only sets how finely the
#: wait is sliced, not how often the GUI is actually touched.
_WAIT_SLICE = 0.05

#: Lines the reader thread has pulled off a non-tty stdin, plus ``None`` once
#: it hits EOF. Created on first use by :func:`_read_line`.
_stdin_lines: "queue.Queue[str | None] | None" = None


def _start_stdin_reader() -> "queue.Queue[str | None]":
    """Read stdin on a daemon thread, so the main thread stays free to pump.

    The split is this way round, not the other, because a GUI toolkit may only
    be touched from the thread that owns it -- Cocoa enforces this outright --
    so the blocking read is what gets banished to a thread, never the pumping.
    """
    lines: "queue.Queue[str | None]" = queue.Queue()

    def reader() -> None:
        while True:
            try:
                line = sys.stdin.readline()
            except (OSError, ValueError):  # stdin closed under us
                line = ""
            if not line:
                lines.put(None)
                return
            lines.put(line.rstrip("\n").rstrip("\r"))

    threading.Thread(target=reader, daemon=True, name="pyrump-stdin").start()
    return lines


def read_prompt(prompt: str) -> str:
    """Write the prompt ourselves, rather than letting readline place it.

    Windows has no readline, so pyRUMP leans on pyreadline3 for history and
    tab completion. Unlike GNU readline or editline, pyreadline3 does not let
    Python write the prompt: it draws it through direct Win32 console calls
    and remembers the coordinates, so it can repaint the line while editing.
    In VS Code's terminal that bookkeeping drifts from what is actually on
    screen, and it strands the cursor mid-way along an earlier line -- the
    text comes out right, only the cursor is wrong.

    Handing it an empty prompt, with the text already written, makes it start
    from wherever the terminal just left the cursor instead of computing a
    position of its own. Editing, history and completion are untouched; the
    cost is that a repaint mid-edit redraws only what was typed, without the
    prompt in front of it.

    GNU readline (Linux) and editline (macOS's usual backend) have no such
    trouble and are given the prompt as before -- they need it to lay the
    line out correctly.
    """
    if sys.platform != "win32":
        return input(prompt)
    sys.stdout.write(prompt)
    sys.stdout.flush()
    return input("")


def _read_line(session, prompt: str) -> str:
    """One line from the user, keeping the plot window alive while we wait.

    On a real terminal this defers to :func:`read_prompt`, so readline (or
    pyreadline3 on Windows) provides history and tab completion, and the GUI
    toolkit's own ``PyOS_InputHook`` -- which CPython invokes only on the
    interactive readline path -- runs the event loop between keystrokes, so
    the figure stays responsive on its own.

    When stdin is *not* a tty that hook is never reached, because :func:`input`
    bypasses ``PyOS_Readline`` entirely. The window then goes unpumped for the
    whole time the prompt waits and the OS declares it hung about five seconds
    in -- the common way to hit this on Windows is running the shell under
    mintty/Git Bash or an IDE's run panel, where stdin is a pipe. There is no
    line editing to preserve in that case, so the read moves to a thread and
    the wait is spent pumping instead of blocking.

    Raises :exc:`EOFError` at end of input, like :func:`input` does, so the
    caller's handling of a closed stdin is the same either way.
    """
    global _stdin_lines

    try:
        interactive_stdin = sys.stdin is not None and sys.stdin.isatty()
    except (AttributeError, ValueError):  # pragma: no cover - exotic stdin
        interactive_stdin = False
    if interactive_stdin:
        return read_prompt(prompt)

    from . import plotting

    if _stdin_lines is None:
        _stdin_lines = _start_stdin_reader()
    sys.stdout.write(prompt)
    sys.stdout.flush()
    while True:
        try:
            line = _stdin_lines.get(timeout=_WAIT_SLICE)
        except queue.Empty:
            plotting.pump(session)
            continue
        if line is None:
            # Keep EOF sticky: a later prompt must see it too, not block.
            _stdin_lines.put(None)
            raise EOFError
        return line


def run_shell(
    data: str | None = None,
    macro: Path | None = None,
    *,
    norc: bool = False,
    batch: bool = False,
    plain_prompt: bool = False,
    faithful: str | None = None,
) -> int:
    """Start the interactive shell. Returns a process exit code."""
    from .commands.rump import Quit

    try:
        session = Session.create(data)
    except SystemExit as error:
        print(error, file=sys.stderr)
        return 1

    # Switch a Windows console into virtual-terminal mode once, so CLS and any
    # other escape sequence behave as they do on Linux and macOS.
    from .commands.system import enable_ansi

    enable_ansi()

    from .. import __version__

    print(f"pyRUMP {__version__} -- interactive shell -- tables from {session.data}")
    print("Type ? for commands, QUIT to leave.")

    stack = ["rump"]

    if not norc and RC_FILE.exists():
        try:
            execute_file(session, RC_FILE, stack)
        except CommandError as error:
            print(f"{RC_FILE}: {error}", file=sys.stderr)

    # An explicit --faithful overrides whatever ~/.pyrumprc set, for this
    # invocation only; a FAITHFUL line inside the macro itself still wins.
    if faithful is not None:
        session.settings.faithful = faithful == "on"

    if macro is not None:
        try:
            execute_file(session, Path(macro), stack)
        except CommandError as error:
            print(error, file=sys.stderr)
            if batch:
                return 1
        except Quit:
            return 0
        if batch:
            return 0
        # The macro ran at XEQ depth, so its lines never redrew: catch the
        # plot up before the first prompt, as an XEQ typed there would.
        from . import plotting

        plotting.follow_simulation(session)

    _setup_readline(session, stack)
    try:
        while True:
            try:
                line = _read_line(session, prompt_for(session, stack, plain_prompt))
            except EOFError:
                print()
                break
            except KeyboardInterrupt:
                print("^C")
                continue
            try:
                execute_line(session, line, stack)
            except Quit:
                break
            except CommandError as error:
                print(error, file=sys.stderr)
            except Exception as error:  # a bug, not a user mistake
                print(f"internal error: {type(error).__name__}: {error}", file=sys.stderr)
    finally:
        _save_history()
        if session.log_file is not None:
            session.log_file.close()
    return 0
