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
import random
import sys
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
    from .commands import pert, rump, sim, system

    return {
        "rump": rump.TABLE,
        "sim": sim.TABLE,
        "pert": pert.TABLE,
        "system": system.TABLE,
    }[mode]


def prompt_for(session, stack: list[str], plain: bool) -> str:
    mode = stack[-1]
    if mode == "sim":
        return "SIM Command: "
    if mode == "pert":
        return "PERT Command: "
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
    commands, mirroring the main loop at rump.c:283-302 where ``LexSystem`` is
    tried after every RUMP table has missed.
    """
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

    # The system tier sits below every mode. Reaching it means no RUMP-level
    # table matched, which in the C also means SIM/PERT have been left behind.
    command = tables_for("system").match(name)
    if command is not None:
        del stack[1:]
        _invoke(session, command, rest, stack)
        return

    raise CommandError(f"unrecognized command: {name}")


def execute_file(session: Session, path: Path, stack: list[str] | None = None) -> None:
    """Run a macro file. Aborts at the first failing line, naming it."""
    path = Path(path)
    if not path.exists() and not path.suffix:
        path = path.with_suffix(".cmd")
    if not path.exists():
        raise CommandError(f"no such command file: {path}")

    depth = session.xeq_depth
    if depth >= MAX_XEQ_DEPTH:
        raise CommandError(f"XEQ nested more than {MAX_XEQ_DEPTH} deep: {path}")

    if stack is None:
        stack = ["rump"]

    frame = XeqFrame(lines=path.read_text().splitlines())
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
            options = _path_completions(text)
        else:
            options = tables_for(stack[-1]).completions(text)
            if stack[-1] != "rump":
                options += tables_for("rump").completions(text)
            options += tables_for("system").completions(text)
        return options[state] if state < len(options) else None

    readline.set_completer(complete)
    readline.set_completer_delims(" \t\n")
    readline.parse_and_bind("tab: complete")


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


def run_shell(
    data: str | None = None,
    macro: Path | None = None,
    *,
    norc: bool = False,
    batch: bool = False,
    plain_prompt: bool = False,
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

    print(f"pyRUMP interactive shell -- tables from {session.data}")
    print("Type ? for commands, QUIT to leave.")

    stack = ["rump"]

    if not norc and RC_FILE.exists():
        try:
            execute_file(session, RC_FILE, stack)
        except CommandError as error:
            print(f"{RC_FILE}: {error}", file=sys.stderr)

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

    _setup_readline(session, stack)
    try:
        while True:
            try:
                line = input(prompt_for(session, stack, plain_prompt))
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
