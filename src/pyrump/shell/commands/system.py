"""General system commands: navigation, file viewing, macros.

This is a port of the ``LexSystem`` tier (``C-code/lexp/system.c:175``), not a
new invention. RUMP's main loop tries its own tables and then **falls through to
this one before reporting an unrecognised command** (rump.c:283-302), so it is
reachable from every level -- which is why it lives at the bottom of the mode
stack rather than inside any one mode.

Names and minimum abbreviations are transcribed from ``cmlist`` at system.c:176.
Left out deliberately:

* the DOS drive-letter forms -- ``DCD``, ``x:``, ``CD -d`` (``CD_Help``,
  system.c:393) -- which mean nothing outside Windows drive-letter semantics;
* the shell escapes ``!`` / ``DOS`` / ``CSH`` / ``SHELL`` (system.c:186-188),
  which would let any ``.cmd`` macro run arbitrary commands.

Everything here is written against :mod:`pathlib` and expands its own globs, so
behaviour does not depend on an OS shell and is the same on Linux, macOS and
Windows.
"""

from __future__ import annotations

import fnmatch
import os
import shutil
import sys
from datetime import datetime
from pathlib import Path

from ..dispatch import ArgReader, CommandError, CommandTable

#: Lines shown before MORE pauses, when there is a terminal to pause for.
_PAGE_MARGIN = 2


def _terminal() -> os.terminal_size:
    # Documented fallback (80x24) when there is no console, e.g. under pytest.
    return shutil.get_terminal_size()


def resolve(token: str) -> Path:
    """Expand ``~`` and make a user-supplied path absolute.

    ``Path`` handles separators and drive letters per-platform, so no branching
    is needed for Windows.
    """
    return Path(token).expanduser()


# ---------------------------------------------------------------------------
# Navigation
# ---------------------------------------------------------------------------


def cmd_pwd(session, args: ArgReader) -> None:
    args.done()
    print(f"  {Path.cwd()}")


def _change_to(session, target: Path) -> None:
    try:
        os.chdir(target)
    except NotADirectoryError:
        raise CommandError(f"not a directory: {target}") from None
    except FileNotFoundError:
        raise CommandError(f"no such directory: {target}") from None
    except PermissionError:
        raise CommandError(f"permission denied: {target}") from None
    print(f"  {Path.cwd()}")


def cmd_cd(session, args: ArgReader) -> None:
    """``CD [directory]`` -- with no argument, go home (system.c:447-453)."""
    token = args.optional()
    args.done()
    _change_to(session, Path.home() if token is None else resolve(token))


def cmd_pushdir(session, args: ArgReader) -> None:
    """Remember where we are, then change directory."""
    token = args.optional()
    args.done()
    here = Path.cwd()
    _change_to(session, Path.home() if token is None else resolve(token))
    session.directory_stack.append(here)


def cmd_popdir(session, args: ArgReader) -> None:
    args.done()
    if not session.directory_stack:
        raise CommandError("nothing on the directory stack; use PUSHDIR first")
    _change_to(session, session.directory_stack.pop())


# ---------------------------------------------------------------------------
# Listing
# ---------------------------------------------------------------------------


def _entries(token: str | None) -> tuple[Path, list[Path]]:
    """Resolve an LS argument into a base directory and its matching entries.

    The argument may be a directory, a glob, or absent. Globs are expanded here
    rather than by a shell, so ``ls *.rbs`` behaves the same on every platform
    (Windows shells do not expand wildcards at all).
    """
    if token is None:
        base = Path.cwd()
        return base, sorted(base.iterdir())

    target = resolve(token)
    if target.is_dir():
        return target, sorted(target.iterdir())

    # A glob: split the pattern off the last component.
    base = target.parent if str(target.parent) else Path.cwd()
    pattern = target.name
    if not base.is_dir():
        raise CommandError(f"no such directory: {base}")
    matched = sorted(
        child for child in base.iterdir() if fnmatch.fnmatch(child.name, pattern)
    )
    if not matched and not any(ch in pattern for ch in "*?["):
        raise CommandError(f"no such file or directory: {target}")
    return base, matched


def _name(path: Path) -> str:
    """Directories are marked with a trailing separator, as ``ls -F`` does."""
    return path.name + "/" if path.is_dir() else path.name


def cmd_ls(session, args: ArgReader) -> None:
    token = args.optional()
    args.done()
    base, entries = _entries(token)
    if not entries:
        print(f"  {base}: empty")
        return

    names = [_name(path) for path in entries]
    width = max(len(name) for name in names) + 2
    columns = max(1, (_terminal().columns - 2) // width)
    for start in range(0, len(names), columns):
        row = names[start : start + columns]
        print("  " + "".join(name.ljust(width) for name in row).rstrip())


def cmd_ll(session, args: ArgReader) -> None:
    """Long listing: size and modification time.

    Deliberately not Unix permission bits -- they carry no meaning on Windows,
    and this listing is the same on every platform.
    """
    token = args.optional()
    args.done()
    base, entries = _entries(token)
    if not entries:
        print(f"  {base}: empty")
        return
    for path in entries:
        try:
            info = path.stat()
        except OSError as error:
            print(f"  {'?':>12}  {'?':>16}  {_name(path)}  ({error.strerror})")
            continue
        when = datetime.fromtimestamp(info.st_mtime).strftime("%Y-%m-%d %H:%M")
        size = "-" if path.is_dir() else f"{info.st_size:,}"
        print(f"  {size:>12}  {when:>16}  {_name(path)}")


# ---------------------------------------------------------------------------
# Viewing
# ---------------------------------------------------------------------------


def cmd_type(session, args: ArgReader) -> None:
    """Show a text file, paging only when there is a terminal to page for."""
    path = resolve(args.token("a file to display"))
    args.done()
    if path.is_dir():
        raise CommandError(f"{path} is a directory")
    try:
        text = path.read_text(errors="replace")
    except FileNotFoundError:
        raise CommandError(f"no such file: {path}") from None
    except OSError as error:
        raise CommandError(f"could not read {path}: {error}") from None

    lines = text.splitlines()
    interactive = sys.stdin.isatty() and sys.stdout.isatty()
    if not interactive:
        # Under XEQ, --batch or pytest, never block for a keypress.
        for line in lines:
            print(line)
        return

    page = max(1, _terminal().lines - _PAGE_MARGIN)
    for start in range(0, len(lines), page):
        for line in lines[start : start + page]:
            print(line)
        if start + page < len(lines):
            try:
                if input("-- more -- (enter to continue, q to stop) ").strip().lower() == "q":
                    return
            except EOFError:
                return


def enable_ansi() -> bool:
    """Make ANSI escapes work on this console. True if they can be used.

    Linux and macOS terminals always can. Windows 10+ can once
    ``ENABLE_VIRTUAL_TERMINAL_PROCESSING`` is set on the console handle, which
    Python does not do for us; older consoles cannot, and fall back to ``cls``.
    Called once at shell start-up, and again lazily by CLS.
    """
    if os.name != "nt":
        return True
    try:
        import ctypes

        kernel32 = ctypes.windll.kernel32
        handle = kernel32.GetStdHandle(-11)  # STD_OUTPUT_HANDLE
        mode = ctypes.c_uint32()
        if not kernel32.GetConsoleMode(handle, ctypes.byref(mode)):
            return False
        return bool(
            kernel32.SetConsoleMode(handle, mode.value | 0x0004)
        )  # ENABLE_VIRTUAL_TERMINAL_PROCESSING
    except Exception:
        return False


def cmd_cls(session, args: ArgReader) -> None:
    """Clear the screen.

    One ANSI sequence everywhere it is supported; ``cls`` only for the older
    Windows consoles that cannot be switched into virtual-terminal mode.
    """
    args.done()
    if enable_ansi():
        print("\x1b[2J\x1b[H", end="")
    else:  # pragma: no cover - pre-Windows-10 console
        os.system("cls")  # noqa: S605 - fixed string, no user input


# ---------------------------------------------------------------------------
# Macros and session logging
# ---------------------------------------------------------------------------


def cmd_xeq(session, args: ArgReader) -> None:
    from ..repl import execute_file

    path = resolve(args.token("a command file"))
    args.done()
    execute_file(session, path)


def cmd_echo(session, args: ArgReader) -> None:
    token = args.optional()
    args.done()
    session.echo = True if token is None else token.lower() not in ("off", "no", "0")
    print(f"echo {'on' if session.echo else 'off'}")


def cmd_quiet(session, args: ArgReader) -> None:
    """``QUIET`` is the C's off-synonym for ECHO (system.c:203)."""
    args.done()
    session.echo = False
    print("echo off")


def cmd_script(session, args: ArgReader) -> None:
    """Record commands to a file, for replay with XEQ.

    RUMP calls this SCRIPT, with LOGFILE as a synonym (system.c:212). Both need
    at least four characters, which is how the original kept them clear of
    ``LOG`` -- the logarithmic yield axis at the RUMP level.
    """
    token = args.optional()
    args.done()
    if session.log_file is not None:
        session.log_file.close()
        session.log_file = None
    if token is None or token.lower() in ("off", "close"):
        print("script session closed")
        return
    path = resolve(token)
    try:
        session.log_file = open(path, "a", encoding="utf-8")
    except OSError as error:
        raise CommandError(f"could not open {path}: {error}") from None
    print(f"recording commands to {path}")


def cmd_help(session, args: ArgReader) -> None:
    args.done()
    print(TABLE.help_text())


TABLE = CommandTable("General System Commands")

_ENTRIES: list[tuple[str, int, object, str]] = [
    # Listing
    ("DIRECTORY", 3, cmd_ls, "list files, optionally matching a pattern"),
    ("LS", 2, cmd_ls, "list files, optionally matching a pattern"),
    ("SL", -2, cmd_ls, "list files (synonym)"),
    ("LL", -2, cmd_ll, "long listing with size and date"),
    # Navigation
    ("CD", 2, cmd_cd, "change directory (no argument: home)"),
    ("CHDIR", -5, cmd_cd, "change directory"),
    ("PUSHDIR", 5, cmd_pushdir, "change directory, remembering this one"),
    ("POPDIR", 4, cmd_popdir, "return to the last PUSHDIR directory"),
    ("WHERE", 3, cmd_pwd, "print the working directory"),
    ("PWD", -3, cmd_pwd, "print the working directory"),
    # Viewing
    ("TYPE", 2, cmd_type, "display a text file"),
    ("CAT", -3, cmd_type, "display a text file"),
    ("MORE", -4, cmd_type, "display a text file"),
    ("CLS", 3, cmd_cls, "clear the screen"),
    # Macros
    ("XEQ", 2, cmd_xeq, "execute a command file"),
    ("CALL", -4, cmd_xeq, "execute a command file"),
    ("EXECUTE", -3, cmd_xeq, "execute a command file"),
    ("ECHO", 4, cmd_echo, "echo commands as they run"),
    ("QUIET", -5, cmd_quiet, "stop echoing commands"),
    ("SCRIPT", 6, cmd_script, "record commands to a file for replay"),
    ("LOGFILE", -4, cmd_script, "record commands to a file"),
    ("RECORD", -3, cmd_script, "record commands to a file"),
]

for _name_, _minlen, _handler, _help in _ENTRIES:
    TABLE.add(_name_, _minlen, _handler, _help)
