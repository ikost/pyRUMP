"""Command tables with RUMP's abbreviation rules.

RUMP matches commands with ``LexCmdl`` (lexp2.c:639) against tables of
``{name, minlen, code}`` (the ``CMTYPE`` at rump.c:98). The rules are:

* ``minlen == 0`` -- the token must equal the name exactly (case-insensitively).
* ``minlen != 0`` -- the token must be at least ``abs(minlen)`` characters long
  and be a case-insensitive prefix of the name.
* ``minlen < 0`` additionally marks a *synonym*: it matches like any other entry
  but is left out of the ``?``/``HELP`` listing.

**The first entry that matches wins.** The C does a linear scan and returns
immediately; despite what its own comment block promises, it never reports
ambiguity. Table order is therefore significant, and we reproduce it rather than
raising on ambiguous prefixes -- old macros depend on which command a short
abbreviation lands on.

The help listing follows ``LexCmdlPrintEx`` (lexp2.c:688), which upper-cases the
required characters and lower-cases the rest, giving the ``REgion``/``PARMS``
convention used throughout the manual.

pyRUMP-only additions with no original counterpart mark ``extension=True``
(:class:`Command`) instead of any of the above -- they're invisible to
matching, listing and completion while the session is FAITHFUL, so a faithful
session's command surface is stock RUMP's, not just its physics.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover - typing only
    from .session import Session

#: A command handler. Receives the session and a reader over the remaining
#: tokens on the line, and mutates the session.
Handler = Callable[["Session", "ArgReader"], None]


class CommandError(Exception):
    """A command was malformed, unknown, or could not be carried out.

    The REPL prints these and continues; ``XEQ`` aborts the macro.
    """


@dataclass(frozen=True, slots=True)
class Command:
    """One entry in a :class:`CommandTable`."""

    name: str
    minlen: int
    handler: Handler
    help: str = ""
    extension: bool = False
    """A pyRUMP-only addition with no original-RUMP counterpart -- invisible
    (unmatched, unlisted, uncompleted) while the session is FAITHFUL, so a
    faithful session's command surface matches stock RUMP's exactly, not just
    its physics. Orthogonal to :attr:`hidden`: a hidden synonym always matches
    but never lists; an extension doesn't even match until FAITHFUL is off."""

    @property
    def hidden(self) -> bool:
        """Synonyms (negative ``minlen``) are matched but not listed."""
        return self.minlen < 0

    @property
    def min_chars(self) -> int:
        """Characters that must be typed. Zero means "the whole name"."""
        return abs(self.minlen)

    def matches(self, token: str) -> bool:
        if self.minlen == 0:
            return token.casefold() == self.name.casefold()
        if len(token) < self.min_chars:
            return False
        return self.name.casefold().startswith(token.casefold())

    @property
    def display(self) -> str:
        """``REgion``: required characters upper-cased, the rest lower."""
        required = self.min_chars or len(self.name)
        return self.name[:required].upper() + self.name[required:].lower()


@dataclass(slots=True)
class CommandTable:
    """An ordered table of commands, matched first-hit-wins."""

    title: str
    commands: list[Command] = field(default_factory=list)

    def add(
        self, name: str, minlen: int, handler: Handler, help: str = "",
        extension: bool = False,
    ) -> None:
        self.commands.append(Command(name, minlen, handler, help, extension))

    def extend(self, entries: Sequence[Command]) -> None:
        self.commands.extend(entries)

    def match(self, token: str, *, faithful: bool = True) -> Command | None:
        """The first command matching ``token``, or None.

        ``faithful=True`` (the default) skips :attr:`Command.extension`
        entries entirely, so an unrecognised pyRUMP-only command behaves
        exactly like stock RUMP's "unrecognized command" rather than working.
        """
        token = token.strip()
        if not token:
            return None
        for command in self.commands:
            if command.extension and faithful:
                continue
            if command.matches(token):
                return command
        return None

    def visible(self, *, faithful: bool = True) -> Iterator[Command]:
        """Commands that appear in the help listing, in table order."""
        return (
            c for c in self.commands
            if not c.hidden and not (c.extension and faithful)
        )

    def completions(self, prefix: str, *, faithful: bool = True) -> list[str]:
        """Command names for tab-completion (visible entries only)."""
        folded = prefix.casefold()
        return [
            c.name.lower()
            for c in self.visible(faithful=faithful)
            if c.name.casefold().startswith(folded)
        ]

    def listing(self, width: int = 80, *, faithful: bool = True) -> str:
        """Render the help listing, RUMP-style."""
        entries = list(self.visible(faithful=faithful))
        if not entries:
            return self.title
        column = max(len(c.display) for c in entries) + 2
        per_line = max(1, width // column)
        lines = [self.title]
        for start in range(0, len(entries), per_line):
            row = entries[start : start + per_line]
            lines.append("  " + "".join(c.display.ljust(column) for c in row).rstrip())
        return "\n".join(lines)

    def help_text(self, *, faithful: bool = True) -> str:
        """One command per line with its description, for ``HELP``."""
        entries = list(self.visible(faithful=faithful))
        if not entries:
            return self.title
        column = max(len(c.display) for c in entries)
        lines = [self.title]
        lines.extend(
            f"  {c.display.ljust(column)}  {c.help}".rstrip() for c in entries
        )
        return "\n".join(lines)


def tokenize(line: str) -> list[str]:
    """Split a command line into tokens.

    Quoted strings are kept whole. ``/`` -- RUMP's "no more arguments"
    terminator, used by ``COMPOSITION`` and friends (script/lcm.py) -- is always
    its own token even when written flush against the previous one.

    An unterminated quote is not an error: RUMP's own lexer just takes
    everything to end of line as the token (lexp.c:506, "If quotes did not
    finish, then just continue on"). ``WRASCII``-written macros rely on this --
    their ``Identifier '...`` line never closes the quote.
    """
    tokens: list[str] = []
    index, length = 0, len(line)
    while index < length:
        while index < length and line[index].isspace():
            index += 1
        if index >= length:
            break
        quote = line[index]
        if quote in ("'", '"'):
            index += 1
            start = index
            while index < length and line[index] != quote:
                index += 1
            token = line[start:index]
            index += 1  # past the closing quote, or harmlessly past end of line
        else:
            start = index
            while index < length and not line[index].isspace():
                index += 1
            token = line[start:index]
        # "3/" -> ["3", "/"], so the terminator need not be spaced off.
        if token != "/" and token.endswith("/"):
            tokens.append(token[:-1])
            tokens.append("/")
        else:
            tokens.append(token)
    return [t for t in tokens if t]


def strip_comment(line: str) -> str:
    """Remove RUMP comments: ``/* ... `` to end of line, or a leading ``#``/``!``."""
    line = line.split("/*", 1)[0].strip()
    if line.startswith("#") or line.startswith("!"):
        return ""
    return line


@dataclass(slots=True)
class ArgReader:
    """Cursor over the tokens following a command name."""

    tokens: list[str]
    command: str = ""
    index: int = 0

    def __bool__(self) -> bool:
        return self.index < len(self.tokens)

    @property
    def remaining(self) -> list[str]:
        return self.tokens[self.index :]

    def _fail(self, what: str) -> CommandError:
        where = f"{self.command}: " if self.command else ""
        return CommandError(f"{where}expected {what}")

    def peek(self) -> str | None:
        return self.tokens[self.index] if self else None

    def token(self, what: str = "an argument") -> str:
        if not self:
            raise self._fail(what)
        value = self.tokens[self.index]
        self.index += 1
        return value

    def optional(self) -> str | None:
        return self.token() if self else None

    def number(self, what: str = "a number") -> float:
        token = self.token(what)
        try:
            return float(token)
        except ValueError:
            raise CommandError(
                f"{self.command}: {token!r} is not a number"
            ) from None

    def integer(self, what: str = "an integer") -> int:
        # RUMP reads everything as a float and truncates, so "3.0" is a valid
        # integer argument (see the int(float(...)) calls in script/lcm.py).
        return int(self.number(what))

    def optional_number(self) -> float | None:
        return self.number() if self else None

    def rest(self) -> str:
        """The remaining tokens rejoined, for free-text arguments."""
        value = " ".join(self.remaining)
        self.index = len(self.tokens)
        return value

    def element_pairs(self) -> dict[str, float]:
        """Parse ``El value El value ... /``.

        Shares its grammar with :func:`pyrump.script.lcm._element_pairs`; the
        trailing ``/`` terminates the list.
        """
        from ..script.lcm import _element_pairs

        try:
            pairs = _element_pairs(self.remaining)
        except ValueError as error:
            raise CommandError(f"{self.command}: {error}") from None
        self.index = len(self.tokens)
        return pairs

    def done(self) -> None:
        """Reject trailing junk, so typos are not silently ignored."""
        if self:
            extra = " ".join(self.remaining)
            raise CommandError(f"{self.command}: unexpected extra argument {extra!r}")
