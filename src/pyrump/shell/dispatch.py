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
"""

from __future__ import annotations

import inspect
import re
from collections.abc import Callable, Collection, Iterator, Sequence
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover - typing only
    from .session import Session

#: A command handler. Receives the session and a reader over the remaining
#: tokens on the line, and mutates the session.
Handler = Callable[["Session", "ArgReader"], None]

#: Many handler docstrings open with their own usage line, e.g.
#: ``"``THICKNESS lo hi element`` -- INTEGRAL plus a thickness conversion."``
#: -- a convention already used throughout commands/*.py. Pull it out for
#: ``HELP <name>`` rather than inventing a second, parallel place to write it.
_USAGE_RE = re.compile(r"^``([^`]+)``")


def _usage(handler: Handler) -> str | None:
    doc = inspect.getdoc(handler)
    if not doc:
        return None
    match = _USAGE_RE.match(doc)
    return match.group(1) if match else None


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
    synonyms: tuple[str, ...] = ()
    """Other names that also reach this command, noted in ``display`` (e.g.
    ``COMPARE / CMP``) via :meth:`CommandTable.note_synonym`. Each synonym is
    normally its own hidden table entry -- kept out of the listing to avoid a
    redundant "synonym for X" row -- so without this it would be invisible to
    anyone reading a bare HELP/``?`` listing rather than already knowing to
    ask ``HELP <synonym>`` by name."""

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
        """``REgion``: required characters upper-cased, the rest lower.

        With ``synonyms`` set, appends them verbatim, e.g. ``COMPARE / CMP``.
        """
        required = self.min_chars or len(self.name)
        rendered = self.name[:required].upper() + self.name[required:].lower()
        if self.synonyms:
            rendered += " / " + " / ".join(self.synonyms)
        return rendered


@dataclass(slots=True)
class CommandTable:
    """An ordered table of commands, matched first-hit-wins."""

    title: str
    commands: list[Command] = field(default_factory=list)

    def add(
        self, name: str, minlen: int, handler: Handler, help: str = ""
    ) -> None:
        self.commands.append(Command(name, minlen, handler, help))

    def extend(self, entries: Sequence[Command]) -> None:
        self.commands.extend(entries)

    def note_synonym(self, name: str, *synonyms: str) -> None:
        """Show ``synonyms`` alongside ``name`` in the listing, e.g.
        ``COMPARE / CMP``. Each synonym is expected to already be its own
        (hidden) entry that actually does the matching -- this only changes
        what ``name`` displays as, so a bare HELP/``?`` listing doesn't leave
        a hidden synonym looking undiscoverable.
        """
        for i, command in enumerate(self.commands):
            if command.name == name:
                self.commands[i] = replace(
                    command, synonyms=command.synonyms + synonyms
                )
                return
        raise KeyError(f"no command named {name!r} to attach a synonym note to")

    def match(self, token: str) -> Command | None:
        """The first command matching ``token``, or None."""
        token = token.strip()
        if not token:
            return None
        for command in self.commands:
            if command.matches(token):
                return command
        return None

    def visible(self, exclude: Collection[str] = ()) -> Iterator[Command]:
        """Commands that appear in the help listing, in table order.

        ``exclude`` leaves a command out of *this* table's own listing without
        touching matching -- e.g. a command cross-referenced into another
        table's grouped display instead (see :meth:`grouped_help_text`).
        """
        return (c for c in self.commands if not c.hidden and c.name not in exclude)

    def completions(self, prefix: str) -> list[str]:
        """Command names for tab-completion (visible entries only)."""
        folded = prefix.casefold()
        return [
            c.name.lower()
            for c in self.visible()
            if c.name.casefold().startswith(folded)
        ]

    def listing(self, width: int = 80) -> str:
        """Render the help listing, RUMP-style."""
        entries = list(self.visible())
        if not entries:
            return self.title
        column = max(len(c.display) for c in entries) + 2
        per_line = max(1, width // column)
        lines = [self.title]
        for start in range(0, len(entries), per_line):
            row = entries[start : start + per_line]
            lines.append("  " + "".join(c.display.ljust(column) for c in row).rstrip())
        return "\n".join(lines)

    def help_text(self, exclude: Collection[str] = ()) -> str:
        """One command per line with its description, for ``HELP``."""
        entries = list(self.visible(exclude))
        if not entries:
            return self.title
        column = max(len(c.display) for c in entries)
        lines = [self.title]
        lines.extend(
            f"  {c.display.ljust(column)}  {c.help}".rstrip() for c in entries
        )
        return "\n".join(lines)

    def uncovered(self, groups: Sequence[tuple[str, Sequence[str | Command]]]) -> list[str]:
        """Visible command names named in none of ``groups``.

        Feed this back in as a trailing ``("Other", ...)`` group so a command
        added later without updating the grouping is never silently dropped
        from the ``?``/``HELP`` listing. A ``Command`` entry (cross-referenced
        from another table, see :meth:`grouped_help_text`) is never one of
        *this* table's own visible commands, so it can't cover anything here
        and is ignored.
        """
        named = {item for _, items in groups for item in items if isinstance(item, str)}
        return [c.name for c in self.visible() if c.name not in named]

    def grouped_help_text(
        self,
        groups: Sequence[tuple[str, Sequence[str | Command]]],
        *,
        show_title: bool = True,
    ) -> str:
        """Like :meth:`help_text`, but under named sections in a given order.

        ``groups`` is a list of ``(heading, [command names...])`` pairs, most
        important section first -- see :meth:`uncovered` for a catch-all
        trailing group. This only changes how the listing is *displayed*;
        matching still runs over ``commands`` in its original (significant)
        order, so abbreviation resolution is untouched.

        An item may also be an actual :class:`Command`, e.g. one looked up
        from a different table with ``match()`` -- lets a caller cross-
        reference one command into a more useful section elsewhere (pair
        this with ``exclude`` on the owning table's own listing, so it isn't
        shown twice).

        ``show_title`` can be turned off to omit the leading table title, for
        a caller making several of these calls in a row to interleave another
        table's own section in between (only the first call needs a title).
        """
        entries = list(self.visible())
        if not entries:
            return self.title if show_title else ""
        by_name = {c.name: c for c in entries}
        column = max(len(c.display) for c in entries)
        lines = [self.title] if show_title else []
        for heading, items in groups:
            rows = [
                item if isinstance(item, Command) else by_name[item]
                for item in items
                if isinstance(item, Command) or item in by_name
            ]
            if not rows:
                continue
            lines.append("")
            lines.append(f"  {heading}")
            lines.extend(
                f"    {c.display.ljust(column)}  {c.help}".rstrip() for c in rows
            )
        return "\n".join(lines)

    def describe(self, token: str) -> str | None:
        """The full entry for a single command, for ``HELP <name>``.

        None if this table has nothing matching ``token``, so callers can fall
        through to the next table the way :func:`execute_line` falls through
        the mode stack. Adds a ``usage:`` line when the handler's own
        docstring states one (see :func:`_usage`) -- the short table
        description alone doesn't say what arguments a command takes.
        """
        command = self.match(token)
        if command is None:
            return None
        lines = [self.title, f"  {command.display}  {command.help}".rstrip()]
        usage = _usage(command.handler)
        if usage:
            lines.append(f"  usage: {usage}")
        return "\n".join(lines)


def _is_number(token: str) -> bool:
    try:
        float(token)
    except ValueError:
        return False
    return True


def tokenize(line: str) -> list[str]:
    """Split a command line into tokens.

    Quoted strings are kept whole. ``/`` -- RUMP's "no more arguments"
    terminator, used by ``COMPOSITION`` and friends (script/lcm.py) -- is always
    its own token even when written flush against a number, e.g. ``3/``. A
    trailing ``/`` on anything else (a directory path from tab completion,
    say) is left alone, so ``CD Documents/`` stays one token.

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
        # "3/" -> ["3", "/"], so the terminator need not be spaced off a
        # number. Anything else ending in "/" (a path) stays one token.
        if token != "/" and token.endswith("/") and _is_number(token[:-1]):
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
