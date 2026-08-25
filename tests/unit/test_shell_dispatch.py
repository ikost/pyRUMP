"""Command matching, tokenizing and argument reading.

The matching rules come from ``LexCmdl`` (lexp2.c:639) and the table format
from ``CMTYPE`` (rump.c:98); the cases below are transcribed from RUMP's own
main-level table so that a regression here shows up as a behaviour change
against the original.
"""

from __future__ import annotations

import pytest

from pyrump.shell.dispatch import (
    ArgReader,
    CommandError,
    CommandTable,
    strip_comment,
    tokenize,
)


def _noop(session, args):
    return None


@pytest.fixture
def table() -> CommandTable:
    # rump.c:147, verbatim in name and minlen.
    entries = [
        ("?", 1),
        ("PARAMETERS", 4),
        ("PARMS", -3),
        ("STATUS", -4),
        ("RESET", 5),
        ("RESTART", 7),
        ("CONFIGURE", 4),
        ("VERSION", 3),
        ("RUMP", -4),
        ("QUIT", 1),
        ("BYE", -2),
    ]
    built = CommandTable("Main Level Commands")
    for name, minlen in entries:
        built.add(name, minlen, _noop, f"{name} help")
    return built


@pytest.mark.parametrize(
    "token, expected",
    [
        ("q", "QUIT"),
        ("quit", "QUIT"),
        ("QUIT", "QUIT"),
        ("by", "BYE"),
        ("bye", "BYE"),
        ("par", "PARMS"),      # PARAMETERS needs 4, so PARMS wins at 3
        ("parm", "PARMS"),     # "parm" is not a prefix of PARAMETERS
        ("para", "PARAMETERS"),
        ("param", "PARAMETERS"),
        ("reset", "RESET"),
        ("restart", "RESTART"),
        ("ver", "VERSION"),
        ("stat", "STATUS"),
    ],
)
def test_abbreviations_resolve_as_in_rump(table, token, expected):
    assert table.match(token).name == expected


@pytest.mark.parametrize("token", ["b", "re", "rese", "ve", "st", "conf1", "nonsense"])
def test_too_short_or_unknown_tokens_do_not_match(table, token):
    assert table.match(token) is None


def test_matching_is_case_insensitive(table):
    assert table.match("ReSeT").name == "RESET"


def test_first_entry_in_table_order_wins():
    """LexCmdl returns the first match and never reports ambiguity."""
    built = CommandTable("t")
    built.add("SIMULATE", 2, _noop)
    built.add("SIMPLE", 2, _noop)
    assert built.match("si").name == "SIMULATE"
    assert built.match("simp").name == "SIMPLE"


def test_zero_minlen_requires_an_exact_match():
    built = CommandTable("t")
    built.add("GO", 0, _noop)
    assert built.match("go").name == "GO"
    assert built.match("g") is None


def test_extension_is_invisible_by_default():
    """A pyRUMP-only addition matches only once the caller opts out of
    ``faithful`` -- the default keeps a faithful session's command surface
    identical to stock RUMP's."""
    built = CommandTable("t")
    built.add("OFFSET", 3, _noop, extension=True)
    assert built.match("offset") is None
    assert built.match("offset", faithful=True) is None
    assert built.match("offset", faithful=False).name == "OFFSET"
    assert "OFFSET" not in [c.name for c in built.visible()]
    assert "OFFSET" in [c.name for c in built.visible(faithful=False)]
    assert built.completions("off") == []
    assert built.completions("off", faithful=False) == ["offset"]


def test_synonyms_are_matched_but_not_listed(table):
    assert table.match("parms").name == "PARMS"
    listed = [command.name for command in table.visible()]
    assert "PARMS" not in listed
    assert "BYE" not in listed
    assert "PARAMETERS" in listed


def test_listing_upper_cases_the_required_characters(table):
    assert "PARAmeters" in table.listing()
    assert "Quit" in table.listing()


def test_completions_cover_visible_commands_only(table):
    assert set(table.completions("r")) == {"reset", "restart"}
    assert "parms" not in table.completions("p")


# -- tokenizing ------------------------------------------------------------


def test_tokenize_splits_the_slash_terminator():
    assert tokenize("comp In 2 O 3 /") == ["comp", "In", "2", "O", "3", "/"]


def test_tokenize_splits_a_slash_flush_against_a_value():
    assert tokenize("comp In 2 O 3/") == ["comp", "In", "2", "O", "3", "/"]


def test_tokenize_keeps_quoted_strings_whole():
    assert tokenize("ident 'a b c'") == ["ident", "a b c"]


def test_tokenize_tolerates_an_unbalanced_quote():
    """RUMP's own lexer just takes the rest of the line (lexp.c:506) --

    WRASCII-written macros rely on it: their ``Identifier '...`` line never
    closes the quote.
    """
    assert tokenize("ident 'unterminated") == ["ident", "unterminated"]


@pytest.mark.parametrize(
    "line, expected",
    [
        ("plot 1 /* draw it", "plot 1"),
        ("  # a comment", ""),
        ("! another", ""),
        ("  ", ""),
        ("region 1 2", "region 1 2"),
    ],
)
def test_strip_comment(line, expected):
    assert strip_comment(line) == expected


# -- argument reading ------------------------------------------------------


def test_arg_reader_numbers_and_integers():
    args = ArgReader(["1.5", "3.0"], command="test")
    assert args.number() == 1.5
    # RUMP reads everything as a float and truncates.
    assert args.integer() == 3


def test_arg_reader_reports_a_missing_argument():
    args = ArgReader([], command="region")
    with pytest.raises(CommandError, match="region: expected the first channel"):
        args.integer("the first channel")


def test_arg_reader_rejects_a_non_number():
    args = ArgReader(["wide"], command="region")
    with pytest.raises(CommandError, match="not a number"):
        args.number()


def test_arg_reader_done_rejects_trailing_junk():
    args = ArgReader(["1", "2", "oops"], command="region")
    args.integer()
    args.integer()
    with pytest.raises(CommandError, match="unexpected extra argument"):
        args.done()


def test_arg_reader_element_pairs_shares_the_lcm_grammar():
    args = ArgReader(["In", "2", "O", "3", "/"], command="composition")
    assert args.element_pairs() == {"In": 2.0, "O": 3.0}
    assert not args


def test_arg_reader_element_pairs_rejects_a_dangling_symbol():
    args = ArgReader(["In", "2", "O"], command="composition")
    with pytest.raises(CommandError, match="no amount"):
        args.element_pairs()
