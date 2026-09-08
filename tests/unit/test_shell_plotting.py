"""The session's persistent plot figure -- creation, reuse, and staleness.

``figure_for`` is the one choke point every drawing command goes through
(PLOT/OVERLAY/REPLOT/SPLOT/AXIS), so it is tested directly rather than through
the REPL -- it needs no buffers, sample data, or atomic tables, only a
``Session`` to hold ``figure``.
"""

from __future__ import annotations

import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

from pyrump.shell import plotting  # noqa: E402
from pyrump.shell.session import Session  # noqa: E402


@pytest.fixture
def session() -> Session:
    """A session with no data tables -- figure_for never touches them."""
    return Session(table=None, registry=None, densities=None, data=None)


def test_figure_for_creates_a_figure_on_first_use(session):
    figure, axes = plotting.figure_for(session)
    assert session.figure is figure
    assert len(figure.axes) == 1
    assert axes is figure.axes[0]


def test_figure_for_reuses_a_live_figure(session):
    figure1, _ = plotting.figure_for(session)
    figure2, _ = plotting.figure_for(session)
    assert figure2 is figure1


def test_figure_for_rebuilds_after_the_window_is_closed(session):
    """Closing the window (the OS close button) deregisters it from pyplot
    without touching ``figure.axes`` or ``session.figure`` -- ``plt.close``
    reproduces exactly that."""
    figure1, _ = plotting.figure_for(session)
    number1 = figure1.number
    plt.close(figure1)
    assert not plt.fignum_exists(number1)

    figure2, axes2 = plotting.figure_for(session)

    assert figure2 is not figure1
    assert plt.fignum_exists(figure2.number)
    assert len(figure2.axes) == 1
    assert axes2 is figure2.axes[0]


def test_figure_for_rebuilds_a_leftover_two_panel_figure_in_place(session):
    """A COMPARE-leftover figure is cleared and rebuilt on the SAME window,
    not closed and reopened -- otherwise switching from a two-panel COMPARE
    back to a single-panel PLOT would jump the window to its default
    position."""
    figure, axes = plt.subplots(2, 1)
    session.figure = figure
    number = figure.number

    new_figure, new_axes = plotting.figure_for(session)

    assert new_figure is figure
    assert new_figure.number == number
    assert len(new_figure.axes) == 1
    assert new_axes is new_figure.axes[0]


def test_compare_figure_for_reuses_a_live_figure_across_layout_changes(session):
    """PLOT (1 panel) then COMPARE (2 panels) then PLOT again must all draw
    into the same window -- only a hand-closed window gets a new one."""
    figure1, _ = plotting.figure_for(session)
    number = figure1.number

    figure2 = plotting.compare_figure_for(session, residuals=True)
    assert figure2 is figure1
    assert figure2.number == number
    assert len(figure2.axes) == 2

    figure3, axes3 = plotting.figure_for(session)
    assert figure3 is figure1
    assert figure3.number == number
    assert len(figure3.axes) == 1
    assert axes3 is figure3.axes[0]


def test_show_skips_the_focus_handoff_on_a_non_interactive_backend(session, monkeypatch):
    """Agg (what the whole test suite runs under) has no real window, so
    there is nothing to steal focus -- capture/restore must not even run."""
    figure, _ = plotting.figure_for(session)

    def boom(*args, **kwargs):
        raise AssertionError("terminal_focus should not run for Agg")

    monkeypatch.setattr(plotting.terminal_focus, "capture", boom)
    monkeypatch.setattr(plotting.terminal_focus, "restore", boom)
    plotting.show(figure)  # no exception means neither was called


def test_show_hands_focus_back_on_an_interactive_backend(session, monkeypatch):
    """A backend that reports a real ``required_interactive_framework`` (every
    GUI backend, unlike Agg) must capture focus before drawing and restore it
    after."""
    figure, _ = plotting.figure_for(session)
    monkeypatch.setattr(
        type(figure.canvas), "required_interactive_framework", "fake", raising=False
    )

    calls = []
    monkeypatch.setattr(plotting.terminal_focus, "capture", lambda: calls.append("capture") or "token")
    monkeypatch.setattr(plotting.terminal_focus, "restore", lambda token: calls.append(("restore", token)))

    plotting.show(figure)

    assert calls == ["capture", ("restore", "token")]
