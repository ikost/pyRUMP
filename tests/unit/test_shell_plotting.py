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


def test_figure_for_rebuilds_a_leftover_two_panel_figure(session):
    """The pre-existing COMPARE-leftover case must still work."""
    figure, axes = plt.subplots(2, 1)
    session.figure = figure

    new_figure, new_axes = plotting.figure_for(session)

    assert new_figure is not figure
    assert len(new_figure.axes) == 1
    assert new_axes is new_figure.axes[0]
