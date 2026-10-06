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


def test_buffer_stem_strips_a_windows_path_stamped_into_the_name():
    """A WRASCII macro's own FILENAME line can stamp a full Windows path
    straight into buffer.name (see cmd_filename) -- the stem must still be
    bare."""
    from types import SimpleNamespace

    buffer = SimpleNamespace(name=r"c:\RBS\data\2026\08\MA8410.RBS", identifier="")
    assert plotting.buffer_stem(buffer) == "MA8410"


def test_buffer_label_sanitizes_a_data_buffers_messy_name(session):
    from types import SimpleNamespace

    buffer = SimpleNamespace(name=r"c:\RBS\data\2026\08\MA8410.RBS", identifier="")
    assert plotting.buffer_label(session, buffer, 1) == "MA8410"


def test_buffer_label_keeps_a_synthetic_splot_label_untouched(session):
    """Buffer 0 is always a simulation (Session.selective_simulation's own
    "SIM(layer 2)" style caption), never data read off disk -- it must not
    run through buffer_stem's path/comment stripping, which would mangle the
    space in "layer 2" into a lost token."""
    from types import SimpleNamespace

    buffer = SimpleNamespace(name="SIM(layer 2)", identifier="")
    assert plotting.buffer_label(session, buffer, 0) == "SIM(layer 2)"


# -- pump: keeping the plot window alive during a long command ---------------
#
# Nothing pumps the GUI event loop while a command runs, so a fit that takes
# tens of seconds leaves the figure unrepainted long enough for the OS to
# declare it hung. ``pump`` is the escape hatch long commands call from their
# own iteration hooks; what matters is that it is cheap enough to call
# indiscriminately and safe when there is no window to pump.


def test_pump_does_nothing_without_a_figure(session):
    """Called from a fit whether or not anything is plotted, so the
    no-figure case has to be free rather than an error."""
    assert session.figure is None
    plotting.pump(session)  # no exception


def test_pump_skips_a_non_interactive_backend(session, monkeypatch):
    """Agg has no event loop; flush_events must not be reached for it."""
    figure, _ = plotting.figure_for(session)

    def boom():
        raise AssertionError("flush_events should not run for Agg")

    monkeypatch.setattr(figure.canvas, "flush_events", boom)
    plotting._last_pump = 0.0
    plotting.pump(session)


def test_pump_reaches_the_event_loop_on_an_interactive_backend(session, monkeypatch):
    figure, _ = plotting.figure_for(session)
    monkeypatch.setattr(
        type(figure.canvas), "required_interactive_framework", "fake", raising=False
    )
    calls = []
    monkeypatch.setattr(figure.canvas, "flush_events", lambda: calls.append("flush"))

    plotting._last_pump = 0.0
    plotting.pump(session)

    assert calls == ["flush"]


def test_pump_throttles_repeated_calls(session, monkeypatch):
    """A fit calls this once per model evaluation -- thousands of times over a
    long run -- so all but one call per PUMP_INTERVAL must be a cheap no-op."""
    figure, _ = plotting.figure_for(session)
    monkeypatch.setattr(
        type(figure.canvas), "required_interactive_framework", "fake", raising=False
    )
    calls = []
    monkeypatch.setattr(figure.canvas, "flush_events", lambda: calls.append("flush"))

    plotting._last_pump = 0.0
    for _ in range(50):
        plotting.pump(session)

    assert calls == ["flush"]


def test_pump_ignores_a_figure_the_user_closed(session, monkeypatch):
    """Pumping also delivers input, so the window can be closed *by* a pump --
    the next one must not touch the dead canvas."""
    figure, _ = plotting.figure_for(session)
    monkeypatch.setattr(
        type(figure.canvas), "required_interactive_framework", "fake", raising=False
    )
    monkeypatch.setattr(
        figure.canvas, "flush_events",
        lambda: (_ for _ in ()).throw(AssertionError("closed figure was pumped")),
    )
    plt.close(figure)

    plotting._last_pump = 0.0
    plotting.pump(session)  # no exception


# -- RUMP's two-scale frame: channel below, energy (MeV) above ----------------
#
# RbsAxdraw (tplot.c:442-483) always labels the bottom axis in channels and
# the top one in MeV, scaled from the reference buffer. ENERGY is pyRUMP's
# own switch to a single keV axis, which has no top axis.


def _buffer(kevch: float = 2.0, kev0: float = 100.0, first: float = 0.0):
    import numpy as np

    from pyrump.model.spectrum import Calibration, Spectrum
    from pyrump.shell.session import Buffer

    calibration = Calibration(kevch=kevch, kev0=kev0, first=first, npt=200)
    counts = np.linspace(10.0, 100.0, 200)
    return Buffer(spectrum=Spectrum(counts=counts, calibration=calibration), name="data")


def _top_limits(ax):
    """The top axis' limits as drawn -- a secondary axis sets them at draw."""
    (top,) = ax.child_axes
    ax.figure.canvas.draw()
    return top.get_xlim(), top.get_xlabel()


def test_default_plot_has_energy_in_mev_on_top(session):
    plotting.add_trace(session, 1, _buffer(first=5.0), clear=True)
    plotting.draw(session)

    ax = session.figure.axes[0]
    assert ax.get_xlabel() == "Channel"
    (low, high), label = _top_limits(ax)
    assert label == "Energy (MeV)"
    # RBSENERGY: 0.001 * (channel * kevch + kev0), with the stored first channel.
    left, right = ax.get_xlim()
    assert low == pytest.approx(0.001 * ((left + 5.0) * 2.0 + 100.0))
    assert high == pytest.approx(0.001 * ((right + 5.0) * 2.0 + 100.0))


def test_energy_on_leaves_one_kev_axis_and_off_restores_the_top_one(session):
    plotting.add_trace(session, 1, _buffer(), clear=True)

    session.plot.energy_axis = True
    plotting.draw(session)
    ax = session.figure.axes[0]
    assert ax.get_xlabel() == "Energy (keV)"
    assert not ax.child_axes

    session.plot.energy_axis = False
    plotting.draw(session)
    assert len(ax.child_axes) == 1


def test_top_axis_follows_the_first_trace(session):
    """Overlays with another calibration don't rescale it -- RUMP scales the
    frame from its reference buffer."""
    plotting.add_trace(session, 1, _buffer(kevch=2.0, kev0=100.0), clear=True)
    plotting.add_trace(session, 2, _buffer(kevch=5.0, kev0=0.0), clear=False)
    plotting.draw(session)

    ax = session.figure.axes[0]
    (low, _), _ = _top_limits(ax)
    assert low == pytest.approx(0.001 * (ax.get_xlim()[0] * 2.0 + 100.0))


def test_uncalibrated_buffer_plots_without_a_top_axis(session):
    plotting.add_trace(session, 1, _buffer(kevch=0.0), clear=True)
    plotting.draw(session)
    assert not session.figure.axes[0].child_axes


def test_compare_puts_the_energy_axis_on_its_upper_panel(session):
    plotting.add_trace(session, 1, _buffer(), clear=True, key=plotting.COMPARE_DATA)
    plotting.add_trace(session, 0, _buffer(), clear=False)
    plotting.draw(session)

    upper, residuals = session.figure.axes
    assert len(upper.child_axes) == 1
    assert not residuals.child_axes
    assert residuals.get_xlabel() == "Channel"

    session.plot.energy_axis = True
    plotting.draw(session)
    upper, _ = session.figure.axes
    assert not upper.child_axes
