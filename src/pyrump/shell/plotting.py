"""The session's persistent plot.

RUMP draws to one graphics device whose state -- region, yield range, yield
scaling -- survives between commands, so ``REGION 100 400`` followed by
``REPLOT`` redraws what is already there. The equivalent here is a single
long-lived matplotlib figure owned by the :class:`~pyrump.shell.session.Session`,
plus a list of traces that ``PLOT`` resets and ``OVERLAY`` appends to.

``PLOT``/``OVERLAY`` draw their traces here rather than through
:mod:`pyrump.plot.spectra`, because :class:`~pyrump.shell.session.PlotState`
owns the axis limits and scaling while ``plot_spectrum`` sets its own. The
whole-figure products -- ``COMPARE`` and ``DISPLAY`` -- do reuse
:func:`~pyrump.plot.spectra.plot_comparison` and
:func:`~pyrump.plot.spectra.plot_depth_profile`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..model.detector import yield_normalisation
from . import terminal_focus
from .dispatch import CommandError

#: Colour cycle for overlaid spectra. The first is RUMP's white-on-black data
#: trace; the rest keep overlays distinguishable in both light and dark themes.
_COLORS = ("0.20", "crimson", "steelblue", "darkgreen", "darkorange", "purple")


def require_matplotlib():
    """Import pyplot, or explain why the shell cannot plot.

    matplotlib ships as a core dependency, but this stays defensive for
    installs where it was stripped out or failed to build.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:  # pragma: no cover - depends on the install
        raise CommandError(
            "plotting needs matplotlib, which is not installed.\n"
            "Install it with:  pip install matplotlib"
        ) from None
    return plt


@dataclass(slots=True)
class Trace:
    """One curve on the plot.

    Holds the buffer rather than its counts, so that switching NORMALIZE/RAW or
    the x axis affects traces that are already on the plot -- and so a trace
    survives its buffer being renumbered.
    """

    buffer: object
    label: str
    index: int


def _values(buffer, state) -> np.ndarray:
    """Counts, normalised if the session is in NORMALIZE mode.

    Normalised units here are counts per msr per uC -- dividing out
    :func:`~pyrump.model.detector.yield_normalisation`, the same factor the
    simulation multiplies in. This makes spectra taken with different doses
    directly comparable; it is not claimed to be bit-identical to RUMP's own
    normalised yield.
    """
    counts = np.asarray(buffer.spectrum.counts, dtype=float)
    if not state.normalized:
        return counts
    factor = yield_normalisation(buffer.measurement)
    return counts / factor if factor else counts


def _x_axis(buffer, state, low: int, high: int):
    """X values and label for the drawn region."""
    if state.energy_axis:
        return buffer.spectrum.energies[low : high + 1], "Energy (keV)"
    return np.arange(low, high + 1, dtype=float), "Channel"


def _apply_scale(ax, state) -> None:
    if state.yscale == "log":
        ax.set_yscale("log")
    elif state.yscale == "sqrt":
        # RUMP's characteristic square-root yield axis: compresses a substrate
        # plateau without losing thin-film peaks the way log does.
        ax.set_yscale(
            "function",
            functions=(
                lambda v: np.sqrt(np.clip(v, 0, None)),
                lambda v: np.square(v),
            ),
        )
    else:
        ax.set_yscale("linear")


def _rebuilt_figure(session, n_axes: int, build, *, default_size: tuple[float, float]):
    """Reuse the session's live figure across layout changes.

    A window the user closed by hand is the one case a new ``Figure`` (and so
    a new OS window, at the backend's default position and size) is
    unavoidable -- ``fignum_exists`` catches that, and only then is
    ``default_size`` applied. Everything else, including switching panel
    counts (a single-panel PLOT after a two-panel COMPARE, or the reverse),
    clears and rebuilds axes on the SAME figure via ``clf()`` instead of
    closing and reopening it -- and leaves its size alone, so a window the
    user resized or moved by hand doesn't snap back on the next COMPARE.
    """
    plt = require_matplotlib()
    figure = session.figure
    if figure is not None and not plt.fignum_exists(figure.number):
        figure = None
    if figure is None:
        plt.ion()
        figure = plt.figure()
        figure.set_size_inches(*default_size)
    if len(figure.axes) != n_axes:
        figure.clf()
        build(figure)
    session.figure = figure
    return figure


def figure_for(session):
    """The session's single-panel figure and axes, created on first use."""
    def build(figure):
        figure.add_subplot(1, 1, 1)

    figure = _rebuilt_figure(session, 1, build, default_size=(9, 5.5))
    return figure, figure.axes[0]


def compare_figure_for(session, *, residuals: bool):
    """The session's COMPARE figure, reused across calls like :func:`figure_for`."""
    def build(figure):
        if residuals:
            figure.subplots(
                2, 1, sharex=True,
                gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
            )
        else:
            figure.add_subplot(1, 1, 1)

    return _rebuilt_figure(session, 2 if residuals else 1, build, default_size=(9, 6))


def draw(session) -> None:
    """Render every trace according to the current :class:`PlotState`."""
    if not session.traces:
        raise CommandError("nothing to plot yet")

    figure, ax = figure_for(session)
    ax.clear()
    state = session.plot

    n_channels = max(t.buffer.n_channels for t in session.traces)
    try:
        low, high = state.region(n_channels)
    except ValueError as error:
        raise CommandError(str(error)) from None

    label_axis = "Channel"
    for position, trace in enumerate(session.traces):
        counts = _values(trace.buffer, state)
        stop = min(high, counts.size - 1)
        if stop <= low:
            continue
        x, label_axis = _x_axis(trace.buffer, state, low, stop)
        ax.step(
            x,
            counts[low : stop + 1],
            where="mid",
            lw=1.0,
            color=_COLORS[position % len(_COLORS)],
            label=trace.label,
        )

    ax.set_xlabel(label_axis)
    ax.set_ylabel("Yield (counts/msr/uC)" if state.normalized else "Counts")
    _apply_scale(ax, state)

    if state.ylow is not None or state.yhigh is not None:
        ax.set_ylim(bottom=state.ylow, top=state.yhigh)
    elif state.yscale == "linear":
        ax.set_ylim(bottom=0)

    if state.labels and any(t.label for t in session.traces):
        ax.legend(frameon=False, fontsize="small")

    figure.tight_layout()
    show(figure)


def show(figure) -> None:
    """Push the figure to the screen without blocking the prompt.

    Raising a GUI window steals focus from wherever the user was typing --
    macOS in particular makes the new window "key". For an interactive
    backend (one with a real ``required_interactive_framework``, unlike Agg
    in tests) the terminal's focus is captured before drawing and restored
    right after, via :mod:`~pyrump.shell.terminal_focus`.
    """
    canvas = figure.canvas
    interactive = getattr(type(canvas), "required_interactive_framework", None) is not None
    token = terminal_focus.capture() if interactive else None
    canvas.draw_idle()
    try:
        canvas.flush_events()
    except (AttributeError, NotImplementedError):  # Agg in tests
        pass
    if interactive:
        terminal_focus.restore(token)


def buffer_label(session, buffer, index: int) -> str:
    """The name PLOT/OVERLAY would show for this buffer in a legend.

    Buffer 0 is always the simulation (Session.simulation()'s convention), so
    when STRUCTLABEL is on and a sample is described, its legend text becomes
    a compact rendering of the layer structure instead of the buffer's own
    name/identifier ("SIM").
    """
    if index == 0 and session.plot.structure_labels:
        from ..script.lcm import structure_label

        label = structure_label(
            session.script, normalize=session.plot.composition_fraction
        )
        if label:
            return label
    return buffer.name or buffer.identifier or f"buffer {index}"


def add_trace(session, index: int, buffer, *, clear: bool, replace: bool = False) -> None:
    """Add a buffer to the plot; ``clear`` makes it a fresh ``PLOT``.

    ``replace`` drops any existing trace for the same buffer index first, so a
    command that re-plots the same buffer on every call (SPLOT re-drawing the
    simulation as the sample changes) updates it in place instead of piling up
    a fresh copy each time. ``OVERLAY`` leaves this off -- stacking distinct
    buffers is the point of it.
    """
    if clear:
        session.traces = []
    elif replace:
        session.traces = [t for t in session.traces if t.index != index]
    label = buffer_label(session, buffer, index)
    session.traces.append(Trace(buffer=buffer, label=label, index=index))
