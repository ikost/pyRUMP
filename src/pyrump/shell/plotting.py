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
from .dispatch import CommandError

#: Colour cycle for overlaid spectra. The first is RUMP's white-on-black data
#: trace; the rest keep overlays distinguishable in both light and dark themes.
_COLORS = ("0.20", "crimson", "steelblue", "darkgreen", "darkorange", "purple")


def require_matplotlib():
    """Import pyplot, or explain why the shell cannot plot.

    matplotlib is an optional dependency -- the library and the batch CLI work
    without it -- but the interactive shell is not much use without a plot.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:  # pragma: no cover - depends on the install
        raise CommandError(
            "plotting needs matplotlib, which is not installed.\n"
            "Install it with:  pip install 'pyrump[plot]'"
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


def figure_for(session):
    """The session's single-panel figure and axes, created on first use.

    COMPARE leaves a two-panel figure behind, so a figure with the wrong number
    of axes is discarded rather than drawn into -- otherwise a PLOT after a
    COMPARE would render into the top panel above an orphaned residual strip.
    Closing the window by hand (the OS close button) is discarded the same
    way: it deregisters the figure from pyplot without clearing its axes or
    touching ``session.figure``, so ``fignum_exists`` is what actually catches
    it -- otherwise PLOT/REPLOT/OVERLAY/SPL would keep drawing into a window
    that no longer exists.
    """
    plt = require_matplotlib()
    stale = session.figure is not None and (
        len(session.figure.axes) != 1 or not plt.fignum_exists(session.figure.number)
    )
    if stale:
        plt.close(session.figure)
        session.figure = None
    if session.figure is None:
        plt.ion()
        session.figure = plt.figure(figsize=(9, 5.5))
        session.figure.add_subplot(1, 1, 1)
    return session.figure, session.figure.axes[0]


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
    """Push the figure to the screen without blocking the prompt."""
    canvas = figure.canvas
    canvas.draw_idle()
    try:
        canvas.flush_events()
    except (AttributeError, NotImplementedError):  # Agg in tests
        pass


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
    label = buffer.name or buffer.identifier or f"buffer {index}"
    session.traces.append(Trace(buffer=buffer, label=label, index=index))
