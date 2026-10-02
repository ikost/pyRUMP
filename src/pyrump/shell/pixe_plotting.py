"""The PIXE window: a second matplotlib figure next to the RBS one.

Kept apart from :mod:`pyrump.shell.plotting` on purpose. The RBS window's
layout, reuse and focus handling stay exactly as they are whether PIXE is on
or not; this module only borrows its yield-axis and show helpers. The PIXE
window shows the ACTIVE buffer's PIXE spectrum against energy, and follows
the RBS window: every RBS redraw redraws it too while PIXE is enabled.
"""

from __future__ import annotations

import numpy as np

from .. import __version__
from .dispatch import CommandError
from .plotting import _apply_limits, _apply_scale, buffer_stem, require_matplotlib, show

#: Default window size in inches, a little smaller than the RBS window's.
DEFAULT_SIZE = (8.0, 5.0)


def figure_for(session):
    """The PIXE figure and its axes, created on first use.

    A window the user closed by hand is reopened at the default size; a live
    one keeps whatever size and position the user gave it.
    """
    plt = require_matplotlib()
    state = session.pixe
    figure = state.figure
    if figure is not None and not plt.fignum_exists(figure.number):
        figure = None
    if figure is None:
        plt.ion()
        figure = plt.figure()
        figure.set_size_inches(*DEFAULT_SIZE)
        figure.canvas.manager.set_window_title(f"pyRUMP {__version__} - PIXE")
    if len(figure.axes) != 1:
        figure.clf()
        figure.add_subplot(1, 1, 1)
    state.figure = figure
    return figure, figure.axes[0]


def close(session) -> None:
    """Close the PIXE window, if it is open."""
    figure = session.pixe.figure
    session.pixe.figure = None
    if figure is None:
        return
    plt = require_matplotlib()
    if plt.fignum_exists(figure.number):
        plt.close(figure)


def pixe_buffer(session):
    """The ACTIVE buffer, if it holds a PIXE spectrum."""
    buffer = session.buffers.active_buffer
    return buffer if buffer is not None and buffer.pixe is not None else None


def draw(session, *, required: bool = True) -> bool:
    """Draw the ACTIVE buffer's PIXE spectrum. Returns whether anything was
    drawn.

    ``required`` makes a missing spectrum an error (the PIXE prompt's own
    PLOT); otherwise -- an RBS redraw refreshing this window -- there is
    simply nothing to do.
    """
    buffer = pixe_buffer(session)
    if buffer is None:
        if required:
            raise CommandError("no PIXE spectrum in the active buffer: PIXE GET <file>")
        return False

    state = session.pixe
    data = buffer.pixe
    counts = np.asarray(data.spectrum.counts, dtype=float)
    energies = data.spectrum.energies
    keep = np.ones(counts.size, dtype=bool)
    if state.emin is not None:
        keep &= energies >= state.emin
    if state.emax is not None:
        keep &= energies <= state.emax
    if keep.sum() < 2:
        raise CommandError(
            f"no PIXE channels between {state.emin} and {state.emax} keV"
        )

    figure, ax = figure_for(session)
    ax.clear()
    ax.step(energies[keep], counts[keep], where="mid", lw=1.0, color="0.20",
            label=data.path.name if data.path else (data.identifier or "PIXE"))
    ax.set_xlim(energies[keep][0], energies[keep][-1])
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel("Counts")
    ax.set_title(buffer_stem(buffer) or "PIXE", fontsize="medium")
    _apply_scale(ax, state.plot)
    _apply_limits(ax, state.plot)
    if state.plot.yscale == "log" and state.plot.ylow is None:
        # Empty channels would drag a log axis down to its clip floor.
        ax.set_ylim(bottom=0.5)
    if state.plot.labels:
        ax.legend(frameon=False, fontsize="small")
    figure.tight_layout()
    show(figure)
    return True


def refresh(session) -> None:
    """Follow an RBS redraw: redraw the PIXE window while PIXE is enabled."""
    if session.pixe.enabled:
        draw(session, required=False)
