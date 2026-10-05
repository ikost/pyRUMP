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
    if state.markers != "off":
        marks = _line_marks(ax, sample_elements(session), major_only=state.markers == "on")
        if marks:
            if state.plot.yhigh is None:
                _make_headroom(ax, counts[keep], state.plot.yscale)
            _draw_marks(ax, marks)
    show(figure)
    return True


def sample_elements(session) -> list[tuple[int, str]]:
    """``(Z, symbol)`` of every element in the SIM sample, in the order SIM
    lists them, leaving out RBS absorber layers (foils in front of the RBS
    detector, which the beam never reaches)."""
    script = session.script
    found: list[tuple[int, str]] = []
    for layer in script.layers[script.absorber_layers:]:
        for name in (*layer.composition, *layer.species):
            try:
                element = session.table.by_z(session.table.parse_ref(name).z)
            except (KeyError, ValueError):
                continue
            if (element.z, element.symbol) not in found:
                found.append((element.z, element.symbol))
    return found


#: Layout of the line markers, in axes fractions from the top edge: the tick
#: at the line's true energy, then a leader to where its label sits.
_TICK_END, _LABEL_TOP = 0.03, 0.07
#: The data are kept below this fraction of the axes height, so the labels
#: hanging from the top edge stay clear of the peaks.
_DATA_CEILING = 0.70
_LABEL_FONT_PT = 7.0
#: Labels this close (in label widths) are one label: the lines can't be
#: told apart on screen anyway ("Mn Lα/Lβ1").
_MERGE_WIDTHS = 0.4


def _line_marks(ax, elements, *, major_only: bool) -> list[tuple[float, str]]:
    """``(energy, label)`` for each line group of ``elements`` inside the
    x range, coincident ones merged."""
    from ..pixe.atomic import atomic_data

    data = atomic_data()
    low, high = ax.get_xlim()
    marks = sorted(
        (group.energy_keV, symbol, group.label)
        for z, symbol in elements
        for group in data.line_groups(z, major_only=major_only)
        if low <= group.energy_keV <= high
    )
    return _merged(marks, _MERGE_WIDTHS * _label_width(ax))


def _label_width(ax) -> float:
    """How much energy one vertical label takes up across the axes."""
    low, high = ax.get_xlim()
    width_pt = ax.get_window_extent().width * 72.0 / ax.figure.dpi
    return (high - low) * 1.3 * _LABEL_FONT_PT / max(width_pt, 1.0)


def _make_headroom(ax, counts, yscale: str) -> None:
    """Raise the top of the y axis so the data stay below the labels."""
    peak = float(np.max(counts)) if counts.size else 0.0
    bottom, _ = ax.get_ylim()
    if peak <= bottom:
        return
    if yscale == "log":
        top = bottom * (peak / bottom) ** (1.0 / _DATA_CEILING)
    elif yscale == "sqrt":
        root = np.sqrt(max(bottom, 0.0))
        top = (root + (np.sqrt(peak) - root) / _DATA_CEILING) ** 2
    else:
        top = bottom + (peak - bottom) / _DATA_CEILING
    ax.set_ylim(top=top)


def _draw_marks(ax, marks: list[tuple[float, str]]) -> None:
    """Ticks hanging from the top edge at each line's energy, each led to a
    vertical label. Labels are spread apart to one label width where lines
    crowd, while the ticks stay at the true energies -- the usual X-ray
    spectrum layout."""
    low, high = ax.get_xlim()
    positions = _spread([energy for energy, _ in marks], _label_width(ax), low, high)
    trans = ax.get_xaxis_transform()
    for (energy, label), x in zip(marks, positions):
        ax.plot([energy, energy, x], [1.0, 1.0 - _TICK_END, 1.0 - _LABEL_TOP],
                transform=trans, color="crimson", lw=0.8)
        ax.annotate(label, (x, 1.0 - _LABEL_TOP), xycoords=trans, xytext=(0, -1),
                    textcoords="offset points", rotation=90, ha="center", va="top",
                    fontsize=_LABEL_FONT_PT, color="crimson")


def _spread(wanted: list[float], gap: float, low: float, high: float) -> list[float]:
    """Positions as close to ``wanted`` (sorted) as possible, at least
    ``gap`` apart: each crowded run is laid out evenly around its own mean,
    and runs that then touch are merged, until none overlap."""
    groups = [[x] for x in wanted]  # each group: the wanted positions it holds
    centres = list(wanted)
    changed = True
    while changed:
        changed = False
        for i in range(len(groups) - 1):
            half_i = (len(groups[i]) - 1) * gap / 2
            half_j = (len(groups[i + 1]) - 1) * gap / 2
            if centres[i + 1] - half_j - (centres[i] + half_i) < gap:
                groups[i] += groups.pop(i + 1)
                centres.pop(i + 1)
                centres[i] = sum(groups[i]) / len(groups[i])
                changed = True
                break
    positions = []
    for group, centre in zip(groups, centres):
        half = (len(group) - 1) * gap / 2
        centre = min(max(centre, low + half), high - half)  # stay inside the axes
        positions += [centre - half + k * gap for k in range(len(group))]
    return positions


def _merged(marks, gap: float) -> list[tuple[float, str]]:
    """Cluster ``(energy, symbol, line)`` marks closer than ``gap`` into one
    ``(mean energy, label)``, naming each element once."""
    clusters: list[list[tuple[float, str, str]]] = []
    for mark in marks:
        if clusters and mark[0] - clusters[-1][-1][0] < gap:
            clusters[-1].append(mark)
        else:
            clusters.append([mark])
    merged = []
    for cluster in clusters:
        by_element: dict[str, list[str]] = {}
        for _, symbol, line in cluster:
            by_element.setdefault(symbol, []).append(line)
        label = " / ".join(f"{symbol} {'/'.join(lines)}" for symbol, lines in by_element.items())
        merged.append((sum(m[0] for m in cluster) / len(cluster), label))
    return merged


def refresh(session) -> None:
    """Follow an RBS redraw: redraw the PIXE window while PIXE is enabled."""
    if session.pixe.enabled:
        draw(session, required=False)
