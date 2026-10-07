"""The PIXE window: a second matplotlib figure next to the RBS one.

Kept apart from :mod:`pyrump.shell.plotting` on purpose: the RBS window's
layout, reuse and focus handling stay exactly as they are whether PIXE is on
or not; this module only borrows its yield-axis and show helpers.

**What the window shows** is a *view* (:attr:`PixeState.view`): a list of
items, each a buffer's PIXE spectrum, the PIXE simulation, or one element's
or one layer's part of it, optionally as a comparison (data, simulation and
residuals). While PIXE is enabled the view mirrors the RBS window -- every
RBS PLOT, OVERLAY, SPLOT, COMPARE or REPLOT rebuilds it from the RBS traces
and draws it (:func:`follow`) -- so the same buffers appear in both. The
PIXE prompt's own PLOT and COMPARE set it directly. Reading data and
changing settings only update a window that is already open
(:func:`refresh`); they never open one.

The x axis is energy (bottom) with channels on top; REGION is in channels,
like the RBS window's.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

from .. import __version__
from ..model.spectrum import Spectrum
from .dispatch import CommandError
from .plotting import (
    COMPARE_DATA,
    _apply_limits,
    _apply_scale,
    buffer_label,
    require_matplotlib,
    show,
)

#: Default window size in inches, a little smaller than the RBS window's.
DEFAULT_SIZE = (8.0, 5.0)


@dataclass(frozen=True, slots=True)
class ViewItem:
    """One thing the PIXE window shows."""

    kind: str
    """``"data"`` (a buffer's PIXE spectrum), ``"sim"`` (the simulation),
    ``"element"`` or ``"layer"`` (part of the simulation)."""

    index: int = 0
    """The buffer, for ``"data"``; the 0-based layer, for ``"layer"``."""

    element: str = ""
    """The symbol, for ``"element"``."""


# -- the window -----------------------------------------------------------------


def figure_for(session, *, residuals: bool = False):
    """The PIXE figure, with one panel or a spectrum + residuals pair.

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
    wanted = 2 if residuals else 1
    if len(figure.axes) != wanted:
        figure.clf()
        if residuals:
            figure.subplots(2, 1, sharex=True,
                            gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
        else:
            figure.add_subplot(1, 1, 1)
    state.figure = figure
    return figure


def is_open(session) -> bool:
    """Whether the PIXE window is on screen (not closed by hand)."""
    figure = session.pixe.figure
    return figure is not None and require_matplotlib().fignum_exists(figure.number)


def close(session) -> None:
    """Close the PIXE window, if it is open."""
    if is_open(session):
        require_matplotlib().close(session.pixe.figure)
    session.pixe.figure = None


# -- when to draw -----------------------------------------------------------------


def view_from_traces(traces) -> tuple[list[ViewItem], bool]:
    """The PIXE view matching the RBS window's traces."""
    items: list[ViewItem] = []
    compare = bool(traces) and traces[0].key == COMPARE_DATA
    for trace in traces:
        key = trace.key
        if isinstance(key, str) and key.startswith("splot:element:"):
            items.append(ViewItem("element", element=key.rsplit(":", 1)[1]))
        elif isinstance(key, str) and key.startswith("splot:layer:"):
            items.append(ViewItem("layer", index=int(key.rsplit(":", 1)[1]) - 1))
        elif trace.index == 0:
            items.append(ViewItem("sim"))
        else:
            items.append(ViewItem("data", index=trace.index))
    return items, compare


def follow(session) -> None:
    """After an RBS redraw: show the same buffers in the PIXE window, while
    PIXE is enabled. This is what opens the PIXE window."""
    if not session.pixe.enabled:
        return
    session.pixe.view, session.pixe.compare = view_from_traces(session.traces)
    draw(session, required=False)


def refresh(session) -> None:
    """After new data or a changed setting: forget the cached simulation and
    redraw the PIXE window -- only if it is already open."""
    session.pixe.cache = None
    if session.pixe.enabled and session.pixe.view and is_open(session):
        draw(session, required=False)


# -- the simulation ---------------------------------------------------------------


def _simulation(session):
    """The PIXE simulation (cached until the sample, the active buffer or a
    PIXE setting changes), or ``None``. A failure is reported once, not on
    every redraw, and never stops the data from being drawn."""
    from .pixe_sim import simulate

    state = session.pixe
    if state.cache is not None:
        return state.cache
    try:
        result = simulate(session)
    except (ValueError, KeyError) as error:
        message = str(error).strip("'")
        if message != state.last_error:
            print(f"  PIXE simulation: {message}")
            state.last_error = message
        return None
    state.last_error = None
    state.cache = result
    return result


@dataclass(slots=True)
class _Curve:
    spectrum: Spectrum | None
    label: str
    simulated: bool = False


def _curve(session, item: ViewItem, simulation) -> _Curve:
    """A view item as a spectrum to draw (``None`` with a note when there is
    nothing to show for it)."""
    if item.kind == "data":
        buffer = session.buffers.get(item.index)
        label = buffer_label(session, buffer, item.index) if buffer is not None else ""
        if buffer is None or buffer.pixe is None:
            return _Curve(None, f"buffer {item.index}: no PIXE spectrum")
        return _Curve(buffer.pixe.spectrum, label or f"buffer {item.index}")
    if simulation is None:
        return _Curve(None, "no PIXE simulation" if session.script.layers else "no SIM sample",
                      simulated=True)
    calibration = simulation.calibration
    if item.kind == "sim":
        counts, label = simulation.counts, "simulation"
    elif item.kind == "element":
        counts = simulation.by_element.get(item.element, np.zeros(calibration.npt))
        label = f"SIM({item.element})"
    else:
        counts, label = simulation.by_layer(item.index), f"SIM(layer {item.index + 1})"
    return _Curve(Spectrum(counts=counts, calibration=calibration), label, simulated=True)


# -- drawing ----------------------------------------------------------------------

_COLORS = ("0.20", "steelblue", "darkgreen", "darkorange", "purple")
_SIM_COLORS = ("crimson", "orchid", "teal", "goldenrod")


def draw(session, *, required: bool = True, elements=()) -> bool:
    """Draw the current view. Returns whether the window was drawn.

    ``required`` makes an empty view an error (the PIXE prompt's PLOT with
    nothing to show); otherwise there is simply nothing to do. ``elements``,
    ``(Z, symbol)`` pairs, get their lines marked as well, whatever MARKERS
    says -- ELEMENT's (:func:`mark_elements`).
    """
    state = session.pixe
    if not state.view:
        if required:
            raise CommandError("nothing to plot: PLOT or COMPARE first")
        return False
    simulation = _simulation(session) if any(i.kind != "data" for i in state.view) else None
    curves = [_curve(session, item, simulation) for item in state.view]
    drawable = [c for c in curves if c.spectrum is not None]
    notes = [c.label for c in curves if c.spectrum is None]

    compare = state.compare and len(curves) >= 2 and all(
        c.spectrum is not None for c in curves[:2]
    )
    reference = (drawable[0].spectrum.calibration if drawable
                 else session.pixe.calibration)
    low, high = region_indices(state, reference)

    if compare:
        figure = _draw_comparison(session, curves, (low, high))
        ax = figure.axes[0]
        shown = [c.spectrum.counts[low:high + 1] for c in drawable]
    else:
        figure = figure_for(session)
        ax = figure.axes[0]
        ax.clear()
        shown = []
        sims = 0
        for position, curve in enumerate(drawable):
            spectrum = curve.spectrum
            stop = min(high, spectrum.counts.size - 1)
            if stop <= low:
                continue
            x = spectrum.energies[low:stop + 1]
            if curve.simulated:
                color = _SIM_COLORS[sims % len(_SIM_COLORS)]
                sims += 1
            else:
                color = _COLORS[position % len(_COLORS)]
            ax.step(x, spectrum.counts[low:stop + 1], where="mid",
                    lw=1.2 if curve.simulated else 1.0, color=color, label=curve.label)
            shown.append(spectrum.counts[low:stop + 1])
        ax.set_xlabel("Energy (keV)")
        ax.set_ylabel("Counts")
        _shade_pixe(session, [ax], reference)
        if drawable and state.plot.labels:
            ax.legend(frameon=False, fontsize="small")

    _channel_axis(ax, reference)
    first, last = reference.edge_energy([low, high + 1])
    ax.set_xlim(first, last)
    _apply_scale(ax, state.plot)
    _apply_limits(ax, state.plot)
    if state.plot.yscale == "log" and state.plot.ylow is None:
        # Empty channels would drag a log axis down to its clip floor.
        ax.set_ylim(bottom=0.5)
    if notes:
        ax.text(0.01, 0.97 if not drawable else 0.80, "\n".join(notes), transform=ax.transAxes,
                va="top", ha="left", fontsize="small", color="0.35")
    with warnings.catch_warnings():
        # The channel axis on top is a secondary axis, which tight_layout
        # can't place in the two-panel comparison; it warns but lays out
        # the panels correctly regardless.
        warnings.filterwarnings("ignore", message=".*not compatible with tight_layout")
        figure.tight_layout()
    marked = sample_elements(session) if state.markers != "off" else []
    marked += [element for element in elements if element not in marked]
    if marked:
        marks = _line_marks(ax, marked, major_only=state.markers != "all")
        if marks:
            if state.plot.yhigh is None and shown:
                _make_headroom(ax, np.concatenate(shown), state.plot.yscale)
            _draw_marks(ax, marks)
    show(figure)
    return True


def mark_elements(session, elements) -> bool:
    """ELEMENT's X-ray lines on the PIXE plot, the way MARKERS labels the
    sample's, if a PIXE plot is showing. Like ELEMENT's ticks on the RBS
    plot, they last until the plot is next drawn. Returns whether there was
    a plot to mark."""
    state = session.pixe
    if state.figure is None or not state.view:
        return False
    return draw(session, required=False, elements=list(elements))


def region_indices(state, calibration) -> tuple[int, int]:
    """REGION (channel numbers, as the file numbers them) as indices into
    a spectrum with ``calibration``; the whole spectrum by default."""
    first = int(round(calibration.first))
    low = 0 if state.low is None else max(0, state.low - first)
    high = calibration.npt - 1 if state.high is None else min(calibration.npt - 1, state.high - first)
    if high <= low:
        raise CommandError(f"no PIXE channels between {state.low} and {state.high}")
    return low, high


def _channel_axis(ax, calibration) -> None:
    """Channel numbers along the top edge, matching the energy axis below."""
    gain, offset = calibration.kevch, calibration.kev0
    top = ax.secondary_xaxis(
        "top", functions=(lambda e: (e - offset) / gain, lambda c: c * gain + offset)
    )
    top.set_xlabel("Channel", fontsize="small")
    top.tick_params(labelsize="small")


def _draw_comparison(session, curves, region):
    """Data, simulation and residuals -- the RBS COMPARE layout, in keV."""
    from ..fit.objective import chi_square
    from ..plot.spectra import plot_comparison

    data, theory, *rest = curves
    observed = np.asarray(data.spectrum.counts, dtype=float)
    expected = np.asarray(theory.spectrum.counts, dtype=float)
    n = min(observed.size, expected.size)
    # Scored, and shaded, over PERT's PIXE windows when there are any -- what
    # GO fits -- otherwise over the REGION shown.
    windows = pert_windows(session, data.spectrum.calibration, n)
    if windows is None:
        mask = np.zeros(n, dtype=bool)
        mask[region[0]:region[1] + 1] = True
    else:
        mask = windows
    summary = chi_square(observed[:n], expected[:n], valid=mask, n_parameters=0)
    gof = f"reduced chi-square {summary.reduced:.4f} ({summary.dof} dof)"
    if windows is not None:
        gof += ", PERT PIXE windows"
    figure = figure_for(session, residuals=True)
    for ax in figure.axes:
        ax.clear()
    figure = plot_comparison(
        data.spectrum, theory.spectrum, energy_axis=True, region=region, figure=figure,
        data_label=data.label, simulation_label=theory.label,
        goodness_of_fit=gof,
        overlays=[(c.spectrum, c.label) for c in rest if c.spectrum is not None],
    )
    _shade_pixe(session, figure.axes[:2], data.spectrum.calibration)
    return figure


def _shade_pixe(session, axes, calibration) -> None:
    """PERT's PIXE windows, under HIGHLIGHT -- channels as the .PIX file
    numbers them, on the energy axis."""
    from .plotting import FIT_TONE, shade_windows

    pert = session.pert
    if pert is None or not pert.highlight or not pert.pixe_windows:
        return
    first = round(calibration.first)
    spans = [(w.low - first, w.high - first) for w in pert.pixe_windows]
    shade_windows(axes, spans, calibration.edge_energy, FIT_TONE)


def pert_windows(session, calibration, n: int) -> np.ndarray | None:
    """PERT's PIXE windows as a channel mask over a spectrum of ``n``
    channels with ``calibration``, or ``None`` with no windows set."""
    pert = session.pert
    if pert is None or not pert.pixe_windows:
        return None
    first = round(calibration.first)
    mask = np.zeros(n, dtype=bool)
    for window in pert.pixe_windows:
        low, high = max(window.low - first, 0), min(window.high - first, n - 1)
        if low <= high:
            mask[low:high + 1] = True
    return mask


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
    peak = float(np.nanmax(counts)) if np.isfinite(counts).any() else 0.0
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
