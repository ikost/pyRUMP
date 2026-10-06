"""The session's persistent plot.

RUMP draws to one graphics device whose state -- region, yield range, yield
scaling -- survives between commands, so ``REGION 100 400`` followed by
``REPLOT`` redraws what is already there. The equivalent here is a single
long-lived matplotlib figure owned by the :class:`~pyrump.shell.session.Session`,
plus a list of traces that ``PLOT`` resets and ``OVERLAY`` appends to.
``COMPARE`` resets it too -- to its data and the simulation -- so that an
OVERLAY/SPLOT afterwards adds to the comparison instead of wiping it, just as
RUMP's own COMPARE is literally ``PLOT NOW ... OV THEORY``.

``PLOT``/``OVERLAY`` draw their traces here rather than through
:mod:`pyrump.plot.spectra`, because :class:`~pyrump.shell.session.PlotState`
owns the axis limits and scaling while ``plot_spectrum`` sets its own. The
whole-figure products -- ``COMPARE`` and ``DISPLAY`` -- do reuse
:func:`~pyrump.plot.spectra.plot_comparison` and
:func:`~pyrump.plot.spectra.plot_depth_profile`.
"""

from __future__ import annotations

import re
import time
from dataclasses import dataclass

import numpy as np

from .. import __version__
from ..model.detector import yield_normalisation
from ..model.spectrum import Spectrum
from . import terminal_focus
from .dispatch import CommandError

#: Colour cycle for overlaid spectra. The first is RUMP's white-on-black data
#: trace; the rest keep overlays distinguishable in both light and dark themes.
_COLORS = ("0.20", "crimson", "steelblue", "darkgreen", "darkorange", "purple")

#: The key COMPARE files its data trace under. As the first trace it switches
#: :func:`draw` to the two-panel comparison; anything that resets the traces
#: (PLOT, DISPLAY, AXIS, ...) drops it, and with it the comparison layout.
COMPARE_DATA = "compare:data"


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
    #: What :func:`add_trace`'s ``replace`` dedupes on -- usually just
    #: ``index``, but a selective SPLOT overlay uses a key of its own (see
    #: ``add_trace``) so it doesn't collide with buffer 0's own trace.
    key: object


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


def _energy_axis(ax, calibration) -> None:
    """RUMP's upper axis: energy in MeV above the channel axis (RbsAxdraw,
    tplot.c:469-483), scaled from ``calibration`` -- the first trace's, as
    RUMP scales it from its reference buffer. A channel's energy is its lower
    edge, RBSENERGY's convention (rumpproto.h:7), so a tick reads what the
    ENERGY axis would put at the same place. Skipped for a buffer with no
    usable calibration."""
    if not calibration.kevch:
        return
    top = ax.secondary_xaxis(
        "top",
        functions=(
            lambda c: calibration.edge_energy(c) / 1000.0,
            lambda e: calibration.channel_of(np.asarray(e) * 1000.0),
        ),
    )
    top.set_xlabel("Energy (MeV)", fontsize="small")
    top.tick_params(labelsize="small")


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


def _apply_limits(ax, state) -> None:
    """COUNTS/YLOW/YHIGH if set, else autoscale -- floored at zero on a
    linear axis.

    A log axis can't show zero or below, so a ``COUNTS 0 <max>``/``BLOWUP``
    floor is left to autoscaling there rather than handed to matplotlib,
    which would ignore it with a warning.
    """
    ax.autoscale(axis="y")
    bottom = state.ylow
    if state.yscale == "log" and bottom is not None and bottom <= 0:
        bottom = None
    if bottom is not None or state.yhigh is not None:
        ax.set_ylim(bottom=bottom, top=state.yhigh)
    elif state.yscale == "linear":
        ax.set_ylim(bottom=0)


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
        figure.canvas.manager.set_window_title(f"pyRUMP {__version__}")
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


def is_comparison(session) -> bool:
    """Whether the traces are a COMPARE (plus anything overlaid on it)."""
    return bool(session.traces) and session.traces[0].key == COMPARE_DATA


def main_axes(session):
    """The figure and the axes the spectra are drawn on -- COMPARE's top
    panel, or PLOT's only one -- without changing the figure's layout."""
    if is_comparison(session):
        figure = compare_figure_for(session, residuals=True)
        return figure, figure.axes[0]
    return figure_for(session)


def draw(session) -> None:
    """Render every trace according to the current :class:`PlotState`."""
    if not session.traces:
        raise CommandError("nothing to plot yet")
    if is_comparison(session):
        _draw_comparison(session)
        return

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
    if not state.energy_axis:
        _energy_axis(ax, session.traces[0].buffer.calibration)
    ax.set_ylabel("Yield (counts/msr/uC)" if state.normalized else "Counts")
    _apply_scale(ax, state)
    _apply_limits(ax, state)

    if state.labels and any(t.label for t in session.traces):
        ax.legend(frameon=False, fontsize="small")

    figure.tight_layout()
    show(figure)
    _refresh_pixe(session)


def _refresh_pixe(session) -> None:
    """While PIXE is enabled, show the same buffers in the PIXE window."""
    from . import pixe_plotting

    pixe_plotting.follow(session)


# -- keeping the simulation curves current ---------------------------------------


def is_open(session) -> bool:
    """Whether the RBS window is on screen (not closed by hand)."""
    figure = session.figure
    return figure is not None and require_matplotlib().fignum_exists(figure.number)


def simulation_for(session, key):
    """The simulation a trace filed under ``key`` shows: the full one for
    ``0`` (buffer 0, cached until the sample changes), or a selective SPLOT's
    ``"splot:element:<symbol>"`` / ``"splot:layer:<n>"``, computed afresh.
    ``None`` when that element or layer is no longer in the sample."""
    if not isinstance(key, str):
        return session.simulation()
    kind, _, which = key.removeprefix("splot:").partition(":")
    if kind == "element":
        z = session.table.by_symbol(which).z
        if z not in {session.table.by_symbol(s).z for s in session.script.elements}:
            return None
        return session.selective_simulation(element_z=z, label=f"SIM({which})")
    layer = int(which)
    if not 1 <= layer <= len(session.script.layers):
        return None
    return session.selective_simulation(layer=layer - 1, label=f"SIM(layer {layer})")


def update_simulations(session) -> bool:
    """Recompute every simulation curve on the plot (the full simulation and
    any selective SPLOTs) against the current sample and settings. Returns
    whether anything changed.

    A SPLOT whose element or layer has gone is dropped, with a note. A
    simulation that fails is reported once and leaves the curves as they
    were -- the plot is never broken by it.
    """
    traces = []
    changed = False
    try:
        for trace in session.traces:
            if trace.index != 0 or trace.key == COMPARE_DATA:
                traces.append(trace)
                continue
            buffer = simulation_for(session, trace.key)
            if buffer is None:
                print(f"  SPLOT {trace.key.rsplit(':', 1)[1]} removed: no longer in the sample")
                changed = True
                continue
            if buffer is not trace.buffer:
                trace = Trace(buffer=buffer, label=buffer_label(session, buffer, 0),
                              index=0, key=trace.key)
                changed = True
            traces.append(trace)
    except (KeyError, ValueError) as error:
        message = str(error).strip("'")
        if message != session.sim_error:
            print(f"  simulation: {message}")
            session.sim_error = message
        return False
    session.sim_error = None
    session.traces = traces
    return changed


def follow_simulation(session) -> None:
    """After a command (a whole XEQ counting as one): if it changed the
    sample or a simulation setting, bring the simulation curves on the RBS
    plot up to date and redraw it -- which brings the PIXE window along --
    while LIVE is on. Neither window is ever opened by this.

    With LIVE off the change is remembered, so turning LIVE back on catches
    the plot up at once.
    """
    if not (session.live and session.sim_changed):
        return
    session.sim_changed = False
    try:
        if session.traces and is_open(session) and update_simulations(session):
            draw(session)
        elif session.pixe.cache is None:
            # Nothing redrawn since the change: the PIXE window catches up on
            # its own (a no-op unless it is open).
            from . import pixe_plotting

            pixe_plotting.refresh(session)
    except CommandError as error:
        print(f"  redraw: {error}")


def goodness_of_fit(
    session, data, theory, n_channels: int, region: tuple[int, int], *, raw: bool = False
) -> str:
    """"reduced chi-square X.XX (Y dof)" -- over PERT's error windows (and
    normalisation window) if any are set, so the number matches what GO
    itself reports; otherwise over the plot's own visible REGION, with no
    normalisation and no parameters subtracted from dof.

    Scored on the yield as plotted, so under NORMALIZE it is no longer a
    count-based chi-square -- the text says so on a second line. ``raw``
    scores the raw counts whatever NORMALIZE says (EXPORTCMP's columns).
    """
    from ..fit.objective import chi_square

    state = session.plot
    normalized = state.normalized and not raw
    if normalized:
        observed = _values(data, state)[:n_channels]
        expected = _values(theory, state)[:n_channels]
    else:
        observed = np.asarray(data.spectrum.counts, dtype=float)[:n_channels]
        expected = np.asarray(theory.spectrum.counts, dtype=float)[:n_channels]
    pert = session.pert
    if pert is not None and pert.windows.error:
        mask = pert.windows.mask(n_channels)
        scale = pert.windows.normalisation_factor(observed, expected)
        n_parameters = len(pert.varying)
    else:
        low, high = region
        mask = np.zeros(n_channels, dtype=bool)
        mask[low : high + 1] = True
        scale = 1.0
        n_parameters = 0
    summary = chi_square(observed * scale, expected, valid=mask, n_parameters=n_parameters)
    text = f"reduced chi-square {summary.reduced:.4f} ({summary.dof} dof)"
    if normalized:
        text += "\nnormalized yield: GOF not from raw counts"
    return text


def _draw_comparison(session) -> None:
    """COMPARE's data and simulation with residuals, plus any OVERLAY/SPLOT
    curves added since, in the top panel.

    The simulation is the most recent buffer-0 trace, so a bare SPLOT after
    COMPARE (which replaces it) updates the residuals and chi-square too; a
    selective SPLOT, keyed apart from buffer 0, is only overlaid. REGION
    applies to both panels, and COUNTS/YLOW/YHIGH and the yield scale to the
    top one, as for PLOT. NORMALIZE applies to every curve, and so to the
    residuals and chi-square too (see :func:`goodness_of_fit`).
    """
    from ..plot.spectra import plot_comparison

    state = session.plot

    def spectrum(buffer) -> Spectrum:
        return Spectrum(counts=_values(buffer, state), calibration=buffer.spectrum.calibration)

    data, *rest = session.traces
    theory = next(t for t in reversed(rest) if t.key == 0)
    overlays = [(spectrum(t.buffer), t.label) for t in rest if t is not theory]

    n_channels = min(data.buffer.n_channels, theory.buffer.n_channels)
    try:
        region = state.region(n_channels)
    except ValueError as error:
        raise CommandError(str(error)) from None

    figure = compare_figure_for(session, residuals=True)
    figure = plot_comparison(
        spectrum(data.buffer), spectrum(theory.buffer),
        energy_axis=state.energy_axis, region=region, figure=figure,
        data_label=data.label, simulation_label=theory.label,
        goodness_of_fit=goodness_of_fit(session, data.buffer, theory.buffer, n_channels, region),
        overlays=overlays,
    )
    top = figure.axes[0]
    if not state.energy_axis:
        _energy_axis(top, data.buffer.calibration)
    if state.normalized:
        top.set_ylabel("Yield (counts/msr/uC)")
    _apply_scale(top, state)
    _apply_limits(top, state)
    session.figure = figure
    show(figure)
    _refresh_pixe(session)


def _interactive(canvas) -> bool:
    """Whether this canvas is backed by a real GUI toolkit.

    Agg -- what the tests and ``--batch`` runs use -- has no
    ``required_interactive_framework``, no event loop, and no window to keep
    alive, so every GUI-facing step here is skipped for it.
    """
    return getattr(type(canvas), "required_interactive_framework", None) is not None


#: Shortest gap between two pumps of the GUI event loop, in seconds. 0.1 s is
#: far below the ~5 s after which Windows declares a window "Not Responding"
#: (and macOS shows the beachball), while costing well under 1% of a fit's
#: runtime -- one ``flush_events`` on an idle window is a fraction of a
#: millisecond.
PUMP_INTERVAL = 0.1

#: When :func:`pump` last reached the event loop. Module state rather than
#: session state because the throttle is about wall-clock GUI latency, not
#: about any one session.
_last_pump = 0.0


def pump(session) -> None:
    """Give the plot window a slice of event-loop time mid-command.

    Nothing pumps the GUI while a command runs, so the figure stops
    repainting for as long as the command takes -- past about five seconds
    the window managers on all three platforms mark it hung (a greyed title
    bar and "Not Responding" on Windows, the beachball on macOS, the WM's own
    force-quit prompt on Linux). A ``PERT GO`` over a fuzzed sample measures
    eighteen seconds, almost all of it in that state.

    Long-running commands call this from whatever per-iteration hook they
    already have; it self-throttles to :data:`PUMP_INTERVAL` so a caller can
    invoke it as often as is convenient. Does nothing when there is no
    figure, when the user has closed it, or on a non-GUI backend.

    Draining the event queue also delivers *input*, so a click on the figure's
    close button during a fit is acted on rather than queued until the fit
    ends. That makes the closed-figure check below load-bearing: a pump can be
    what closes the window that the next pump would otherwise reach into.
    """
    global _last_pump
    figures = [f for f in (session.figure, session.pixe.figure) if f is not None]
    if not figures:
        return
    now = time.monotonic()
    if now - _last_pump < PUMP_INTERVAL:
        return
    _last_pump = now

    plt = require_matplotlib()
    for figure in figures:
        if not plt.fignum_exists(figure.number) or not _interactive(figure.canvas):
            continue
        try:
            figure.canvas.flush_events()
        except (AttributeError, NotImplementedError, RuntimeError):  # pragma: no cover
            # A backend without an event loop, or a window torn down mid-pump.
            pass


def show(figure) -> None:
    """Push the figure to the screen without blocking the prompt.

    Raising a GUI window steals focus from wherever the user was typing --
    macOS in particular makes the new window "key". For an interactive
    backend (one with a real ``required_interactive_framework``, unlike Agg
    in tests) the terminal's focus is captured before drawing and restored
    right after, via :mod:`~pyrump.shell.terminal_focus`.
    """
    canvas = figure.canvas
    interactive = _interactive(canvas)
    token = terminal_focus.capture() if interactive else None
    canvas.draw_idle()
    try:
        canvas.flush_events()
    except (AttributeError, NotImplementedError):  # Agg in tests
        pass
    if interactive:
        terminal_focus.restore(token)


def mark_whatisit(session, candidates, target_keV) -> bool:
    """Overlay WHATISIT's candidates as ticks along the plot's bottom edge.

    ``RbsMark``'s ``MK_WH1``/``MK_WHL``/``MK_WHR`` (tplot.c:260-292): the best
    match gets a solid tick, its neighbors dashed ones, each labelled with the
    element symbol. Skips the C's pixel-exact anti-collision jog -- readability
    over imitation, the same call this module's ``spectra.py`` sibling already
    makes -- and only marks a plot that's already on screen, mirroring
    ``RbsLocate``'s own gate: with no plot device it reports that and draws
    nothing, rather than opening one.

    Returns ``False`` (having drawn nothing) if there is no active
    PLOT/OVERLAY/COMPARE to mark. On a COMPARE the ticks go on its top panel.
    """
    if not session.traces:
        return False
    require_matplotlib()
    figure, ax = main_axes(session)
    trans = ax.get_xaxis_transform()
    low, high = ax.get_xlim()
    best = min(candidates, key=lambda c: abs(c.energy_keV - target_keV))
    for row, candidate in enumerate(candidates):
        x = candidate.energy_keV if session.plot.energy_axis else candidate.channel
        if x < low or x > high:
            continue
        is_best = candidate is best
        y = 0.05 if is_best else 0.05 + 0.05 * (row % 2)
        ax.plot(
            [x, x], [0.0, y], transform=trans, clip_on=False, color="0.15",
            lw=1.4 if is_best else 1.0, ls="-" if is_best else "--",
        )
        ax.annotate(
            candidate.symbol, (x, y), xycoords=trans,
            xytext=(0, 3), textcoords="offset points",
            ha="center", va="bottom", fontsize="small",
            fontweight="bold" if is_best else "normal",
        )
    show(figure)
    return True


def mark_element(session, marks) -> bool:
    """Tick ELEMENT's surface edges along the plot's bottom edge.

    ``RbsMark(MK_TKL, ...)`` (anlytc.c:250): one solid tick per element at its
    predicted edge, labelled with the symbol, drawn the same way as
    :func:`mark_whatisit`'s best match. ``marks`` is a sequence of
    ``(energy_keV, channel, label)``. Ticks alternate between two heights in
    order along the axis, whatever order the elements were typed in, so
    labels of close edges don't print on top of each other. The
    C's optional marker height (a typed value or a cursor pick) is not
    carried over: the ticks always sit on the axis.

    Draws nothing, silently, when there is no active PLOT/OVERLAY/COMPARE
    (anlytc.c:250's ``if (PlotSystem...)``), and skips any edge outside the
    current x range. Returns whether a plot was there to mark. On a COMPARE
    the ticks go on its top panel.
    """
    if not session.traces:
        return False
    require_matplotlib()
    figure, ax = main_axes(session)
    trans = ax.get_xaxis_transform()
    low, high = ax.get_xlim()
    xs = sorted(
        (energy_keV if session.plot.energy_axis else channel, label)
        for energy_keV, channel, label in marks
    )
    visible = [(x, label) for x, label in xs if low <= x <= high]
    for row, (x, label) in enumerate(visible):
        y = 0.05 + 0.05 * (row % 2)
        ax.plot([x, x], [0.0, y], transform=trans, clip_on=False, color="0.15", lw=1.4)
        ax.annotate(
            label, (x, y), xycoords=trans,
            xytext=(0, 3), textcoords="offset points",
            ha="center", va="bottom", fontsize="small", fontweight="bold",
        )
    show(figure)
    return True


def mark_matrix(session, buffer, energy_keV: float, channel: float, height: float, symbol: str) -> bool:
    """Crosshair MATRIX's predicted point directly onto the plot.

    ``RbsMark(MK_XHR, ...)`` (anlytc.c:274): a small "+" at the predicted
    ``(energy, height)``, labelled with the element symbol -- unlike
    WHATISIT's bottom ticks, this sits right on the spectrum where the step
    should be. ``height`` arrives in normalized yield units
    (counts/uC/keV/msr, :attr:`MatrixResult.height`'s own convention) and is
    converted to raw counts here when the plot isn't in NORMALIZE mode --
    the same factor :func:`_values` already applies the other way.

    Draws nothing (and RUMP's own C prints no warning for it either,
    anlytc.c:274's silent ``if (PlotSystem...)``) when there is no active
    PLOT/OVERLAY/COMPARE, or when the point falls outside the current x range. The
    y range, unlike x, is expanded to fit the marker when it's taller than
    the current view -- ``draw()`` freezes it via ``set_ylim`` on every
    PLOT/OVERLAY, so without this a tall prediction would be silently
    clipped -- unless YLOW/YHIGH pinned it explicitly, which still wins.
    """
    if not session.traces:
        return False
    require_matplotlib()
    figure, ax = main_axes(session)
    x = energy_keV if session.plot.energy_axis else channel
    low, high = ax.get_xlim()
    if x < low or x > high:
        return False
    y = height if session.plot.normalized else height * yield_normalisation(buffer.measurement)
    ax.plot([x], [y], marker="+", markersize=12, markeredgewidth=1.6, color="0.15", linestyle="none")
    ax.annotate(
        symbol, (x, y), xytext=(0, 8), textcoords="offset points",
        ha="center", va="bottom", fontsize="small", fontweight="bold",
    )
    state = session.plot
    ylow, yhigh = ax.get_ylim()
    if state.yhigh is None and y > yhigh:
        ax.set_ylim(top=y * 1.1)
    if state.ylow is None and y < ylow:
        ax.set_ylim(bottom=y * 1.1 if y < 0 else 0)
    show(figure)
    return True


#: Extensions buffer_stem drops -- the spectrum formats GET reads. Anything
#: else after a dot is part of the sample name (``Ta2.5``), not a suffix.
_SPECTRUM_SUFFIXES = {".rbs", ".dat", ".asc", ".ascii", ".txt", ".xls"}


def buffer_stem(buffer) -> str:
    """A short, filesystem- and legend-safe sample ID for this buffer.

    The spectrum's own ``IDENTIFIER`` names the sample, so it wins; the
    buffer's name (usually the file it came from) is only the fallback for a
    spectrum without one. Either may hold a full path rather than a bare
    filename -- e.g. a WRASCII macro's own ``FILENAME`` line stamping a
    Windows path (``C:\\RBS\\data\\...\\MA8410.RBS``) straight into
    ``buffer.name``. Splits on both slash conventions regardless of host
    OS (``Path.stem`` alone only understands the platform's own separator),
    then keeps just the first whitespace-separated token, so a descriptive
    trailing comment (``"MA8410.RBS  170 Degree RBS LT = ..."``) doesn't
    leak in either, and drops a spectrum-file extension. ``""`` if neither
    yields anything usable, leaving the fallback to the caller.
    """

    def sanitize(text: str) -> str:
        if not text:
            return ""
        tail = re.split(r"[\\/]", text.strip())[-1]
        tail = tail.split()[0] if tail.split() else tail
        stem, dot, suffix = tail.rpartition(".")
        if dot and stem and f".{suffix.lower()}" in _SPECTRUM_SUFFIXES:
            return stem
        return tail

    return sanitize(buffer.identifier) or sanitize(buffer.name)


def buffer_label(session, buffer, index: int) -> str:
    """The name PLOT/OVERLAY would show for this buffer in a legend.

    Buffer 0 is usually the full simulation (Session.simulation()'s
    convention), so when STRUCTLABEL is on and a sample is described, its
    legend text becomes a compact rendering of the layer structure instead of
    the buffer's own name/identifier ("SIM"). But a selective SPLOT overlay
    also lands at index 0 without ever being stored into buffer 0 (see
    Session.selective_simulation) precisely so it keeps its own "SIM(Ru)" /
    "SIM(layer 2)" name here instead of being mistaken for the whole sample --
    checked by identity, not position, against what buffer 0 actually holds.
    """
    if index == 0 and session.plot.structure_labels and buffer is session.buffers.get(0):
        from ..script.lcm import structure_label

        label = structure_label(
            session.script, normalize=session.plot.composition_fraction
        )
        if label:
            return label
    if index == 0:
        # Buffer 0 is always a simulation (full or a selective SPLOT), never
        # data read off disk, so its name is already a clean, purpose-built
        # caption ("SIM", "SIM(Ru)", "SIM(layer 2)") -- not a filename in
        # need of buffer_stem's path/comment stripping.
        return buffer.name or buffer.identifier or f"buffer {index}"
    return buffer_stem(buffer) or f"buffer {index}"


def add_trace(
    session, index: int, buffer, *, clear: bool, replace: bool = False, key: object = None
) -> None:
    """Add a buffer to the plot; ``clear`` makes it a fresh ``PLOT``.

    ``replace`` drops any existing trace with the same ``key`` first (default
    ``index``, i.e. the same buffer slot), so a command that re-plots the same
    thing on every call updates it in place instead of piling up a fresh copy
    each time. ``OVERLAY`` leaves this off -- stacking distinct buffers is the
    point of it.

    A plain re-plot of buffer 0 (PLOT/OVERLAY 0, bare SPLOT) all key on
    ``index`` and so dedupe against each other, since they show the identical
    theory spectrum. A *selective* SPLOT (one element or layer's contribution)
    passes its own ``key`` instead -- distinct from plain ``0`` and from each
    other -- so it neither collides with a full-simulation OVERLAY 0 already
    on the plot nor with a different SPLOT selection, while a repeated call
    for the *same* selection still updates in place.
    """
    if key is None:
        key = index
    if clear:
        session.traces = []
    elif replace:
        session.traces = [t for t in session.traces if t.key != key]
    label = buffer_label(session, buffer, index)
    session.traces.append(Trace(buffer=buffer, label=label, index=index, key=key))
