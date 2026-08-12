"""The top-level RUMP command table.

Command names, minimum abbreviations and grouping follow
``C-code/html/RUMP/cmds.htm`` and the ``cmlist`` table at rump.c:147. Handlers
are deliberately thin: parse arguments, call into ``pyrump.*``, mutate the
session.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ...model.geometry import GeometryKind
from ...model.spectrum import Calibration, Spectrum
from ..dispatch import ArgReader, CommandError, CommandTable
from ..session import Buffer
from .. import plotting

#: RUMP's own maximum channel count (rumpdata.h: CMAX), used to bound SWALLOW.
_CMAX = 16384


class Quit(Exception):
    """Raised by QUIT/BYE to unwind out of the REPL."""


class EnterMode(Exception):
    """Raised by SIM/PERT to push a sub-level onto the mode stack."""

    def __init__(self, name: str):
        super().__init__(name)
        self.name = name


class Return(Exception):
    """Raised by RETURN to pop a sub-level."""


# ---------------------------------------------------------------------------
# Buffer loading
# ---------------------------------------------------------------------------


def read_into_buffer(session, path: Path, index: int | None = None) -> int:
    """Read a spectrum file into a buffer and make it ACTIVE.

    If a buffer already holds this file, it is selected instead of re-read --
    the behaviour cmds.htm describes for PLOT.
    """
    from ...cli._common import read_spectrum
    from ...io.rbs import RbsSpectrum

    path = Path(path).expanduser()
    if not path.exists() and not path.suffix:
        # RUMP defaults the extension to .RBS.
        candidate = path.with_suffix(".rbs")
        if candidate.exists():
            path = candidate
    if not path.exists():
        raise CommandError(f"no such file: {path}")
    # Resolve once, here: buffers must hold absolute paths so that a CD between
    # two GETs of the same file does not load it into a second buffer.
    path = path.resolve()

    existing = session.buffers.find_path(path)
    if existing is not None:
        session.buffers.active = existing
        session.touch()
        return existing

    try:
        source = read_spectrum(path)
    except Exception as error:
        raise CommandError(f"could not read {path}: {error}") from None

    if isinstance(source, RbsSpectrum):
        buffer = Buffer.from_rbs(source, path)
    else:
        buffer = Buffer.from_ascii(source, path)

    slot = session.buffers.load(buffer, index)
    session.buffers.active = slot
    session.touch()
    return slot


def _buffer_argument(session, args: ArgReader, *, read: bool = True) -> tuple[int, object]:
    """Resolve a buffer token, optionally reading a file if it names one.

    With no argument at all, the ACTIVE buffer is used.
    """
    token = args.optional()
    if token is None:
        index = session.buffers.active
        buffer = session.buffers.get(index)
        if buffer is None and index == 0:
            buffer = session.simulation()
        if buffer is None:
            raise CommandError("no active buffer: read a spectrum first (GET <file>)")
        return index, buffer

    try:
        index = session.resolve(token)
    except KeyError:
        if not read:
            raise CommandError(f"no buffer named {token!r}") from None
        index = read_into_buffer(session, Path(token))
        return index, session.buffers[index]

    if index == 0:
        return 0, session.simulation()
    buffer = session.buffers.get(index)
    if buffer is None:
        raise CommandError(f"buffer {index} is empty")
    return index, buffer


# ---------------------------------------------------------------------------
# Session control
# ---------------------------------------------------------------------------


def cmd_help(session, args: ArgReader) -> None:
    """List every command, one section per table -- as rump.c:351-360 does."""
    args.done()
    from .system import TABLE as SYSTEM_TABLE

    print(TABLE.help_text())
    print()
    print(SYSTEM_TABLE.help_text())


def cmd_quit(session, args: ArgReader) -> None:
    raise Quit()


def cmd_sim(session, args: ArgReader) -> None:
    if args:
        # "SIM Reset", "SIM thick 500" -- one-shot form, handled by SIM's table.
        from .sim import execute_in_sim

        execute_in_sim(session, args)
        return
    raise EnterMode("sim")


def cmd_pert(session, args: ArgReader) -> None:
    args.done()
    raise EnterMode("pert")


def cmd_return(session, args: ArgReader) -> None:
    raise Return()


def cmd_data(session, args: ArgReader) -> None:
    """Report, or switch, the directory the atomic tables came from."""
    token = args.optional()
    args.done()
    if token is None:
        print(session.data)
        return
    from ...cli._common import data_dir, load_tables

    directory = data_dir(token)
    session.table, session.registry, session.densities = load_tables(directory)
    session.data = directory
    session.touch()
    print(f"tables reloaded from {directory}")


# ---------------------------------------------------------------------------
# Buffer control
# ---------------------------------------------------------------------------


def cmd_buffers(session, args: ArgReader) -> None:
    args.done()
    print(session.buffers.listing())


def cmd_get(session, args: ArgReader) -> None:
    token = args.token("a file or buffer number")
    args.done()
    try:
        index = session.resolve(token)
    except KeyError:
        index = read_into_buffer(session, Path(token))
    else:
        if index and session.buffers.get(index) is None:
            raise CommandError(f"buffer {index} is empty")
        session.buffers.active = index
        session.touch()
    print(f"active buffer is now {index}")


def cmd_pointat(session, args: ArgReader) -> None:
    index = args.integer("a buffer number")
    args.done()
    if index and session.buffers.get(index) is None:
        raise CommandError(f"buffer {index} is empty")
    session.buffers.active = index
    session.touch()
    print(f"active buffer is now {index}")


def cmd_release(session, args: ArgReader) -> None:
    token = args.optional()
    args.done()
    index = session.buffers.active if token is None else session.resolve(token)
    try:
        session.buffers.release(index)
    except KeyError as error:
        raise CommandError(str(error)) from None
    print(f"released buffer {index}; active is now {session.buffers.active}")


def cmd_newall(session, args: ArgReader) -> None:
    args.done()
    session.buffers.clear()
    session.traces = []
    session.touch()
    print("all buffers cleared")


def cmd_copy(session, args: ArgReader) -> None:
    source = session.resolve(args.token("a source buffer"))
    target = args.integer("a target buffer")
    args.done()
    buffer = session.simulation() if source == 0 else session.buffers.get(source)
    if buffer is None:
        raise CommandError(f"buffer {source} is empty")
    if target == 0:
        raise CommandError("buffer 0 is the simulation and cannot be written to")
    session.buffers.set(target, buffer.copy())
    print(f"copied buffer {source} to {target}")


def cmd_move(session, args: ArgReader) -> None:
    first = session.resolve(args.token("a buffer"))
    second = session.resolve(args.token("a buffer"))
    args.done()
    if 0 in (first, second):
        raise CommandError("buffer 0 is the simulation and cannot be exchanged")
    a, b = session.buffers.get(first), session.buffers.get(second)
    session.buffers.set(first, b)
    session.buffers.set(second, a)
    print(f"exchanged buffers {first} and {second}")


def cmd_write(session, args: ArgReader) -> None:
    from ...io.rbs import write_rbs

    target = Path(args.token("an output file"))
    args.done()
    buffer = session.buffers.require_active()
    write_rbs(target, buffer.to_rbs())
    print(f"wrote {target}")


def cmd_wrascii(session, args: ArgReader) -> None:
    from ...io.ascii import write_ascii

    target = Path(args.token("an output file"))
    args.done()
    buffer = session.buffers.require_active()
    write_ascii(target, buffer.spectrum.counts, identifier=buffer.identifier)
    print(f"wrote {target}")


def cmd_recalculate(session, args: ArgReader) -> None:
    args.done()
    session.touch()
    print("simulation marked for recalculation")


# ---------------------------------------------------------------------------
# Spectrum parameters
# ---------------------------------------------------------------------------


def cmd_active(session, args: ArgReader) -> None:
    args.done()
    index = session.buffers.active
    buffer = session.buffers.get(index)
    if buffer is None:
        raise CommandError("no active buffer")
    print(f"Buffer {index}:")
    print(buffer.describe())


def _chain(session, args: ArgReader) -> None:
    """Dispatch any tokens left on the line as a further RUMP command.

    RUMP's buffer-parameter setters only ever consume as many tokens as they
    need and hand the rest onward -- cmds.htm's own SWALLOW example chains
    ``empty swallow``, and real WRASCII output chains ``Choff 0 FWHM 15`` --
    unlike most of the shell's other commands, which reject any leftover.
    """
    if not args:
        return
    name = args.token()
    command = TABLE.match(name)
    if command is None:
        raise CommandError(f"unrecognized command: {name}")
    command.handler(session, ArgReader(args.remaining, command=command.name.lower()))


def _setter(apply, describe):
    """Build a handler that reports the value with no argument and sets it with one."""

    def handler(session, args: ArgReader) -> None:
        buffer = session.buffers.require_active()
        if not args:
            print(describe(buffer))
            return
        apply(session, buffer, args)
        session.touch()
        print(describe(buffer))
        _chain(session, args)

    return handler


def _numeric(field, target, label, unit=""):
    """A setter for one numeric field of a frozen parameter dataclass.

    ``target`` names the ``Buffer`` attribute holding it -- ``calibration``,
    ``geometry`` or ``measurement`` -- and the matching ``set_*`` method
    rebuilds it, since all three are frozen.
    """

    def apply(session, buffer, args: ArgReader) -> None:
        getattr(buffer, f"set_{target}")(**{field: args.number(label)})

    def describe(buffer) -> str:
        return f"  {label} = {getattr(getattr(buffer, target), field):g}{unit}"

    return _setter(apply, describe)


def cmd_beam(session, args: ArgReader) -> None:
    """``BEAM 4He++`` -- species and charge state."""
    from ...cli._common import resolve_beam

    buffer = session.buffers.require_active()
    if not args:
        print(
            f"  beam Z={buffer.beam.z} mass={buffer.beam.mass:.4f}"
            f" charge state {buffer.measurement.charge_state}"
        )
        return
    token = args.token("a beam species, e.g. 4He++")
    charge = token.count("+") or None
    symbol = token.rstrip("+")
    try:
        z, mass = resolve_beam(session.table, symbol)
    except Exception as error:
        raise CommandError(f"unknown beam {symbol!r}: {error}") from None
    buffer.beam.z = z
    buffer.beam.mass = mass
    if charge is not None:
        buffer.beam.charge_state = charge
        buffer.set_measurement(charge_state=charge)
    session.touch()
    print(f"  beam Z={z} mass={mass:.4f} charge state {buffer.measurement.charge_state}")
    _chain(session, args)


def cmd_mev(session, args: ArgReader) -> None:
    buffer = session.buffers.require_active()
    if not args:
        print(f"  MeV = {buffer.beam.e0_MeV:g}")
        return
    buffer.beam.e0_MeV = args.number("the beam energy in MeV")
    session.touch()
    print(f"  MeV = {buffer.beam.e0_MeV:g}")
    _chain(session, args)


def cmd_geometry(session, args: ArgReader) -> None:
    buffer = session.buffers.require_active()
    if not args:
        print(f"  geometry = {buffer.geometry.kind.name.lower()}")
        return
    token = args.token("cornell, ibm or general")
    try:
        kind = GeometryKind[token.upper()]
    except KeyError:
        raise CommandError(
            f"unknown geometry {token!r}: expected cornell, ibm or general"
        ) from None
    buffer.set_geometry(kind=kind)
    session.touch()
    print(f"  geometry = {kind.name.lower()}")
    _chain(session, args)


def cmd_conversion(session, args: ArgReader) -> None:
    """``CONVERSION <keV/ch> [keV(0)]``."""
    buffer = session.buffers.require_active()
    if not args:
        c = buffer.calibration
        print(f"  {c.kevch:g} keV/channel, offset {c.kev0:g} keV")
        return
    kevch = args.number("keV per channel")
    kev0 = args.optional_number()
    changes = {"kevch": kevch}
    if kev0 is not None:
        changes["kev0"] = kev0
    buffer.set_calibration(**changes)
    session.touch()
    c = buffer.calibration
    print(f"  {c.kevch:g} keV/channel, offset {c.kev0:g} keV")
    _chain(session, args)


def cmd_identifier(session, args: ArgReader) -> None:
    buffer = session.buffers.require_active()
    if not args:
        print(f"  {buffer.identifier}")
        return
    buffer.identifier = args.rest()
    print(f"  {buffer.identifier}")


def cmd_date(session, args: ArgReader) -> None:
    """``DATE <text>`` -- when the spectrum was measured (bmanip.c: B1_DATE)."""
    buffer = session.buffers.require_active()
    if not args:
        print(f"  {buffer.date}")
        return
    buffer.date = args.token("a date")
    session.touch()
    print(f"  {buffer.date}")
    _chain(session, args)


def cmd_filename(session, args: ArgReader) -> None:
    """``FILENAME <name>`` -- record a buffer's source name (bmanip.c: B1_FILE)."""
    buffer = session.buffers.require_active()
    if not args:
        print(f"  {buffer.name}")
        return
    buffer.name = args.token("a filename")
    print(f"  {buffer.name}")
    _chain(session, args)


def cmd_empty(session, args: ArgReader) -> None:
    """``EMPTY [n]`` -- reset buffer *n* to blank (bmanip.c: B_EMPTY).

    RUMP always has a "current" buffer to reset; pyRUMP buffers do not exist
    until something is loaded, so with no buffer active this opens a fresh
    one instead. That is how WRASCII macros bootstrap themselves -- they
    always open with an ``Empty ...`` line (cmds.htm's SWALLOW example is
    literally ``empty swallow``).
    """
    token = args.peek()
    index: int | None = None
    if token is not None:
        try:
            index = int(float(token))
        except ValueError:
            index = None
    if index is not None:
        args.token()
        if index == 0:
            raise CommandError("buffer 0 is the simulation and cannot be emptied")
    elif session.buffers.active_buffer is not None:
        index = session.buffers.active
    else:
        index = session.buffers.first_free()
    session.buffers.set(index, Buffer(spectrum=Spectrum.zeros(Calibration())))
    session.buffers.active = index
    session.touch()
    print(f"buffer {index} emptied")
    _chain(session, args)


def cmd_swallow(session, args: ArgReader) -> None:
    """``SWALLOW`` -- read the macro's following lines as channel data.

    Consumes real-number tokens straight off the running ``XEQ`` file, one
    channel per token (or channel/value pairs with ``-twocolumn``), stopping
    at a blank line -- the same mechanism RUMP's own lexer uses to fill
    ``ibf->counts`` (bmanip.c: B1_SWALLOW, cmds.htm: "Swallow").
    """
    twocol = False
    token = args.optional()
    if token is not None:
        if token.lower() in ("-twocolumn", "-2column"):
            twocol = True
        elif token.lower() in ("-onecolumn", "-1column"):
            twocol = False
        else:
            raise CommandError(f"swallow: unrecognized option {token!r}")
    args.done()

    frame = session.xeq_frame()
    if frame is None:
        raise CommandError("swallow: only valid inside a command file (XEQ)")

    values: dict[int, float] = {}
    channel = 0
    while frame.index < len(frame.lines):
        line = frame.lines[frame.index].strip()
        frame.index += 1
        if not line:
            break
        parts = line.replace(",", " ").split()
        try:
            numbers = [float(part) for part in parts]
        except ValueError:
            frame.index -= 1  # not data after all -- leave it for the next command
            break
        if twocol:
            for column, value in zip(numbers[0::2], numbers[1::2]):
                values[max(0, int(column))] = value
        else:
            for value in numbers:
                values[channel] = value
                channel += 1
        if len(values) >= _CMAX:
            break

    npt = (max(values) + 1) if values else 0
    counts = np.zeros(npt)
    for channel, value in values.items():
        counts[channel] = value

    buffer = session.buffers.active_buffer
    if buffer is None:
        buffer = Buffer(spectrum=Spectrum.zeros(Calibration(npt=npt)))
        index = session.buffers.load(buffer)
        session.buffers.active = index
    else:
        index = session.buffers.active
    buffer.spectrum.counts = counts
    buffer.set_calibration(npt=npt)
    session.touch()
    print(f"{npt} points entered into buffer {index}")


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def cmd_plot(session, args: ArgReader) -> None:
    index, buffer = _buffer_argument(session, args)
    args.done()
    if index != 0:
        session.buffers.active = index
    plotting.add_trace(session, index, buffer, clear=True)
    plotting.draw(session)


def cmd_overlay(session, args: ArgReader) -> None:
    index, buffer = _buffer_argument(session, args)
    args.done()
    plotting.add_trace(session, index, buffer, clear=False)
    plotting.draw(session)


def cmd_replot(session, args: ArgReader) -> None:
    args.done()
    plotting.draw(session)


def cmd_compare(session, args: ArgReader) -> None:
    """Plot the active buffer against the simulation, with residuals."""
    args.done()
    from ...plot.spectra import plot_comparison

    plotting.require_matplotlib()
    data = session.buffers.require_active()
    theory = session.simulation()

    if session.figure is not None:
        import matplotlib.pyplot as plt

        plt.close(session.figure)
    figure = plot_comparison(
        data.spectrum, theory.spectrum, energy_axis=session.plot.energy_axis
    )
    session.figure = figure
    session.traces = []
    plotting.show(figure)


def cmd_display(session, args: ArgReader) -> None:
    """Plot the composition of the SIM sample against depth."""
    args.done()
    from ...plot.spectra import plot_depth_profile
    from ...script.lcm import to_sample
    from ...sim.engine import build_sample_grid

    if not session.script.layers:
        raise CommandError("no sample described: use SIM to build one")

    plt = plotting.require_matplotlib()
    sample = to_sample(session.script, session.table, session.densities)
    reference = session.buffers.require_active()
    grid = build_sample_grid(sample, reference.geometry, session.table)

    if session.figure is not None:
        plt.close(session.figure)
    figure, ax = plt.subplots(figsize=(9, 5.5))
    plot_depth_profile(grid, list(session.script.elements), ax=ax)
    ax.set_title(session.script.description or "Sample depth profile")
    session.figure = figure
    session.traces = []
    plotting.show(figure)


def cmd_region(session, args: ArgReader) -> None:
    if not args:
        print(session.plot.describe())
        return
    low = args.integer("the first channel")
    high = args.integer("the last channel")
    args.done()
    if high <= low:
        raise CommandError(f"empty region: {low} to {high}")
    session.plot.low, session.plot.high = low, high
    if session.traces:
        plotting.draw(session)


def cmd_expand(session, args: ArgReader) -> None:
    """Narrow the region to a subset of the current one and replot."""
    low = args.integer("the first channel")
    high = args.integer("the last channel")
    args.done()
    current_low = session.plot.low or 0
    current_high = session.plot.high
    if low < current_low or (current_high is not None and high > current_high):
        raise CommandError(
            f"EXPAND takes a subset of the current region "
            f"({current_low} to {current_high if current_high is not None else 'end'})"
        )
    session.plot.low, session.plot.high = low, high
    plotting.draw(session)


def cmd_counts(session, args: ArgReader) -> None:
    if not args:
        print(session.plot.describe())
        return
    low = args.number("the lowest yield")
    high = args.optional_number()
    args.done()
    if high is None:
        low, high = 0.0, low
    session.plot.ylow, session.plot.yhigh = low, high
    if session.traces:
        plotting.draw(session)


def cmd_blowup(session, args: ArgReader) -> None:
    """Expand the vertical scale: ``BLOWUP <max counts>``."""
    ceiling = args.number("the maximum yield")
    args.done()
    session.plot.ylow, session.plot.yhigh = 0.0, ceiling
    plotting.draw(session)


def _scale(name):
    def handler(session, args: ArgReader) -> None:
        args.done()
        session.plot.yscale = name
        if session.traces:
            plotting.draw(session)
        else:
            print(f"yield axis is {name}")

    return handler


def _flag(attribute, value, message):
    def handler(session, args: ArgReader) -> None:
        args.done()
        setattr(session.plot, attribute, value)
        if session.traces:
            plotting.draw(session)
        else:
            print(message)

    return handler


def cmd_labels(session, args: ArgReader) -> None:
    token = args.optional()
    args.done()
    session.plot.labels = token is None or token.lower() not in ("off", "no", "none")
    if session.traces:
        plotting.draw(session)
    print(f"labels {'on' if session.plot.labels else 'off'}")


def cmd_parms(session, args: ArgReader) -> None:
    args.done()
    print(session.plot.describe())
    # LexSystem's U_PARM case prints the working directory too (system.c:312).
    print(f"  Directory  {Path.cwd()}")
    if session.log_file is not None:
        print(f"  Script     {getattr(session.log_file, 'name', '(open)')}")


def cmd_axis(session, args: ArgReader) -> None:
    """Draw the axes with no data."""
    args.done()
    figure, ax = plotting.figure_for(session)
    ax.clear()
    session.traces = []
    ax.set_xlabel("Energy (keV)" if session.plot.energy_axis else "Channel")
    ax.set_ylabel("Counts")
    plotting.show(figure)


def cmd_energy(session, args: ArgReader) -> None:
    """Switch the x axis between channel (RUMP's default) and energy."""
    token = args.optional()
    args.done()
    session.plot.energy_axis = token is None or token.lower() not in ("off", "no")
    if session.traces:
        plotting.draw(session)
    print(f"x axis is {'energy' if session.plot.energy_axis else 'channel'}")


def cmd_integral(session, args: ArgReader) -> None:
    """Sum counts over a channel range in the active buffer."""
    low = args.integer("the first channel")
    high = args.integer("the last channel")
    args.done()
    buffer = session.buffers.require_active()
    counts = buffer.spectrum.counts
    if not 0 <= low <= high < counts.size:
        raise CommandError(f"channels {low}-{high} outside 0-{counts.size - 1}")
    window = counts[low : high + 1]
    total = float(window.sum())
    energies = buffer.spectrum.energies
    print(f"  channels {low}-{high} ({energies[low]:.1f}-{energies[high]:.1f} keV)")
    if total <= 0:
        print("  integral 0 counts")
        return
    centroid = float(np.dot(window, np.arange(low, high + 1)) / total)
    print(f"  integral {total:.1f} counts, centroid channel {centroid:.2f}")


# ---------------------------------------------------------------------------
# The table. Order matters: LexCmdl returns the first match (lexp2.c:639).
# ---------------------------------------------------------------------------

TABLE = CommandTable("Main Level Commands")

_ENTRIES: list[tuple[str, int, object, str]] = [
    # Session
    ("?", 1, cmd_help, "list the commands"),
    ("HELP", 4, cmd_help, "list the commands"),
    ("QUIT", 1, cmd_quit, "leave pyRUMP"),
    ("BYE", -2, cmd_quit, "leave pyRUMP"),
    ("DATA", 4, cmd_data, "show or change the atomic data directory"),
    # Sub-processors
    ("SIM", 3, cmd_sim, "enter the sample-description editor"),
    ("PERT", 3, cmd_pert, "enter the fitting sub-processor"),
    ("RETURN", 3, cmd_return, "return from a sub-level"),
    # Plotting
    ("PLOT", 2, cmd_plot, "erase and plot a buffer or file"),
    ("OVERLAY", 2, cmd_overlay, "overlay a buffer or file on the current plot"),
    ("REPLOT", 3, cmd_replot, "redraw the current plot"),
    ("COMPARE", 4, cmd_compare, "plot the active buffer against the simulation"),
    ("AXIS", 2, cmd_axis, "draw axes only"),
    ("BLOWUP", 2, cmd_blowup, "expand the vertical scale"),
    ("EXPAND", 2, cmd_expand, "narrow the region and replot"),
    # Plot scaling
    ("PARAMETERS", 4, cmd_parms, "display the plot parameters"),
    ("PARMS", -3, cmd_parms, "display the plot parameters"),
    ("REGION", 3, cmd_region, "set the channel range"),
    ("COUNTS", 2, cmd_counts, "set the yield range"),
    ("LINEAR", 2, _scale("linear"), "linear yield axis"),
    ("SQRT", 2, _scale("sqrt"), "square-root yield axis"),
    ("LOG", 2, _scale("log"), "logarithmic yield axis"),
    ("NORMALIZE", 2, _flag("normalized", True, "yield is normalized"),
     "plot yield in normalized units"),
    ("RAW", 2, _flag("normalized", False, "yield is raw"), "plot raw yield"),
    ("LABELS", 2, cmd_labels, "label the axes (LABELS OFF to suppress)"),
    ("ENERGY", 4, cmd_energy, "x axis in energy rather than channel"),
    # Buffers
    ("BUFFERS", 3, cmd_buffers, "display the buffer status"),
    ("GET", 3, cmd_get, "point at a buffer, or read a file into one"),
    ("READ", 2, cmd_get, "read a file into a buffer"),
    ("POINTAT", 2, cmd_pointat, "point at a buffer by number"),
    ("RELEASE", 7, cmd_release, "release the active buffer"),
    ("EMPTY", 3, cmd_empty, "reset a buffer to blank, or open a new one"),
    ("NEWALL", 3, cmd_newall, "clear every buffer"),
    ("COPY", 4, cmd_copy, "copy one buffer to another"),
    ("MOVE", 4, cmd_move, "exchange two buffers"),
    ("WRITE", 5, cmd_write, "write the active buffer to a .rbs file"),
    ("WRASCII", 3, cmd_wrascii, "write the active buffer as text"),
    ("RECALCULATE", 6, cmd_recalculate, "force the simulation to recompute"),
    # Parameters
    ("ACTIVE", 2, cmd_active, "list the active buffer's parameters"),
    ("BEAM", 4, cmd_beam, "incident beam species, e.g. 4He++"),
    ("MEV", 3, cmd_mev, "beam energy"),
    ("THETA", 3, _numeric("theta", "geometry", "theta", " deg"), "sample tilt"),
    ("PHI", 3, _numeric("phi", "geometry", "phi", " deg"),
     "supplement of the scattering angle"),
    ("PSI", 3, _numeric("psi", "geometry", "psi", " deg"), "exit angle"),
    ("GEOMETRY", 4, cmd_geometry, "cornell, ibm or general"),
    ("CONVERSION", 4, cmd_conversion, "keV per channel and offset"),
    ("CORRECTION", 3, _numeric("correction", "measurement", "corr"),
     "normalization fudge factor"),
    ("CHARGE", 2, _numeric("charge_uC", "measurement", "charge", " uC"), "beam dose"),
    ("CURRENT", 4, _numeric("current_nA", "measurement", "current", " nA"),
     "average beam current, for pileup"),
    ("CHOFF", 3, _numeric("first", "calibration", "first"),
     "channel number of the first data point"),
    ("FWHM", 4, _numeric("fwhm_keV", "measurement", "FWHM", " keV"),
     "detector resolution"),
    ("OMEGA", 5, _numeric("omega_msr", "measurement", "omega", " msr"),
     "detector solid angle"),
    ("TAU", 3, _numeric("tau_us", "measurement", "tau", " us"), "MCA shaping time"),
    ("IDENTIFIER", 3, cmd_identifier, "description of the spectrum"),
    ("DATE", 4, cmd_date, "when the spectrum was measured"),
    ("FILENAME", 4, cmd_filename, "record the buffer's source filename"),
    ("SWALLOW", 7, cmd_swallow, "read the following macro lines as channel data"),
    # Analysis
    ("INTEGRAL", 3, cmd_integral, "sum counts over a channel range"),
    ("DISPLAY", 3, cmd_display, "plot the sample composition against depth"),
]

for _name, _minlen, _handler, _help in _ENTRIES:
    TABLE.add(_name, _minlen, _handler, _help)
