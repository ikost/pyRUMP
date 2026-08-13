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


def _integral_report(session, result) -> None:
    label = "Interpolated" if session.integration_interp else "Discrete"
    print(f"  {label} integration on buffer {session.buffers.active}")
    print(
        f"  Region: {result.lo_channel:6.1f} to {result.hi_channel:6.1f}"
        f"  Gross: {result.gross:11.2f}  Net: {result.net:11.2f}  (#/uC/msr)"
    )


def cmd_integral(session, args: ArgReader) -> None:
    """``INTEGRAL lo hi`` -- gross/net counts over a channel range.

    ``RbsThickn``'s ``TH_INT`` (anlytc.c:349-351): reports gross and
    background-corrected net counts in normalized-yield units, honoring
    INTSET's rounding mode. Region bounds are plain 0-based channel indices,
    not RUMP's own ``first``-relative channel numbers.
    """
    from ...analysis.integration import integrate_region

    low = args.number("the first channel")
    high = args.number("the last channel")
    args.done()
    buffer = session.buffers.require_active()
    try:
        result = integrate_region(buffer, low, high, interp=session.integration_interp)
    except ValueError as error:
        raise CommandError(f"integral: {error}") from None
    _integral_report(session, result)


def cmd_thickness_analysis(session, args: ArgReader) -> None:
    """``THICKNESS lo hi element`` -- INTEGRAL plus a thickness conversion.

    ``RbsThickn``'s ``TH_THK`` (anlytc.c:353-356): surface-approximation
    thickness in atoms/cm^2 and Angstroms, plus (if ``INTSET`` is in
    ESTIMATED or QUERY mode) a second, energy-loss-ratio "compensated" pass.
    """
    from ...analysis.integration import integrate_region

    low = args.number("the first channel")
    high = args.number("the last channel")
    token = args.token("an element")
    alpha = args.number("a value for alpha") if session.integration_qmode == 2 else None
    args.done()

    buffer = session.buffers.require_active()
    try:
        result = integrate_region(
            buffer, low, high, interp=session.integration_interp,
            registry=session.registry, table=session.table, target_token=token,
            qmode=session.integration_qmode, alpha_override=alpha,
        )
    except (KeyError, ValueError) as error:
        raise CommandError(f"thickness: {error}") from None

    _integral_report(session, result)
    t = result.thickness
    print(f"  {token} surface approximation, density {t.density_g_cc:5.2f} g/cc")
    print(f"   (Gross) {t.gross_atoms_cm2:11.4e} Atoms/cm**2 ({t.gross_angstrom:7.1f} Angstroms)")
    print(f"   ( Net ) {t.net_atoms_cm2:11.4e} Atoms/cm**2 ({t.net_angstrom:7.1f} Angstroms)")
    if t.compensated is not None:
        c = t.compensated
        print("  Compensated calculation (Chu et al. page 65)")
        print(
            f"   (Gross) {c.gross_atoms_cm2:11.4e} Atoms/cm**2 ({c.gross_angstrom:7.1f} Angstroms)"
        )
        print(
            f"   ( Net ) {c.net_atoms_cm2:11.4e} Atoms/cm**2 ({c.net_angstrom:7.1f} Angstroms)"
        )


def cmd_intset(session, args: ArgReader) -> None:
    """``INTSET [Round|Interp|Surface|Estimated|Query|?]``.

    The two mode flags ``RbsThickn`` shares across INTEGRAL/THICKNESS
    (anlytc.c:1513-1541) -- rounding (discrete vs interpolated) and how
    THICKNESS's compensated pass picks its alpha.
    """
    token = args.optional()
    args.done()
    choice = (token or "?").strip().lower()
    if choice.startswith("i"):
        session.integration_interp = True
    elif choice.startswith("r"):
        session.integration_interp = False
    elif choice.startswith("s"):
        session.integration_qmode = 0
    elif choice.startswith("e"):
        session.integration_qmode = 1
    elif choice.startswith("q"):
        session.integration_qmode = 2
    elif choice == "?":
        interp, qmode = session.integration_interp, session.integration_qmode
        print(
            "  current mode - single letter command - explanation of mode\n"
            f"  {'*' if not interp else ' '} R - integrals are rounded to nearest channel\n"
            f"  {'*' if interp else ' '} I - integrals are interpolated between channels\n"
            f"  {'*' if qmode == 0 else ' '} S - thickness calculation based on surface "
            "approximation only\n"
            f"  {'*' if qmode == 1 else ' '} E - compensated thickness calculation with an "
            "estimated alpha\n"
            f"  {'*' if qmode == 2 else ' '} Q - compensated thickness calculation with query "
            "for alpha"
        )
    else:
        raise CommandError(f"intset: unrecognized mode {token!r}; use INTSET ? for help")


def cmd_element(session, args: ArgReader) -> None:
    """``ELEMENT el [el ...]`` -- expected energy/channel of each surface edge.

    ``RbsQueryElement``/``RbsKappa`` per element (anlytc.c:208-254); RUMP's
    optional cursor-driven height marker is dropped -- headless, this command
    is purely a report.
    """
    from ...analysis.elements import matrix_result

    tokens: list[str] = []
    while args:
        tokens.append(args.token())
    if not tokens:
        raise CommandError("element: expected one or more element names")

    buffer = session.buffers.require_active()
    for token in tokens:
        try:
            result = matrix_result(buffer, session.table, session.registry, token)
        except (KeyError, ValueError) as error:
            raise CommandError(f"element: {error}") from None
        if result.k == 0.0:
            print(f"  {result.symbol:2s}  Z={result.z:2d}  Mass={result.mass:7.3f}"
                  "  scattering event cannot occur")
            continue
        print(
            f"  {result.symbol:2s}  Z={result.z:2d}  Mass={result.mass:7.3f}"
            f"  K(ion)={result.k:6.4f}  Energy={result.energy_keV:8.1f} eV"
            f"  Channel={result.channel:8.3f}"
        )


def cmd_matrix(session, args: ArgReader) -> None:
    """``MATRIX el`` -- expected energy, channel and matrix height.

    ``RbsSigma``/``RbsEpsilon`` combined into a predicted step height
    (anlytc.c:258-275).
    """
    from ...analysis.elements import matrix_result

    token = args.token("an element")
    args.done()
    buffer = session.buffers.require_active()
    try:
        result = matrix_result(buffer, session.table, session.registry, token)
    except (KeyError, ValueError) as error:
        raise CommandError(f"matrix: {error}") from None
    if result.height is None:
        raise CommandError(f"matrix: scattering event cannot occur for {result.symbol}")
    print(
        f"  {result.symbol} expected at {result.energy_keV:8.1f} eV"
        f" ({result.channel:6.1f}) and height {result.height:8.3f}"
    )


def cmd_whatisit(session, args: ArgReader) -> None:
    """``WHATISIT <channel>`` -- identify elements near a channel.

    ``RbsLocate`` (anlytc.c:1398): the best-matching element by predicted
    surface-edge energy, plus its 2 neighbors on each side by Z.
    """
    from ...analysis.elements import locate_candidates

    channel = args.number("a channel number")
    args.done()
    buffer = session.buffers.require_active()
    target_keV = float(buffer.calibration.edge_energy(channel))
    try:
        candidates = locate_candidates(buffer, session.table, target_keV)
    except ValueError as error:
        raise CommandError(f"whatisit: {error}") from None
    print(f"  near channel {channel:g} ({target_keV:.1f} keV):")
    for candidate in candidates:
        print(
            f"    {candidate.symbol:2s} (Z={candidate.z:2d})"
            f"  {candidate.energy_keV:8.1f} eV  channel {candidate.channel:7.2f}"
        )


def cmd_info(session, args: ArgReader) -> None:
    """``INFO el`` -- density, K, cross section, stopping factors, isotopes.

    ``AN_INFO`` (anlytc.c:290-345), the fullest of the per-element reports.
    """
    from ...analysis.elements import matrix_result, resolve_element

    token = args.token("an element")
    args.done()
    buffer = session.buffers.require_active()
    try:
        element, mass, _ = resolve_element(session.table, token)
        result = matrix_result(buffer, session.table, session.registry, token)
    except (KeyError, ValueError) as error:
        raise CommandError(f"info: {error}") from None

    density_g_cc = element.atomic_density / 6.022e23 * mass
    print("-----------------------------------------------------")
    print(
        f"{element.symbol:2s}  Z: {element.z:2d}  Mass: {mass:6.2f}"
        f"  Density: {element.atomic_density:11.4e} at/cc ({density_g_cc:5.2f} g/cc)"
    )
    print(
        f"Parameters: Energy {buffer.beam.e0_MeV * 1000.0:8.1f} eV"
        f"  Theta{buffer.geometry.theta:6.2f}     Phi{buffer.geometry.phi:7.2f}"
    )
    c = buffer.calibration
    print(f"            keV/channel{c.kevch:6.3f}     keV(0){c.kev0:8.3f}")
    if result.k == 0.0:
        print("Scattering event cannot occur")
        return
    print(
        f"Surface Scattering:          {result.k:6.4f} at {result.energy_keV:8.1f} eV"
        f" (Channel: {result.channel:5.1f})"
    )
    if result.height is not None:
        print(f"Matrix scattering height:    {result.height:6.2f} Counts/uC/keV/msr")
        print(f"Scattering cross section:    {result.cross_section_barns:7.3f} (1E-24 cm2/ster)")
        # eps in "1e-15 eV-cm2" display units, then to eV/A -- anlytc.c:341-344.
        eps_display = result.epsilon_eVcm2 * 1e15
        eps_eVA = eps_display * element.atomic_density * 1e-23
        print(
            f"Stopping Factors:            [e] = {eps_display:5.1f} (1E-15 eV-cm2)"
            f"   [S] = {eps_eVA:5.1f} eV/A"
        )
    print()
    for isotope in element.isotopes:
        if isotope.mass <= 0:
            break
        print(f"  Mass: {isotope.mass:6.2f}  Abundance: {isotope.fraction:7.5f}")


def cmd_width_thick(session, args: ArgReader) -> None:
    """``WIDTH_THICK ch1 ch2 element`` -- thickness from a peak's half-height width.

    The energy-loss-ratio method (Chu et al. p.65), from two half-height
    channel positions instead of a whole region (anlytc.c:513-552).
    """
    from ...analysis.elements import cosines, cross_section_barns, kappa, resolve_element, stopper

    ch1 = args.number("the first half-height channel")
    ch2 = args.number("the second half-height channel")
    token = args.token("an element")
    args.done()

    buffer = session.buffers.require_active()
    e1_MeV = float(buffer.calibration.edge_energy(ch1)) / 1000.0
    e2_MeV = float(buffer.calibration.edge_energy(ch2)) / 1000.0
    width_MeV = abs(e2_MeV - e1_MeV)

    try:
        element, mass, _ = resolve_element(session.table, token)
    except (KeyError, ValueError) as error:
        raise CommandError(f"width_thick: {error}") from None

    z1, m1, e0 = buffer.beam.z, buffer.beam.mass, buffer.beam.e0_MeV
    angle = buffer.geometry.scattering_angle
    k = kappa(m1, mass, angle)
    if k == 0.0:
        raise CommandError(f"width_thick: scattering event cannot occur for {element.symbol}")
    cosin, cosout = cosines(buffer.geometry)
    registry = session.registry
    alpha = (
        stopper(registry, z1, m1, element.z, k * e0)
        / stopper(registry, z1, m1, element.z, e0)
        * cosin / cosout
    )

    e_in_mean = e0 - 0.5 * width_MeV / (k + alpha)
    e_out_mean = k * e0 - 0.5 * width_MeV * alpha / (k + alpha)
    mean_stopping = (
        k * stopper(registry, z1, m1, element.z, e_in_mean) / cosin
        + stopper(registry, z1, m1, element.z, e_out_mean) / cosout
    )
    mean_stopping *= 1e15  # "1e-15 eV-cm2" display units, anlytc.c:539

    density = element.atomic_density if element.atomic_density > 0 else 1.0e22
    stopping_eVA = mean_stopping * density * 1e-23

    width_eV = width_MeV * 1.0e6
    areal_density = width_eV / mean_stopping * 1e15
    thickness_A = width_eV / stopping_eVA

    sigma = cross_section_barns(z1, m1, element.z, mass, angle, e0 * 1000.0)
    # anlytc.c prints mean_stopping under this header -- a mislabel in the
    # original (it is epsilon, not a cross section); reproduced as-is.
    print(f"  Scattering cross section: {mean_stopping:7.3f}  (sigma = {sigma:.3f} barns/sr)")
    print(f"  Width: {width_eV:9.1f} eV")
    print(f"  Areal density: {areal_density:11.4e} Atoms/cm**2")
    print(f"  Thickness: {thickness_A:9.1f} Angstroms")


def cmd_calibrate(session, args: ArgReader) -> None:
    """``CALIBRATE ch1 el1 ch2 el2 [energy_eV marker_channel]``.

    Two-point energy calibration (anlytc.c:564-612): two (channel, element)
    pairs fix the K-vs-channel line; an optional third, precisely-known
    marker energy also fixes the absolute beam energy.
    """
    from ...analysis.elements import kappa, resolve_element

    ch1 = args.number("the first peak's channel")
    el1 = args.token("the first element")
    ch2 = args.number("the second peak's channel")
    el2 = args.token("the second element")
    marker_energy = args.optional_number()
    marker_channel = args.number("the marker's channel") if marker_energy is not None else None
    args.done()

    if abs(ch2 - ch1) < 2:
        raise CommandError("calibrate: needs two DIFFERENT channel positions")

    buffer = session.buffers.require_active()
    try:
        _, mass1, _ = resolve_element(session.table, el1)
        _, mass2, _ = resolve_element(session.table, el2)
    except (KeyError, ValueError) as error:
        raise CommandError(f"calibrate: {error}") from None

    angle = buffer.geometry.scattering_angle
    kh1 = kappa(buffer.beam.mass, mass1, angle)
    kh2 = kappa(buffer.beam.mass, mass2, angle)
    if abs(kh1 - kh2) < 0.001:
        raise CommandError("calibrate: needs two DIFFERENT elements")

    slope = (kh1 - kh2) / (ch1 - ch2)
    e0 = buffer.beam.e0_MeV
    if marker_energy is not None:
        e0 = 0.001 * marker_energy / (kh2 + slope * (marker_channel - ch2))

    kevch = 1000.0 * slope * e0
    kev0 = 1000.0 * kh1 * e0 - ch1 * kevch

    buffer.beam.e0_MeV = e0
    buffer.set_calibration(kevch=kevch, kev0=kev0)
    session.touch()
    print(f"  Energy={e0:.4f} MeV    Conversion:{kevch:.4f} keV/ch   {kev0:.4f} keV(0)")


def cmd_background(session, args: ArgReader) -> None:
    """``BACKGROUND lo1 hi1 lo2 hi2 order [-inplace] [-noplot]``.

    ``RbsBackground`` (anlytc.c:1011): weighted polynomial fit of the two
    flanking regions, subtracted across the whole span including the peak.
    Defaults to a new, cropped buffer (RUMP's ``.cut`` convention);
    ``-inplace`` subtracts on the active buffer instead.
    """
    from ...analysis.background import fit_background

    i0 = args.integer("the start of the lower flanking region")
    i1 = args.integer("the end of the lower flanking region")
    i2 = args.integer("the start of the upper flanking region")
    i3 = args.integer("the end of the upper flanking region")
    order = args.integer("the polynomial order")
    inplace = False
    noplot = False
    while args:
        flag = args.token().lower()
        if flag == "-inplace":
            inplace = True
        elif flag == "-noplot":
            noplot = True
        else:
            raise CommandError(f"background: unrecognized option {flag!r}")

    buffer = session.buffers.require_active()
    try:
        fit = fit_background(buffer.spectrum.counts, i0, i1, i2, i3, order)
    except ValueError as error:
        raise CommandError(f"background: {error}") from None

    if inplace:
        buffer.spectrum.counts[i0 : i3 + 1] = fit.stripped
        session.touch()
        index = session.buffers.active
        print(f"  background subtracted in place on buffer {index}")
    else:
        new_calibration = Calibration(
            kevch=buffer.calibration.kevch,
            kev0=buffer.calibration.kev0,
            first=buffer.calibration.first + i0,
            npt=fit.stripped.size,
        )
        new_buffer = Buffer(
            spectrum=Spectrum(counts=fit.stripped.copy(), calibration=new_calibration),
            beam=buffer.beam, geometry=buffer.geometry, measurement=buffer.measurement,
            name=(Path(buffer.name).stem + ".cut") if buffer.name else "background.cut",
            identifier=buffer.identifier,
        )
        index = session.buffers.first_free()
        session.buffers.set(index, new_buffer)
        session.buffers.active = index
        session.touch()
        print(f"  background subtracted; result in buffer {index} ({fit.stripped.size} channels)")

    if not noplot:
        plt = plotting.require_matplotlib()
        if session.figure is not None:
            plt.close(session.figure)
        figure, ax = plt.subplots(figsize=(9, 5.5))
        ax.plot(fit.channels, buffer.spectrum.counts[i0 : i3 + 1], color="0.5", lw=1.0,
                label="data")
        ax.plot(fit.channels, fit.fit, lw=1.5, label="fit")
        ax.plot(fit.channels, fit.stripped, lw=1.0, label="background-subtracted")
        ax.set_xlabel("Channel")
        ax.set_ylabel("Counts")
        ax.legend(frameon=False, fontsize="small")
        session.figure = figure
        session.traces = []
        plotting.show(figure)


def cmd_smooth(session, args: ArgReader) -> None:
    """``SMOOTH [-sv|-conv|-fft] [-range lo hi] [n]``.

    ``RbsSmooth_SV``/``_Conv``/``_FFT`` (anlytc.c:454-509). Defaults to
    ``-sv`` over the whole buffer; ``n`` sets ``-conv``'s iteration count
    (default 2), or ``-fft``'s smoothing width (required).
    """
    from ...analysis import smoothing

    mode = "sv"
    low: int | None = None
    high: int | None = None
    n_iterations = 2
    width: float | None = None

    while args:
        token = args.token()
        lowered = token.lower()
        if lowered == "-sv":
            mode = "sv"
        elif lowered in ("-conv", "-convolution", "-convolute"):
            mode = "conv"
        elif lowered == "-fft":
            mode = "fft"
        elif lowered == "-range":
            low = args.integer("the first channel")
            high = args.integer("the last channel")
        else:
            try:
                value = float(token)
            except ValueError:
                raise CommandError(f"smooth: unrecognized option {token!r}") from None
            if mode == "fft":
                width = value
            else:
                n_iterations = int(value)

    buffer = session.buffers.require_active()
    counts = buffer.spectrum.counts
    if low is None or high is None:
        low, high = 0, counts.size - 1

    try:
        if mode == "sv":
            new_counts = smoothing.smooth_sv(counts, low, high)
        elif mode == "conv":
            new_counts, rms_history = smoothing.smooth_conv(
                counts, low, high, buffer.measurement.fwhm_keV, buffer.calibration.kevch,
                n_iterations,
            )
        else:
            if width is None:
                raise CommandError("smooth -fft: needs a smoothing width")
            new_counts = smoothing.smooth_fft(counts, low, high, width)
    except ValueError as error:
        raise CommandError(f"smooth: {error}") from None

    buffer.spectrum.counts = new_counts
    session.touch()
    if mode == "conv":
        print("  RMS: " + " ".join(f"{value:.4f}" for value in rms_history))
    print(f"  buffer {session.buffers.active} smoothed ({mode}), channels {low}-{high}")


def cmd_fft(session, args: ArgReader) -> None:
    """``FFT lo hi width`` -- same as ``SMOOTH -fft -range lo hi width``
    (anlytc.c:454-456: the command is literally an alias in the original).
    """
    from ...analysis import smoothing

    low = args.integer("the first channel")
    high = args.integer("the last channel")
    width = args.number("the smoothing width")
    args.done()
    buffer = session.buffers.require_active()
    try:
        new_counts = smoothing.smooth_fft(buffer.spectrum.counts, low, high, width)
    except ValueError as error:
        raise CommandError(f"fft: {error}") from None
    buffer.spectrum.counts = new_counts
    session.touch()
    print(f"  buffer {session.buffers.active} smoothed (fft), channels {low}-{high}")


def cmd_profile(session, args: ArgReader) -> None:
    """``PROFILE`` -- not implemented, never was.

    ``RbsNewprf`` (null.c): dead code even in the original C, which prints
    this exact message and does nothing. Reproduced verbatim, the same
    precedent as THICKFILM's author confession in profiles/equations.py.
    """
    args.done()
    print("OOPS: Didn't think anyone used this routine anymore - sorry not implemented")


def cmd_cursor(session, args: ArgReader) -> None:
    """``CURSOR`` -- not available in this shell.

    RUMP's own behavior with no interactive graphics device (anlytc.c:187):
    print this message and do nothing. There is no keyboard fallback for
    CURSOR specifically (unlike every other command in this family).
    """
    args.done()
    print("Cursor not enabled or illegal device")


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
    # Analysis (anlytc.c's own cmlist order; abbreviation lengths from there
    # too, except DISPLAY -- kept at its already-shipped minlen 3 ("DIS"),
    # not the C's 4, to avoid changing already-tested behavior)
    ("CURSOR", 3, cmd_cursor, "graphics cursor (not available in this shell)"),
    ("ELEMENT", 2, cmd_element, "expected energy/channel of an element's surface peak"),
    ("MATRIX", 3, cmd_matrix, "expected energy, channel and matrix height"),
    ("WHATISIT", 4, cmd_whatisit, "identify elements near a channel"),
    ("INFO", 3, cmd_info, "detailed report on an element"),
    ("INTEGRAL", 3, cmd_integral, "sum counts over a channel range"),
    ("THICKNESS", 4, cmd_thickness_analysis, "integral plus thickness conversion"),
    ("BACKGROUND", 4, cmd_background, "fit and subtract a polynomial background"),
    ("SMOOTH", 3, cmd_smooth, "smooth the active buffer (-sv, -conv, -fft)"),
    ("WIDTH_THICK", 3, cmd_width_thick, "thickness from a peak's half-height width"),
    ("PROFILE", 3, cmd_profile, "not implemented -- never was, in the original"),
    ("INTSET", 6, cmd_intset, "change INTEGRAL/THICKNESS rounding and alpha mode"),
    ("CALIBRATE", 3, cmd_calibrate, "energy-calibrate from two known peaks"),
    ("DISPLAY", 3, cmd_display, "plot the sample composition against depth"),
    ("FFT", 3, cmd_fft, "FFT smooth (same as SMOOTH -FFT -RANGE)"),
]

for _name, _minlen, _handler, _help in _ENTRIES:
    TABLE.add(_name, _minlen, _handler, _help)
