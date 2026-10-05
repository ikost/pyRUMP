"""PIXE: the X-ray sub-processor.

Entering ``PIXE`` turns PIXE on: from then on the PIXE window follows the
RBS one, showing the ACTIVE buffer's PIXE spectrum. ``DISABLE`` turns it off
and leaves the prompt in one step; ``RETURN`` leaves PIXE on. One-shot
``PIXE <command>`` from the RUMP level runs a single command without
entering, and does not turn PIXE on -- so ``~/.pyrumprc`` can set up the
detector with one-shots, or with a ``pixe`` ... ``disable`` block.

A PIXE spectrum lives on a buffer (:attr:`~pyrump.shell.session.Buffer.pixe`)
next to that buffer's RBS spectrum, and shares its beam, geometry and
charge. The detector, calibration default and plot settings are the
session's (:class:`~pyrump.shell.session.PixeState`).
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from ...model.spectrum import Spectrum
from ...pixe.data import PixeData
from ...pixe.detector import COMPOUNDS, Absorber, disc_solid_angle_msr
from ..dispatch import ArgReader, CommandError, CommandTable
from ..session import Buffer
from .. import pixe_plotting
from .rump import FIGSAVE_DPI, EnterMode, Return, _input_path, describe_topic

#: Extensions a bare ``GET name`` tries, in order.
PIXE_SUFFIXES = (".PIX", ".pix")


# ---------------------------------------------------------------------------
# Entering and leaving
# ---------------------------------------------------------------------------


def enter(session, args: ArgReader) -> None:
    """``PIXE [command]`` at the RUMP level: enter the PIXE prompt (turning
    PIXE on), or run one PIXE command without entering it."""
    if args:
        execute_in_pixe(session, args)
        return
    state = session.pixe
    if not state.enabled:
        state.enabled = True
        print("PIXE enabled (DISABLE turns it off)")
        pixe_plotting.draw(session, required=False)
    raise EnterMode("pixe")


def execute_in_pixe(session, args: ArgReader) -> None:
    """Run a one-shot ``PIXE <command>`` from the RUMP level."""
    name = args.token("a PIXE command")
    command = TABLE.match(name)
    if command is None:
        raise CommandError(f"unrecognized PIXE command: {name}")
    try:
        command.handler(session, ArgReader(args.remaining, command=command.name.lower()))
    except Return:
        pass  # RETURN/DISABLE one-shots: there is no PIXE level to leave
    args.index = len(args.tokens)


def cmd_return(session, args: ArgReader) -> None:
    """``RETURN`` -- back to the RUMP level, leaving PIXE on."""
    raise Return()


def cmd_disable(session, args: ArgReader) -> None:
    """``DISABLE`` -- turn PIXE off and return to the RUMP level.

    The PIXE window closes; the settings and the buffers' PIXE spectra are
    kept, and entering ``PIXE`` again turns it back on."""
    args.done()
    if session.pixe.enabled:
        session.pixe.enabled = False
        print("PIXE disabled")
    pixe_plotting.close(session)
    raise Return()


def cmd_help(session, args: ArgReader) -> None:
    """``HELP`` lists the PIXE commands; ``HELP <name>`` describes one."""
    topic = args.optional()
    args.done()
    if topic is None:
        leftover = TABLE.uncovered(_HELP_GROUPS)
        groups = [*_HELP_GROUPS, ("Other", leftover)] if leftover else _HELP_GROUPS
        print(TABLE.grouped_help_text(groups))
        return
    from .rump import TABLE as RUMP_TABLE
    from .system import TABLE as SYSTEM_TABLE

    print(describe_topic(session, topic, (TABLE, RUMP_TABLE, SYSTEM_TABLE)))


def setup_lines(session) -> list[str]:
    """The settings as PIXE commands -- what SHOW prints, ready for
    ``~/.pyrumprc``."""
    state = session.pixe
    d, c = state.detector, state.calibration
    lines = [
        f"angle {d.angle_deg:g}",
        f"tiltsign {d.tilt_sign:d}",
        f"solid {d.solid_angle_msr:.6g}",
        f"window {d.window.material} {d.window.thickness_um:g}",
        f"crystal {d.crystal.material} {d.crystal.thickness_um:g}",
        f"fwhm {d.fwhm_eV:g}",
        f"fano {d.fano:g}",
        # Replayed, FILTER adds to what is there: start from none.
        "filter clear",
    ]
    lines += [
        f"filter {f.material} {f.thickness_um:g} {f.hole_percent:g}" for f in d.filters
    ]
    lines += [
        f"calib {c.kevch:.8g} {c.kev0:.8g}",
        f"h {state.h[0]:g} {state.h[1]:g} {state.h[2]:g}",
        f"escape {'on' if state.escape else 'off'}",
        f"pair {'on' if state.pair else 'off'}",
        f"markers {state.markers}",
    ]
    return lines


def cmd_show(session, args: ArgReader) -> None:
    """``SHOW`` -- the PIXE settings, as commands for ``~/.pyrumprc``, and
    the active buffer's PIXE spectrum."""
    args.done()
    state = session.pixe
    print(f"  ! PIXE {'enabled' if state.enabled else 'disabled'}")
    for line in setup_lines(session):
        print(f"  {line}")
    from ..pixe_sim import reference_buffer

    reference = reference_buffer(session)
    theta = reference.geometry.theta
    print(
        f"  ! with THETA {theta:g}: beam {abs(theta):g} deg, X-rays "
        f"{state.detector.exit_angle(theta):g} deg to the sample normal"
    )
    buffer = session.buffers.active_buffer
    if buffer is not None and buffer.pixe is not None:
        print(f"  ! buffer {session.buffers.active}:")
        print(f"  !{buffer.pixe.describe()[1:]}")
    else:
        print("  ! no PIXE spectrum in the active buffer")


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------


def cmd_get(session, args: ArgReader) -> None:
    """``GET <file>`` -- read a PIXE spectrum into the active buffer.

    With no data buffer active, the spectrum gets a buffer of its own, whose
    beam, geometry and charge are the defaults (``MEV``, ``BEAM``, ``CHARGE``
    ... as ``~/.pyrumprc`` sets them)."""
    path = _pixe_path(Path(args.token("a PIXE spectrum file")))
    args.done()
    index = session.buffers.active
    buffer = session.buffers.get(index) if index else None
    if buffer is None:
        buffer = _pixe_only_buffer(session, path)
        index = session.buffers.scroll_in(buffer)
        session.buffers.active = index
    attach(session, buffer, path)
    print(f"  {path.name} -> PIXE of buffer {index}")
    if session.pixe.enabled:
        pixe_plotting.draw(session, required=False)


def attach(session, buffer: Buffer, path: Path) -> None:
    """Read ``path`` as ``buffer``'s PIXE spectrum, replacing any before."""
    try:
        data, notices = PixeData.read(path, session.pixe.calibration)
    except (ValueError, OSError) as error:
        raise CommandError(f"could not read {path}: {error}") from None
    for notice in notices:
        print(f"  WARNING: {notice}")
    buffer.pixe = data


def pair(session, buffer: Buffer | None, rbs_path: Path) -> None:
    """With ``PAIR ON``, give a freshly read RBS buffer the ``.PIX`` file of
    the same name from the same folder. Called by GET and XEQ."""
    if not session.pixe.pair or buffer is None or buffer.pixe is not None:
        return
    companion = _companion(Path(rbs_path))
    if companion is None:
        print(f"  PAIR: no {Path(rbs_path).stem}.PIX next to {Path(rbs_path).name}")
        return
    attach(session, buffer, companion)
    print(f"  PAIR: {companion.name} -> PIXE of this buffer")
    if session.pixe.enabled:
        pixe_plotting.draw(session, required=False)


def _companion(rbs_path: Path) -> Path | None:
    """``x.PIX`` (any case) next to ``x.RBS``."""
    folder = rbs_path.parent if str(rbs_path.parent) else Path(".")
    try:
        entries = sorted(folder.iterdir())
    except OSError:
        return None
    stem = rbs_path.stem.lower()
    for entry in entries:
        if entry.is_file() and entry.stem.lower() == stem and entry.suffix.lower() == ".pix":
            return entry.resolve()
    return None


def _pixe_path(path: Path) -> Path:
    """The PIXE file, trying ``.PIX`` and ``.pix`` for a bare name."""
    path = path.expanduser()
    if not path.exists() and not path.suffix:
        for suffix in PIXE_SUFFIXES:
            if path.with_suffix(suffix).exists():
                return path.with_suffix(suffix).resolve()
    return _input_path(path, PIXE_SUFFIXES[0])


def _pixe_only_buffer(session, path: Path) -> Buffer:
    """A buffer with no RBS data, its parameters from the session defaults."""
    defaults = session.settings.experiment_defaults
    return Buffer(
        spectrum=Spectrum.zeros(defaults.calibration),
        beam=replace(defaults.beam),
        geometry=defaults.geometry,
        measurement=defaults.measurement,
        name=path.name,
        identifier=f"{path.stem} (PIXE only)",
    )


def cmd_pair(session, args: ArgReader) -> None:
    """``PAIR ON|OFF`` -- whether reading ``x.RBS`` (GET or XEQ) also loads
    ``x.PIX`` from the same folder."""
    token = args.optional()
    args.done()
    if token is not None:
        if token.lower() not in ("on", "off"):
            raise CommandError("PAIR: expected ON or OFF")
        session.pixe.pair = token.lower() == "on"
    print(f"  pair {'on' if session.pixe.pair else 'off'}")


# ---------------------------------------------------------------------------
# Detector and calibration
# ---------------------------------------------------------------------------


def _detector_number(field: str, label: str, unit: str, *, low: float = 0.0,
                     high: float = float("inf")):
    """A setter for a number that must lie in ``(low, high)``."""
    def handler(session, args: ArgReader) -> None:
        value = args.optional_number()
        args.done()
        detector = session.pixe.detector
        if value is not None:
            if not low < value < high:
                bound = f"between {low:g} and {high:g}" if high < float("inf") else "positive"
                raise CommandError(f"{label.upper()}: must be {bound}")
            detector = replace(detector, **{field: value})
            session.pixe.detector = detector
            pixe_plotting.refresh(session)
        print(f"  {label} {getattr(detector, field):g}{unit}")

    handler.__doc__ = f"``{label.upper()} [<value>]`` -- show or set the detector's {label}."
    return handler


cmd_angle = _detector_number(
    "angle_deg", "angle", "  ! deg, detector axis to the untilted sample's normal",
    low=-90.0, high=90.0,
)


def cmd_tiltsign(session, args: ArgReader) -> None:
    """``TILTSIGN 1|-1|0`` -- how the sample's tilt THETA moves the X-rays'
    exit angle, ``|ANGLE + TILTSIGN x THETA|``: 1 when a negative THETA
    turns the sample towards the PIXE detector, -1 when away, 0 when the
    tilt leaves the detector direction alone."""
    token = args.optional()
    args.done()
    detector = session.pixe.detector
    if token is not None:
        if token not in ("1", "-1", "0", "+1"):
            raise CommandError("TILTSIGN: expected 1, -1 or 0")
        detector = replace(detector, tilt_sign=int(token))
        session.pixe.detector = detector
        pixe_plotting.refresh(session)
    print(f"  tiltsign {detector.tilt_sign:d}")


def cmd_solid(session, args: ArgReader) -> None:
    """``SOLID <msr>``, or ``SOLID <area mm^2> <distance> [MM|IN]`` -- the
    detector's solid angle, given directly or from its active area and its
    distance to the sample (a round detector seen on axis)."""
    detector = session.pixe.detector
    if args:
        first = args.number("a solid angle in msr, or an area in mm^2")
        second = args.optional_number()
        unit = (args.optional() or "mm").lower()
        args.done()
        if unit not in ("mm", "in"):
            raise CommandError("SOLID: the distance unit is MM or IN")
        if second is None:
            value, source = first, ""
        else:
            distance = second * (25.4 if unit == "in" else 1.0)
            if distance <= 0:
                raise CommandError("SOLID: the distance must be positive")
            value = disc_solid_angle_msr(first, distance)
            source = f" ({first:g} mm^2 at {distance:g} mm)"
        if value <= 0:
            raise CommandError("SOLID: must be positive")
        detector = replace(detector, solid_angle_msr=value)
        session.pixe.detector = detector
        pixe_plotting.refresh(session)
        print(f"  solid {value:.6g}  ! msr{source}")
        return
    print(f"  solid {detector.solid_angle_msr:.6g}  ! msr")
cmd_fwhm = _detector_number("fwhm_eV", "fwhm", "  ! eV at Mn Ka")
cmd_fano = _detector_number("fano", "fano", "")


def _absorber(session, args: ArgReader, what: str, *, hole: bool = False,
              compounds: bool = True) -> Absorber:
    name = args.token("an element symbol" + (" or MYLAR/KAPTON" if compounds else ""))
    compound = COMPOUNDS.get(name.upper()) if compounds else None
    if compound is not None:
        material = compound.name
    else:
        try:
            material = session.table.by_symbol(name).symbol
        except KeyError as error:
            known = f" (or {', '.join(c.name for c in COMPOUNDS.values())})" if compounds else ""
            raise CommandError(f"{what}: {str(error).strip(chr(39))}{known}") from None
    thickness = args.number("a thickness in µm")
    if thickness <= 0:
        raise CommandError(f"{what}: thickness must be positive")
    open_area = args.optional_number() if hole else None
    if open_area is not None and not 0 <= open_area < 100:
        raise CommandError(f"{what}: hole area must be 0 to 100 %")
    args.done()
    return Absorber(material, thickness, open_area or 0.0)


def _absorber_command(field: str, label: str, *, compounds: bool):
    def handler(session, args: ArgReader) -> None:
        if args:
            absorber = _absorber(session, args, label.upper(), compounds=compounds)
            session.pixe.detector = replace(session.pixe.detector, **{field: absorber})
            pixe_plotting.refresh(session)
        absorber = getattr(session.pixe.detector, field)
        print(f"  {label} {absorber.material} {absorber.thickness_um:g}  ! µm")

    what = "<element|MYLAR|KAPTON>" if compounds else "<element>"
    handler.__doc__ = f"``{label.upper()} [{what} <µm>]`` -- show or set the detector's {label}."
    return handler


cmd_window = _absorber_command("window", "window", compounds=True)
cmd_crystal = _absorber_command("crystal", "crystal", compounds=False)


def cmd_filter(session, args: ArgReader) -> None:
    """``FILTER <element|MYLAR|KAPTON> <µm> [<hole %>]`` adds an absorber in
    front of the window; ``FILTER CLEAR`` removes them all."""
    if args and args.peek().lower() == "clear":
        args.token()
        args.done()
        session.pixe.detector = replace(session.pixe.detector, filters=())
        pixe_plotting.refresh(session)
    elif args:
        absorber = _absorber(session, args, "FILTER", hole=True)
        detector = session.pixe.detector
        session.pixe.detector = replace(detector, filters=(*detector.filters, absorber))
        pixe_plotting.refresh(session)
    filters = session.pixe.detector.filters
    if not filters:
        print("  no filters")
    for f in filters:
        print(f"  filter {f.material} {f.thickness_um:g} {f.hole_percent:g}  ! µm, hole %")


def cmd_calib(session, args: ArgReader) -> None:
    """``CALIB [<keV/ch> <offset keV>]`` -- the PIXE energy calibration.

    Sets the default for spectra whose header has none, and the active
    buffer's PIXE spectrum if it has one."""
    if args:
        gain = args.number("a gain in keV/channel")
        offset = args.number("an offset in keV")
        args.done()
        if gain <= 0:
            raise CommandError("CALIB: the gain must be positive")
        state = session.pixe
        state.calibration = replace(state.calibration, kevch=gain, kev0=offset)
        buffer = session.buffers.active_buffer
        if buffer is not None and buffer.pixe is not None:
            spectrum = buffer.pixe.spectrum
            spectrum.calibration = replace(spectrum.calibration, kevch=gain, kev0=offset)
            if state.enabled:
                pixe_plotting.draw(session, required=False)
    buffer = session.buffers.active_buffer
    calibration = (
        buffer.pixe.calibration
        if buffer is not None and buffer.pixe is not None
        else session.pixe.calibration
    )
    print(f"  calib {calibration.kevch:.8g} {calibration.kev0:.8g}  ! keV/ch, keV")


def cmd_h(session, args: ArgReader) -> None:
    """``H <K> <L> <M>`` or ``H K|L|M <value>`` -- the instrumental constant
    for each shell's lines: measured yield = H x calculated yield. It takes
    up the solid angle, charge and database errors; set it from standards."""
    state = session.pixe
    if args:
        first = args.token("H values, or K, L or M")
        if first.upper() in ("K", "L", "M"):
            value = args.number("an H value")
            values = list(state.h)
            values["KLM".index(first.upper())] = value
        else:
            try:
                values = [float(first), args.number("H for L lines"), args.number("H for M lines")]
            except ValueError:
                raise CommandError("H: expected K L M values, or K|L|M <value>") from None
        args.done()
        if any(v <= 0 for v in values):
            raise CommandError("H: values must be positive")
        state.h = tuple(values)
        pixe_plotting.refresh(session)
    print(f"  h {state.h[0]:g} {state.h[1]:g} {state.h[2]:g}  ! K, L, M")


def cmd_escape(session, args: ArgReader) -> None:
    """``ESCAPE ON|OFF`` -- whether the simulation includes Si escape peaks."""
    token = args.optional()
    args.done()
    if token is not None:
        if token.lower() not in ("on", "off"):
            raise CommandError("ESCAPE: expected ON or OFF")
        session.pixe.escape = token.lower() == "on"
        pixe_plotting.refresh(session)
    print(f"  escape {'on' if session.pixe.escape else 'off'}")


def cmd_lines(session, args: ArgReader) -> None:
    """``LINES [ALL]`` -- the simulated lines: energy, cross section at the
    beam energy, detector efficiency and counts, strongest first. Lines
    under 0.1 % of the strongest are left out unless ALL is given."""
    from ..pixe_sim import simulate

    token = args.optional()
    args.done()
    show_all = token is not None and token.lower() == "all"
    try:
        result = simulate(session)
    except (ValueError, KeyError) as error:
        raise CommandError(f"LINES: {str(error).strip(chr(39))}") from None
    if result is None:
        raise CommandError("LINES: no SIM sample to simulate")
    if not result.lines:
        print("  no film lines: the sample is a substrate only")
        return
    strongest = result.lines[0].counts
    print("  element  line    E (keV)   sigma (b)   efficiency      counts")
    for line in result.lines:
        if not show_all and line.counts < 1e-3 * strongest:
            continue
        print(
            f"  {line.symbol:<7}  {line.line.line:<6}  {line.energy_keV:7.4f}"
            f"   {line.sigma_barn:9.4g}   {line.efficiency:10.4f}  {line.counts:10.4g}"
        )
    totals: dict[str, float] = {}
    for line in result.lines:
        key = f"{line.symbol} {line.family}"
        totals[key] = totals.get(key, 0.0) + line.counts
    print("  totals:  " + ", ".join(f"{k} {v:.4g}" for k, v in totals.items()))


# ---------------------------------------------------------------------------
# The PIXE window
# ---------------------------------------------------------------------------


def _redraw_or(session, message: str) -> None:
    if not pixe_plotting.draw(session, required=False):
        print(message)


def cmd_plot(session, args: ArgReader) -> None:
    """``PLOT`` -- draw the active buffer's PIXE spectrum."""
    args.done()
    pixe_plotting.draw(session)


def cmd_region(session, args: ArgReader) -> None:
    """``REGION <keV> <keV>`` -- the energy range shown; ``REGION ALL`` the
    whole spectrum; no argument shows it."""
    state = session.pixe
    if args and args.peek().lower() == "all":
        args.token()
        state.emin = state.emax = None
    elif args:
        low = args.number("the lowest energy, keV")
        high = args.number("the highest energy, keV")
        if low == high:
            raise CommandError(f"empty region: {low} to {high} keV")
        state.emin, state.emax = min(low, high), max(low, high)
    args.done()
    span = (
        "all"
        if state.emin is None and state.emax is None
        else f"{state.emin:g} to {state.emax:g} keV"
    )
    _redraw_or(session, f"  region {span}")


def cmd_counts(session, args: ArgReader) -> None:
    """``COUNTS <low> <high>`` -- the yield range; ``COUNTS ALL`` autoscales."""
    plot = session.pixe.plot
    if args and args.peek().lower() == "all":
        args.token()
        plot.ylow = plot.yhigh = None
    elif args:
        low = args.number("the lowest count")
        high = args.number("the highest count")
        plot.ylow, plot.yhigh = min(low, high), max(low, high)
    args.done()
    span = "auto" if plot.ylow is None else f"{plot.ylow:g} to {plot.yhigh:g}"
    _redraw_or(session, f"  counts {span}")


def _scale(name: str):
    def handler(session, args: ArgReader) -> None:
        args.done()
        session.pixe.plot.yscale = name
        _redraw_or(session, f"  PIXE yield axis is {name}")

    handler.__doc__ = f"``{name.upper()}`` -- {name} yield axis in the PIXE window."
    return handler


def cmd_markers(session, args: ArgReader) -> None:
    """``MARKERS ON|ALL|OFF`` -- label the lines of the SIM sample's
    elements: ON the main ones (Kα, Kβ, Lα, Lβ1, Lβ2, Lγ1, Mα, Mβ), ALL
    also Ll, Mζ and Mγ."""
    token = args.optional()
    args.done()
    if token is not None:
        if token.lower() not in ("on", "all", "off"):
            raise CommandError("MARKERS: expected ON, ALL or OFF")
        session.pixe.markers = token.lower()
    _redraw_or(session, f"  markers {session.pixe.markers}")


def cmd_figsave(session, args: ArgReader) -> None:
    """``FIGSAVE <file>`` -- save the PIXE window to an image file."""
    path = Path(args.token("an output image file"))
    args.done()
    figure = session.pixe.figure
    if figure is None:
        raise CommandError("no PIXE plot yet -- PLOT first")
    if not path.suffix:
        path = path.with_suffix(".png")
    figure.savefig(path, dpi=FIGSAVE_DPI)
    print(f"wrote {path}")


# ---------------------------------------------------------------------------
# The table
# ---------------------------------------------------------------------------

TABLE = CommandTable("PIXE Commands")

_ENTRIES: list[tuple[str, int, object, str]] = [
    ("?", -1, cmd_help, "synonym for HELP"),
    ("HELP", 2, cmd_help, "list the PIXE commands"),
    ("RETURN", 3, cmd_return, "return to the RUMP level, PIXE stays on"),
    ("QUIT", -1, cmd_return, "synonym for RETURN (not exit pyRUMP)"),
    ("Q", -1, cmd_return, "synonym for RETURN"),
    ("DISABLE", 3, cmd_disable, "turn PIXE off and return to the RUMP level"),
    ("SHOW", 2, cmd_show, "the settings, as commands for ~/.pyrumprc"),
    # Data
    ("GET", 3, cmd_get, "read a PIXE spectrum (.PIX) into the active buffer"),
    ("PAIR", 3, cmd_pair, "ON: reading x.RBS also reads x.PIX"),
    # Detector
    ("ANGLE", 2, cmd_angle, "detector axis to the untilted sample's normal, degrees"),
    ("TILTSIGN", 2, cmd_tiltsign, "how THETA moves the exit angle: 1, -1 or 0"),
    ("SOLID", 2, cmd_solid, "solid angle: msr, or area mm^2 and distance [MM|IN]"),
    ("WINDOW", 2, cmd_window, "detector window: element and thickness in µm"),
    ("CRYSTAL", 2, cmd_crystal, "detector crystal: element and thickness in µm"),
    ("FWHM", 2, cmd_fwhm, "resolution at Mn Ka, eV"),
    ("FANO", 2, cmd_fano, "Fano factor"),
    ("FILTER", 3, cmd_filter, "add an absorber (element, µm, hole %), or FILTER CLEAR"),
    ("CALIB", 3, cmd_calib, "energy calibration: keV/channel and offset in keV"),
    ("H", 1, cmd_h, "instrumental constant for K, L, M lines (H K L M, or H K|L|M v)"),
    ("ESCAPE", 2, cmd_escape, "Si escape peaks in the simulation: ON or OFF"),
    ("LINES", 2, cmd_lines, "table of the simulated lines (ALL: weak ones too)"),
    # The PIXE window
    ("PLOT", 2, cmd_plot, "draw the active buffer's PIXE spectrum"),
    ("REGION", 3, cmd_region, "energy range shown, keV (ALL for everything)"),
    ("COUNTS", 2, cmd_counts, "yield range (ALL to autoscale)"),
    ("LINEAR", 2, _scale("linear"), "linear yield axis"),
    ("SQRT", 2, _scale("sqrt"), "square-root yield axis"),
    ("LOG", 2, _scale("log"), "logarithmic yield axis (the default)"),
    ("MARKERS", 2, cmd_markers, "label the sample's lines: ON, ALL or OFF"),
    ("FIGSAVE", 3, cmd_figsave, "save the PIXE window to an image file"),
    ("HCOPY", -5, cmd_figsave, "synonym for FIGSAVE"),
]

for _name, _minlen, _handler, _help in _ENTRIES:
    TABLE.add(_name, _minlen, _handler, _help)
TABLE.note_synonym("HELP", "?")
TABLE.note_synonym("RETURN", "QUIT", "Q")
TABLE.note_synonym("FIGSAVE", "HCOPY")

_HELP_GROUPS: list[tuple[str, list[str]]] = [
    ("Getting around", ["HELP", "RETURN", "DISABLE", "SHOW"]),
    ("Data", ["GET", "PAIR"]),
    ("Detector",
     ["ANGLE", "TILTSIGN", "SOLID", "WINDOW", "CRYSTAL", "FWHM", "FANO", "FILTER", "CALIB",
      "ESCAPE"]),
    ("Simulation", ["H", "LINES"]),
    ("PIXE window",
     ["PLOT", "REGION", "COUNTS", "LINEAR", "SQRT", "LOG", "MARKERS", "FIGSAVE"]),
]
