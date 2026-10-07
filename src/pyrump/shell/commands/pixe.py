"""PIXE: the X-ray sub-processor.

Entering ``PIXE`` turns PIXE on: from then on every RBS PLOT, OVERLAY,
SPLOT and COMPARE shows the same buffers in the PIXE window
(:mod:`pyrump.shell.pixe_plotting`). ``DISABLE`` turns it off
and leaves the prompt in one step; ``RETURN`` leaves PIXE on. One-shot
``PIXE <command>`` from the RUMP level runs a single command without
entering, and does not turn PIXE on -- so ``~/.pyrumprc`` can set up the
detector with one-shots, or with a ``pixe`` ... ``disable`` block.

A PIXE spectrum lives on a buffer (:attr:`~pyrump.shell.session.Buffer.pixe`)
next to that buffer's RBS spectrum, and shares its beam and charge. The
detector, its geometry (GEOMETRY, THETA, PHI, PSI -- the tilt the buffer's
under THETA RBS), the calibration default and plot settings are the
session's (:class:`~pyrump.shell.session.PixeState`).
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np

from ...model.geometry import GeometryKind
from ...model.spectrum import Spectrum
from ...pixe.data import PixeData
from ...pixe.detector import COMPOUNDS, Absorber, angles, disc_solid_angle_msr
from ...pixe.yields import Exposure
from ..dispatch import ArgReader, CommandError, CommandTable
from ..session import Buffer
from .. import pixe_plotting
from .rump import FIGSAVE_DPI, EnterMode, Return, describe_topic

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
        print("PIXE enabled: PLOT and COMPARE now show PIXE too (DISABLE turns it off)")
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
    d, c, g = state.detector, state.calibration, state.geometry
    lines = [
        f"geometry {g.kind.name.lower()}",
        f"theta {'rbs' if state.theta_rbs else f'{g.theta:g}'}",
        f"phi {g.phi:g}",
        f"psi {g.psi:g}",
        f"solid {d.solid_angle_msr:.6g}",
        f"window {d.window.material} {d.window.thickness_um:g}",
        f"crystal {d.crystal.material} {d.crystal.thickness_um:g}",
        f"fwhm {d.fwhm_eV:g}",
        f"fano {d.fano:g}",
        # Replayed, FILTER adds to what is there: start from none.
        "filter clear",
    ]
    lines += [_filter_line(number, f) for number, f in enumerate(d.filters, start=1)]
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
    print(f"  {_angles_note(session)}")
    # The dose is the RBS buffer's, CORRECTION included: tuning CORR for
    # RBS rescales the PIXE simulation too.
    m = reference.measurement
    dose = f"  ! dose, shared with RBS: CHARGE {m.charge_uC:g} uC / CORRECTION {m.correction:g}"
    if m.correction and m.charge_state:
        exposure = Exposure(
            charge_uC=m.charge_uC, charge_state=m.charge_state, correction=m.correction
        )
        dose += (
            f" = {m.charge_uC / m.correction:.6g} uC,"
            f" {exposure.ions:.4g} ions (charge state {m.charge_state})"
        )
    print(dose)
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
    path = _pixe_path(Path(args.token("a PIXE spectrum file")), session)
    args.done()
    index = session.buffers.active
    buffer = session.buffers.get(index) if index else None
    if buffer is None:
        buffer = _pixe_only_buffer(session, path)
        index = session.buffers.scroll_in(buffer)
        session.buffers.active = index
    attach(session, buffer, path)
    print(f"  {path.name} -> PIXE of buffer {index}")
    pixe_plotting.refresh(session)


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
    pixe_plotting.refresh(session)


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


def _pixe_path(path: Path, session) -> Path:
    """The PIXE file, trying ``.PIX`` and ``.pix`` for a bare name (and,
    inside a macro, the macro's own folder first)."""
    found = session.locate(path, PIXE_SUFFIXES)
    if found is None:
        raise CommandError(f"no such file: {path}")
    return found.resolve()


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
    ``x.PIX`` from the same folder. ON also makes the PIXE tilt the RBS
    buffer's THETA: the two spectra of one run, one sample."""
    token = args.optional()
    args.done()
    if token is not None:
        if token.lower() not in ("on", "off"):
            raise CommandError("PAIR: expected ON or OFF")
        session.pixe.pair = token.lower() == "on"
        pixe_plotting.refresh(session)
    print(f"  pair {'on' if session.pixe.pair else 'off'}")
    if token is not None:
        print(f"  {_angles_note(session)}")


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


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------


def _angles_note(session) -> str:
    """What the PIXE geometry gives next to the ACTIVE buffer: the beam's
    and the X-rays' angles to the sample normal."""
    from ..pixe_sim import reference_buffer

    state = session.pixe
    geometry = state.geometry_for(reference_buffer(session).geometry)
    beam, out = angles(geometry)
    tilt = f"THETA {geometry.theta:g}"
    if state.pair:
        tilt += " from the RBS buffer (PAIR ON)"
    elif state.theta_rbs:
        tilt += " from the RBS buffer"
    return f"! {tilt}: beam {beam:g} deg, X-rays {out:.4g} deg to the sample normal"


def _set_geometry(session, **changes) -> None:
    session.pixe.geometry = replace(session.pixe.geometry, **changes)
    pixe_plotting.refresh(session)


def cmd_geometry(session, args: ArgReader) -> None:
    """``GEOMETRY [IBM|CORNELL|GENERAL]`` -- how the X-rays' exit angle PSI
    follows from THETA and PHI, as RBS's GEOMETRY does for its detector."""
    token = args.optional()
    args.done()
    if token is not None:
        try:
            kind = GeometryKind[token.upper()]
        except KeyError:
            raise CommandError("GEOMETRY: expected IBM, CORNELL or GENERAL") from None
        _set_geometry(session, kind=kind)
    print(f"  geometry {session.pixe.geometry.kind.name.lower()}")
    print(f"  {_angles_note(session)}")


cmd_geometry.details = """\
Where the PIXE detector stands, in the same terms as RBS's GEOMETRY (the
manual's Experimental geometry page draws all three):
  IBM      beside the beam, in the plane the sample turns in:
           PSI = |THETA + PHI|, a negative THETA turning it towards it
  CORNELL  above or below the beam: cos PSI = cos THETA x cos PHI
  GENERAL  anywhere: PSI as typed
Default: IBM, PHI 45, THETA 0 -- 45 deg from the beam, normal incidence."""


def cmd_theta(session, args: ArgReader) -> None:
    """``THETA [<deg>|RBS]`` -- the sample tilt for PIXE: its own, or
    the RBS buffer's (RBS), which then follows PERT THETA too."""
    token = args.optional()
    args.done()
    state = session.pixe
    if token is not None:
        if token.upper() == "RBS":
            state.theta_rbs = True
        else:
            try:
                value = float(token)
            except ValueError:
                raise CommandError("THETA: expected degrees or RBS") from None
            if state.pair:
                raise CommandError(
                    "THETA: PAIR ON takes the tilt from the RBS buffer; PAIR OFF first"
                )
            if not -90.0 < value < 90.0:
                raise CommandError("THETA: must be between -90 and 90")
            state.theta_rbs = False
            state.geometry = replace(state.geometry, theta=value)
        pixe_plotting.refresh(session)
    print(f"  theta {'rbs' if state.theta_rbs else f'{state.geometry.theta:g}'}")
    print(f"  {_angles_note(session)}")


cmd_theta.details = """\
The beam's angle to the sample normal, in degrees: 0 (the default) is
normal incidence. THETA RBS takes the RBS buffer's THETA instead -- the
same sample, tilted once for both detectors -- and follows it while PERT
fits THETA; a number unlinks it again. PAIR ON always takes the RBS
buffer's (the two spectra of one run) and refuses a number. The sign says which way the sample
turns, as for RBS: with GEOMETRY IBM a negative THETA turns it towards the
PIXE detector."""


def _geometry_angle(field: str, low: float, high: float, note: str):
    """A setter for PHI or PSI, in degrees, ``low <= value < high``."""
    def handler(session, args: ArgReader) -> None:
        value = args.optional_number()
        args.done()
        if value is not None:
            if not low <= value < high:
                raise CommandError(f"{field.upper()}: must be from {low:g} up to {high:g}")
            _set_geometry(session, **{field: value})
        print(f"  {field} {getattr(session.pixe.geometry, field):g}  ! {note}")
        print(f"  {_angles_note(session)}")

    handler.__doc__ = f"``{field.upper()} [<deg>]`` -- show or set the PIXE {note}."
    return handler


cmd_phi = _geometry_angle("phi", 0.0, 180.0, "deg, the detector to the beam (looking back)")
cmd_psi = _geometry_angle("psi", 0.0, 90.0, "deg, the detector to the normal: GENERAL only")


# ---------------------------------------------------------------------------
# Detector
# ---------------------------------------------------------------------------


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
cmd_fwhm.details = """\
The detector's resolution: the full width at half maximum of a line at
Mn Ka (5.9 keV), in eV -- what a detector's specification quotes. Lines
at other energies get their width from it and FANO.
RC43: 122 eV, fitted on the Mn Ka doublet of examples/MnPt."""

cmd_fano = _detector_number("fano", "fano", "")
cmd_fano.details = """\
How the peak width grows with X-ray energy. The width has two parts,
electronic noise (the same for every line) and charge statistics
(growing with energy E):
    FWHM(E)^2 = noise^2 + 2.355^2 x 3.64 eV x F x E
FWHM fixes the total at Mn Ka; F splits it between the two parts, so a
larger F makes lines below 5.9 keV (Si K) narrower and lines above it
(Pt L) wider. Si detectors: about 0.1. RC43: 0.104, from the Si Ka
(78.5 eV) and Mn Ka (122 eV) widths. Rarely needs changing."""


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
    open_area = _hole(args, what) if hole else None
    args.done()
    return Absorber(material, thickness, open_area or 0.0)


def _hole(args: ArgReader, what: str) -> float | None:
    """An optional open area: ``60``, ``60%``, ``HOLE 60`` or ``HOLE 60%``."""
    token = args.optional()
    if token is not None and token.lower() == "hole":
        token = args.token("the hole area, %")
    if token is None:
        return None
    try:
        value = float(token.rstrip("%"))
    except ValueError:
        raise CommandError(f"{what}: expected a hole area in %, not {token!r}") from None
    if not 0 <= value < 100:
        raise CommandError(f"{what}: hole area must be 0 to 100 %")
    return value


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


def _filter_number(token: str, count: int, *, existing: bool) -> int:
    """A filter's number, 1 facing the sample: an existing one, or (to set
    one) up to one past the last, so the list has no gaps."""
    try:
        number = int(token)
    except ValueError:
        raise CommandError(
            "FILTER: give the filter's number first, e.g. FILTER 1 MYLAR 50 "
            "(1 faces the sample)"
        ) from None
    if existing and not 1 <= number <= count:
        raise CommandError(f"FILTER: there is no filter {number}")
    if not existing and not 1 <= number <= count + 1:
        raise CommandError(f"FILTER: set filter {count + 1} first (filters are numbered without gaps)")
    return number


def _filter_line(number: int, f: Absorber) -> str:
    hole = f" hole {f.hole_percent:g}%" if f.hole_percent else ""
    return f"filter {number} {f.material} {f.thickness_um:g}{hole}"


#: Energies at which FILTER reports the filters' total transmission, keV.
_TRANSMISSION_AT = (1.5, 2.5, 5.0, 10.0)


def cmd_filter(session, args: ArgReader) -> None:
    """``FILTER`` lists the absorbers between sample and window, numbered
    from the sample outwards (1 faces the sample), with their total
    transmission. ``FILTER <n> <element|MYLAR|KAPTON> <µm> [[HOLE] <%>]``
    sets filter *n* -- replacing it, or adding it after the last; ``FILTER
    CLEAR`` removes them all, ``FILTER CLEAR <n>`` one of them (the ones
    after it move up)."""
    from ...pixe.atomic import atomic_data
    from ...pixe.yields import transmission

    detector = session.pixe.detector
    filters = list(detector.filters)
    if args:
        first = args.token()
        if first.lower() == "clear":
            token = args.optional()
            args.done()
            if token is None:
                filters = []
            else:
                del filters[_filter_number(token, len(filters), existing=True) - 1]
        else:
            number = _filter_number(first, len(filters), existing=False)
            absorber = _absorber(session, args, f"FILTER {number}", hole=True)
            if number > len(filters):
                filters.append(absorber)
            else:
                filters[number - 1] = absorber
        session.pixe.detector = replace(detector, filters=tuple(filters))
        pixe_plotting.refresh(session)
    if not filters:
        print("  no filters")
        return
    for number, f in enumerate(filters, start=1):
        print(f"  {_filter_line(number, f)}  ! µm")
    atomic = atomic_data()
    total = [
        float(np.prod([transmission(f, e, session.table, atomic) for f in filters]))
        for e in _TRANSMISSION_AT
    ]
    print("  ! transmission: " + ", ".join(
        f"{e:g} keV {value:.3g}" for e, value in zip(_TRANSMISSION_AT, total)))


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
        pixe_plotting.refresh(session)
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
        print("  no lines: the beam reaches no element of the sample")
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
# Exports
# ---------------------------------------------------------------------------


def pixe_target(target: Path) -> Path:
    """The PIXE export's file next to ``target``: ``fit.txt`` -> ``fit_pixe.txt``."""
    return target.with_name(f"{target.stem}_pixe{target.suffix}")


def _pixe_export_buffer(session) -> tuple[int, Buffer]:
    index = session.buffers.active
    buffer = session.buffers.get(index) if index else None
    if buffer is None or buffer.pixe is None:
        raise CommandError("no PIXE spectrum in the active buffer: PIXE GET <file>")
    return index, buffer


def _export_angles(session, rbs) -> str:
    """The PIXE geometry for an export header, next to the RBS geometry
    ``rbs`` of the buffer exported."""
    geometry = session.pixe.geometry_for(rbs)
    beam, out = angles(geometry)
    return (
        f"geometry {geometry.kind.name.lower()}  theta {geometry.theta:g}"
        f"{' (RBS)' if session.pixe.follows_rbs else ''}  phi {geometry.phi:g}"
        f"  psi {geometry.psi:g}: beam {beam:g} deg, X-rays {out:.4g} deg to the sample normal"
    )


def _pixe_header(session, index: int, buffer: Buffer, what: str) -> list[str]:
    """What is needed to read a PIXE export: the spectrum, the run it shares
    with the RBS spectrum, and the PIXE setup."""
    import datetime

    from ... import __version__
    from .rump import _wrascii_beam_code

    data, b, g, m = buffer.pixe, buffer.beam, buffer.geometry, buffer.measurement
    c = data.calibration
    times = (f"live {data.live_time_s:g} s  real {data.real_time_s:g} s"
             f"  (live fraction {data.live_fraction:.6f})"
             if data.live_time_s and data.real_time_s else "live/real time unknown")
    return [
        f"pyRUMP {__version__} PIXE {what} export, {datetime.date.today().isoformat()}",
        f"PIXE data      {data.path.name if data.path else data.identifier}"
        f"  ({data.path if data.path else '-'}, buffer {index})",
        f"Date           {data.date}" if data.date else f"Date           {buffer.date}",
        f"Times          {times}",
        f"Calibration    {c.kevch:.8g} keV/ch  offset {c.kev0:.8g} keV"
        f"  first channel {c.first:g}  {c.npt} channels",
        f"Beam           {_wrascii_beam_code(session, b)}  {b.e0_MeV:.6f} MeV",
        f"Angles         {_export_angles(session, g)}",
        f"Dose           charge {m.charge_uC:.6f} uC  charge state {m.charge_state}"
        f"  corr {m.correction:.6f}",
        f"Setup          {'; '.join(setup_lines(session))}",
    ]


def _pixe_columns(buffer: Buffer) -> dict[str, np.ndarray]:
    spectrum = buffer.pixe.spectrum
    return {
        "channel": np.arange(spectrum.counts.size) + int(round(spectrum.calibration.first)),
        "energy_keV": spectrum.energies,
    }


def write_pixe_export(session, target: Path) -> None:
    """PIXE ``EXPORT``: the active buffer's PIXE spectrum as columns."""
    from ...io.ascii import write_columns
    from .rump import _export_delimiter

    index, buffer = _pixe_export_buffer(session)
    counts = np.asarray(buffer.pixe.spectrum.counts, dtype=np.float64)
    columns = _pixe_columns(buffer)
    columns["counts"] = counts
    columns["error"] = np.sqrt(np.clip(counts, 0.0, None))
    write_columns(
        target, columns,
        comments=[*_pixe_header(session, index, buffer, "spectrum"),
                  "channel is numbered as in the spectrum file; energy_keV is the lower edge "
                  "of each channel; error is sqrt(counts)"],
        delimiter=_export_delimiter(target), formats={"channel": "d", "energy_keV": ".4f"},
    )
    print(f"wrote {target}: {counts.size} PIXE channels from buffer {index}")


def write_pixe_exportcmp(session, target: Path) -> None:
    """PIXE ``EXPORTCMP``: data, simulation (total and per element),
    difference and Poisson residual."""
    from ...fit.objective import chi_square, poisson_residuals
    from ...io.ascii import write_columns
    from ..pixe_sim import simulate
    from .rump import _export_delimiter, _export_sample_lines

    index, buffer = _pixe_export_buffer(session)
    try:
        result = simulate(session)
    except (ValueError, KeyError) as error:
        raise CommandError(f"EXPORTCMP: {str(error).strip(chr(39))}") from None
    if result is None:
        raise CommandError("EXPORTCMP: no SIM sample to simulate")
    spectrum = buffer.pixe.spectrum
    counts = np.asarray(spectrum.counts, dtype=np.float64)
    simulation = np.asarray(result.counts, dtype=np.float64)
    residual, _ = poisson_residuals(counts, simulation)
    residual[simulation <= 0] = np.nan

    low, high = pixe_plotting.region_indices(session.pixe, spectrum.calibration)
    mask = np.zeros(counts.size, dtype=bool)
    mask[low:high + 1] = True
    summary = chi_square(counts, simulation, valid=mask)
    span = (f"channels {spectrum.calibration.first + low:g}-{spectrum.calibration.first + high:g}")
    gof = [f"GOF            reduced chi-square {summary.reduced:.4f} ({summary.dof} dof) over {span}"]

    columns = _pixe_columns(buffer)
    columns.update(counts=counts, simulation=simulation, diff=counts - simulation,
                   residual=residual)
    for symbol, part in result.by_element.items():
        columns[f"sim_{symbol}"] = np.asarray(part, dtype=np.float64)
    notes = ["channel is numbered as in the spectrum file; energy_keV is the lower edge of "
             "each channel; residual is the Poisson residual in sigma (nan where the "
             "simulation is zero); sim_<element> are each element's simulated lines"]
    write_columns(
        target, columns,
        comments=[*_pixe_header(session, index, buffer, "COMPARE"),
                  *_export_sample_lines(session), *gof, *notes],
        delimiter=_export_delimiter(target), formats={"channel": "d", "energy_keV": ".4f"},
    )
    print(f"wrote {target}: {counts.size} PIXE channels, reduced chi-square "
          f"{summary.reduced:.4f} ({summary.dof} dof)")


def paired_export(session, target: Path, *, compare: bool) -> None:
    """After a RUMP-level EXPORT/EXPORTCMP: with PAIR ON and a PIXE spectrum
    in the active buffer, the PIXE file too, as ``<name>_pixe``."""
    buffer = session.buffers.active_buffer
    if not session.pixe.pair or not session.buffers.active or buffer is None or buffer.pixe is None:
        return
    (write_pixe_exportcmp if compare else write_pixe_export)(session, pixe_target(target))


def cmd_export(session, args: ArgReader) -> None:
    """``EXPORT <file>`` -- write the active buffer's PIXE spectrum as plain
    columns (channel, energy, counts, error) to ``<file>_pixe``, under a
    ``#``-commented header with the PIXE setup. At the RUMP level, EXPORT
    writes the RBS spectrum, and with PAIR ON this file as well."""
    from .rump import _export_target

    write_pixe_export(session, pixe_target(_export_target(args)))


def cmd_exportcmp(session, args: ArgReader) -> None:
    """``EXPORTCMP <file>`` (synonym ``EC``) -- write the PIXE comparison as
    columns to ``<file>_pixe``: channel, energy, counts, simulation, diff,
    Poisson residual and each element's simulated lines; chi-square in the
    header."""
    from .rump import _export_target

    write_pixe_exportcmp(session, pixe_target(_export_target(args)))


# ---------------------------------------------------------------------------
# The PIXE window
# ---------------------------------------------------------------------------


def _redraw_or(session, message: str) -> None:
    """Redraw an open PIXE window; otherwise just report the setting."""
    if pixe_plotting.is_open(session) and session.pixe.view:
        pixe_plotting.draw(session, required=False)
    else:
        print(message)


def cmd_plot(session, args: ArgReader) -> None:
    """``PLOT [buffer]`` -- draw a buffer's PIXE spectrum (default: the
    active one); ``PLOT 0`` the simulation. Leaves the RBS window alone."""
    token = args.optional()
    args.done()
    index = session.buffers.active if token is None else _resolve(session, token)
    item = pixe_plotting.ViewItem("sim") if index == 0 else pixe_plotting.ViewItem("data", index)
    session.pixe.view, session.pixe.compare = [item], False
    pixe_plotting.draw(session)


def cmd_compare(session, args: ArgReader) -> None:
    """``COMPARE`` -- the active buffer's PIXE spectrum against the
    simulation, with residuals. Leaves the RBS window alone."""
    args.done()
    index = session.buffers.active
    if session.buffers.get(index) is None or index == 0:
        raise CommandError("no active data buffer: GET a spectrum first")
    session.pixe.view = [pixe_plotting.ViewItem("data", index), pixe_plotting.ViewItem("sim")]
    session.pixe.compare = True
    pixe_plotting.draw(session)


def _resolve(session, token: str) -> int:
    try:
        index = session.resolve(token)
    except KeyError as error:
        raise CommandError(str(error).strip("'")) from None
    if index != 0 and session.buffers.get(index) is None:
        raise CommandError(f"buffer {index} is empty")
    return index


def cmd_region(session, args: ArgReader) -> None:
    """``REGION <channel> <channel>`` -- the channel range shown (as the
    spectrum file numbers its channels); ``REGION ALL`` the whole spectrum;
    no argument shows it. The energy axis below and the channel axis on top
    follow."""
    state = session.pixe
    if args and args.peek().lower() == "all":
        args.token()
        state.low = state.high = None
    elif args:
        low = args.integer("the first channel")
        high = args.integer("the last channel")
        if low == high:
            raise CommandError(f"empty region: {low} to {high}")
        state.low, state.high = min(low, high), max(low, high)
    args.done()
    span = "all" if state.low is None and state.high is None else f"{state.low} to {state.high}"
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
    ("PAIR", 3, cmd_pair, "ON: reading x.RBS also reads x.PIX, THETA follows RBS"),
    # Detector
    # Geometry
    ("GEOMETRY", 4, cmd_geometry, "where the detector is: IBM, CORNELL or GENERAL"),
    ("THETA", 3, cmd_theta, "sample tilt, degrees (0: normal incidence), or RBS"),
    ("PHI", 3, cmd_phi, "detector to the beam (looking back), degrees"),
    ("PSI", 3, cmd_psi, "detector to the sample normal, degrees: GENERAL only"),
    ("SOLID", 2, cmd_solid, "solid angle: msr, or area mm^2 and distance [MM|IN]"),
    ("WINDOW", 2, cmd_window, "detector window: element and thickness in µm"),
    ("CRYSTAL", 2, cmd_crystal, "detector crystal: element and thickness in µm"),
    ("FWHM", 2, cmd_fwhm, "resolution: peak width (FWHM) at Mn Ka, eV"),
    ("FANO", 2, cmd_fano, "Fano factor: how the peak width grows with energy"),
    ("FILTER", 3, cmd_filter, "list filters; FILTER n material µm [hole %]; FILTER CLEAR [n]"),
    ("CALIB", 3, cmd_calib, "energy calibration: keV/channel and offset in keV"),
    ("H", 1, cmd_h, "instrumental constant for K, L, M lines (H K L M, or H K|L|M v)"),
    ("ESCAPE", 2, cmd_escape, "Si escape peaks in the simulation: ON or OFF"),
    # Typed in full only: LIN, LINE stay LINEAR's, as at the RUMP level.
    ("LINES", 0, cmd_lines, "table of the simulated lines (ALL: weak ones too)"),
    # The PIXE window
    ("PLOT", 2, cmd_plot, "draw a buffer's PIXE spectrum (0: the simulation)"),
    ("COMPARE", 0, cmd_compare, "PIXE data against the simulation, with residuals"),
    ("CMP", -3, cmd_compare, "synonym for COMPARE"),
    ("REGION", 3, cmd_region, "channel range shown (ALL for everything)"),
    ("COUNTS", 2, cmd_counts, "yield range (ALL to autoscale)"),
    ("LINEAR", 2, _scale("linear"), "linear yield axis"),
    ("SQRT", 2, _scale("sqrt"), "square-root yield axis"),
    ("LOG", 2, _scale("log"), "logarithmic yield axis (the default)"),
    ("MARKERS", 2, cmd_markers, "label the sample's lines: ON, ALL or OFF"),
    ("FIGSAVE", 3, cmd_figsave, "save the PIXE window to an image file"),
    ("EXPORTCMP", 7, cmd_exportcmp, "write the PIXE comparison as columns to <file>_pixe"),
    ("EC", -2, cmd_exportcmp, "synonym for EXPORTCMP"),
    ("EXPORT", 4, cmd_export, "write the PIXE spectrum as columns to <file>_pixe"),
    ("HCOPY", -5, cmd_figsave, "synonym for FIGSAVE"),
]

for _name, _minlen, _handler, _help in _ENTRIES:
    TABLE.add(_name, _minlen, _handler, _help)
TABLE.note_synonym("HELP", "?")
TABLE.note_synonym("RETURN", "QUIT", "Q")
TABLE.note_synonym("FIGSAVE", "HCOPY")
TABLE.note_synonym("COMPARE", "CMP")
TABLE.note_synonym("EXPORTCMP", "EC")

_HELP_GROUPS: list[tuple[str, list[str]]] = [
    ("Getting around", ["HELP", "RETURN", "DISABLE", "SHOW"]),
    ("Data", ["GET", "PAIR"]),
    ("Geometry", ["GEOMETRY", "THETA", "PHI", "PSI"]),
    ("Detector",
     ["SOLID", "WINDOW", "CRYSTAL", "FWHM", "FANO", "FILTER", "CALIB",
      "ESCAPE"]),
    ("Simulation", ["H", "LINES"]),
    ("PIXE window",
     ["PLOT", "COMPARE", "REGION", "COUNTS", "LINEAR", "SQRT", "LOG", "MARKERS", "FIGSAVE"]),
    ("Exports", ["EXPORT", "EXPORTCMP"]),
]
