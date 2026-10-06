"""SNAPSHOT: save the session as it stands, to pick it up again later.

A pyRUMP addition. Whatever brought the session here -- a PERT ``GO``, hand
tuning in SIM, a CORRECTION typed at the prompt -- ``SNAPSHOT`` writes it out
under the sample's name:

- ``<sample>_fit.xeq`` -- a restore macro: it reloads the spectrum from the
  file it was read from, sets every parameter of the buffer (IDENTIFIER and
  DATE too, in case they were corrected), the simulation settings, the PIXE
  setup (with ``PAIR ON``), the sample and the PERT setup, and ends on
  ``COMPARE``. ``pyrump <sample>_fit.xeq`` (or ``XEQ <sample>_fit`` in
  pyRUMP) is the session again. ``.xeq`` is RUMP's own macro extension, and
  unlike ``.cmd`` not one Windows runs as a batch file on a double-click.
- ``<sample>.lcm``, ``<sample>.pert`` -- the sample and the PERT selection,
  which the restore macro reads back
- ``<sample>.png`` -- the ``COMPARE`` plot (and ``<sample>_pixe.png``)
- ``<sample>.report`` -- one block appended per snapshot: the fit with its
  uncertainties if nothing has changed since ``GO``, otherwise the current
  values and chi-square, marked as set by hand

``REPORT ON`` is ``GO`` followed by this.

The restore macro is plain commands, so it can be read and edited. Paths in
it are relative to its own folder (:meth:`Session.locate`), so the folder can
be moved, along with the data, as a whole.
"""

from __future__ import annotations

import os
from datetime import datetime
from pathlib import Path, PurePath

from .dispatch import ArgReader, CommandError
from .session import Buffer

#: What the restore macro's name adds to the sample's: ``XEQ <sample>``
#: still finds the raw ``<sample>.RBS``, and a changed spectrum saved as
#: ``<sample>_fit.rbs`` can't overwrite it.
RESTORE_SUFFIX = "_fit"

#: The restore macro's extension.
RESTORE_EXTENSION = ".xeq"


# ---------------------------------------------------------------------------
# What is snapshotted
# ---------------------------------------------------------------------------


def data_buffer(session) -> Buffer | None:
    """The active data buffer -- not the simulation -- or ``None``."""
    if session.buffers.active < 1:
        return None
    return session.buffers.active_buffer


def stem(buffer) -> str:
    """The sample name every snapshot file starts with: the spectrum's
    IDENTIFIER, else its file name (:func:`~pyrump.shell.plotting.buffer_stem`).
    A spectrum reloaded from a snapshot's own ``<sample>_fit.rbs`` keeps the
    sample's name, not that file's."""
    from . import plotting

    name = plotting.buffer_stem(buffer) or "buffer"
    if name.endswith(RESTORE_SUFFIX) and len(name) > len(RESTORE_SUFFIX):
        name = name[: -len(RESTORE_SUFFIX)]
    return name


def _missing(session, buffer: Buffer | None) -> list[str]:
    """What a snapshot needs that the session doesn't have."""
    missing = []
    if buffer is None:
        missing.append(
            "no measured spectrum is active -- load it first "
            "(XEQ <file>.RBS for an acquisition file, GET <file> otherwise), "
            "or POINTAT its buffer"
        )
    if not session.script.layers:
        missing.append(
            "no SIM sample -- read one (SIM GET <file>.lcm) or build one in SIM"
        )
    return missing


def _no_source(session, buffer: Buffer) -> str | None:
    """Why the restore macro could not reload ``buffer``, or ``None``."""
    if buffer.path is not None or buffer.macro is not None:
        return None
    return (
        f"buffer {session.buffers.active} was not read from a file "
        "(EMPTY, COPY or similar), so there is nothing to reload it from -- "
        "WRITE it to a .rbs file and GET that first"
    )


def _wants_pixe(session, buffer: Buffer) -> bool:
    """PIXE goes into the snapshot with PAIR ON -- or whenever the spectrum
    has a PIXE spectrum, or PIXE is on, however it got there."""
    state = session.pixe
    return state.pair or state.enabled or buffer.pixe is not None


# ---------------------------------------------------------------------------
# The restore macro's lines
# ---------------------------------------------------------------------------


def _number(value) -> str:
    """A float exactly as stored: ``repr`` round-trips to the last bit."""
    return repr(float(value))


def _quoted(path: Path, folder: Path) -> str:
    """``path`` for a macro in ``folder``: relative to it when possible, with
    forward slashes (Windows reads those too), quoted against spaces."""
    try:
        text = PurePath(os.path.relpath(path, folder)).as_posix()
    except ValueError:  # another drive, on Windows
        text = str(path)
    return f'"{text}"'


def _settings_lines(session) -> list[str]:
    settings = session.settings
    return [
        f"FAITHFUL {'on' if settings.faithful else 'off'}",
        f"SCREENING {settings.screening.name}",
        f"MODE {session.thickness_mode}",
    ]


def _beam_line(session, buffer: Buffer) -> str:
    """``BEAM He+`` or ``BEAM 4He+`` -- whichever gives back exactly this
    mass (the element's average mass, or the isotope's) -- or a comment, if
    neither does (a mass read from the data file, which reloading it
    restores)."""
    from ..cli._common import resolve_beam

    beam = buffer.beam
    charge = "+" * max(buffer.measurement.charge_state, 0)
    try:
        symbol = session.table.by_z(beam.z).symbol
    except Exception:  # noqa: BLE001 - an unknown Z: fall back to the comment
        symbol = None
    for spec in (symbol, f"{round(beam.mass)}{symbol}") if symbol else ():
        try:
            exact = resolve_beam(session.table, spec) == (beam.z, beam.mass)
        except Exception:  # noqa: BLE001 - not a spec the table reads
            exact = False
        if exact:
            return f"BEAM {spec}{charge}"
    return f"! beam Z={beam.z} mass={beam.mass!r}, as the data file gives it"


def _buffer_lines(session, buffer: Buffer) -> list[str]:
    """Every parameter of the buffer, as the commands that set it."""
    c, g, m = buffer.calibration, buffer.geometry, buffer.measurement
    return [
        _beam_line(session, buffer),
        f"MEV {_number(buffer.beam.e0_MeV)}",
        f"GEOMETRY {g.kind.name.lower()}",
        f"THETA {_number(g.theta)}",
        f"PHI {_number(g.phi)}",
        f"PSI {_number(g.psi)}",
        f"CONVERSION {_number(c.kevch)} {_number(c.kev0)}",
        f"CHOFF {_number(c.first)}",
        f"FWHM {_number(m.fwhm_keV)}",
        f"OMEGA {_number(m.omega_msr)}",
        f"TAU {_number(m.tau_us)}",
        f"CHARGE {_number(m.charge_uC)}",
        f"CORRECTION {_number(m.correction)}",
        f"CURRENT {_number(m.current_nA)}",
    ]


def _quote_text(text: str) -> str | None:
    """``text`` as one token, spaces and all: in whichever quotes it doesn't
    contain, or ``None`` if it contains both (the tokenizer has no escapes)."""
    for quote in ("'", '"'):
        if quote not in text:
            return f"{quote}{text}{quote}"
    return None


def _label_lines(buffer: Buffer) -> list[str]:
    """IDENTIFIER and DATE, which reloading the file would otherwise put back
    as the file has them, even after a typo was corrected."""
    lines = []
    for command, text in (("IDENTIFIER", buffer.identifier), ("DATE", buffer.date)):
        if not text:
            continue
        quoted = _quote_text(text)
        if quoted is None:
            lines.append(f"! {command} has both kinds of quote, kept as the file has it")
        else:
            lines.append(f"{command} {quoted}")
    return lines


def _pixe_setup_lines(session) -> list[str]:
    """The PIXE settings, PAIR left out (it goes after the data is read, so
    reading the data doesn't pair it a second time)."""
    from .commands.pixe import setup_lines

    return [f"PIXE {line}" for line in setup_lines(session) if not line.startswith("pair ")]


def _view_lines(session) -> list[str]:
    """The RBS plot as it was drawn."""
    plot = session.plot
    lines = []
    if plot.low is not None and plot.high is not None:
        lines.append(f"REGION {plot.low} {plot.high}")
    if plot.ylow is not None and plot.yhigh is not None:
        lines.append(f"COUNTS {plot.ylow:g} {plot.yhigh:g}")
    lines += [
        plot.yscale.upper(),
        "NORMALIZE" if plot.normalized else "RAW",
        "ENERGY" if plot.energy_axis else "ENERGY OFF",
        "LABELS" if plot.labels else "LABELS OFF",
        "STRUCTLABEL" if plot.structure_labels else "STRUCTLABEL OFF",
        "COMPFRAC" if plot.composition_fraction else "COMPFRAC OFF",
    ]
    return lines


def _pixe_view_lines(session) -> list[str]:
    state = session.pixe
    plot = state.plot
    region = "ALL" if state.low is None else f"{state.low} {state.high}"
    counts = "ALL" if plot.ylow is None else f"{plot.ylow:g} {plot.yhigh:g}"
    return [f"PIXE REGION {region}", f"PIXE COUNTS {counts}", f"PIXE {plot.yscale.upper()}"]


def _pert_lines(session) -> list[str]:
    from .commands.pert import _to_lines

    return _to_lines(session.pert) if session.pert is not None else []


def fingerprint(session, *, labels: bool = True) -> str | None:
    """Everything a snapshot saves that shapes the simulation or the fit --
    the data, the buffer's parameters, the settings, the sample, the PERT and
    PIXE setup -- as one string. Equal strings: nothing to save. ``None``
    when there is nothing to snapshot. Plot scaling is left out: a REGION
    changes how the result looks, not the result.

    ``labels`` adds IDENTIFIER and DATE: worth saving, but a corrected
    IDENTIFIER doesn't make a fit any less the fit (GO's own record leaves
    them out)."""
    from ..script.lcm import write_lcm
    from .session import counts_digest

    buffer = data_buffer(session)
    if _missing(session, buffer):
        return None
    parts = [
        f"source {buffer.path or buffer.macro}",
        f"counts {counts_digest(buffer.spectrum.counts)}",
        *_buffer_lines(session, buffer),
        *(_label_lines(buffer) if labels else []),
        *_settings_lines(session),
        write_lcm(session.script),
        *_pert_lines(session),
    ]
    if _wants_pixe(session, buffer):
        from .commands.pixe import setup_lines

        pixe = buffer.pixe
        parts += setup_lines(session)
        parts.append(f"pixe {pixe.path if pixe is not None else None} {session.pixe.enabled}")
        if pixe is not None:
            parts.append(f"pixe counts {counts_digest(pixe.spectrum.counts)}")
    return "\n".join(parts)


def unsaved(session) -> bool:
    """Whether the session holds something no SNAPSHOT has saved yet."""
    state = fingerprint(session)
    return state is not None and state != session.saved_state


# ---------------------------------------------------------------------------
# The report block
# ---------------------------------------------------------------------------


def _current_values(session, buffer: Buffer) -> list[str]:
    """The PERT selection's current values, laid out as GO lays out a fit's
    -- without uncertainties, which only a fit has."""
    from ..fit.parameters import FitInputs
    from ..script.lcm import to_sample
    from .commands.pert import _display, _units_per_areal

    state = session.pert
    if state is None or not state.varying:
        return []
    inputs = FitInputs(
        sample=to_sample(session.script, session.table, session.densities),
        beam=buffer.beam,
        geometry=buffer.geometry,
        calibration=buffer.calibration,
        measurement=buffer.measurement,
    )
    lines = []
    for entry in state.varying:
        per_areal = (
            _units_per_areal(session, inputs, entry.layer) if entry.kind == "thickness" else None
        )
        factor, fmt, unit = _display(session, entry, per_areal)
        value = entry.parameter.get(inputs) * factor
        lines.append(f"  {entry.name:26s} {fmt(value) + (f' {unit}' if unit else ''):>14s}")
    return lines


def _manual_lines(session, buffer: Buffer) -> list[str]:
    """The report block for a state no GO produced."""
    from ..script.lcm import structure_label
    from . import plotting
    from .commands.rump import _fit_window_description

    label = structure_label(session.script, normalize=session.plot.composition_fraction)
    lines = [f"  {stem(buffer)}: {label}"]
    if session.last_fit is None:
        lines.append("  set by hand -- not by GO, so there are no uncertainties")
    else:
        lines.append(
            "  changed since the last GO -- that fit's uncertainties no longer apply"
        )
    theory = session.simulation()
    n_channels = min(buffer.n_channels, theory.n_channels)
    try:
        region = session.plot.region(n_channels)
    except ValueError:
        region = (0, n_channels - 1)
    gof = plotting.goodness_of_fit(session, buffer, theory, n_channels, region, raw=True)
    lines.append(f"\n  {gof}")
    lines.append(f"  over {_fit_window_description(session, region)}")
    values = _current_values(session, buffer)
    if values:
        lines.append("")
        lines += values
    return lines


def _context_lines(session, buffer: Buffer) -> list[str]:
    """The sample and the experiment, in full, for the record."""
    from .commands.sim import describe, editor_for

    return [
        "\n  sample:",
        *(f"  {line}" for line in describe(session, editor_for(session)).splitlines()),
        "\n  data:",
        *(f"  {line}" for line in buffer.describe().splitlines()),
    ]


# ---------------------------------------------------------------------------
# Taking one
# ---------------------------------------------------------------------------


def _target(session, name: str | None, buffer: Buffer) -> tuple[Path, str]:
    """The folder to write in and the sample name, from ``SNAPSHOT [name]``."""
    if name is None:
        return Path("."), stem(buffer)
    path = Path(name).expanduser()
    for suffix in (RESTORE_EXTENSION, ".cmd", ".report", ".lcm", ".pert", ".png"):
        if path.name.lower().endswith(suffix):
            path = path.with_name(path.name[: -len(suffix)])
    if path.name.endswith(RESTORE_SUFFIX) and len(path.name) > len(RESTORE_SUFFIX):
        path = path.with_name(path.name[: -len(RESTORE_SUFFIX)])
    folder = path.parent
    if not folder.is_dir():
        raise CommandError(f"SNAPSHOT: no such folder: {folder}")
    return folder, path.name


def _restore_macro(
    session, buffer: Buffer, folder: Path, base: str, fitted: bool, wrote_pert: bool,
    timestamp: str,
) -> list[str]:
    """The restore macro's lines, writing ``<base>_fit.rbs`` first if the
    counts have changed since they were read."""
    from ..io.rbs import write_rbs

    restore = f"{base}{RESTORE_SUFFIX}"
    how = "as GO fitted it" if fitted else "set by hand"
    lines = [
        f"! pyRUMP SNAPSHOT of {base}, {timestamp} -- {how}",
        f"! Restore with:  pyrump {restore}{RESTORE_EXTENSION}"
        f"    (in pyRUMP:  XEQ {restore})",
        f"! The next SNAPSHOT of {base} here replaces this file.",
        "",
        "! Simulation settings",
        *_settings_lines(session),
    ]
    pixe = _wants_pixe(session, buffer)
    if pixe:
        lines += ["", "! PIXE detector and simulation", *_pixe_setup_lines(session)]

    lines += [""]
    source = buffer.path or buffer.macro
    if buffer.counts_changed():
        target = folder / f"{restore}.rbs"
        write_rbs(target, buffer.to_rbs())
        print(f"  wrote {target}")
        lines += [
            f"! The spectrum, changed since it was read from {source.name}",
            "! (BACKGROUND, SMOOTH ...): reloaded as it is now",
            f"GET {_quoted(target.resolve(), folder.resolve())}",
        ]
    elif buffer.path is not None:
        lines += [
            "! The spectrum, from the file it was read from",
            f"GET {_quoted(buffer.path, folder.resolve())}",
        ]
    else:
        lines += [
            "! The spectrum, from the acquisition file it was read from",
            f"XEQ {_quoted(buffer.macro, folder.resolve())}",
        ]
    lines += [
        "! Its parameters as they were (a fit may have changed some)",
        *_buffer_lines(session, buffer),
        *_label_lines(buffer),
    ]
    if pixe:
        spectrum = buffer.pixe
        if spectrum is not None and spectrum.path is not None:
            lines.append(f"PIXE GET {_quoted(spectrum.path, folder.resolve())}")
        elif spectrum is not None:
            lines.append("! its PIXE spectrum was not read from a file: not restored")
        lines.append(f"PIXE pair {'on' if session.pixe.pair else 'off'}")

    lines += [
        "",
        "! The sample and the fit setup",
        f"SIM GET {_quoted((folder / f'{base}.lcm').resolve(), folder.resolve())}",
    ]
    if wrote_pert:
        lines.append(f"PERT GET {_quoted((folder / f'{base}.pert').resolve(), folder.resolve())}")
    else:
        lines.append("PERT CLEAR")

    lines += ["", "! The plot", *_view_lines(session)]
    if pixe:
        lines += _pixe_view_lines(session)
        if session.pixe.enabled:
            lines += ["PIXE", "RETURN"]
    lines += ["COMPARE", "", "! Marks the session as saved, for QUIT", "SNAPSHOT -restored"]
    return lines


def _save_plots(session, folder: Path, base: str, buffer: Buffer) -> None:
    """``<base>.png`` from a fresh COMPARE, and ``<base>_pixe.png`` when the
    PIXE window shows this snapshot's PIXE too."""
    from .commands.rump import FIGSAVE_DPI, cmd_compare

    try:
        cmd_compare(session, ArgReader([], command="compare"))
    except CommandError as error:
        print(f"  no {base}.png: {error}")
        return
    png = folder / f"{base}.png"
    session.figure.savefig(png, dpi=FIGSAVE_DPI)
    print(f"  wrote {png}")
    figure = session.pixe.figure
    if _wants_pixe(session, buffer) and session.pixe.enabled and figure is not None:
        png = folder / f"{base}_pixe.png"
        figure.savefig(png, dpi=FIGSAVE_DPI)
        print(f"  wrote {png}")


def take(session, name: str | None = None, *, after_go: bool = False) -> None:
    """Write a snapshot of the session (see the module docstring).

    Refuses, naming what is missing, without a data buffer and a SIM
    sample. A data buffer not read from any file can't be reloaded: an
    explicit SNAPSHOT refuses that too, while REPORT (``after_go``) still
    records the fit and only says the restore macro is missing.
    """
    from ..script.lcm import write_lcm

    buffer = data_buffer(session)
    missing = _missing(session, buffer)
    if missing:
        raise CommandError(
            "SNAPSHOT not possible:\n" + "\n".join(f"  - {line}" for line in missing)
        )
    no_source = _no_source(session, buffer)
    if no_source is not None and not after_go:
        raise CommandError(f"SNAPSHOT not possible:\n  - {no_source}")

    folder, base = _target(session, name, buffer)
    state = fingerprint(session)
    fitted = session.last_fit is not None and (
        session.last_fit[0] == fingerprint(session, labels=False)
    )
    timestamp = datetime.now().isoformat(timespec="seconds")

    block = list(session.last_fit[1]) if fitted else _manual_lines(session, buffer)
    block += _context_lines(session, buffer)
    restore = folder / f"{base}{RESTORE_SUFFIX}{RESTORE_EXTENSION}"
    if no_source is None:
        block.append(f"\n  restore: pyrump {restore}")
    else:
        block.append(f"\n  no restore macro: {no_source}")

    report = folder / f"{base}.report"
    with report.open("a", encoding="utf-8") as handle:
        handle.write(f"--- {timestamp}  {'fit (GO)' if fitted else 'set by hand'} ---\n")
        handle.write("\n".join(block) + "\n\n")
    print(f"\n  updated {report}")

    pert_lines = _pert_lines(session)
    if pert_lines:
        pert = folder / f"{base}.pert"
        pert.write_text("\n".join(pert_lines) + "\n")
        print(f"  wrote {pert}")

    lcm = folder / f"{base}.lcm"
    lcm.write_text(write_lcm(session.script))
    print(f"  wrote {lcm}")

    _save_plots(session, folder, base, buffer)

    if no_source is None:
        lines = _restore_macro(
            session, buffer, folder, base, fitted, bool(pert_lines), timestamp
        )
        restore.write_text("\n".join(lines) + "\n")
        print(f"  wrote {restore}")
        session.saved_state = state
    else:
        print(f"  no restore macro: {no_source}")


def cmd_snapshot(session, args: ArgReader) -> None:
    """``SNAPSHOT [name]`` (``SNAP``) -- save the session as it stands, to
    pick it up later with ``pyrump <name>_fit.xeq`` or ``XEQ <name>_fit``.

    Writes ``<name>_fit.xeq`` (the restore macro), ``<name>.lcm``,
    ``<name>.pert``, ``<name>.png`` and appends to ``<name>.report``. The name
    defaults to the sample's, as REPORT names it, and may carry a folder
    (``SNAP results/MA8410``). Needs the measured spectrum, read from a file,
    and a SIM sample. Works the same after a GO or after tuning by hand: the
    report says which, with uncertainties only for an unchanged fit. A
    pyRUMP addition; can be typed at the SIM, PERT and PIXE prompts too.

    ``SNAPSHOT -restored`` ends each restore macro: it writes nothing, and
    only marks the session as saved, so QUIT doesn't warn.
    """
    token = args.optional()
    args.done()
    if token is not None and len(token) > 1 and "-restored".startswith(token.lower()):
        session.saved_state = fingerprint(session)
        print("  session restored from its snapshot")
        return
    take(session, token)
