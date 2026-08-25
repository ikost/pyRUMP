"""SIM: the sample-description editor.

SIM is not the simulator -- it edits the theoretical sample, and the spectrum is
recomputed on demand by :meth:`~pyrump.shell.session.Session.simulation`
(sim.htm: "The SIM command enters a sub-process that is not the spectrum
simulator, but rather primarily an editor of the theoretical sample description
table").

The sample-definition commands are applied by
:class:`~pyrump.script.lcm.SampleEditor`, the same object that replays ``.lcm``
files, so the file format and the interactive language cannot drift apart.
"""

from __future__ import annotations

from pathlib import Path

from ...script.lcm import EQUATION_NAMES, SampleEditor, read_lcm, write_lcm
from ..dispatch import ArgReader, CommandError, CommandTable
from .rump import Return, cmd_compare


def editor_for(session) -> SampleEditor:
    """The session's editor, rebuilt if the script was replaced under it."""
    editor = session.editor
    if editor is None or editor.script is not session.script:
        editor = SampleEditor(session.script)
        session.editor = editor
    return editor


def _apply(session, verb: str, args: ArgReader) -> None:
    """Feed one command to the SampleEditor and mark the simulation stale."""
    editor = editor_for(session)
    line = " ".join([verb, *args.remaining])
    args.index = len(args.tokens)
    try:
        editor.execute(line)
    except (ValueError, IndexError) as error:
        raise CommandError(f"{verb}: {error}") from None
    session.script = editor.script
    session.touch()


def _editor_command(verb: str):
    def handler(session, args: ArgReader) -> None:
        _apply(session, verb, args)

    return handler


def cmd_layer(session, args: ArgReader) -> None:
    editor = editor_for(session)
    if args:
        number = args.integer("a layer number")
        args.done()
        if number < 1:
            raise CommandError("layers are numbered from 1")
        if number > len(editor.script.layers) + 1:
            raise CommandError(
                f"layer {number} is not defined "
                f"({len(editor.script.layers)} layers, plus one blank)"
            )
        editor.select(number - 1)
    session.touch()
    print(_where(editor))


def cmd_next(session, args: ArgReader) -> None:
    args.done()
    editor = editor_for(session)
    editor.next_layer()
    session.touch()
    print(_where(editor))


def cmd_open(session, args: ArgReader) -> None:
    args.done()
    editor = editor_for(session)
    editor.open_layer()
    session.touch()
    print("new layer opened up")
    print(_where(editor))


def _where(editor: SampleEditor) -> str:
    total = len(editor.script.layers)
    number = editor.current + 1
    if editor.layer is None:
        return f"  you are now working on a fresh layer # {number} (of {total})"
    return f"  you are now working on layer # {number} of {total}"


def cmd_reset(session, args: ArgReader) -> None:
    args.done()
    editor = editor_for(session)
    editor.reset()
    session.script = editor.script
    session.touch()
    print("sample reset to empty space")


def cmd_show(session, args: ArgReader) -> None:
    args.done()
    print(describe(session, editor_for(session)))


def describe(session, editor: SampleEditor) -> str:
    """The ``SHOW`` listing: the sample as RUMP prints it."""
    script = editor.script
    lines = []
    if script.description:
        lines.append(f"  {script.description}")
    if not script.layers:
        lines.append("  (empty space)")
    for index, layer in enumerate(script.layers):
        mark = ">" if index == editor.current else " "
        composition = " ".join(
            f"{symbol} {value:g}" for symbol, value in layer.composition.items()
        )
        lines.append(
            f" {mark}{index + 1:3d}  {layer.thickness:12g} {layer.unit:<8s} {composition}"
        )
        if layer.species:
            species = " ".join(
                f"{symbol} {value:g}" for symbol, value in layer.species.items()
            )
            lines.append(f"          species: {species}")
        if layer.profile is not None:
            values = " ".join(f"{v:g}" for v in layer.profile.parameters)
            lines.append(
                f"          equation: {layer.profile.kind.name.lower()} {values}"
            )
        if layer.sublayers:
            lines.append(f"          sublayers: {layer.sublayers}")
        if layer.fuzz_amount:
            lines.append(
                f"          fuzz: {layer.fuzz_amount:g} in {layer.fuzz_steps} steps"
            )
    if editor.layer is None and script.layers:
        lines.append(f" >{len(script.layers) + 1:3d}  (blank)")
    lines.append(
        f"  maxpth {script.maxpth:g}   straggle {script.straggle:g}"
        f"   multiple {script.multiple:g}"
        f"   absorber {script.absorber_layers}"
    )
    return "\n".join(lines)


def cmd_status(session, args: ArgReader) -> None:
    args.done()
    script = session.script
    print(
        f"  {len(script.layers)} layers, maxpth {script.maxpth:g},"
        f" straggle {script.straggle:g}, multiple {script.multiple:g}\n"
        f"  elements: {' '.join(script.elements) or '(none)'}\n"
        f"  simulation is {'stale' if session.dirty else 'current'}"
    )


def cmd_get(session, args: ArgReader) -> None:
    """Read a sample description from a ``.lcm`` file."""
    path = Path(args.token("a .lcm file"))
    args.done()
    if not path.exists() and not path.suffix:
        path = path.with_suffix(".lcm")
    if not path.exists():
        raise CommandError(f"no such file: {path}")
    try:
        session.script = read_lcm(path)
    except (ValueError, OSError) as error:
        raise CommandError(f"could not read {path}: {error}") from None
    session.editor = SampleEditor(session.script)
    session.touch()
    print(f"read {path}: {len(session.script.layers)} layers")


def cmd_save(session, args: ArgReader) -> None:
    """Write the sample description, in RUMP's own format."""
    path = Path(args.token("an output .lcm file"))
    args.done()
    if not path.suffix:
        path = path.with_suffix(".lcm")
    path.write_text(write_lcm(session.script))
    print(f"wrote {path}")


def cmd_density(session, args: ArgReader) -> None:
    """List the thickness units and compound densities RUMP knows.

    Both share one namespace, as they do in the C (``DensityTable``): a Thick
    unit may be a length, an areal density, or any compound in density.tab.
    """
    pattern = args.optional()
    args.done()
    densities = session.densities

    units = [unit.name for unit in densities.units]
    print("  units:     " + "  ".join(units))

    names = sorted(densities.compounds)
    if pattern:
        needle = pattern.upper()
        names = [name for name in names if needle in name]
        if not names:
            # Worth saying plainly: an unknown Thick unit is not an error, it
            # silently uses silicon's density (atomic/density.py, SILICON_DENSITY).
            # That is why "Thick 151 ITO" gives 75.5e15 at/cm^2 in RUMP too.
            print(
                f"  no compound matching {pattern!r} -- as a Thick unit it would"
                f" fall back to silicon, 0.4997e23 atoms/cm^3"
            )
            return
    print(f"  compounds: {len(names)} in density.tab (1e23 atoms/cm^3)")
    for start in range(0, len(names), 4):
        row = names[start : start + 4]
        print(
            "    "
            + "  ".join(
                f"{name:<12s}{densities.density(name):7.4f}" for name in row
            ).rstrip()
        )


def cmd_splot(session, args: ArgReader) -> None:
    """Overlay the simulation on the current plot.

    Replaces any simulation trace already on the plot rather than stacking a
    new one, so repeated SPL calls while tweaking the sample update the
    overlay in place instead of piling up copies.
    """
    args.done()
    from .. import plotting

    buffer = session.simulation()
    plotting.add_trace(session, 0, buffer, clear=False, replace=True)
    plotting.draw(session)


def cmd_equation_help(session, args: ArgReader) -> None:
    args.done()
    print("  known equations:")
    names = sorted(set(EQUATION_NAMES))
    for start in range(0, len(names), 5):
        print("    " + "  ".join(n.ljust(14) for n in names[start : start + 5]).rstrip())


def cmd_help(session, args: ArgReader) -> None:
    args.done()
    print(TABLE.help_text(faithful=session.settings.faithful))


def cmd_return(session, args: ArgReader) -> None:
    raise Return()


def cmd_abort(session, args: ArgReader) -> None:
    raise Return()


def execute_in_sim(session, args: ArgReader) -> None:
    """Run a one-shot ``SIM <command>`` from the RUMP level."""
    name = args.token("a SIM command")
    command = TABLE.match(name, faithful=session.settings.faithful)
    if command is None:
        raise CommandError(f"unrecognized SIM command: {name}")
    command.handler(session, ArgReader(args.remaining, command=command.name.lower()))
    args.index = len(args.tokens)


TABLE = CommandTable("SIM Commands")

_ENTRIES: list[tuple[str, int, object, str]] = [
    ("?", 1, cmd_help, "list the SIM commands"),
    ("HELP", 2, cmd_help, "list the SIM commands"),
    ("RETURN", 3, cmd_return, "return to the RUMP level"),
    ("ABORT", 5, cmd_abort, "return to the RUMP level"),
    ("QUIT", 1, cmd_return, "return to the RUMP level (not exit pyRUMP)"),
    # Layer navigation
    ("LAYER", 2, cmd_layer, "move to a layer by number"),
    ("NEXT", 2, cmd_next, "move to the next layer"),
    ("OPEN", 2, cmd_open, "insert a blank layer above this one"),
    ("RESET", 5, cmd_reset, "reset the sample to empty space"),
    ("SHOW", 2, cmd_show, "display the sample"),
    ("STATUS", 2, cmd_status, "summarise the SIM parameters"),
    # Layer contents
    ("THICKNESS", 2, _editor_command("thick"), "set this layer's thickness"),
    ("COMPOSITION", 1, _editor_command("composition"),
     "set this layer's composition, e.g. In 2 O 3 /"),
    ("SPECIES", 2, _editor_command("species"), "impurity species for EQUATION"),
    ("EQUATION", 2, _editor_command("equation"), "impurity distribution equation"),
    ("EQLIST", 3, cmd_equation_help, "list the known equations"),
    ("FUZZ", 2, _editor_command("fuzz"), "fuzz an interface"),
    ("SUBLAYER", 3, _editor_command("sublayer"), "sublayers in this layer"),
    ("STHICKNESS", 3, _editor_command("sthick"), "thickness of each sublayer"),
    # Global sample parameters
    ("MAXPTH", 3, _editor_command("maxpth"), "maximum internal layer thickness"),
    ("STRAGGLE", 4, _editor_command("straggle"), "Bohr straggling multiplier"),
    ("ABSORBER", 3, _editor_command("absorber"), "stopper-foil layer count"),
    ("MULTIPLE", 3, _editor_command("multiple"), "multiple-scattering amount"),
    # Files and plotting
    ("GET", 3, cmd_get, "read a sample description from a file"),
    ("SAVE", 2, cmd_save, "write the sample description to a file"),
    ("DENSITY", 2, cmd_density, "list the known thickness units"),
    ("SPLOT", 3, cmd_splot, "overlay the simulation on the plot"),
    ("COMPARE", 4, cmd_compare, "plot the active buffer against the simulation"),
]

for _name, _minlen, _handler, _help in _ENTRIES:
    TABLE.add(_name, _minlen, _handler, _help)
