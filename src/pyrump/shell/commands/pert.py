"""PERT: the fitting sub-processor.

PERT selects what may vary and over which channels, then ``GO`` runs the search
(pert.htm). Everything it needs already exists in :mod:`pyrump.fit`: the
parameter constructors, the Poisson objective, the error and normalisation
windows, and the Levenberg-Marquardt driver. This module is the interactive
front end to those, plus the write-back that puts fitted values into the SIM
sample description so they survive into ``SIM SAVE``.

Divergences from the original, both noted in the plan:

* RUMP requires the data in buffer 1; here the ACTIVE buffer is used.
* ``MULTI`` is the default, because :func:`pyrump.fit.lm.fit` is a
  simultaneous least-squares solve. ``SINGLE`` is emulated by fitting one
  parameter at a time, in the order they were selected.
* ``GET``/``SAVE`` round-trip a ``.pert`` file, as in the original -- but as
  plain PERT commands (``WINDOW``, ``THICKNESS``, ...) rather than RUMP's own
  bounded-``VARY`` format (``pert.c``'s ``PertWriteParms``): pyRUMP has no
  generic ``VARY`` verb, so a search bound is instead an optional trailing
  ``<min> <max>`` on the selecting command itself (``THICKNESS 1 100 500``),
  which round-trips the same way GET/SAVE always have.
"""

from __future__ import annotations

import time
import warnings
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

from ...fit.parameters import (
    FitInputs,
    composition,
    equation_parameter,
    parameter,
    thickness,
)
from ...fit.windows import MAX_ERROR_WINDOWS, Window, WindowSet
from ..dispatch import ArgReader, CommandError, CommandTable
from .rump import Return, cmd_compare
from .sim import describe as _describe_sample
from .sim import editor_for


def _format_windows(windows: list[Window]) -> str:
    """One-line, 1-based ``[n] lo-hi`` rendering of the error windows."""
    if not windows:
        return "(none -- the whole spectrum)"
    return "  ".join(f"[{i}] {w.low}-{w.high}" for i, w in enumerate(windows, start=1))


def _format_varying(varying: list[Vary]) -> list[str]:
    """Lines for the 'varying:' block, one per parameter, 1-based ``[n]``."""
    if not varying:
        return ["  varying:    (nothing selected)"]
    lines = ["  varying:"]
    for i, v in enumerate(varying, start=1):
        bound = f"  bounds {v.bounds[0]:g}-{v.bounds[1]:g}" if v.bounds else ""
        lines.append(f"    [{i}] {v.name}{bound}")
    return lines


#: The few simple-parameter rump-names whose command word differs (pert.py's
#: OFFSET adds ``kev(0)``, the only one where the two aren't the same string).
_SIMPLE_PARAMETER_COMMANDS = {"kev(0)": "offset"}


def _vary_command(entry: Vary) -> str:
    """The PERT command line that would recreate this one selection."""
    if entry.kind == "thickness":
        base = f"thickness {entry.layer + 1}"
    elif entry.kind == "composition":
        base = f"composition {entry.layer + 1} {entry.symbol}"
    elif entry.kind == "species":
        base = f"species {entry.layer + 1} {entry.symbol}"
    elif entry.kind == "equation":
        base = f"equation {entry.layer + 1} {entry.index + 1}"
    else:
        base = _SIMPLE_PARAMETER_COMMANDS.get(entry.name, entry.name)
    if entry.bounds is not None:
        base += f" {entry.bounds[0]:g} {entry.bounds[1]:g}"
    return base


def _to_lines(state: PertState) -> list[str]:
    """Command lines that recreate ``state``, for PERT SAVE/GET."""
    lines = [f"window {w.low} {w.high}" for w in state.windows.error]
    norm = state.windows.normalisation
    if norm is not None:
        lines.append(f"normalize {norm.low} {norm.high}")
    if not state.multi:
        lines.append("single")
    lines.extend(_vary_command(v) for v in state.varying)
    return lines


@dataclass(slots=True)
class Vary:
    """One selected parameter, with what it maps back onto.

    ``name`` is what PERT prints and echoes back to the user, so it uses the
    same 1-based layer numbering the user typed (``THICK 1`` selects layer 1).
    ``parameter.name`` is the fit engine's own identity for the same
    parameter (:mod:`pyrump.fit.parameters`), which indexes layers the
    Python-native 0-based way and is what ``FitResult.parameters``/
    ``uncertainties`` are keyed by -- the two are deliberately not the same
    string.

    ``bounds``, if given, is the user's own ``<min> <max>`` from the
    selecting command -- kept separately from ``parameter.lower``/``.upper``
    (which the solver actually reads, and which already default to a
    physically sensible range for some parameters, e.g. thickness >= 0) so
    that ``PARMS``/``SAVE`` only echo a bound the user actually typed, not
    every parameter's built-in default range.
    """

    parameter: object
    kind: str                 # thickness | composition | equation | simple | sample
    layer: int = -1
    index: int = -1
    symbol: str = ""
    name: str = ""
    bounds: tuple[float, float] | None = None


@dataclass(slots=True)
class PertState:
    """What PERT has been told so far."""

    varying: list[Vary] = field(default_factory=list)
    windows: WindowSet = field(default_factory=WindowSet)
    multi: bool = True
    verbose: bool = False

    def describe(self) -> str:
        lines = [f"  mode        {'multiple' if self.multi else 'single'} variable"]
        lines.append(f"  error win   {_format_windows(self.windows.error)}")
        norm = self.windows.normalisation
        lines.append(
            f"  norm win    {f'{norm.low}-{norm.high}' if norm else '(none)'}"
        )
        lines.extend(_format_varying(self.varying))
        return "\n".join(lines)


def state_for(session) -> PertState:
    if session.pert is None:
        session.pert = PertState()
    return session.pert


def layers_at_risk(session, from_index: int) -> list[str]:
    """Names of PERT selections whose layer is at or after ``from_index``.

    Inserting or removing a SIM layer at ``from_index`` shifts every later
    layer's position, which silently misdirects any such selection: the fit
    parameters it built (:mod:`pyrump.fit.parameters`'s ``thickness()``/
    ``composition()``/``equation_parameter()``) close over a plain integer
    layer index, not the layer's identity, so they now read/write whatever
    layer ended up at that index. Used by SIM's DELETE/OPEN to warn, not
    block -- callers decide what to do with the names.
    """
    state = session.pert
    if state is None:
        return []
    return [v.name for v in state.varying if v.layer >= from_index]


def _layer_argument(session, args: ArgReader) -> int:
    """A 1-based layer number, validated against the sample."""
    number = args.integer("a layer number")
    layers = len(session.script.layers)
    if not 1 <= number <= layers:
        raise CommandError(f"layer {number} is outside 1-{layers}")
    return number - 1


def _optional_bounds(args: ArgReader) -> tuple[float, float] | None:
    """Parse an optional trailing ``<min> <max>`` search bound.

    Matches the original's own ``VARY`` grammar (pert.c's
    ``PertGetVariable``): give neither to leave the parameter at its default
    range, or both to constrain the search to it. There is no way to give
    just one -- as in the original, min and max travel together.
    """
    if not args:
        return None
    low = args.number("a minimum bound")
    high = args.number("a maximum bound")
    if high <= low:
        raise CommandError(f"empty bound: {low} to {high}")
    return low, high


def _add(session, entry: Vary) -> None:
    state = state_for(session)
    if any(v.name == entry.name for v in state.varying):
        raise CommandError(f"{entry.name} is already being varied")
    state.varying.append(entry)
    print(f"  varying {entry.name}")


# ---------------------------------------------------------------------------
# Selecting parameters
# ---------------------------------------------------------------------------


def cmd_thickness(session, args: ArgReader) -> None:
    layer = _layer_argument(session, args)
    bound = _optional_bounds(args)
    args.done()
    param = thickness(layer)
    if bound is not None:
        param = replace(param, lower=bound[0], upper=bound[1])
    _add(
        session,
        Vary(
            parameter=param,
            kind="thickness",
            layer=layer,
            name=f"layer {layer + 1} thickness",
            bounds=bound,
        ),
    )


def _element_index(session, symbol: str) -> int:
    elements = session.script.elements
    for index, name in enumerate(elements):
        if name.lower() == symbol.lower():
            return index
    raise CommandError(
        f"{symbol} is not in the sample; it has: {' '.join(elements) or '(nothing)'}"
    )


def _layer_element(session, layer: int, symbol: str, kind: str) -> str:
    """Confirm ``symbol`` is declared in *this* layer's composition/species,
    not just somewhere in the sample, and return it in its stored casing.
    """
    layer_obj = session.script.layers[layer]
    table = layer_obj.species if kind == "species" else layer_obj.composition
    for name in table:
        if name.lower() == symbol.lower():
            return name
    raise CommandError(
        f"{symbol} is not part of layer {layer + 1}'s {kind} "
        f"(it has: {', '.join(table) or '(nothing)'}) -- check your input"
    )


def cmd_composition(session, args: ArgReader) -> None:
    layer = _layer_argument(session, args)
    symbol = args.token("an element symbol")
    bound = _optional_bounds(args)
    args.done()
    index = _element_index(session, symbol)
    canonical = _layer_element(session, layer, symbol, "composition")
    param = composition(layer, index)
    if bound is not None:
        param = replace(param, lower=bound[0], upper=bound[1])
    _add(
        session,
        Vary(
            parameter=param,
            kind="composition",
            layer=layer,
            index=index,
            symbol=canonical,
            name=f"layer {layer + 1} composition {canonical}",
            bounds=bound,
        ),
    )


def cmd_species(session, args: ArgReader) -> None:
    layer = _layer_argument(session, args)
    symbol = args.token("an element symbol")
    bound = _optional_bounds(args)
    args.done()
    index = _element_index(session, symbol)
    canonical = _layer_element(session, layer, symbol, "species")
    param = composition(layer, index)
    if bound is not None:
        param = replace(param, lower=bound[0], upper=bound[1])
    _add(
        session,
        Vary(
            parameter=param,
            kind="species",
            layer=layer,
            index=index,
            symbol=canonical,
            name=f"layer {layer + 1} species {canonical}",
            bounds=bound,
        ),
    )


def cmd_equation(session, args: ArgReader) -> None:
    layer = _layer_argument(session, args)
    index = args.integer("a parameter number") - 1
    bound = _optional_bounds(args)
    args.done()
    profile = session.script.layers[layer].profile
    if profile is None:
        raise CommandError(f"layer {layer + 1} has no equation")
    if not 0 <= index < len(profile.parameters):
        raise CommandError(
            f"equation parameter {index + 1} is outside 1-{len(profile.parameters)}"
        )
    param = equation_parameter(layer, index)
    if bound is not None:
        param = replace(param, lower=bound[0], upper=bound[1])
    _add(
        session,
        Vary(
            parameter=param,
            kind="equation",
            layer=layer,
            index=index,
            name=f"layer {layer + 1} equation parameter {index + 1}",
            bounds=bound,
        ),
    )


def _simple(rump_name: str, kind: str = "simple"):
    """A parameter from SIMPLE_PARAMETERS, with an optional trailing bound."""

    def handler(session, args: ArgReader) -> None:
        bound = _optional_bounds(args)
        args.done()
        try:
            param = parameter(rump_name)
        except KeyError as error:
            raise CommandError(str(error)) from None
        if bound is not None:
            param = replace(param, lower=bound[0], upper=bound[1])
        _add(session, Vary(parameter=param, kind=kind, name=rump_name, bounds=bound))

    return handler


def cmd_fuzz(session, args: ArgReader) -> None:
    args.done()
    raise CommandError(
        "varying FUZZ is not implemented: pyrump.fit.parameters has no fuzz "
        "parameter yet, and fuzz changes the number of simulated replicas "
        "rather than a continuous value"
    )


# ---------------------------------------------------------------------------
# Windows and modes
# ---------------------------------------------------------------------------


def cmd_window(session, args: ArgReader) -> None:
    state = state_for(session)
    if not args:
        print(state.describe())
        return
    token = args.peek()
    if token is not None and token.lower() in ("clear", "none", "reset"):
        args.token()
        args.done()
        state.windows.error = []
        print("  error windows cleared")
        return
    if token is not None and token.lower() == "remove":
        args.token()
        n = args.integer("a window number")
        args.done()
        windows = state.windows.error
        if not windows:
            raise CommandError("no error windows are set")
        if not 1 <= n <= len(windows):
            raise CommandError(f"window {n} is outside 1-{len(windows)}")
        del windows[n - 1]
        print(f"  error windows {_format_windows(windows)}")
        return
    low = args.integer("the first channel")
    high = args.integer("the last channel")
    args.done()
    if high <= low:
        raise CommandError(f"empty window: {low} to {high}")
    if len(state.windows.error) >= MAX_ERROR_WINDOWS:
        raise CommandError(f"at most {MAX_ERROR_WINDOWS} error windows")
    state.windows.error.append(Window(low, high))
    print(f"  error windows {_format_windows(state.windows.error)}")


def cmd_normalize(session, args: ArgReader) -> None:
    state = state_for(session)
    if not args:
        norm = state.windows.normalisation
        print(f"  normalisation window {f'{norm.low}-{norm.high}' if norm else '(none)'}")
        return
    token = args.peek()
    if token is not None and token.lower() in ("clear", "none", "off"):
        args.token()
        args.done()
        state.windows.normalisation = None
        print("  normalisation window cleared")
        return
    low = args.integer("the first channel")
    high = args.integer("the last channel")
    args.done()
    if high <= low:
        raise CommandError(f"empty window: {low} to {high}")
    state.windows.normalisation = Window(low, high)
    print(f"  normalisation window {low}-{high}")


def cmd_single(session, args: ArgReader) -> None:
    args.done()
    state_for(session).multi = False
    print("  single-variable mode")


def cmd_multi(session, args: ArgReader) -> None:
    args.done()
    state_for(session).multi = True
    print("  multiple-variable mode")


def cmd_parms(session, args: ArgReader) -> None:
    args.done()
    print(state_for(session).describe())


def cmd_show(session, args: ArgReader) -> None:
    """Display the sample description, as SIM SHOW does.

    RUMP's own ``pert.c`` has this too (``PE_SHOW``, calling the same
    ``SimShowSample``) -- undocumented there; documented here since it's the
    natural way to check what a layer actually contains before COMPOSITION
    or SPECIES.
    """
    args.done()
    print(_describe_sample(session, editor_for(session)))


def cmd_get(session, args: ArgReader) -> None:
    """Replay a saved PERT selection from a ``.pert`` file, replacing this one.

    A trailing ``GO`` runs the fit immediately after loading, so ``PERT GET
    usual.pert GO`` works as a single line (from the RUMP prompt or in a
    macro) -- the common "load my usual setup and fit" idiom.
    """
    from ..repl import execute_file

    path = Path(args.token("a .pert file"))
    run_go = False
    token = args.peek()
    if token is not None and token.lower() == "go":
        args.token()
        run_go = True
    args.done()
    if not path.suffix:
        path = path.with_suffix(".pert")
    if not path.exists():
        raise CommandError(f"no such file: {path}")
    session.pert = PertState()
    execute_file(session, path, stack=["rump", "pert"])
    print(f"read {path}")
    print(state_for(session).describe())
    if run_go:
        cmd_go(session, ArgReader([], command="go"))


def cmd_save(session, args: ArgReader) -> None:
    """Write the current selection to a ``.pert`` file, for GET to replay later."""
    path = Path(args.token("an output .pert file"))
    args.done()
    if not path.suffix:
        path = path.with_suffix(".pert")
    lines = _to_lines(state_for(session))
    path.write_text("\n".join(lines) + ("\n" if lines else ""))
    print(f"wrote {path}")


def cmd_clear(session, args: ArgReader) -> None:
    if not args:
        session.pert = PertState()
        print("  PERT settings cleared")
        return
    n = args.integer("a parameter number")
    args.done()
    state = state_for(session)
    varying = state.varying
    if not varying:
        raise CommandError("nothing selected to clear")
    if not 1 <= n <= len(varying):
        raise CommandError(f"parameter {n} is outside 1-{len(varying)}")
    del varying[n - 1]
    print("\n".join(_format_varying(varying)))


def cmd_volume(session, args: ArgReader) -> None:
    """``VOLUME [off]`` -- print a line for every model evaluation during
    ``GO`` (default off).

    Each evaluation is a full simulation, so a fit with several varying
    parameters can run for seconds with nothing else printed in between --
    this is the difference between that and an apparently hung prompt.
    """
    token = args.optional()
    args.done()
    state = state_for(session)
    state.verbose = token is None or token.lower() not in ("off", "no", "0")
    print(f"  messages {'on' if state.verbose else 'off'}")


def cmd_help(session, args: ArgReader) -> None:
    """``HELP`` lists the PERT commands; ``HELP <name>`` describes one.

    A name not in PERT's own table falls through to RUMP and then the system
    tier, mirroring how an unrecognised command escapes this mode
    (repl.py's ``execute_line``).
    """
    topic = args.optional()
    args.done()
    if topic is None:
        print(TABLE.help_text())
        return
    from .rump import TABLE as RUMP_TABLE
    from .system import TABLE as SYSTEM_TABLE

    text = TABLE.describe(topic) or RUMP_TABLE.describe(topic) or SYSTEM_TABLE.describe(topic)
    if text is None:
        raise CommandError(f"no help for {topic!r} -- try HELP with no argument")
    print(text)


def cmd_return(session, args: ArgReader) -> None:
    raise Return()


def execute_in_pert(session, args: ArgReader) -> None:
    """Run a one-shot ``PERT <command>`` from the RUMP level, as SIM does.

    ``PERT GET usual.pert GO`` reaches here as ``args`` = ``["GET",
    "usual.pert", "GO"]``: this dispatches only the first word (``GET``,
    with the rest as its own args) -- GET's own handler is what understands
    the trailing ``GO`` and runs the fit after loading.
    """
    name = args.token("a PERT command")
    command = TABLE.match(name)
    if command is None:
        raise CommandError(f"unrecognized PERT command: {name}")
    command.handler(session, ArgReader(args.remaining, command=command.name.lower()))
    args.index = len(args.tokens)


# ---------------------------------------------------------------------------
# GO
# ---------------------------------------------------------------------------


def _write_back(session, entry: Vary, inputs: FitInputs, before: float) -> None:
    """Put one fitted value back where the user can see and save it."""
    layers = session.script.layers
    value = entry.parameter.get(inputs)

    if entry.kind == "thickness":
        # The script keeps a magnitude and a unit ("151 ITO"); the sample keeps
        # areal density. SimThickConvert is linear in the magnitude, so scaling
        # by the ratio is exact whatever the unit.
        if before > 0:
            layers[entry.layer].thickness *= value / before
    elif entry.kind in ("composition", "species"):
        target = (
            layers[entry.layer].species
            if entry.kind == "species"
            else layers[entry.layer].composition
        )
        symbol = entry.symbol or session.script.elements[entry.index]
        target[symbol] = value
    elif entry.kind == "equation":
        profile = layers[entry.layer].profile
        params = list(profile.parameters)
        params[entry.index] = value
        layers[entry.layer].profile = replace(profile, parameters=tuple(params))
    elif entry.kind == "sample":
        setattr(session.script, entry.name, value)
    else:
        # A buffer parameter. fit() mutates ``inputs`` in place and leaves it
        # holding the best-fit objects, so the buffer just adopts them.
        buffer = session.buffers.require_active()
        buffer.beam = inputs.beam
        buffer.geometry = inputs.geometry
        buffer.measurement = inputs.measurement
        buffer.spectrum.calibration = inputs.calibration


def cmd_go(session, args: ArgReader) -> None:
    args.done()
    from ...fit.lm import fit
    from ...script.lcm import thickness_label, to_sample
    from ...sim.engine import simulate

    state = state_for(session)
    if not state.varying:
        raise CommandError("nothing selected to vary")
    if not session.script.layers:
        raise CommandError("no sample described: use SIM to build one")

    data_buffer = session.buffers.require_active()
    observed = np.asarray(data_buffer.spectrum.counts, dtype=float)

    sample = to_sample(session.script, session.table, session.densities)
    inputs = FitInputs(
        sample=sample,
        beam=data_buffer.beam,
        geometry=data_buffer.geometry,
        calibration=data_buffer.calibration,
        measurement=data_buffer.measurement,
    )

    def run(current: FitInputs) -> np.ndarray:
        return simulate(
            current.sample,
            current.beam,
            current.geometry,
            session.registry,
            session.table,
            current.calibration,
            current.measurement,
        ).counts

    starting = {v.name: v.parameter.get(inputs) for v in state.varying}
    groups = (
        [state.varying] if state.multi else [[v] for v in state.varying]
    )

    def _progress(evaluation: int, reduced: float) -> None:
        print(f"    eval {evaluation:3d}   chi2/dof {reduced:.4f}")

    result = None
    started = time.perf_counter()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("always", RuntimeWarning)
            for group in groups:
                result = fit(
                    run,
                    observed,
                    inputs,
                    [v.parameter for v in group],
                    windows=state.windows,
                    progress=_progress if state.verbose else None,
                )
                if not state.multi:
                    print(
                        f"  {group[0].name}: {result.parameters[group[0].parameter.name]:.6g}"
                        f"   chi2/dof {result.reduced_chi_square:.4f}"
                    )
    except ValueError as error:
        raise CommandError(f"go: {error}") from None
    elapsed = time.perf_counter() - started

    for entry in state.varying:
        _write_back(session, entry, inputs, starting[entry.name])
    session.editor = None
    session.touch()

    print(f"  fit took {elapsed:.2f} s")
    print(f"\n  reduced chi-square {result.reduced_chi_square:.4f} on {result.dof} dof")
    print(f"  {result.n_evaluations} evaluations, {result.message}")
    if result.normalisation != 1.0:
        print(f"  data scaled by {result.normalisation:.5f} over the norm window")
    if result.n_invalid:
        print(
            f"  warning: {result.n_invalid} windowed channels had "
            "no predicted counts"
        )
    for entry in state.varying:
        name = entry.parameter.name
        value = entry.parameter.get(inputs)
        sigma = result.uncertainties.get(name)
        before = starting[entry.name]
        if entry.kind == "thickness":
            value_text, before_text = thickness_label(value), thickness_label(before)
        else:
            value_text, before_text = f"{value:.6g}", f"{before:.6g}"
        line = f"  {entry.name:26s} {value_text:>14s}"
        if sigma:
            line += f"  +/- {sigma:.4g}"
        print(line + f"   (was {before_text})")
    sample_id = data_buffer.identifier or data_buffer.name or "(unnamed)"
    print(f"\n  sample id      {sample_id}")


TABLE = CommandTable("PERT Commands")

_ENTRIES: list[tuple[str, int, object, str]] = [
    ("?", 1, cmd_help, "list the PERT commands"),
    ("HELP", 1, cmd_help, "list the PERT commands"),
    ("RETURN", 1, cmd_return, "return to the RUMP level"),
    ("QUIT", 1, cmd_return, "return to the RUMP level (not exit pyRUMP)"),
    ("GO", 2, cmd_go, "run the search"),
    ("PARMS", 2, cmd_parms, "display the current settings"),
    ("SHOW", 2, cmd_show, "display the sample description (same as SIM SHOW)"),
    ("GET", 2, cmd_get, "replay a saved PERT selection from a .pert file, or GET <file> GO"),
    ("SAVE", 2, cmd_save, "save the current PERT selection to a .pert file"),
    ("CLEAR", 2, cmd_clear, "forget every selected parameter and window, or CLEAR <n> one parameter"),
    # Windows and mode
    ("WINDOW", 2, cmd_window, "set an error window in channels, or WINDOW REMOVE <n>"),
    ("NORMALIZE", 2, cmd_normalize, "set the normalisation window"),
    ("SINGLE", 2, cmd_single, "vary one parameter at a time"),
    ("MULTI", 3, cmd_multi, "vary all parameters together (default)"),
    ("VOLUME", 3, cmd_volume, "verbose progress messages"),
    # Parameters -- all take an optional trailing "<min> <max>" search bound
    ("THICKNESS", 2, cmd_thickness, "vary a layer thickness, e.g. THICKNESS <layer> [<min> <max>]"),
    ("COMPOSITION", 3, cmd_composition, "vary an element in a layer [<min> <max>]"),
    ("SPECIES", 2, cmd_species, "vary the species composition [<min> <max>]"),
    ("EQUATION", 2, cmd_equation, "vary an equation parameter [<min> <max>]"),
    ("MEV", 3, _simple("mev"), "vary the beam energy [<min> <max>]"),
    ("FWHM", 2, _simple("fwhm"), "vary the detector resolution [<min> <max>]"),
    ("STRAGGLE", 5, _simple("straggle", "sample"), "vary the straggling constant [<min> <max>]"),
    ("FUZZ", 2, cmd_fuzz, "vary the fuzz parameter (not implemented)"),
    ("CORRECTION", 3, _simple("correction"), "vary the normalization correction [<min> <max>]"),
    ("THETA", 4, _simple("theta"), "vary the sample tilt [<min> <max>]"),
    ("OFFSET", 3, _simple("kev(0)"),
     "vary the calibration energy offset (e.g. a sample-charging shift) [<min> <max>]"),
    ("COMPARE", 0, cmd_compare, "plot the active buffer against the simulation"),
    ("CMP", -3, cmd_compare, "synonym for COMPARE"),
]

for _name, _minlen, _handler, _help in _ENTRIES:
    TABLE.add(_name, _minlen, _handler, _help)
