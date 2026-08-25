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
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

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
from .rump import Return


@dataclass(slots=True)
class Vary:
    """One selected parameter, with what it maps back onto."""

    parameter: object
    kind: str                 # thickness | composition | equation | simple | sample
    layer: int = -1
    index: int = -1
    symbol: str = ""
    name: str = ""


@dataclass(slots=True)
class PertState:
    """What PERT has been told so far."""

    varying: list[Vary] = field(default_factory=list)
    windows: WindowSet = field(default_factory=WindowSet)
    multi: bool = True
    verbose: bool = False

    def describe(self) -> str:
        lines = [f"  mode        {'multiple' if self.multi else 'single'} variable"]
        if self.windows.error:
            spans = ", ".join(f"{w.low}-{w.high}" for w in self.windows.error)
        else:
            spans = "(none -- the whole spectrum)"
        lines.append(f"  error win   {spans}")
        norm = self.windows.normalisation
        lines.append(
            f"  norm win    {f'{norm.low}-{norm.high}' if norm else '(none)'}"
        )
        if self.varying:
            lines.append("  varying:")
            lines.extend(f"    {v.name}" for v in self.varying)
        else:
            lines.append("  varying:    (nothing selected)")
        return "\n".join(lines)


def state_for(session) -> PertState:
    if session.pert is None:
        session.pert = PertState()
    return session.pert


def _layer_argument(session, args: ArgReader) -> int:
    """A 1-based layer number, validated against the sample."""
    number = args.integer("a layer number")
    layers = len(session.script.layers)
    if not 1 <= number <= layers:
        raise CommandError(f"layer {number} is outside 1-{layers}")
    return number - 1


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
    args.done()
    _add(
        session,
        Vary(
            parameter=thickness(layer),
            kind="thickness",
            layer=layer,
            name=f"thickness[{layer}]",
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


def cmd_composition(session, args: ArgReader) -> None:
    layer = _layer_argument(session, args)
    symbol = args.token("an element symbol")
    args.done()
    index = _element_index(session, symbol)
    _add(
        session,
        Vary(
            parameter=composition(layer, index),
            kind="composition",
            layer=layer,
            index=index,
            symbol=symbol,
            name=f"composition[{layer},{index}]",
        ),
    )


def cmd_species(session, args: ArgReader) -> None:
    layer = _layer_argument(session, args)
    symbol = args.token("an element symbol")
    args.done()
    index = _element_index(session, symbol)
    _add(
        session,
        Vary(
            parameter=composition(layer, index),
            kind="species",
            layer=layer,
            index=index,
            symbol=symbol,
            name=f"species[{layer},{symbol}]",
        ),
    )


def cmd_equation(session, args: ArgReader) -> None:
    layer = _layer_argument(session, args)
    index = args.integer("a parameter number") - 1
    args.done()
    profile = session.script.layers[layer].profile
    if profile is None:
        raise CommandError(f"layer {layer + 1} has no equation")
    if not 0 <= index < len(profile.parameters):
        raise CommandError(
            f"equation parameter {index + 1} is outside 1-{len(profile.parameters)}"
        )
    _add(
        session,
        Vary(
            parameter=equation_parameter(layer, index),
            kind="equation",
            layer=layer,
            index=index,
            name=f"equation[{layer},{index}]",
        ),
    )


def _simple(rump_name: str, kind: str = "simple"):
    """A no-argument parameter from SIMPLE_PARAMETERS."""

    def handler(session, args: ArgReader) -> None:
        args.done()
        try:
            param = parameter(rump_name)
        except KeyError as error:
            raise CommandError(str(error)) from None
        _add(session, Vary(parameter=param, kind=kind, name=rump_name))

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
    low = args.integer("the first channel")
    high = args.integer("the last channel")
    args.done()
    if high <= low:
        raise CommandError(f"empty window: {low} to {high}")
    if len(state.windows.error) >= MAX_ERROR_WINDOWS:
        raise CommandError(f"at most {MAX_ERROR_WINDOWS} error windows")
    state.windows.error.append(Window(low, high))
    print(f"  error window {low}-{high}")


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


def cmd_clear(session, args: ArgReader) -> None:
    args.done()
    session.pert = PertState()
    print("  PERT settings cleared")


def cmd_volume(session, args: ArgReader) -> None:
    token = args.optional()
    args.done()
    state = state_for(session)
    state.verbose = token is None or token.lower() not in ("off", "no", "0")
    print(f"  messages {'on' if state.verbose else 'off'}")


def cmd_help(session, args: ArgReader) -> None:
    args.done()
    print(TABLE.help_text())


def cmd_return(session, args: ArgReader) -> None:
    raise Return()


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
        from dataclasses import replace as _replace

        profile = layers[entry.layer].profile
        params = list(profile.parameters)
        params[entry.index] = value
        layers[entry.layer].profile = _replace(profile, parameters=tuple(params))
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
    from ...script.lcm import to_sample
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

    result = None
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
                )
                if not state.multi:
                    print(
                        f"  {group[0].name}: {result.parameters[group[0].parameter.name]:.6g}"
                        f"   chi2/dof {result.reduced_chi_square:.4f}"
                    )
    except ValueError as error:
        raise CommandError(f"go: {error}") from None

    for entry in state.varying:
        _write_back(session, entry, inputs, starting[entry.name])
    session.editor = None
    session.touch()

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
        line = f"  {entry.name:26s} {value:14.6g}"
        if sigma:
            line += f"  +/- {sigma:.4g}"
        print(line + f"   (was {starting[entry.name]:.6g})")


TABLE = CommandTable("PERT Commands")

_ENTRIES: list[tuple[str, int, object, str]] = [
    ("?", 1, cmd_help, "list the PERT commands"),
    ("HELP", 1, cmd_help, "list the PERT commands"),
    ("RETURN", 1, cmd_return, "return to the RUMP level"),
    ("QUIT", 1, cmd_return, "return to the RUMP level (not exit pyRUMP)"),
    ("GO", 2, cmd_go, "run the search"),
    ("PARMS", 2, cmd_parms, "display the current settings"),
    ("CLEAR", 2, cmd_clear, "forget every selected parameter and window"),
    # Windows and mode
    ("WINDOW", 2, cmd_window, "set an error window in channels"),
    ("NORMALIZE", 2, cmd_normalize, "set the normalisation window"),
    ("SINGLE", 2, cmd_single, "vary one parameter at a time"),
    ("MULTI", 3, cmd_multi, "vary all parameters together (default)"),
    ("VOLUME", 3, cmd_volume, "verbose progress messages"),
    # Parameters
    ("THICKNESS", 2, cmd_thickness, "vary a layer thickness"),
    ("COMPOSITION", 3, cmd_composition, "vary an element in a layer"),
    ("SPECIES", 2, cmd_species, "vary the species composition"),
    ("EQUATION", 2, cmd_equation, "vary an equation parameter"),
    ("MEV", 3, _simple("mev"), "vary the beam energy"),
    ("FWHM", 2, _simple("fwhm"), "vary the detector resolution"),
    ("STRAGGLE", 5, _simple("straggle", "sample"), "vary the straggling constant"),
    ("FUZZ", 2, cmd_fuzz, "vary the fuzz parameter (not implemented)"),
    ("CORRECTION", 3, _simple("correction"), "vary the normalization correction"),
    ("THETA", 4, _simple("theta"), "vary the sample tilt"),
]

for _name, _minlen, _handler, _help in _ENTRIES:
    TABLE.add(_name, _minlen, _handler, _help)
