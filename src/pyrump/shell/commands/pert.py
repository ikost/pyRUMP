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
    atoms,
    composition,
    equation_parameter,
    parameter,
    pixe_h,
    thickness,
)
from ...fit.windows import MAX_ERROR_WINDOWS, Window, WindowSet
from ...script.lcm import thickness_label, thickness_mode_views
from ..dispatch import ArgReader, CommandError, CommandTable
from .. import plotting, snapshot
from .rump import (
    Return,
    cmd_compare,
    cmd_export,
    cmd_exportcmp,
    cmd_mode,
    describe_topic,
    needs_mode,
)
from .sim import describe as _describe_sample
from .sim import editor_for


def _format_windows(windows: list[Window], empty: str = "(none -- the whole spectrum)") -> str:
    """One-line, 1-based ``[n] lo-hi`` rendering of the error windows."""
    if not windows:
        return empty
    return "  ".join(f"[{i}] {w.low}-{w.high}" for i, w in enumerate(windows, start=1))


def _format_pixe_windows(windows: list[Window]) -> str:
    return _format_windows(windows, "(none -- RBS only)")


#: RUMP's areal-density unit, 1e15 at/cm^2, spelt the way SIM SHOW's brackets
#: spell it. The fit engine's THICKNESS and ATOMS parameters both work in it
#: (:mod:`pyrump.fit.parameters`).
_AREAL = "/CM2"


def _bound_unit(entry: Vary, script) -> str:
    """The unit ``entry``'s bound is typed in, or "" if it has none worth naming."""
    if entry.kind == "thickness" and script is not None and entry.layer < len(script.layers):
        return script.layers[entry.layer].unit
    if entry.kind == "atoms":
        return _AREAL
    return ""


def _format_varying(varying: list[Vary], script=None) -> list[str]:
    """Lines for the 'varying:' block, one per parameter, 1-based ``[n]``.

    Given the ``script``, a THICKNESS bound is labelled with its layer's unit.
    """
    if not varying:
        return ["  varying:    (nothing selected)"]
    lines = ["  varying:"]
    for i, v in enumerate(varying, start=1):
        bound = ""
        if v.bounds:
            unit = _bound_unit(v, script)
            bound = f"  bounds {v.bounds[0]:g}-{v.bounds[1]:g}{f' {unit}' if unit else ''}"
        lines.append(f"    [{i}] {v.name}{bound}")
    return lines


def _vary_command(entry: Vary) -> str:
    """The PERT command line that would recreate this one selection."""
    if entry.kind == "thickness":
        base = f"thickness {entry.layer + 1}"
    elif entry.kind == "composition":
        base = f"composition {entry.layer + 1} {entry.symbol}"
    elif entry.kind == "atoms":
        base = f"atoms {entry.layer + 1} {entry.symbol}"
    elif entry.kind == "species":
        base = f"species {entry.layer + 1} {entry.symbol}"
    elif entry.kind == "equation":
        base = f"equation {entry.layer + 1} {entry.index + 1}"
    elif entry.kind == "pixe_h":
        base = f"pixe h {entry.symbol}"
    else:
        base = entry.name
    if entry.bounds is not None:
        base += f" {entry.bounds[0]:g} {entry.bounds[1]:g}"
    return base


def _to_lines(state: PertState) -> list[str]:
    """Command lines that recreate ``state``, for PERT SAVE/GET."""
    lines = [f"window {w.low} {w.high}" for w in state.windows.error]
    lines += [f"pixe {w.low} {w.high}" for w in state.pixe_windows]
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
    every parameter's built-in default range. A THICKNESS bound is the one
    exception that is *not* copied onto ``parameter`` at selection time: it
    is in the layer's own unit, and GO converts it (:func:`_fit_parameter`).
    """

    parameter: object
    kind: str                 # thickness | composition | atoms | equation | simple | sample | pixe_h
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
    #: PIXE channels (numbered as the .PIX file numbers them) fitted
    #: together with the RBS spectrum -- ``PIXE <lo> <hi>``.
    pixe_windows: list[Window] = field(default_factory=list)
    multi: bool = True
    verbose: bool = False
    autocmp: bool = False
    report: bool = False

    def describe(self, script=None) -> str:
        lines = [f"  mode        {'multiple' if self.multi else 'single'} variable"]
        lines.append(f"  autocmp     {'on' if self.autocmp else 'off'}")
        lines.append(f"  report      {'on' if self.report else 'off'}")
        lines.append(f"  error win   {_format_windows(self.windows.error)}")
        norm = self.windows.normalisation
        lines.append(
            f"  norm win    {f'{norm.low}-{norm.high}' if norm else '(none)'}"
        )
        lines.append(f"  PIXE win    {_format_pixe_windows(self.pixe_windows)}")
        lines.extend(_format_varying(self.varying, script))
        return "\n".join(lines)


def state_for(session) -> PertState:
    if session.pert is None:
        session.pert = PertState()
    return session.pert


#: The MODE each mode-specific selection kind needs (see ``MODE``); every
#: other kind works in either.
_KIND_MODE = {"thickness": "comp", "composition": "comp", "atoms": "atoms"}


def _wrong_mode(entry: Vary, mode: str) -> bool:
    return _KIND_MODE.get(entry.kind, mode) != mode


def drop_other_mode(session, mode: str) -> str | None:
    """Drop every PERT selection that needs a MODE other than ``mode``.

    MODE converts each layer to the other convention, so a selection made
    for the old one no longer means what it did: an ATOMS variable would
    read a COMP fraction as an areal density, and vice versa. Called by
    MODE; returns the one line to tell the user, or None if nothing was
    selected that way.
    """
    state = session.pert
    if state is None:
        return None
    dropped = [v for v in state.varying if _wrong_mode(v, mode)]
    if not dropped:
        return None
    state.varying = [v for v in state.varying if not _wrong_mode(v, mode)]

    by_layer: dict[int, list[str]] = {}
    for v in dropped:
        by_layer.setdefault(v.layer, []).append(f"{v.kind} {v.symbol}".rstrip())
    groups = "; ".join(
        f"{', '.join(names)} (layer {layer + 1})" for layer, names in by_layer.items()
    )
    needs = "ATOMS needs MODE ATOMS" if mode == "comp" else (
        "THICKNESS/COMPOSITION need MODE COMP"
    )
    return f"  PERT: dropped {groups} -- {needs}"


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
    """Add ``entry`` to what PERT varies, or replace it in place if the same
    parameter is already selected -- RUMP's own semantics (pert.c's
    ``PertAddVar``, "Gets name and adds (or modifies) a variable"): reselecting
    a parameter resets its search bound to whatever this call gave (or the
    default range, if none), rather than being rejected as a duplicate.
    """
    state = state_for(session)
    state.varying = [v for v in state.varying if v.name != entry.name]
    state.varying.append(entry)
    print(f"  varying {entry.name}")


# ---------------------------------------------------------------------------
# Selecting parameters
# ---------------------------------------------------------------------------


@needs_mode("comp")
def cmd_thickness(session, args: ArgReader) -> None:
    """``THICKNESS <layer> [<min> <max>]`` -- vary a layer's thickness.

    The bound is in the layer's own unit (whatever SIM SHOW lists it in, e.g.
    Angstroms), not the 1e15 at/cm^2 the fit parameter itself works in, so it
    is left off the parameter here and converted at GO (:func:`_fit_parameter`).
    """
    layer = _layer_argument(session, args)
    bound = _optional_bounds(args)
    args.done()
    _add(
        session,
        Vary(
            parameter=thickness(layer),
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


@needs_mode("comp")
def cmd_composition(session, args: ArgReader) -> None:
    """``COMPOSITION <layer> <element> [<min> <max>]`` -- vary one element's
    stoichiometry in a layer.
    """
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


@needs_mode("atoms")
def cmd_atoms(session, args: ArgReader) -> None:
    """``ATOMS <layer> <element> [<min> <max>]`` -- vary one element's own
    areal density.

    Unlike ``COMPOSITION``, the layer's total thickness is re-derived as the
    sum of its composition row on every trial, so the other elements' own
    areal density stays fixed instead of being reallocated (see
    :func:`pyrump.fit.parameters.atoms`). Only meaningful for a layer whose
    composition values are already each element's 1e15 at/cm^2 -- set up
    with ``SIM ATOMS``.
    """
    layer = _layer_argument(session, args)
    symbol = args.token("an element symbol")
    bound = _optional_bounds(args)
    args.done()
    index = _element_index(session, symbol)
    canonical = _layer_element(session, layer, symbol, "composition")
    param = atoms(layer, index)
    if bound is not None:
        param = replace(param, lower=bound[0], upper=bound[1])
    _add(
        session,
        Vary(
            parameter=param,
            kind="atoms",
            layer=layer,
            index=index,
            symbol=canonical,
            name=f"layer {layer + 1} atoms {canonical}",
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
        print(state.describe(session.script))
        return
    token = args.peek()
    if token is not None and token.lower() in ("clear", "none", "reset"):
        args.token()
        if not args:
            args.done()
            state.windows.error = []
            print("  error windows cleared")
            return
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


def cmd_pixe(session, args: ArgReader) -> None:
    """``PIXE <lo> <hi>`` -- fit the PIXE spectrum over these channels too,
    together with the RBS one: the same sample drives both, and each
    spectrum weighs in by its own counting statistics. Channels are
    numbered as the ``.PIX`` file numbers them, as PIXE ``REGION`` takes
    them. Put the windows on clear peaks: the PIXE continuum is not
    simulated. Meant for elements whose RBS signals overlap but whose X-ray
    lines don't (Ta-W, Fe-Ni, Ni-Co) -- best with lines of the same shell,
    so H cancels in their ratio.

    ``PIXE CLEAR [<n>]`` removes window *n*, or all of them.
    ``PIXE H K|L|M [<min> <max>]`` varies that shell's instrumental
    constant H, so the PIXE spectrum decides the ratio and RBS the amount.
    ``PIXE`` alone lists the windows. A pyRUMP addition. (In PERT, ``PIXE``
    means this: RETURN first for the PIXE prompt.)
    """
    state = state_for(session)
    if not args:
        print(f"  PIXE windows {_format_pixe_windows(state.pixe_windows)}")
        return
    token = args.peek().lower()
    if token in ("clear", "none", "reset"):
        args.token()
        if not args:
            state.pixe_windows = []
            print("  PIXE windows cleared")
            return
        n = args.integer("a PIXE window number")
        args.done()
        windows = state.pixe_windows
        if not windows:
            raise CommandError("no PIXE windows are set")
        if not 1 <= n <= len(windows):
            raise CommandError(f"PIXE window {n} is outside 1-{len(windows)}")
        del windows[n - 1]
        print(f"  PIXE windows {_format_pixe_windows(windows)}")
        return
    if token == "h":
        args.token()
        family = args.token("K, L or M").upper()
        if family not in ("K", "L", "M"):
            raise CommandError(f"PIXE H: expected K, L or M, not {family!r}")
        bound = _optional_bounds(args)
        args.done()
        param = pixe_h(family)
        if bound is not None:
            param = replace(param, lower=bound[0], upper=bound[1])
        _add(session, Vary(
            parameter=param, kind="pixe_h", symbol=family, name=f"PIXE H {family}", bounds=bound,
        ))
        return
    low = args.integer("the first PIXE channel")
    high = args.integer("the last PIXE channel")
    args.done()
    if high <= low:
        raise CommandError(f"empty window: {low} to {high}")
    if len(state.pixe_windows) >= MAX_ERROR_WINDOWS:
        raise CommandError(f"at most {MAX_ERROR_WINDOWS} PIXE windows")
    state.pixe_windows.append(Window(low, high))
    print(f"  PIXE windows {_format_pixe_windows(state.pixe_windows)}")


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
    print(state_for(session).describe(session.script))


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

    typed = Path(args.token("a .pert file"))
    run_go = False
    token = args.peek()
    if token is not None and token.lower() == "go":
        args.token()
        run_go = True
    args.done()
    if not typed.suffix:
        typed = typed.with_suffix(".pert")
    path = session.locate(typed)
    if path is None:
        raise CommandError(f"no such file: {typed}")
    # autocmp/report are standing preferences (typically set once from
    # .pyrumprc), not part of the file-specific selection GET replaces --
    # carry them over so a fresh GET doesn't silently turn them back off.
    state = state_for(session)
    session.pert = PertState(autocmp=state.autocmp, report=state.report)
    execute_file(session, path, stack=["rump", "pert"])
    print(f"read {path}")
    print(state_for(session).describe(session.script))
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
        # Same standing-preference carve-out as GET (see cmd_get).
        state = state_for(session)
        session.pert = PertState(autocmp=state.autocmp, report=state.report)
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
    print("\n".join(_format_varying(varying, session.script)))


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


def cmd_autocmp(session, args: ArgReader) -> None:
    """``AUTOCMP [off]`` -- run COMPARE automatically at the end of GO
    (default off).

    Saves the separate COMPARE call after every fit, for anyone who always
    wants to eyeball the residuals right away.
    """
    token = args.optional()
    args.done()
    state = state_for(session)
    state.autocmp = token is None or token.lower() not in ("off", "no", "0")
    print(f"  autocmp {'on' if state.autocmp else 'off'}")


def cmd_report(session, args: ArgReader) -> None:
    """``REPORT [off]`` -- after every ``GO``, take a ``SNAPSHOT`` of the fit
    under the sample's own name (default off).

    A pyRUMP-only addition, not part of legacy RUMP. Once on, no further
    action is needed per sample: each ``GO`` writes the files named after
    ``<sample>``, the active data buffer's own file stem (e.g. data loaded
    from ``MA8410.RBS`` writes files starting ``MA8410``) -- so switching
    samples with ``XEQ``/``GET`` naturally routes later fits to different
    files with no filename to remember:

    - ``<sample>.report`` -- fitted values, uncertainties and chi-square,
      appended (one block per ``GO``, so refitting the same sample keeps a
      history)
    - ``<sample>.pert`` -- the PERT selection, as ``PERT SAVE`` would write
    - ``<sample>.lcm`` -- the sample description, as ``SIM SAVE`` would write
    - ``<sample>.png`` -- the ``COMPARE`` plot, as ``FIGSAVE``/``HCOPY``
      would write (drawn fresh for this, even if ``AUTOCMP`` is off)
    - ``<sample>_fit.xeq`` -- the restore macro: ``pyrump <sample>_fit.xeq``
      is this session again (see :mod:`pyrump.shell.snapshot`)

    ``GO`` echoes each one back ("updated .../wrote ...") at the end of its
    own output, so none of the writes are silent. Tuning by hand after the
    fit is saved by a ``SNAPSHOT`` of its own. A standing preference like
    ``AUTOCMP``: survives ``GET``/``CLEAR``, not saved by ``SAVE``.
    """
    token = args.optional()
    args.done()
    state = state_for(session)
    state.report = token is None or token.lower() not in ("off", "no", "0")
    print(f"  report {'on' if state.report else 'off'}")


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

    print(describe_topic(session, topic, (TABLE, RUMP_TABLE, SYSTEM_TABLE)))


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
    elif entry.kind in ("composition", "species", "atoms"):
        target = (
            layers[entry.layer].species
            if entry.kind == "species"
            else layers[entry.layer].composition
        )
        symbol = entry.symbol or session.script.elements[entry.index]
        target[symbol] = value
        if entry.kind == "atoms":
            # composition values are each element's own 1e15 at/cm^2 here
            # (SIM ATOMS's convention), so the layer's total is their sum.
            layer = layers[entry.layer]
            layer.thickness = sum(layer.composition.values())
            layer.unit = "/CM2"
    elif entry.kind == "equation":
        profile = layers[entry.layer].profile
        params = list(profile.parameters)
        params[entry.index] = value
        layers[entry.layer].profile = replace(profile, parameters=tuple(params))
    elif entry.kind == "sample":
        setattr(session.script, entry.name, value)
    elif entry.kind == "pixe_h":
        session.pixe.h = tuple(inputs.pixe_h)
    else:
        # A buffer parameter. fit() mutates ``inputs`` in place and leaves it
        # holding the best-fit objects, so the buffer just adopts them.
        buffer = session.buffers.require_active()
        buffer.beam = inputs.beam
        buffer.geometry = inputs.geometry
        buffer.measurement = inputs.measurement
        buffer.spectrum.calibration = inputs.calibration


def _report_stem(buffer) -> str:
    """A short, filesystem-safe name for this buffer's spectrum, shared by
    every file ``REPORT`` writes and by the plot legend -- see
    :func:`~pyrump.shell.snapshot.stem`."""
    return snapshot.stem(buffer)


def _report_path(buffer) -> Path:
    """``<sample>.report``, where ``<sample>`` is :func:`_report_stem`."""
    return Path(f"{_report_stem(buffer)}.report")


def _units_per_areal(session, inputs: FitInputs, layer: int) -> float | None:
    """How many of a layer's own thickness unit make 1e15 at/cm^2.

    The ratio of the layer's script thickness to the areal density
    ``to_sample`` made of it: SimThickConvert is linear in the magnitude (see
    :func:`_write_back`), so this one factor converts a fitted value, its
    uncertainty and a bound alike, whatever the unit. None for a layer with
    no thickness to take the ratio from.
    """
    areal = float(inputs.sample.thicknesses[layer])
    if areal <= 0:
        return None
    return session.script.layers[layer].thickness / areal


def _fit_parameter(entry: Vary, per_areal: float | None):
    """The parameter the solver gets for ``entry``: a THICKNESS bound, typed
    in the layer's own unit, converted to the 1e15 at/cm^2 the parameter
    works in. Every other kind already carries its bound (see :class:`Vary`).
    """
    if entry.kind != "thickness" or entry.bounds is None:
        return entry.parameter
    if per_areal is None:
        raise CommandError(
            f"go: {entry.name} is zero, so its bound has no unit to convert from"
        )
    low, high = entry.bounds
    return replace(entry.parameter, lower=low / per_areal, upper=high / per_areal)


def _g6(value: float) -> str:
    return f"{value:.6g}"


def _display(session, entry: Vary, per_areal: float | None):
    """How GO shows ``entry``'s values, as ``(factor, format, unit)``.

    A THICKNESS is shown in its layer's own unit, the way SIM SHOW lists it
    (``factor`` converts from the fit's 1e15 at/cm^2), and an ATOMS amount in
    /CM2; anything else as the fit holds it, with no unit.
    """
    if entry.kind == "thickness":
        if per_areal is None:
            return 1.0, thickness_label, _AREAL
        return per_areal, thickness_label, session.script.layers[entry.layer].unit
    if entry.kind == "atoms":
        return 1.0, _g6, _AREAL
    return 1.0, _g6, ""


def _pixe_extra(session, state: PertState, data_buffer) -> list:
    """The PIXE spectrum as GO's second spectrum, over the PIXE windows --
    or nothing, without windows. Checked here, before the fit starts."""
    from ...fit.lm import ExtraSpectrum
    from ..pixe_sim import simulate_with

    if not state.pixe_windows:
        if any(v.kind == "pixe_h" for v in state.varying):
            raise CommandError(
                "go: PIXE H sets the PIXE spectrum's scale, but no PIXE windows "
                "are set -- PIXE <lo> <hi>"
            )
        return []
    pixe = data_buffer.pixe
    if pixe is None:
        raise CommandError(
            "go: PIXE windows are set, but the active buffer has no PIXE spectrum "
            "-- PIXE GET <file>, or PIXE PAIR ON and read the .RBS again"
        )
    counts = np.asarray(pixe.spectrum.counts, dtype=float)
    first = round(pixe.calibration.first)
    mask = np.zeros(counts.size, dtype=bool)
    for window in state.pixe_windows:
        low, high = window.low - first, window.high - first
        if high < 0 or low >= counts.size:
            raise CommandError(
                f"go: PIXE window {window.low}-{window.high} is outside the PIXE "
                f"spectrum's channels {first}-{first + counts.size - 1}"
            )
        mask[max(low, 0):min(high, counts.size - 1) + 1] = True

    def run_pixe(current: FitInputs) -> np.ndarray:
        return simulate_with(
            session, current.sample, current.beam, current.geometry,
            current.measurement, current.pixe_h, pixe,
        ).counts

    return [ExtraSpectrum(simulate=run_pixe, data=counts, mask=mask)]


def cmd_go(session, args: ArgReader) -> None:
    args.done()
    from ...fit.lm import fit
    from ...script.lcm import structure_label, to_sample
    from ...sim.engine import simulate

    state = state_for(session)
    if not state.varying:
        raise CommandError("nothing selected to vary")
    stale = [v.name for v in state.varying if _wrong_mode(v, session.thickness_mode)]
    if stale:
        raise CommandError(
            f"go: {', '.join(stale)} do not match MODE "
            f"{session.thickness_mode.upper()} -- CLEAR them or switch MODE"
        )
    if not session.script.layers:
        raise CommandError("no sample described: use SIM to build one")

    data_buffer = session.buffers.require_active()
    observed = np.asarray(data_buffer.spectrum.counts, dtype=float)
    extra = _pixe_extra(session, state, data_buffer)

    sample_id = _report_stem(data_buffer)
    initial_structure = structure_label(
        session.script, normalize=session.plot.composition_fraction
    )
    header_lines = [f"  Fitting {sample_id}: {initial_structure}"]
    if extra:
        header_lines.append(
            f"  with the PIXE spectrum over channels "
            f"{', '.join(f'{w.low}-{w.high}' for w in state.pixe_windows)}"
        )
    for line in header_lines:
        print(line)
    report_lines = list(header_lines)

    sample = to_sample(session.script, session.table, session.densities)
    inputs = FitInputs(
        sample=sample,
        beam=data_buffer.beam,
        geometry=data_buffer.geometry,
        calibration=data_buffer.calibration,
        measurement=data_buffer.measurement,
        pixe_h=tuple(session.pixe.h),
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
            screening=session.settings.screening,
            faithful=session.settings.faithful,
        ).counts

    starting = {v.name: v.parameter.get(inputs) for v in state.varying}
    per_areal = {
        v.layer: _units_per_areal(session, inputs, v.layer)
        for v in state.varying
        if v.kind == "thickness"
    }
    fit_parameters = {v.name: _fit_parameter(v, per_areal.get(v.layer)) for v in state.varying}
    groups = (
        [state.varying] if state.multi else [[v] for v in state.varying]
    )

    initial_reduced = None

    def _progress(evaluation: int, reduced: float) -> None:
        nonlocal initial_reduced
        if initial_reduced is None:
            initial_reduced = reduced
        if state.verbose:
            print(f"    eval {evaluation:3d}   chi2/dof {reduced:.4f}")
        # Each evaluation is a full simulation, so a fit blocks the prompt --
        # and the plot window's event loop -- for seconds to minutes. This is
        # the only per-iteration hook the fit has; pumping here keeps the
        # figure repainting instead of being declared hung.
        plotting.pump(session)

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
                    [fit_parameters[v.name] for v in group],
                    windows=state.windows,
                    progress=_progress,
                    extra=extra,
                )
                if not state.multi:
                    entry = group[0]
                    factor, fmt, unit = _display(session, entry, per_areal.get(entry.layer))
                    value = result.parameters[entry.parameter.name] * factor
                    print(
                        f"  {entry.name}: {fmt(value)}{f' {unit}' if unit else ''}"
                        f"   chi2/dof {result.reduced_chi_square:.4f}"
                    )
    except ValueError as error:
        raise CommandError(f"go: {error}") from None
    elapsed = time.perf_counter() - started

    for entry in state.varying:
        _write_back(session, entry, inputs, starting[entry.name])
    session.editor = None
    session.touch()

    report_lines.append(f"\n  fit took {elapsed:.2f} s")
    chi_line = f"\n  reduced chi-square {result.reduced_chi_square:.4f} on {result.dof} dof"
    if initial_reduced is not None:
        chi_line += f"   (was {initial_reduced:.4f})"
    report_lines.append(chi_line)
    if result.parts:
        report_lines.append("  " + ",   ".join(
            f"{label} chi-square {part.total:.1f} over {part.n_used} channels"
            for label, part in result.parts.items()
        ))
    status = "converged" if result.success else "did not converge"
    report_lines.append(f"  {result.n_evaluations} evaluations, {status}")
    if result.normalisation != 1.0:
        # RUMP writes the fitted scale back into the buffer's CORR factor
        # (pert.c:1402-1403, "Estimated correction factor set for buffer"),
        # rather than only reporting it -- CORRECTION and a normalisation
        # window can't be varied together (WindowSet.validate_against), so
        # this never fights with a PERT CORRECTION selection.
        new_correction = data_buffer.measurement.correction * result.normalisation
        data_buffer.measurement = replace(data_buffer.measurement, correction=new_correction)
        report_lines.append(
            f"  data scaled by {result.normalisation:.5f} over the norm window"
            f"   (correction factor set to {new_correction:.5g})"
        )
    if result.n_invalid:
        report_lines.append(
            f"  warning: {result.n_invalid} windowed channels had "
            "no predicted counts"
        )
    # An ATOMS line brackets its layer's resulting thickness in Angstroms --
    # the same view SIM SHOW brackets in MODE ATOMS.
    angstroms = [None] * len(session.script.layers)
    if any(v.kind == "atoms" for v in state.varying):
        try:
            angstroms = thickness_mode_views(
                session.script, session.table, session.densities, to_atoms=False
            )
        except KeyError:
            pass  # an element the table doesn't know: no brackets, as in SIM SHOW
    for entry in state.varying:
        name = entry.parameter.name
        value = entry.parameter.get(inputs)
        sigma = result.uncertainties.get(name)
        before = starting[entry.name]
        factor, fmt, unit = _display(session, entry, per_areal.get(entry.layer))
        suffix = f" {unit}" if unit else ""
        line = f"  {entry.name:26s} {fmt(value * factor) + suffix:>14s}"
        if sigma:
            line += f"  +/- {sigma * factor:.4g}{suffix}"
        line += f"   (was {fmt(before * factor)}{suffix})"
        if entry.kind == "thickness" and unit.upper() != _AREAL:
            line += f"   [{thickness_label(value)} {_AREAL}]"
        elif entry.kind == "atoms" and angstroms[entry.layer] is not None:
            line += f"   [{thickness_label(angstroms[entry.layer][0])} A]"
        report_lines.append(line)
    new_structure = structure_label(session.script, normalize=session.plot.composition_fraction)
    report_lines.append(f"\n  {sample_id}: {new_structure}")

    for line in report_lines[len(header_lines):]:
        print(line)

    # What this fit left behind, so a SNAPSHOT of exactly this state records
    # it as fitted, uncertainties and all.
    session.last_fit = (snapshot.fingerprint(session, labels=False), report_lines)

    # REPORT's snapshot draws its own COMPARE, for its .png.
    if state.report:
        snapshot.take(session, after_go=True)
    elif state.autocmp:
        cmd_compare(session, ArgReader([], command="compare"))


TABLE = CommandTable("PERT Commands")

_ENTRIES: list[tuple[str, int, object, str]] = [
    ("?", -1, cmd_help, "synonym for HELP"),
    ("HELP", 1, cmd_help, "list the PERT commands"),
    ("RETURN", 1, cmd_return, "return to the RUMP level"),
    ("QUIT", -1, cmd_return, "synonym for RETURN (not exit pyRUMP)"),
    ("Q", -1, cmd_return, "synonym for RETURN"),
    ("GO", 2, cmd_go, "run the search"),
    ("PARMS", 2, cmd_parms, "display the current settings"),
    ("SHOW", 2, cmd_show, "display the sample description (same as SIM SHOW)"),
    ("GET", 2, cmd_get, "replay a saved PERT selection from a .pert file, or GET <file> GO"),
    ("SAVE", 2, cmd_save, "save the current PERT selection to a .pert file"),
    ("CLEAR", 2, cmd_clear, "forget every selected parameter and window, or CLEAR <n> one parameter"),
    # Windows and mode
    ("WINDOW", 2, cmd_window, "set an error window in channels, or WINDOW CLEAR [<n>]"),
    ("NORMALIZE", 2, cmd_normalize, "set the normalisation window, or NORMALIZE CLEAR"),
    ("PIXE", 2, cmd_pixe,
     "fit PIXE too: PIXE <lo> <hi> adds a window, PIXE CLEAR [<n>], PIXE H K|L|M [<min> <max>]"),
    ("SINGLE", 2, cmd_single, "vary one parameter at a time"),
    ("MULTI", 3, cmd_multi, "vary all parameters together (default)"),
    ("VOLUME", 3, cmd_volume, "verbose progress messages"),
    ("AUTOCMP", 4, cmd_autocmp, "run COMPARE automatically at the end of GO"),
    ("REPORT", 3, cmd_report,
     "SNAPSHOT after every GO: <sample>.report, _fit.xeq ... (REPORT OFF to stop)"),
    # Registered here too, not only at the RUMP level, so switching doesn't
    # fall through and leave PERT (repl.py's execute_line).
    ("MODE", 4, cmd_mode,
     "switch between COMP and ATOMS (converts every layer, drops mismatched selections)"),
    # Parameters -- all take an optional trailing "<min> <max>" search bound
    ("THICKNESS", 2, cmd_thickness,
     "vary a layer thickness [<min> <max> in the layer's unit, e.g. A] (needs MODE COMP)"),
    ("COMPOSITION", 3, cmd_composition,
     "vary an element in a layer [<min> <max>] (needs MODE COMP)"),
    ("ATOMS", 3, cmd_atoms,
     "vary one element's areal density, others held fixed [<min> <max> in /CM2] "
     "(needs MODE ATOMS)"),
    ("SPECIES", 2, cmd_species, "vary the species composition [<min> <max>]"),
    ("EQUATION", 2, cmd_equation, "vary an equation parameter [<min> <max>]"),
    ("MEV", 3, _simple("mev"), "vary the beam energy [<min> <max>]"),
    ("FWHM", 2, _simple("fwhm"), "vary the detector resolution [<min> <max>]"),
    ("STRAGGLE", 5, _simple("straggle", "sample"), "vary the straggling constant [<min> <max>]"),
    ("FUZZ", 2, cmd_fuzz, "vary the fuzz parameter (not implemented)"),
    ("CORRECTION", 3, _simple("correction"), "vary the normalization correction [<min> <max>]"),
    ("THETA", 4, _simple("theta"), "vary the sample tilt [<min> <max>]"),
    ("SLOPE", 5, _simple("kev/ch"),
     "vary the calibration slope, keV per channel [<min> <max>]"),
    ("KEV/CH", -4, _simple("kev/ch"), "synonym for SLOPE"),
    ("OFFSET", 3, _simple("kev(0)"),
     "vary the calibration energy offset (e.g. a sample-charging shift) [<min> <max>]"),
    ("KEV(0)", -6, _simple("kev(0)"), "synonym for OFFSET"),
    ("COMPARE", 0, cmd_compare, "plot the active buffer against the simulation"),
    ("CMP", -3, cmd_compare, "synonym for COMPARE"),
    ("EXPORTCMP", 7, cmd_exportcmp,
     "write the active buffer, simulation, difference and GOF as columns"),
    ("EC", -2, cmd_exportcmp, "synonym for EXPORTCMP"),
    ("EXPORT", 4, cmd_export, "write the active buffer as columns: channel, energy, counts, error"),
]

for _name, _minlen, _handler, _help in _ENTRIES:
    TABLE.add(_name, _minlen, _handler, _help)
TABLE.note_synonym("HELP", "?")
TABLE.note_synonym("RETURN", "QUIT", "Q")
TABLE.note_synonym("SLOPE", "KEV/CH")
TABLE.note_synonym("OFFSET", "KEV(0)")
TABLE.note_synonym("COMPARE", "CMP")
TABLE.note_synonym("EXPORTCMP", "EC")
