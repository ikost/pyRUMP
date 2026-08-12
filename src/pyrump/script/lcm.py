"""RUMP ``.lcm`` / ``.sim`` sample-description files.

These are plain-text SIM command scripts. RUMP writes them with
``SimWriteSample`` (sim2.c:2070-2160) and reads them back by feeding each line
to the SIM command interpreter, so the format is simply a subset of the
interactive command language::

    Sim Reset
    Layer 1
     Thick 151 ITO
     Composition In 2 O 3 Sn 0.1 /
    Next
     Thick 10 um
     Composition O 4 C 14 H 10 /
    Maxpth 200
    Foil disable

Scope
-----

The **sample-definition** commands only — enough to round-trip the files RUMP
itself writes. Deliberately excluded: the ``G_*`` global-profile subsystem, the
``gvparse`` expression language, and everything to do with plotting, buffers or
fitting. Those are thousands of lines of UI, and the Python API is pyRUMP's
primary interface.

Unrecognised commands are collected in :attr:`Script.ignored` rather than
raising, so a real-world file with plotting commands in it still yields its
sample.

Thickness units
---------------

``Thick`` takes a magnitude and a unit, where the "unit" may be a length
(``A``, ``nm``, ``um``), an explicit areal density (``/CM2``, ``M/CM2``), or a
**compound name** from ``density.tab``. Conversion to areal density follows
``SimThickConvert`` (sim2.c:2349):

* length units multiply by the layer's own atomic density
* ``/CM2`` is already 1e15 at/cm^2
* ``M/CM2`` multiplies by the composition sum (molecules to atoms)
* a compound name uses that compound's tabulated density
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from ..atomic.density import DensityTable, ThicknessKind
from ..atomic.tables import PeriodicTable
from ..profiles.equations import EquationType, Profile

#: Equation names RUMP accepts, mapped to pyRUMP's types. Several are aliases
#: (sim2.c:101-131): "Error" for ERFC, "Thicfilm" for ThickFilm, "Implant" for
#: Gaussian.
EQUATION_NAMES: dict[str, EquationType] = {
    "none": EquationType.NONE,
    "constant": EquationType.CONSTANT,
    "linear": EquationType.LINEAR,
    "erfc": EquationType.ERFC,
    "error": EquationType.ERFC,
    "exponential": EquationType.EXPONENTIAL,
    "semi-infinite": EquationType.SEMI_INFINITE,
    "thinfilm": EquationType.THINFILM,
    "buriedthinfilm": EquationType.BURIEDTHINFILM,
    "thickfilm": EquationType.THICKFILM,
    "thicfilm": EquationType.THICKFILM,
    "timedependent": EquationType.TIMEDEPENDENT,
    "gaussian": EquationType.GAUSSIAN,
    "implant": EquationType.GAUSSIAN,
    "edgeworth": EquationType.EDGEWORTH,
}


@dataclass(slots=True)
class LcmLayer:
    """One layer as written in the file, before unit conversion."""

    thickness: float = 0.0
    unit: str = "A"
    composition: dict[str, float] = field(default_factory=dict)
    sublayers: int = 0
    sub_thickness: float = 0.0
    sub_unit: str = "A"
    profile: Profile | None = None
    species: dict[str, float] = field(default_factory=dict)
    fuzz_amount: float = 0.0
    fuzz_steps: int = 0


@dataclass(slots=True)
class Script:
    """A parsed sample description."""

    layers: list[LcmLayer] = field(default_factory=list)
    description: str = ""
    maxpth: float = 200.0
    absorber_layers: int = 0
    straggle: float = 0.0
    multiple: float = 0.0
    foil: bool = False
    ignored: list[str] = field(default_factory=list)

    @property
    def elements(self) -> list[str]:
        """Element symbols in first-seen order, matrix before species."""
        seen: list[str] = []
        for layer in self.layers:
            for source in (layer.composition, layer.species):
                for name in source:
                    if name not in seen:
                        seen.append(name)
        return seen


def _element_pairs(tokens: list[str]) -> dict[str, float]:
    """Parse ``El value El value ... /`` into a mapping.

    The trailing ``/`` terminates the list; RUMP's command reader uses it to
    mean "no more arguments".
    """
    out: dict[str, float] = {}
    index = 0
    while index < len(tokens):
        token = tokens[index]
        if token == "/":
            break
        if index + 1 >= len(tokens) or tokens[index + 1] == "/":
            raise ValueError(f"element {token!r} has no amount")
        out[token] = float(tokens[index + 1])
        index += 2
    return out


class SampleEditor:
    """Applies SIM sample-definition commands one line at a time.

    RUMP's ``.lcm`` files *are* SIM command scripts, replayed through the same
    interpreter the user types at (see the module docstring), so file parsing
    and interactive editing share this one implementation. :func:`parse_lcm` is
    a loop over :meth:`execute`; :mod:`pyrump.shell.commands.sim` drives the
    same object from the prompt.

    Layer bookkeeping follows sim.htm: there is always one blank layer at the
    bottom of the structure, and a layer left at zero thickness "disappears
    automatically if the active layer changes". :meth:`_prune` is that rule, and
    it is also why a file ending in a bare ``Next`` does not gain an empty layer.
    """

    __slots__ = ("script", "current")

    def __init__(self, script: Script | None = None):
        self.script = script if script is not None else Script()
        #: Index into ``script.layers``; may equal ``len(layers)``, meaning the
        #: implicit blank layer at the bottom.
        self.current = 0

    # -- layer pointer ----------------------------------------------------

    @property
    def layer(self) -> LcmLayer | None:
        """The layer under the pointer, or None when it is on the blank one."""
        if 0 <= self.current < len(self.script.layers):
            return self.script.layers[self.current]
        return None

    def _writable(self) -> LcmLayer:
        """The layer under the pointer, materialising the blank one on write."""
        if self.layer is None:
            self.script.layers.append(LcmLayer())
            self.current = len(self.script.layers) - 1
        return self.script.layers[self.current]

    def _prune(self) -> None:
        """Drop layers that were never given substance (sim.htm, OPEN)."""
        kept = [
            (index, layer)
            for index, layer in enumerate(self.script.layers)
            if layer.thickness > 0 or layer.composition
        ]
        if len(kept) == len(self.script.layers):
            return
        surviving = {index for index, _ in kept}
        before = sum(1 for index in range(self.current) if index in surviving)
        self.script.layers = [layer for _, layer in kept]
        self.current = min(before, len(self.script.layers))

    def select(self, index: int) -> None:
        """Move the pointer, pruning what the old layer left behind."""
        self._prune()
        self.current = max(0, min(index, len(self.script.layers)))

    def next_layer(self) -> None:
        self.select(self.current + 1)

    def open_layer(self) -> None:
        """Insert a blank layer above the pointer and select it."""
        self._prune()
        index = min(self.current, len(self.script.layers))
        self.script.layers.insert(index, LcmLayer())
        self.current = index

    def reset(self) -> None:
        self.script = Script()
        self.current = 0

    # -- command application ----------------------------------------------

    def execute(self, line: str) -> None:
        """Apply one command line. Unknown commands land in ``script.ignored``."""
        line = line.split("/*", 1)[0].strip()
        if not line or line.startswith("#") or line.startswith("!"):
            return
        tokens = line.split()
        command = tokens[0].lower()
        rest = tokens[1:]
        script = self.script

        if command == "sim" and rest and rest[0].lower().startswith("res"):
            self.reset()
        elif command == "layer":
            # "Layer n" navigates (sim.htm); with no argument, or past the end,
            # it lands on the blank layer at the bottom, which writing fills in.
            self.select(int(float(rest[0])) - 1 if rest else len(script.layers))
            self._writable()
        elif command == "next":
            self.next_layer()
            self._writable()
        elif command == "open":
            self.open_layer()
        elif command.startswith("desc"):
            script.description = line.split(None, 1)[1].strip().strip("'\"")
        elif command.startswith("thick"):
            layer = self._writable()
            layer.thickness = float(rest[0])
            layer.unit = rest[1] if len(rest) > 1 else "A"
        elif command.startswith("sublay"):
            self._writable().sublayers = int(float(rest[0]))
        elif command.startswith("sthick"):
            layer = self._writable()
            layer.sub_thickness = float(rest[0])
            layer.sub_unit = rest[1] if len(rest) > 1 else "A"
        elif command.startswith("comp"):
            self._writable().composition = _element_pairs(rest)
        elif command.startswith("species"):
            self._writable().species = _element_pairs(rest)
        elif command.startswith("eq"):
            name = rest[0].lower()
            if name not in EQUATION_NAMES:
                raise ValueError(f"unknown equation {rest[0]!r}")
            # Parameters may be followed by a unit token, which we drop --
            # doses are carried in the values themselves.
            values = [float(t) for t in rest[1:] if _is_number(t)]
            self._writable().profile = Profile(EQUATION_NAMES[name], tuple(values))
        elif command.startswith("fuzz"):
            layer = self._writable()
            layer.fuzz_amount = float(rest[0])
            layer.fuzz_steps = int(float(rest[1])) if len(rest) > 1 else 0
        elif command.startswith("maxp"):
            script.maxpth = float(rest[0])
        elif command.startswith("absorb"):
            script.absorber_layers = int(float(rest[0]))
        elif command.startswith("strag"):
            script.straggle = float(rest[0])
        elif command.startswith("multi"):
            script.multiple = float(rest[0])
        elif command == "foil":
            script.foil = bool(rest) and rest[0].lower() not in ("disable", "off", "no")
        else:
            script.ignored.append(line)

    def finish(self) -> Script:
        """Prune the trailing blank layer and hand back the script."""
        self._prune()
        return self.script


def parse_lcm(text: str) -> Script:
    """Parse a ``.lcm`` sample description."""
    editor = SampleEditor()
    for raw in text.splitlines():
        editor.execute(raw)
    return editor.finish()


def _is_number(token: str) -> bool:
    try:
        float(token)
    except ValueError:
        return False
    return True


def read_lcm(path: str | Path) -> Script:
    return parse_lcm(Path(path).read_text())


def to_sample(
    script: Script,
    periodic_table: PeriodicTable,
    densities: DensityTable | None = None,
):
    """Convert a parsed script into a :class:`UniformSample`.

    Thickness units are resolved here, which needs the periodic table (for
    atomic densities) and optionally the compound table.
    """
    from ..atomic.density import layer_atomic_density
    from ..sim.engine import UniformSample

    import numpy as np

    symbols = script.elements
    if not symbols:
        raise ValueError("script defines no elements")
    element_z = [periodic_table.by_symbol(s).z for s in symbols]
    atomic_densities = np.array(
        [periodic_table.by_z(z).atomic_density for z in element_z], dtype=np.float64
    )

    thicknesses: list[float] = []
    compositions: list[list[float]] = []
    species: list[list[float]] = []
    profiles: list[Profile | None] = []
    sublayers: list[int] = []
    fuzz_amounts: list[float] = []
    fuzz_steps: list[int] = []

    for layer in script.layers:
        row = [layer.composition.get(s, 0.0) for s in symbols]
        density = layer_atomic_density(row, atomic_densities)
        areal = _to_areal(layer.thickness, layer.unit, density, sum(row), densities)

        thicknesses.append(areal)
        compositions.append(row)
        species.append([layer.species.get(s, 0.0) for s in symbols])
        profiles.append(layer.profile)
        sublayers.append(layer.sublayers)
        fuzz_amounts.append(
            _to_areal(layer.fuzz_amount, layer.unit, density, sum(row), densities)
            if layer.fuzz_amount
            else 0.0
        )
        fuzz_steps.append(layer.fuzz_steps)

    return UniformSample(
        thicknesses=thicknesses,
        element_z=element_z,
        compositions=compositions,
        sublayers=sublayers if any(sublayers) else None,
        maxpth=script.maxpth,
        straggle=script.straggle,
        multiple=script.multiple,
        absorber_layers=script.absorber_layers,
        profiles=profiles if any(p is not None for p in profiles) else None,
        species=species if any(any(r) for r in species) else None,
        fuzz_amounts=fuzz_amounts if any(fuzz_steps) else None,
        fuzz_steps=fuzz_steps if any(fuzz_steps) else None,
    )


def _to_areal(
    magnitude: float,
    unit: str,
    density: float,
    composition_sum: float,
    densities: DensityTable | None,
) -> float:
    """Convert a thickness to areal density in 1e15 at/cm^2 (sim2.c:2349)."""
    densities = densities or DensityTable(compounds={})

    known = densities.unit(unit)
    if known is not None:
        if known.kind is ThicknessKind.ANGSTROMS:
            return magnitude * known.scale * density
        if known.kind is ThicknessKind.ATOMIC:
            return magnitude * known.scale
        return magnitude * known.scale * composition_sum

    # Otherwise it names a compound, whose tabulated density replaces the
    # layer's own (SimThickConvert's ABSOLUTE branch writes back *pdensity).
    return magnitude * densities.density(unit)


def write_lcm(script: Script) -> str:
    """Render a script in RUMP's own format (``SimWriteSample``)."""
    lines = ["Sim Reset"]
    if script.description:
        lines.append(f"Description '{script.description}'")
    lines.append("Layer 1")

    for index, layer in enumerate(script.layers):
        lines.append(f" Thick {_g(layer.thickness)} {layer.unit}")
        if layer.sublayers:
            lines.append(f"  Sublayer {layer.sublayers}")
        elif layer.sub_thickness:
            lines.append(f"  Sthickness {_g(layer.sub_thickness)} {layer.sub_unit}")

        pairs = " ".join(f"{k} {_g(v)}" for k, v in layer.composition.items())
        lines.append(f" Composition {pairs} /")

        if layer.profile is not None:
            name = next(
                key for key, value in EQUATION_NAMES.items()
                if value is layer.profile.type and key not in ("error", "thicfilm", "implant")
            )
            values = " ".join(_g(v) for v in layer.profile.parameters)
            lines.append(f" Equation {name.capitalize()} {values}")
            if layer.species:
                pairs = " ".join(f"{k} {_g(v)}" for k, v in layer.species.items())
                lines.append(f" Species {pairs} /")

        if layer.fuzz_steps:
            lines.append(f" Fuzzy {layer.fuzz_amount:f} {layer.fuzz_steps}")

        if index != len(script.layers) - 1:
            lines.append("Next")

    lines.append(f"Maxpth {_g(script.maxpth)}")
    if script.absorber_layers:
        lines.append(f"Absorber {script.absorber_layers}")
    lines.append("Foil disable")
    return "\n".join(lines) + "\n"


def _g(value: float) -> str:
    """Format like C's ``%g``, which is what RUMP writes."""
    return f"{value:g}"
