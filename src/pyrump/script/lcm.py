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


def parse_lcm(text: str) -> Script:
    """Parse a ``.lcm`` sample description."""
    script = Script()
    current: LcmLayer | None = None

    for raw in text.splitlines():
        line = raw.split("/*", 1)[0].strip()
        if not line or line.startswith("#") or line.startswith("!"):
            continue
        tokens = line.split()
        command = tokens[0].lower()
        rest = tokens[1:]

        if command == "sim" and rest and rest[0].lower().startswith("res"):
            script = Script()
            current = None
        elif command == "layer":
            current = LcmLayer()
            script.layers.append(current)
        elif command == "next":
            current = LcmLayer()
            script.layers.append(current)
        elif command.startswith("desc"):
            script.description = line.split(None, 1)[1].strip().strip("'\"")
        elif command.startswith("thick") and current is not None:
            current.thickness = float(rest[0])
            current.unit = rest[1] if len(rest) > 1 else "A"
        elif command.startswith("sublay") and current is not None:
            current.sublayers = int(float(rest[0]))
        elif command.startswith("sthick") and current is not None:
            current.sub_thickness = float(rest[0])
            current.sub_unit = rest[1] if len(rest) > 1 else "A"
        elif command.startswith("comp") and current is not None:
            current.composition = _element_pairs(rest)
        elif command.startswith("species") and current is not None:
            current.species = _element_pairs(rest)
        elif command.startswith("eq") and current is not None:
            name = rest[0].lower()
            if name not in EQUATION_NAMES:
                raise ValueError(f"unknown equation {rest[0]!r}")
            # Parameters may be followed by a unit token, which we drop --
            # doses are carried in the values themselves.
            values = [float(t) for t in rest[1:] if _is_number(t)]
            script.layers[-1].profile = Profile(
                EQUATION_NAMES[name], tuple(values)
            )
        elif command.startswith("fuzz") and current is not None:
            current.fuzz_amount = float(rest[0])
            current.fuzz_steps = int(float(rest[1])) if len(rest) > 1 else 0
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

    # A trailing "Next" with nothing after it leaves an empty layer.
    script.layers = [
        layer for layer in script.layers if layer.thickness > 0 or layer.composition
    ]
    return script


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
