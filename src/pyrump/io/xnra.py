"""SIMNRA's ``.xnra`` file format.

An ``.xnra`` file is an IDF (IBA Data Format) XML document -- default namespace
``http://idf.schemas.itn.pt`` -- extended with SIMNRA's own elements in the
``http://www.simnra.com/simnra`` namespace (SIMNRA User's Guide, section 3.22).
It carries one measured spectrum, the instrument parameters, the sample's
layer structure and SIMNRA's own simulated spectra.

Real files differ from the published description in several places, and the
code follows the real files (SIMNRA 7.03 and 7.04):

* ``simnra:xnraversionnr`` is ``1.1``, not ``1.0``;
* simulated spectra are several ``process/simulations/simulation`` elements,
  one ``total`` plus ``pileup`` and ``partialelement`` curves, with
  ``simpledata`` directly under each one (no ``calculateddatas`` wrapper);
* calibration parameters carry no count, only units -- ``keV``,
  ``keV/channel``, ``keV/channel^2`` -- in that order;
* numbers are Delphi-formatted, ``5.10000000000000E+0001``.

Two quantities need converting rather than copying:

* **Dose.** SIMNRA has no separate solid angle: ``beamfluence`` is
  particles x solid angle (msr), and ``solidangle`` is written as 1 msr. RUMP's
  yield scales as omega * charge / (charge state * CORR)
  (:func:`~pyrump.model.detector.yield_normalisation`), so reading keeps that
  product: OMEGA = ``solidangle``, CHARGE = N * e, charge state 1, CORR 1.
* **Geometry.** SIMNRA stores alpha (incidence), the scattering angle and beta
  (exit); RUMP stores THETA, PHI = 180 - scattering angle, and PSI. The file is
  taken as it stands. The one wrinkle is IBM, where RUMP's exit angle is
  theta + phi and so depends on theta's sign, which SIMNRA does not store: the
  sign is taken from the file's own beta.

Writing patches a tree rather than generating one: the buffer's original file
when it came from ``.xnra`` (so SIMNRA-only settings pyRUMP does not model
survive), else a skeleton modelled on the files SIMNRA itself writes.
"""

from __future__ import annotations

import copy
import math
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import numpy as np

from ..model.detector import Measurement
from ..model.geometry import Geometry, GeometryKind
from ..model.spectrum import Calibration
from ..physics.xsec.rutherford import ScreeningModel
from ..script.lcm import LcmLayer, Script
from .rbs import RbsSpectrum

IDF_NS = "http://idf.schemas.itn.pt"
SIMNRA_NS = "http://www.simnra.com/simnra"
_NS = {"i": IDF_NS, "s": SIMNRA_NS}

ET.register_namespace("", IDF_NS)
ET.register_namespace("simnra", SIMNRA_NS)

#: Written into ``simnra:xnraversionnr`` -- what SIMNRA 7.03/7.04 write.
XNRA_VERSION = "1.1"

#: Exact, SI 2019.
ELEMENTARY_CHARGE = 1.602176634e-19

#: Unit conversions to the canonical unit each quantity is read in, keyed by
#: (unit as written, lower-cased; canonical unit).
_CONVERSIONS: dict[tuple[str, str], float] = {
    ("kev", "keV"): 1.0,
    ("mev", "keV"): 1e3,
    ("ev", "keV"): 1e-3,
    ("degree", "degree"): 1.0,
    ("deg", "degree"): 1.0,
    ("rad", "degree"): 180.0 / math.pi,
    ("msr", "msr"): 1.0,
    ("sr", "msr"): 1e3,
    ("s", "s"): 1.0,
    ("ms", "s"): 1e-3,
    ("us", "s"): 1e-6,
    ("us", "us"): 1.0,
    ("s", "us"): 1e6,
    ("ns", "us"): 1e-3,
    ("amu", "amu"): 1.0,
    ("1e15at/cm2", "1e15at/cm2"): 1.0,
    ("at/cm2", "1e15at/cm2"): 1e-15,
    ("kev/channel", "keV/channel"): 1.0,
    ("kev/channel^2", "keV/channel^2"): 1.0,
    ("#particles", "#particles"): 1.0,
    ("fraction", "fraction"): 1.0,
}

#: SIMNRA's ``geometrytype`` strings.
_GEOMETRY_NAMES = {
    "cornell": GeometryKind.CORNELL,
    "ibm": GeometryKind.IBM,
    "general": GeometryKind.GENERAL,
}
_GEOMETRY_WRITTEN = {
    GeometryKind.CORNELL: "Cornell",
    GeometryKind.IBM: "IBM",
    GeometryKind.GENERAL: "general",
}

#: ``28Si``, ``Si``, ``2D`` -- an optional mass number, then a symbol.
_NUCLIDE = re.compile(r"^\s*(\d*)\s*([A-Z][a-z]?)\s*$")

#: Hydrogen's isotopes under their own symbols, which RUMP's table lacks.
_HYDROGEN_ALIASES = {"D": "H", "T": "H"}


class XnraFormatError(ValueError):
    """The file is not a SIMNRA ``.xnra`` file pyRUMP can read."""


@dataclass(slots=True)
class XnraFile:
    """What :func:`read_xnra` found."""

    spectrum: RbsSpectrum | None
    """The measured spectrum with its metadata, or ``None`` if the file holds
    only a simulation."""

    simulation: RbsSpectrum | None
    """SIMNRA's ``total`` simulated spectrum, with the same metadata, or
    ``None`` if there is none (SIMNRA writes a single zero point when nothing
    has been calculated)."""

    script: Script | None
    """The sample as SIM layers, or ``None`` if the file has no sample."""

    physics: dict[str, str] = field(default_factory=dict)
    """SIMNRA's calculation settings, for telling the user what differs."""

    notices: list[str] = field(default_factory=list)
    """Everything that was ignored or converted, one line each, for the user."""

    xml: bytes = b""
    """The original document, kept so a later write can patch it."""

    @property
    def geometry(self) -> Geometry:
        source = self.spectrum or self.simulation
        return source.geometry


@dataclass(slots=True)
class XnraLayer:
    """One layer as written: areal density and atomic fractions."""

    thickness: float
    """1e15 atoms/cm^2."""

    composition: dict[str, float]
    """Element symbol -> atomic fraction, summing to 1."""


# ------------------------------------------------------------------ helpers


def _find(parent: ET.Element, path: str) -> ET.Element | None:
    return parent.find(path, _NS)


def _text(parent: ET.Element, path: str, default: str = "") -> str:
    element = _find(parent, path)
    if element is None or element.text is None:
        return default
    return element.text.strip()


def _number(parent: ET.Element, path: str, unit: str, default: float | None = None) -> float | None:
    """A quantity at ``path`` in ``unit``, converting from its ``@units``."""
    element = _find(parent, path)
    if element is None or not (element.text or "").strip():
        return default
    return _value(element, unit)


def _value(element: ET.Element, unit: str) -> float:
    tag = element.tag.split("}")[-1]
    try:
        value = float(element.text)
    except (TypeError, ValueError):
        raise XnraFormatError(f"<{tag}>: {element.text!r} is not a number") from None
    written = element.get("units")
    if written is None or written == unit:
        return value
    try:
        return value * _CONVERSIONS[(written.lower(), unit)]
    except KeyError:
        raise XnraFormatError(
            f"<{tag}>: unit {written!r} not understood (expected {unit})"
        ) from None


def _boolean(text: str) -> bool:
    return text.strip().lower() in ("true", "1", "yes", "on")


def _numbers(text: str) -> np.ndarray:
    return np.array(text.split(), dtype=np.float64) if text.strip() else np.zeros(0)


def format_number(value: float) -> str:
    """Delphi's ``%.14E`` with a four-digit exponent, as SIMNRA writes it."""
    mantissa, exponent = f"{float(value):.14E}".split("E")
    return f"{mantissa}E{int(exponent):+05d}"


# ------------------------------------------------------------------ reading


def read_xnra(path: str | Path) -> XnraFile:
    """Read a SIMNRA ``.xnra`` file."""
    path = Path(path)
    blob = path.read_bytes()
    try:
        root = ET.fromstring(blob)
    except ET.ParseError as error:
        raise XnraFormatError(f"{path.name} is not valid XML: {error}") from None
    result = _parse(root, path)
    result.xml = blob
    return result


def _parse(root: ET.Element, path: Path) -> XnraFile:
    # The validity rules of the SIMNRA User's Guide, section 3.22.
    if root.tag != f"{{{IDF_NS}}}idf":
        raise XnraFormatError(f"{path.name} is not an IDF file (root is not <idf>)")
    if _text(root, "i:attributes/s:filetype").lower() != "xnra":
        raise XnraFormatError(
            f"{path.name} is a generic IDF file, not one written by SIMNRA "
            "(no simnra:filetype = xnra)"
        )
    version = _text(root, "i:attributes/s:xnraversionnr")
    notices: list[str] = []
    if version and not version.startswith("1."):
        notices.append(f"xnra version {version}: only 1.x is known, reading anyway")

    samples = root.findall("i:sample", _NS)
    if len(samples) != 1:
        raise XnraFormatError(f"{path.name} has {len(samples)} samples, expected exactly one")
    sample = samples[0]
    spectra = sample.findall("i:spectra/i:spectrum", _NS)
    if not spectra:
        raise XnraFormatError(f"{path.name} contains no spectrum")
    if len(spectra) > 1:
        notices.append(f"the file holds {len(spectra)} spectra: only the first is loaded")
    spectrum_node = spectra[0]

    totals = [
        s for s in spectrum_node.findall("i:process/i:simulations/i:simulation", _NS)
        if _text(s, "i:simulationtype", "total").lower() == "total"
    ]
    if len(totals) > 1:
        raise XnraFormatError(f"{path.name} has {len(totals)} total simulations, expected one")

    metadata = _metadata(root, spectrum_node, notices)
    measured = _channel_data(_find(spectrum_node, "i:data"), "measured spectrum", notices)
    simulated = _channel_data(totals[0], "SIMNRA simulation", notices) if totals else None
    if measured is None and simulated is None:
        raise XnraFormatError(f"{path.name} contains neither measured nor simulated data")

    identifier = (
        _text(spectrum_node, "i:data/s:graphics/s:legend")
        or _text(sample, "i:description")
        or path.stem
    )
    measured_spectrum = (
        _spectrum(metadata, *measured, identifier) if measured is not None else None
    )
    simulated_spectrum = (
        _spectrum(metadata, *simulated, f"SIMNRA simulation: {path.name}")
        if simulated is not None
        else None
    )

    return XnraFile(
        spectrum=measured_spectrum,
        simulation=simulated_spectrum,
        script=_sample(sample, notices),
        physics=_physics(spectrum_node),
        notices=notices,
    )


@dataclass(slots=True)
class _Metadata:
    geometry: Geometry
    measurement: Measurement
    kevch: float
    kev0: float
    e0_MeV: float
    zbeam: int
    mbeam: float
    date: str
    livetime: str
    comments: list[str]


def _metadata(root: ET.Element, spectrum: ET.Element, notices: list[str]) -> _Metadata:
    technique = _text(spectrum, "i:reactions/i:technique", "RBS")
    if technique.upper() != "RBS":
        notices.append(
            f"technique {technique}: pyRUMP simulates RBS only -- the counts are "
            "loaded, but buffer 0 will not describe them"
        )

    # Beam.
    zbeam = int(_number(spectrum, "i:beam/i:beamZ", "", 2.0) or 2)
    mbeam = _number(spectrum, "i:beam/i:beammass", "amu", 4.0026) or 4.0026
    energy_keV = _number(spectrum, "i:beam/i:beamenergy", "keV", 0.0) or 0.0
    spread = _number(spectrum, "i:beam/i:beamenergyspread", "keV", 0.0) or 0.0
    if spread:
        notices.append(f"beam energy spread {spread:g} keV FWHM ignored")
    if _layers(spectrum, "i:beam/i:beamfoil/i:foillayers"):
        notices.append("beam foil ignored")

    # Geometry: taken as it stands, see the module docstring.
    kind = _GEOMETRY_NAMES.get(
        _text(spectrum, "i:geometry/i:geometrytype", "general").lower(), GeometryKind.GENERAL
    )
    alpha = _number(spectrum, "i:geometry/i:incidenceangle", "degree", 0.0) or 0.0
    scattering = _number(spectrum, "i:geometry/i:scatteringangle", "degree", 170.0) or 170.0
    beta = _number(spectrum, "i:geometry/i:exitangle", "degree", 0.0) or 0.0
    phi = 180.0 - scattering
    theta = alpha
    if kind is GeometryKind.IBM:
        # RUMP's IBM exit angle is |theta + phi|: pick the sign of theta that
        # gives the file's own beta.
        plus = abs(abs(abs(alpha) + phi) - beta)
        minus = abs(abs(-abs(alpha) + phi) - beta)
        theta = abs(alpha) if plus <= minus else -abs(alpha)
    geometry = Geometry(theta=theta, phi=phi, psi=beta, kind=kind)

    # Detector and dose.
    omega = _number(spectrum, "i:detection/i:detector/i:solidangle", "msr", 1.0) or 1.0
    fluence = _number(spectrum, "i:beam/i:beamfluence", "#particles", 0.0) or 0.0
    charge_uC = fluence * ELEMENTARY_CHARGE * 1e6
    if fluence:
        notices.append(
            "SIMNRA stores charge x solid angle together: loaded as "
            f"OMEGA {omega:g} msr, CHARGE {charge_uC:.4g} uC. If you know them "
            "separately, set them with OMEGA and CHARGE"
        )
    if _layers(spectrum, "i:detection/i:stoppingfoil/i:foillayers"):
        notices.append("stopping foil in front of the detector ignored")
    shaping_s = _number(spectrum, "i:detection/i:electronics/i:amplifier/i:shapingtime", "us", 0.0)

    resolutions = spectrum.findall(
        "i:calibrations/i:detectorresolutions/i:detectorresolution"
        "/i:resolutionparameters/i:resolutionparameter", _NS,
    )
    fwhm = _value(resolutions[0], "keV") if resolutions else 15.0
    if len(resolutions) > 1:
        notices.append(
            f"detector resolution has {len(resolutions)} parameters: only the "
            f"first (FWHM {fwhm:g} keV) is used"
        )

    measurement = Measurement(
        omega_msr=omega,
        charge_uC=charge_uC or Measurement().charge_uC,
        correction=1.0,
        charge_state=1,
        fwhm_keV=fwhm,
        tau_us=shaping_s or Measurement().tau_us,
    )

    # Energy calibration: E(ch) = offset + gain*ch + quad*ch^2.
    calibration = _find(spectrum, "i:calibrations/i:energycalibrations/i:energycalibration")
    offset, gain, quad = 0.0, 1.0, 0.0
    if calibration is not None:
        mode = _text(calibration, "i:calibrationmode", "energy")
        if mode.lower() != "energy":
            notices.append(f"calibration mode {mode!r} read as an energy calibration")
        parameters = calibration.findall("i:calibrationparameters/i:calibrationparameter", _NS)
        by_unit = {p.get("units", "").lower(): p for p in parameters}
        if {"kev", "kev/channel"} <= by_unit.keys():
            offset = _value(by_unit["kev"], "keV")
            gain = _value(by_unit["kev/channel"], "keV/channel")
            if "kev/channel^2" in by_unit:
                quad = _value(by_unit["kev/channel^2"], "keV/channel^2")
        else:  # no units: positional, as IDF orders them
            values = [float(p.text) for p in parameters]
            offset, gain, quad = (values + [0.0, 1.0, 0.0][len(values):])[:3]
    if quad:
        notices.append(
            f"quadratic calibration term {quad:g} keV/channel^2 ignored -- RUMP's "
            "calibration is linear"
        )

    livetime = _number(spectrum, "i:log/i:livetime", "s")
    comments = [
        note.text.strip()
        for note in root.findall("i:notes/i:note", _NS) + spectrum.findall("i:notes/i:note", _NS)
        if note.text and note.text.strip()
    ]
    return _Metadata(
        geometry=geometry,
        measurement=measurement,
        kevch=gain,
        kev0=offset,
        e0_MeV=energy_keV / 1000.0,
        zbeam=zbeam,
        mbeam=mbeam,
        date=_text(spectrum, "i:log/i:starttime"),
        livetime=f"{livetime:g}" if livetime else "",
        comments=comments,
    )


def _channel_data(
    node: ET.Element | None, what: str, notices: list[str]
) -> tuple[np.ndarray, float] | None:
    """Counts and the first channel number from ``node/simpledata``.

    A lone zero point -- what SIMNRA writes for a simulation never
    calculated -- counts as no data.
    """
    if node is None:
        return None
    x = _numbers(_text(node, "i:simpledata/i:x"))
    y = _numbers(_text(node, "i:simpledata/i:y"))
    if y.size == 0 or (y.size == 1 and y[0] == 0.0):
        return None
    if x.size != y.size:
        raise XnraFormatError(f"{what}: {x.size} channel numbers but {y.size} counts")
    if np.any(np.diff(x) != 1.0):
        raise XnraFormatError(f"{what}: channel numbers are not consecutive")
    mode = _text(node, "i:channelmode", "left")
    if mode.lower() != "left":
        notices.append(
            f"{what}: channel mode {mode!r} read as 'left' (energy at the channel's lower edge)"
        )
    return y, float(x[0])


def _spectrum(meta: _Metadata, counts: np.ndarray, first: float, identifier: str) -> RbsSpectrum:
    return RbsSpectrum(
        counts=counts,
        calibration=Calibration(kevch=meta.kevch, kev0=meta.kev0, first=first, npt=counts.size),
        geometry=meta.geometry,
        measurement=meta.measurement,
        e0_MeV=meta.e0_MeV,
        zbeam=meta.zbeam,
        mbeam=meta.mbeam,
        identifier=identifier,
        date=meta.date,
        livetime=meta.livetime,
        comments=list(meta.comments),
    )


def _layers(parent: ET.Element, path: str) -> list[ET.Element]:
    """Layers under ``path`` that actually hold material."""
    container = _find(parent, path)
    if container is None:
        return []
    return [
        layer for layer in container.findall("i:layer", _NS)
        if (_number(layer, "i:layerthickness", "1e15at/cm2", 0.0) or 0.0) > 0
    ]


def _sample(sample: ET.Element, notices: list[str]) -> Script | None:
    structure = _find(sample, "i:structure/i:layeredstructure")
    if structure is None:
        return None
    layers: list[LcmLayer] = []
    for number, node in enumerate(structure.findall("i:layers/i:layer", _NS), start=1):
        thickness = _number(node, "i:layerthickness", "1e15at/cm2", 0.0) or 0.0
        composition: dict[str, float] = {}
        for element in node.findall("i:layerelements/i:layerelement", _NS):
            fraction = _number(element, "i:concentration", "fraction", 0.0) or 0.0
            if fraction <= 0:
                continue
            name = _text(element, "i:name")
            symbol = _symbol(name, number, notices)
            composition[symbol] = composition.get(symbol, 0.0) + fraction
        if thickness <= 0 or not composition:
            if composition:
                notices.append(f"layer {number}: zero thickness, skipped")
            continue
        _layer_extras(node, number, notices)
        layers.append(LcmLayer(thickness=thickness, unit="/CM2", composition=composition))

    distribution = _text(structure, "s:substraterougness/s:distribution", "none")
    if distribution.lower() != "none":
        notices.append(f"substrate roughness ({distribution}) ignored")
    if not layers:
        return None
    return Script(layers=layers, description=_text(sample, "i:description"))


def _symbol(name: str, layer: int, notices: list[str]) -> str:
    """The element symbol of ``Si``, ``28Si`` or ``2D``.

    pyRUMP's samples are natural elements, so an isotope folds into its
    element, with a notice.
    """
    match = _NUCLIDE.match(name)
    if match is None:
        raise XnraFormatError(f"layer {layer}: element name {name!r} not understood")
    mass_number, symbol = match.groups()
    element = _HYDROGEN_ALIASES.get(symbol, symbol)
    if mass_number or element != symbol:
        notices.append(
            f"layer {layer}: isotope {name} read as natural {element} "
            "(pyRUMP samples use natural isotopic abundance)"
        )
    return element


def _layer_extras(node: ET.Element, number: int, notices: list[str]) -> None:
    """Notices for SIMNRA layer features pyRUMP has no equivalent of."""
    if _boolean(_text(node, "s:hasroughness", "false")):
        notices.append(f"layer {number}: roughness ignored")
    if _boolean(_text(node, "s:porosity/s:hasporosity", "false")):
        notices.append(f"layer {number}: porosity ignored")
    for child in node.iter():
        tag = child.tag.split("}")[-1].lower()
        if "correctionfactor" in tag and (child.text or "").strip():
            try:
                factor = float(child.text)
            except ValueError:
                continue
            if factor != 1.0:
                notices.append(f"layer {number}: correction factor {tag} = {factor:g} ignored")


def _physics(spectrum: ET.Element) -> dict[str, str]:
    """SIMNRA's calculation settings that decide what its simulation shows."""
    defaults = _find(spectrum, "i:process/i:physicsdefaults")
    if defaults is None:
        return {}
    physics = {
        "stopping": _text(defaults, "i:stoppingpowerdefault/i:computercode/i:name"),
        "straggling": _text(defaults, "i:energyspreaddefault/i:energylossstraggling"),
        "screening": _text(defaults, "i:crosssectiondefault/i:screening"),
        "rutherford": _text(defaults, "i:crosssectiondefault/i:Rutherford"),
        "multiple scattering": _text(defaults, "i:energyspreaddefault/i:multiplescattering"),
        "dual scattering": _text(defaults, "s:dualscatteringdefault/s:dualscattering"),
        "pileup": _text(defaults, "s:pileupcalculationdefault/s:pileup"),
    }
    return {key: value for key, value in physics.items() if value}


def describe_physics(physics: dict[str, str], screening: ScreeningModel) -> list[str]:
    """Lines telling the user which SIMNRA settings pyRUMP does not share."""
    lines = []
    if physics.get("stopping"):
        lines.append(f"stopping {physics['stopping']!r} (pyRUMP: its own tables)")
    if physics.get("straggling"):
        lines.append(f"straggling {physics['straggling']!r} (pyRUMP: RUMP's own model)")
    if _boolean(physics.get("rutherford", "false")):
        simnra_screening, same = "none (pure Rutherford)", screening is ScreeningModel.NONE
    else:
        simnra_screening = physics.get("screening", "")
        same = simnra_screening.lower() == screening.name.lower()
    if simnra_screening and not same:
        lines.append(f"screening {simnra_screening!r} (pyRUMP: {screening.name})")
    for key in ("multiple scattering", "dual scattering", "pileup"):
        if _boolean(physics.get(key, "false")):
            lines.append(f"{key} on in SIMNRA")
    return lines


# ------------------------------------------------------------------ writing


def layers_from_script(script: Script, periodic_table, densities, geometry: Geometry):
    """The SIM sample as :class:`XnraLayer` s, plus notices.

    Uniform layers are written as they are. A graded layer (profile or
    species) has no SIMNRA equivalent, so it is written as the uniform
    sublayers pyRUMP's own simulation slices it into.
    """
    from ..script.lcm import to_sample
    from ..sim.engine import build_sample_grid

    notices: list[str] = []
    sample = to_sample(script, periodic_table, densities)
    symbols = [periodic_table.by_z(z).symbol for z in sample.element_z]
    grid = None
    if sample.profiles is not None or sample.species is not None:
        grid = build_sample_grid(sample, geometry, periodic_table)

    layers: list[XnraLayer] = []
    for index, thickness in enumerate(sample.thicknesses):
        number = index + 1
        if index < sample.absorber_layers:
            notices.append(f"layer {number} is an absorber layer: not written")
            continue
        graded = grid is not None and (
            (sample.profiles and sample.profiles[index] is not None)
            or (sample.species and any(sample.species[index]))
        )
        if graded:
            rows = np.flatnonzero(grid.layer_index == index)
            notices.append(f"layer {number} is graded: written as {rows.size} uniform sublayers")
            for row in rows:
                layers.append(
                    XnraLayer(float(grid.areal_density[row]), _fractions(symbols, grid.composition[row]))
                )
            continue
        if sample.fuzz_steps and sample.fuzz_steps[index]:
            notices.append(f"layer {number}: interface fuzz not written")
        layers.append(XnraLayer(float(thickness), _fractions(symbols, sample.compositions[index])))
    return layers, notices


def _fractions(symbols: list[str], row) -> dict[str, float]:
    total = float(sum(row))
    if total <= 0:
        return {}
    return {s: float(v) / total for s, v in zip(symbols, row) if v}


def write_xnra(
    path: str | Path,
    buffer,
    *,
    layers: list[XnraLayer] | None = None,
    simulation: np.ndarray | None = None,
    screening: ScreeningModel | None = None,
    periodic_table=None,
    base: bytes | None = None,
    description: str = "",
    version: str = "",
) -> list[str]:
    """Write ``buffer`` (a :class:`~pyrump.shell.session.Buffer` or anything
    with its fields) as an ``.xnra`` file. Returns notices for the user.

    ``layers`` of ``None`` leaves the sample as the base file has it (an empty
    one for a fresh file). ``simulation`` is pyRUMP's own simulated spectrum,
    on the buffer's channels.
    """
    path = Path(path)
    notices: list[str] = []
    root = ET.fromstring(base) if base else ET.fromstring(_SKELETON)
    spectrum = _find(root, "i:sample/i:spectra/i:spectrum")
    if spectrum is None:
        raise XnraFormatError("base document has no sample/spectra/spectrum")

    # File attributes.
    _set(root, "i:attributes/i:idfversion", "1.01")
    _set(root, "i:attributes/i:filename", path.name)
    _set(root, "i:attributes/i:createtime", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    _set(root, "i:attributes/s:filetype", "xnra")
    _set(root, "i:attributes/s:xnraversionnr", XNRA_VERSION)
    if not base:
        note = "file created by pyRUMP" + (f" {version}" if version else "")
        _set(root, "i:notes/i:note", note)

    # Beam.
    beam, g, m, c = buffer.beam, buffer.geometry, buffer.measurement, buffer.calibration
    if periodic_table is not None:
        symbol = periodic_table.by_z(beam.z).symbol
        _set(spectrum, "i:beam/i:beamparticle", f"{round(beam.mass)}{symbol}")
    _set(spectrum, "i:beam/i:beamZ", str(beam.z))
    _set(spectrum, "i:beam/i:beammass", format_number(beam.mass), units="amu")
    _set(spectrum, "i:beam/i:beamenergy", format_number(beam.e0_MeV * 1000.0), units="keV")
    _set(spectrum, "i:beam/i:beamfluence", format_number(fluence(m)), units="#particles")

    # Geometry: beta is the exit angle pyRUMP itself simulates with.
    _set(spectrum, "i:geometry/i:geometrytype", _GEOMETRY_WRITTEN[g.kind])
    _set(spectrum, "i:geometry/i:incidenceangle", format_number(abs(g.theta)), units="degree")
    _set(spectrum, "i:geometry/i:scatteringangle", format_number(g.scattering_angle), units="degree")
    _set(spectrum, "i:geometry/i:exitangle", format_number(exit_angle(g)), units="degree")

    # Detector: the solid angle already lives in the fluence.
    _set(spectrum, "i:detection/i:detector/i:solidangle", format_number(1.0), units="msr")
    resolution = _set(
        spectrum,
        "i:calibrations/i:detectorresolutions/i:detectorresolution"
        "/i:resolutionparameters/i:resolutionparameter",
        format_number(m.fwhm_keV), units="keV",
    )
    resolution.set("mode", "FWHM")
    _calibration(spectrum, c)

    livetime = _float_or_none(buffer.livetime)
    if livetime:
        _set(spectrum, "i:log/i:livetime", format_number(livetime), units="s")

    # Measured data.
    first = int(round(c.first))
    counts = np.asarray(buffer.spectrum.counts, dtype=np.float64)
    data = _ensure(spectrum, "i:data")
    _simpledata(data, first, counts)
    legend = _ensure(data, "s:graphics/s:legend")
    legend.text = buffer.identifier or path.stem
    processed = _find(spectrum, "s:processeddata")
    if processed is not None:  # smoothed copies of the old counts are stale
        for child in list(processed):
            processed.remove(child)

    # Simulations: pyRUMP's own replaces the total's data. The total element
    # itself stays, since it also holds SIMNRA's cross-section choices; the
    # pileup and per-element curves no longer match and are dropped.
    simulations = _ensure(spectrum, "i:process/i:simulations")
    total = None
    for child in simulations.findall("i:simulation", _NS):
        if total is None and _text(child, "i:simulationtype", "total").lower() == "total":
            total = child
        else:
            simulations.remove(child)
    if total is None:
        total = ET.SubElement(simulations, f"{{{IDF_NS}}}simulation")
        _set(total, "i:simulationtype", "total")
    if simulation is not None:
        _simpledata(total, first, np.asarray(simulation, dtype=np.float64))
        _ensure(total, "s:graphics/s:legend").text = "pyRUMP simulation"
    else:
        _simpledata(total, 0, np.zeros(1))

    if screening is not None:
        notices.extend(_screening(spectrum, screening))

    if layers is not None:
        _sample_layers(root, layers, description)

    tree = ET.ElementTree(root)
    ET.indent(tree, space="\t")
    tree.write(path, encoding="utf-8", xml_declaration=True)
    return notices


def fluence(measurement: Measurement) -> float:
    """SIMNRA's particles x solid angle (msr) for a RUMP dose."""
    denominator = measurement.charge_state * measurement.correction * ELEMENTARY_CHARGE
    if denominator == 0:
        return 0.0
    return measurement.omega_msr * measurement.charge_uC * 1e-6 / denominator


def exit_angle(geometry: Geometry) -> float:
    """The exit angle, in degrees from the normal, that pyRUMP simulates with."""
    if geometry.kind is GeometryKind.GENERAL:
        return abs(geometry.psi)
    sec_out = geometry.sec_out
    if sec_out == 0.0:
        return 90.0
    return math.degrees(math.acos(max(-1.0, min(1.0, 1.0 / sec_out))))


def _float_or_none(text: str) -> float | None:
    try:
        return float(str(text).split()[0])
    except (ValueError, IndexError):
        return None


def _ensure(parent: ET.Element, path: str) -> ET.Element:
    """The element at ``path``, creating any missing steps."""
    node = parent
    for step in path.split("/"):
        found = node.find(step, _NS)
        if found is None:
            prefix, tag = step.split(":")
            found = ET.SubElement(node, f"{{{_NS[prefix]}}}{tag}")
        node = found
    return node


def _set(parent: ET.Element, path: str, text: str, *, units: str | None = None) -> ET.Element:
    element = _ensure(parent, path)
    element.text = text
    if units is not None:
        element.set("units", units)
    return element


def _calibration(spectrum: ET.Element, calibration: Calibration) -> None:
    node = _ensure(spectrum, "i:calibrations/i:energycalibrations/i:energycalibration")
    _set(node, "i:calibrationmode", "energy")
    parameters = _ensure(node, "i:calibrationparameters")
    for child in list(parameters):
        parameters.remove(child)
    for value, unit in (
        (calibration.kev0, "keV"),
        (calibration.kevch, "keV/channel"),
        (0.0, "keV/channel^2"),
    ):
        element = ET.SubElement(parameters, f"{{{IDF_NS}}}calibrationparameter", units=unit)
        element.text = format_number(value)


def _simpledata(node: ET.Element, first: int, counts: np.ndarray) -> None:
    _set(node, "i:datamode", "simple")
    _set(node, "i:channelmode", "left")
    data = _ensure(node, "i:simpledata")
    _set(data, "i:xaxis/i:axisname", "channel")
    _set(data, "i:xaxis/i:axisunit", "#")
    _set(data, "i:yaxis/i:axisname", "yield")
    _set(data, "i:yaxis/i:axisunit", "counts")
    _set(data, "i:x", " ".join(str(first + i) for i in range(counts.size)))
    _set(data, "i:y", " ".join(format_number(v) for v in counts))


def _screening(spectrum: ET.Element, screening: ScreeningModel) -> list[str]:
    defaults = _ensure(spectrum, "i:process/i:physicsdefaults/i:crosssectiondefault")
    if screening is ScreeningModel.NONE:
        _set(defaults, "i:Rutherford", "true")
        return []
    if screening is ScreeningModel.ANDERSEN:
        _set(defaults, "i:Rutherford", "false")
        _set(defaults, "i:screening", "Andersen")
        return []
    current = _text(defaults, "i:screening") or "(unset)"
    return [
        f"screening: SIMNRA's name for L'Ecuyer is not known, so the file keeps "
        f"{current!r}"
    ]


def _sample_layers(root: ET.Element, layers: list[XnraLayer], description: str) -> None:
    sample = _find(root, "i:sample")
    if description:
        _set(sample, "i:description", description)

    symbols: list[str] = []
    for layer in layers:
        for symbol in layer.composition:
            if symbol not in symbols:
                symbols.append(symbol)
    elements = _ensure(sample, "i:elementsandmolecules/i:elements")
    for child in list(elements):
        elements.remove(child)
    _set(elements, "i:nelements", str(len(symbols)))
    for symbol in symbols:
        _set(ET.SubElement(elements, f"{{{IDF_NS}}}element"), "i:name", symbol)

    structure = _ensure(sample, "i:structure/i:layeredstructure")
    _set(structure, "i:nlayers", str(max(1, len(layers))))
    container = _ensure(structure, "i:layers")
    previous = container.findall("i:layer", _NS)
    for child in list(container):
        container.remove(child)
    for index, layer in enumerate(layers or [XnraLayer(0.0, {})]):  # SIMNRA's own empty sample
        node = ET.SubElement(container, f"{{{IDF_NS}}}layer")
        _set(node, "i:layerthickness", format_number(layer.thickness), units="1e15at/cm2")
        members = _ensure(node, "i:layerelements")
        for symbol, fraction in layer.composition.items():
            member = ET.SubElement(members, f"{{{IDF_NS}}}layerelement")
            _set(member, "i:name", symbol)
            _set(member, "i:concentration", format_number(fraction), units="fraction")
        # SIMNRA's own per-layer settings (roughness, porosity, ...) carry
        # over from the layer that stood in the same place.
        extras = [
            copy.deepcopy(child)
            for child in (previous[index] if index < len(previous) else [])
            if child.tag.startswith(f"{{{SIMNRA_NS}}}")
        ]
        if extras:
            node.extend(extras)
        else:
            _set(node, "s:hasroughness", "false")
            _set(node, "s:roughness/s:distribution", "none")
            _set(node, "s:porosity/s:hasporosity", "false")


#: The document a fresh file starts from, following the structure SIMNRA 7.04
#: itself writes. ``simnraversionnr`` names the SIMNRA release whose layout this
#: mirrors. Settings pyRUMP has no opinion on use the values SIMNRA saved in
#: tests/data/xnra/MnPt.xnra.
_SKELETON = f"""<?xml version="1.0" encoding="utf-8"?>
<idf xmlns="{IDF_NS}" xmlns:simnra="{SIMNRA_NS}" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
<users><user/></users>
<notes><note/></notes>
<attributes>
<idfversion>1.01</idfversion>
<filename/>
<createtime/>
<simnra:simnraversionnr>7.04</simnra:simnraversionnr>
<simnra:filetype>xnra</simnra:filetype>
<simnra:xnraversionnr>{XNRA_VERSION}</simnra:xnraversionnr>
</attributes>
<sample>
<users><user/></users>
<notes><note/></notes>
<description/>
<elementsandmolecules><elements><nelements>0</nelements></elements></elementsandmolecules>
<structure>
<notes><note/></notes>
<layeredstructure>
<nlayers>1</nlayers>
<layers>
<layer>
<layerthickness units="1e15at/cm2">0.00000000000000E+0000</layerthickness>
<layerelements/>
<simnra:hasroughness>false</simnra:hasroughness>
<simnra:roughness><simnra:distribution>none</simnra:distribution></simnra:roughness>
<simnra:porosity><simnra:hasporosity>false</simnra:hasporosity></simnra:porosity>
</layer>
</layers>
<simnra:substraterougness>
<simnra:distribution>none</simnra:distribution>
<simnra:width units="degree" mode="FWHM">0.00000000000000E+0000</simnra:width>
</simnra:substraterougness>
</layeredstructure>
</structure>
<spectra>
<spectrum>
<notes><note/></notes>
<log/>
<beam>
<beamparticle>4He</beamparticle>
<beamZ>2</beamZ>
<beammass units="amu">4.00260325413000E+0000</beammass>
<beamenergy units="keV">2.00000000000000E+0003</beamenergy>
<beamenergyspread mode="FWHM" units="keV">0.00000000000000E+0000</beamenergyspread>
<beamfluence units="#particles">1.00000000000000E+0000</beamfluence>
<beamangularspread units="degree" mode="FWHM">0</beamangularspread>
</beam>
<geometry>
<geometrytype>general</geometrytype>
<incidenceangle units="degree">0.00000000000000E+0000</incidenceangle>
<scatteringangle units="degree">1.70000000000000E+0002</scatteringangle>
<exitangle units="degree">1.00000000000000E+0001</exitangle>
</geometry>
<detection>
<detector>
<detectortype>SSB</detectortype>
<solidangle units="msr">1.00000000000000E+0000</solidangle>
</detector>
</detection>
<calibrations>
<detectorresolutions/>
<energycalibrations/>
</calibrations>
<reactions><technique>RBS</technique></reactions>
<data>
<datamode>simple</datamode>
<channelmode>left</channelmode>
<simpledata/>
<simnra:graphics>
<simnra:legend/>
<simnra:linecolor>red</simnra:linecolor>
<simnra:active>true</simnra:active>
<simnra:symbol><simnra:shape>circle</simnra:shape><simnra:visible>true</simnra:visible></simnra:symbol>
</simnra:graphics>
</data>
<simnra:processeddata/>
<process>
<physicsdefaults>
<crosssectiondefault>
<Rutherford>false</Rutherford>
<screening>Andersen</screening>
<simnra:isotopes>true</simnra:isotopes>
</crosssectiondefault>
</physicsdefaults>
<simulations/>
</process>
<simnra:subspectra>
<simnra:elementspectra>false</simnra:elementspectra>
<simnra:isotopespectra>false</simnra:isotopespectra>
</simnra:subspectra>
</spectrum>
</spectra>
</sample>
</idf>
"""
