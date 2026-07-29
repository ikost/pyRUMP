"""Binding named fit parameters onto the simulation inputs.

Port of the ``VARY`` parameter table (``pert.h:1-7``, ``pert.c:156-168``).

A fit varies quantities scattered across several immutable dataclasses — beam
energy lives on :class:`Beam`, detector resolution on :class:`Measurement`,
thickness inside a :class:`UniformSample`. Rather than mutate those, each
parameter is a getter/setter pair over a *bundle* of inputs, so the optimiser
sees a flat vector and the simulation always receives well-formed objects.

Parameter names follow RUMP's, so a ``.pert`` file's ``VARY`` list maps across
directly.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Callable

import numpy as np


@dataclass(slots=True)
class FitInputs:
    """Everything a simulation call needs, as one mutable bundle."""

    sample: object
    beam: object
    geometry: object
    calibration: object
    measurement: object


@dataclass(slots=True)
class Parameter:
    """One fittable quantity."""

    name: str
    get: Callable[[FitInputs], float]
    set: Callable[[FitInputs, float], None]
    lower: float = -np.inf
    upper: float = np.inf
    scale: float = 1.0
    """Typical magnitude, used to size the numerical-derivative step."""


def _beam(field: str) -> tuple[Callable, Callable]:
    return (
        lambda inp: getattr(inp.beam, field),
        lambda inp, v: setattr(inp.beam, field, v),
    )


def _measurement(field: str) -> tuple[Callable, Callable]:
    # Measurement is frozen, so setting replaces the whole object.
    return (
        lambda inp: getattr(inp.measurement, field),
        lambda inp, v: setattr(inp, "measurement", replace(inp.measurement, **{field: v})),
    )


def _calibration(field: str) -> tuple[Callable, Callable]:
    return (
        lambda inp: getattr(inp.calibration, field),
        lambda inp, v: setattr(inp, "calibration", replace(inp.calibration, **{field: v})),
    )


def _geometry(field: str) -> tuple[Callable, Callable]:
    return (
        lambda inp: getattr(inp.geometry, field),
        lambda inp, v: setattr(inp, "geometry", replace(inp.geometry, **{field: v})),
    )


def _sample(field: str) -> tuple[Callable, Callable]:
    return (
        lambda inp: getattr(inp.sample, field),
        lambda inp, v: setattr(inp.sample, field, v),
    )


def thickness(layer: int) -> Parameter:
    """Areal thickness of one layer, in 1e15 at/cm^2."""

    def get(inp: FitInputs) -> float:
        return float(inp.sample.thicknesses[layer])

    def set_(inp: FitInputs, value: float) -> None:
        thicknesses = list(inp.sample.thicknesses)
        thicknesses[layer] = value
        inp.sample.thicknesses = thicknesses

    return Parameter(f"thickness[{layer}]", get, set_, lower=0.0, scale=100.0)


def composition(layer: int, element: int) -> Parameter:
    """Stoichiometric coefficient of one element in one layer."""

    def get(inp: FitInputs) -> float:
        return float(inp.sample.compositions[layer][element])

    def set_(inp: FitInputs, value: float) -> None:
        rows = [list(r) for r in inp.sample.compositions]
        rows[layer][element] = value
        inp.sample.compositions = rows

    return Parameter(f"composition[{layer},{element}]", get, set_, lower=0.0, scale=1.0)


def equation_parameter(layer: int, index: int) -> Parameter:
    """One parameter of a layer's depth-profile equation."""

    def get(inp: FitInputs) -> float:
        return float(inp.sample.profiles[layer].parameters[index])

    def set_(inp: FitInputs, value: float) -> None:
        profile = inp.sample.profiles[layer]
        params = list(profile.parameters)
        params[index] = value
        profiles = list(inp.sample.profiles)
        profiles[layer] = replace(profile, parameters=tuple(params))
        inp.sample.profiles = profiles

    return Parameter(f"equation[{layer},{index}]", get, set_)


#: Parameters that need no arguments, by RUMP's own names.
SIMPLE_PARAMETERS: dict[str, tuple[tuple[Callable, Callable], float, float, float]] = {
    "mev": (_beam("e0_MeV"), 0.0, np.inf, 1.0),
    "theta": (_geometry("theta"), -89.9, 89.9, 10.0),
    "phi": (_geometry("phi"), 0.0, 179.9, 10.0),
    "psi": (_geometry("psi"), -89.9, 89.9, 10.0),
    "fwhm": (_measurement("fwhm_keV"), 0.0, np.inf, 10.0),
    "tau": (_measurement("tau_us"), 0.0, np.inf, 1.0),
    "current": (_measurement("current_nA"), 0.0, np.inf, 10.0),
    "correction": (_measurement("correction"), 0.0, np.inf, 1.0),
    "kev/ch": (_calibration("kevch"), 0.0, np.inf, 1.0),
    "kev(0)": (_calibration("kev0"), -np.inf, np.inf, 10.0),
    "straggle": (_sample("straggle"), 0.0, np.inf, 1.0),
    "multiple_scatter": (_sample("multiple"), 0.0, np.inf, 1.0),
}


def parameter(name: str) -> Parameter:
    """Look up a no-argument parameter by RUMP's name."""
    key = name.lower()
    if key not in SIMPLE_PARAMETERS:
        raise KeyError(
            f"unknown fit parameter {name!r}; known: "
            + ", ".join(sorted(SIMPLE_PARAMETERS))
            + ", plus thickness(), composition() and equation_parameter()"
        )
    (get, set_), lower, upper, scale = SIMPLE_PARAMETERS[key]
    return Parameter(key, get, set_, lower=lower, upper=upper, scale=scale)


def pack(parameters: list[Parameter], inputs: FitInputs) -> np.ndarray:
    return np.array([p.get(inputs) for p in parameters], dtype=np.float64)


def unpack(parameters: list[Parameter], inputs: FitInputs, values: np.ndarray) -> None:
    for p, value in zip(parameters, values):
        p.set(inputs, float(value))


def bounds(parameters: list[Parameter]) -> tuple[np.ndarray, np.ndarray]:
    return (
        np.array([p.lower for p in parameters], dtype=np.float64),
        np.array([p.upper for p in parameters], dtype=np.float64),
    )
