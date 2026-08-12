r"""Depth-profile equations: the ``EQUATION`` forms.

Port of the two ``switch`` blocks in ``FillSimStructure`` -- setup at
creatr.c:705-810, per-sublayer evaluation at creatr.c:843-908.

A layer normally has one composition throughout. An equation replaces that with
a depth-dependent *species fraction* :math:`f(x)`, blended against the layer's
matrix (creatr.c:948):

.. math::
    c_i(x) = f(x)\,\frac{\text{species}_i}{\sum_j \text{species}_j}
           + \big(1-f(x)\big)\,\frac{\text{matrix}_i}{\sum_j \text{matrix}_j}

**Both compositions are normalised to 1 before mixing** -- so ``Si 1 O 2`` and
``Si 0.5 O 1`` describe the same material, and the equation controls only the
blend.

Two families
------------

*Point forms* evaluate :math:`f` at the **sublayer centre**,
``x = (i + 0.5)/n``. Accuracy therefore depends on having enough sublayers.

*Integral forms* (THINFILM, BURIEDTHINFILM, GAUSSIAN) instead evaluate the
**cumulative** distribution at each sublayer's **back edge** and difference it,
giving the exact dose in each sublayer regardless of how coarse the grid is.
That is stateful -- each value depends on the previous edge -- so those must be
evaluated in order.

Not implemented
---------------

``SPLINE`` and ``USEREQN`` require GENPLOT's spline fitter and expression
evaluator respectively. They raise rather than silently returning zero.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum

import numpy as np
from scipy.special import erfc, ndtr

#: Angstroms -> cm.
_ANGSTROM_CM = 1e-8

#: The C clamps |k| beyond this in the Gaussian-family forms.
_K_LIMIT = 5.0


class EquationType(Enum):
    """``EQNTYPE`` (sample.h:30-33), by the names the SIM command uses."""

    NONE = "none"
    CONSTANT = "constant"
    LINEAR = "linear"
    ERFC = "erfc"
    EXPONENTIAL = "exponential"
    SEMI_INFINITE = "semi-infinite"
    THINFILM = "thinfilm"
    BURIEDTHINFILM = "buriedthinfilm"
    GAUSSIAN = "gaussian"
    EDGEWORTH = "edgeworth"
    THICKFILM = "thickfilm"
    TIMEDEPENDENT = "timedependent"
    SPLINE = "spline"
    USEREQN = "usereqn"


#: Forms whose value is the integral over a sublayer rather than a point sample.
INTEGRAL_FORMS = frozenset(
    {EquationType.THINFILM, EquationType.BURIEDTHINFILM, EquationType.GAUSSIAN}
)

#: Forms pyRUMP cannot evaluate without GENPLOT.
UNSUPPORTED = frozenset({EquationType.SPLINE, EquationType.USEREQN})


@dataclass(slots=True)
class Profile:
    """An equation plus its parameters, as a SIM ``EQUATION`` command supplies."""

    type: EquationType
    parameters: tuple[float, ...] = ()
    tags: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.type in UNSUPPORTED:
            raise NotImplementedError(
                f"{self.type.value} needs GENPLOT's spline fitter / expression "
                "evaluator, which pyRUMP does not implement"
            )

    def parameter(self, index: int) -> float:
        return self.parameters[index] if index < len(self.parameters) else 0.0


def _resolve_length(value: float) -> float:
    """Interpret a length parameter as Angstroms, or cm if implausibly small.

    creatr.c:721-722 (and :746, :773): ``x0 = p*1E-8; if (x0 < 1E-10) x0 *= 1E8;``
    -- i.e. if converting from Angstroms yields under 1e-10 cm, assume the user
    meant centimetres and undo it. A heuristic, reproduced as written.
    """
    centimetres = value * _ANGSTROM_CM
    return value if centimetres < 1e-10 else centimetres


def species_fraction(
    profile: Profile,
    n_sublayers: int,
    layer_thickness_cm: float,
    *,
    areal_thickness: float = 0.0,
    matrix_density: float = 1.0,
    species_density: float = 1.0,
) -> np.ndarray:
    """Species fraction for each sublayer of one layer.

    Returns ``(n_sublayers,)``. Values may exceed 1 or fall below 0 for some
    parameter choices; :func:`mix_composition` applies the C's clamping.
    """
    if profile.type is EquationType.NONE:
        return np.zeros(n_sublayers, dtype=np.float64)

    p = profile.parameters
    index = np.arange(n_sublayers, dtype=np.float64)

    # Point forms sample the sublayer centre; REVERSE mirrors it (creatr.c:845).
    centre = (index + 0.5) / n_sublayers
    reverse = profile.type in {
        EquationType.ERFC,
        EquationType.SEMI_INFINITE,
        EquationType.EXPONENTIAL,
    } and profile.parameter(0) < 0
    if reverse:
        centre = 1.0 - centre
    depth_cm = centre * layer_thickness_cm

    if profile.type is EquationType.CONSTANT:
        return np.full(n_sublayers, profile.parameter(0))

    elif profile.type is EquationType.LINEAR:
        c0 = profile.parameter(0)
        return c0 + (profile.parameter(1) - c0) * centre

    elif profile.type is EquationType.ERFC:
        c0 = abs(profile.parameter(0))
        four_dt = math.sqrt(4.0 * profile.parameter(1) * profile.parameter(2))
        return c0 * erfc(depth_cm / four_dt)

    elif profile.type is EquationType.SEMI_INFINITE:
        c0 = abs(profile.parameter(0)) / 2.0
        x0 = _resolve_length(profile.parameter(3))
        four_dt = math.sqrt(4.0 * profile.parameter(1) * profile.parameter(2))
        return c0 * erfc((depth_cm - x0) / four_dt)

    elif profile.type is EquationType.EXPONENTIAL:
        c0 = abs(profile.parameter(0))
        d_over_v = profile.parameter(1) / profile.parameter(2)
        return c0 * np.exp(-depth_cm / d_over_v)

    elif profile.type is EquationType.THICKFILM:
        # NOTE: sim.htm carries the author's own disclaimer on this one --
        # "this equation makes no sense to me right now but this is how it
        # is currently implemented". Reproduced verbatim.
        dt = profile.parameter(2) * profile.parameter(3)
        sigma = math.sqrt(2.0 * dt)
        c0 = profile.parameter(0)
        c1 = profile.parameter(1) - c0
        x0 = _resolve_length(profile.parameter(4))
        return c0 + c1 * (
            2.0 - ndtr((x0 - depth_cm) / sigma) - ndtr((x0 + depth_cm) / sigma)
        )

    elif profile.type is EquationType.EDGEWORTH:
        dose = profile.parameter(0) * dose_scale(profile.type, species_density)
        sigma = profile.parameter(2) * _ANGSTROM_CM
        c0 = dose / profile.parameter(2) / math.sqrt(2.0 * math.pi)
        c0 = c0 / matrix_density
        x0 = profile.parameter(1) * _ANGSTROM_CM
        k = (depth_cm - x0) / sigma
        skew, kurtosis = profile.parameter(3), profile.parameter(4)
        value = (
            c0
            * np.exp(-k * k / 2.0)
            * (
                1.0
                + skew / 6.0 * (k * (k * k - 3.0))
                + kurtosis / 24.0 * ((k * k - 6.0) * k * k + 3.0)
                + skew * skew / 72.0 * (((k * k - 15.0) * k * k + 45.0) * k * k - 15.0)
            )
        )
        value = np.where(np.abs(k) > _K_LIMIT, 0.0, value)
        return np.maximum(value, 0.0)

    elif profile.type is EquationType.TIMEDEPENDENT:
        # creatr.c's own comment: "Don't understand this one but keep for
        # those that do."
        c0 = profile.parameter(0)
        x0 = (
            profile.parameter(3) ** 2
            * profile.parameter(2)
            / profile.parameter(1)
        )
        sigma = math.sqrt(2.0 * x0)
        d_over_v = profile.parameter(1) / profile.parameter(3)
        k = depth_cm / d_over_v
        return c0 * (
            np.exp(-k) * ndtr(-(k - x0) / sigma) + ndtr(-(k + x0) / sigma)
        )

    elif profile.type in INTEGRAL_FORMS:
        return _integral_form(
            profile,
            n_sublayers,
            layer_thickness_cm,
            areal_thickness,
            species_density,
        )

    raise ValueError(f"unhandled equation type {profile.type}")


#: Dose units per equation, from ``eqlist``'s ``dflt_units`` (sim2.c:101-131).
#: THINFILM/BURIEDTHINFILM take a dose in Angstroms of species; the Gaussian
#: family takes it directly in 1e15 at/cm^2.
_DOSE_IN_ANGSTROMS = frozenset(
    {EquationType.THINFILM, EquationType.BURIEDTHINFILM}
)


def dose_scale(equation: EquationType, species_density: float) -> float:
    """``sthick_to_cm2`` (creatr.c:651, via ``SimThickConvert``).

    For Angstrom units the dose is multiplied by the species' atomic density,
    converting A x 1e23/cm^3 into 1e15/cm^2. For ATOMIC units it is already in
    1e15/cm^2 and the factor is 1.
    """
    return species_density if equation in _DOSE_IN_ANGSTROMS else 1.0


def _integral_form(
    profile: Profile,
    n_sublayers: int,
    layer_thickness_cm: float,
    areal_thickness: float,
    species_density: float,
) -> np.ndarray:
    """THINFILM / BURIEDTHINFILM / GAUSSIAN (creatr.c:872-880).

    These evaluate the cumulative normal at each sublayer's **back edge** and
    difference successive values, so each sublayer receives its exact share of
    the dose however coarse the grid. The running ``integ_frac`` makes this
    order-dependent.
    """
    scale = dose_scale(profile.type, species_density)
    per_sublayer = areal_thickness / n_sublayers

    if profile.type is EquationType.THINFILM:
        x0 = 0.0
        sigma = math.sqrt(2.0 * profile.parameter(1) * profile.parameter(2))
        # Doubled because only half the Gaussian lies inside the layer.
        dose = 2.0 * profile.parameter(0) * scale / per_sublayer
        running = 0.5
    elif profile.type is EquationType.BURIEDTHINFILM:
        x0 = _resolve_length(profile.parameter(3))
        sigma = math.sqrt(2.0 * profile.parameter(1) * profile.parameter(2))
        dose = profile.parameter(0) * scale / per_sublayer
        running = float(ndtr(-x0 / sigma))
    else:  # GAUSSIAN -- parameter 2 is FWHM, not sigma
        x0 = profile.parameter(1) * _ANGSTROM_CM
        sigma = profile.parameter(2) / math.sqrt(8.0 * math.log(2.0)) * _ANGSTROM_CM
        dose = profile.parameter(0) * scale / per_sublayer
        running = float(ndtr(-x0 / sigma))

    out = np.zeros(n_sublayers, dtype=np.float64)
    for i in range(n_sublayers):
        back_edge = (i + 1.0) / n_sublayers * layer_thickness_cm
        k = (back_edge - x0) / sigma
        if k < -_K_LIMIT or (k > _K_LIMIT and running > 0.999):
            out[i] = 0.0
            continue
        cumulative = float(ndtr(k))
        out[i] = dose * (cumulative - running)
        running = cumulative
    return out


def mix_composition(
    fraction: np.ndarray,
    matrix: np.ndarray,
    species: np.ndarray,
) -> np.ndarray:
    """Blend matrix and species by the species fraction (creatr.c:927-951).

    ``matrix`` and ``species`` are ``(n_element,)`` and are normalised here;
    ``fraction`` is ``(n_slab,)``. Returns ``(n_slab, n_element)`` fractional
    composition.

    Clamping follows the C: the species fraction is capped at 1 but *negatives
    are allowed through*, with any resulting negative element density clipped to
    zero afterwards.
    """
    fraction = np.atleast_1d(np.asarray(fraction, dtype=np.float64))
    matrix = np.asarray(matrix, dtype=np.float64)
    species = np.asarray(species, dtype=np.float64)

    matrix_sum = matrix.sum()
    species_sum = species.sum()
    matrix_norm = matrix / matrix_sum if matrix_sum else matrix
    species_norm = species / species_sum if species_sum else species

    capped = np.minimum(fraction, 1.0)
    composition = (
        capped[:, None] * species_norm[None, :]
        + (1.0 - capped)[:, None] * matrix_norm[None, :]
    )
    return np.maximum(composition, 0.0)


#: Per-equation sublayer counts from ``eqlist`` (sim2.c:101-131). When a layer
#: carries an equation these override ``maxpath`` entirely (creatr.c:700).
RECOMMENDED_SUBLAYERS: dict[EquationType, int] = {
    EquationType.NONE: 0,
    EquationType.CONSTANT: 5,
    EquationType.LINEAR: 10,
    EquationType.ERFC: 10,
    EquationType.EXPONENTIAL: 20,
    EquationType.SEMI_INFINITE: 20,
    EquationType.THINFILM: 30,
    EquationType.BURIEDTHINFILM: 30,
    EquationType.THICKFILM: 20,
    EquationType.TIMEDEPENDENT: 20,
    EquationType.GAUSSIAN: 20,
    EquationType.EDGEWORTH: 30,
}


def recommended_sublayers(equation: EquationType) -> int:
    return RECOMMENDED_SUBLAYERS.get(equation, 0)
