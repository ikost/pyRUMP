r"""Lateral thickness non-uniformity ("fuzz").

Port of the FUZZ block in ``FillSimStructure`` (creatr.c:662-686) and the
iteration count at creatr.c:433-436.

A real sample is not perfectly flat. RUMP models roughness by simulating the
whole sample several times with the layer thickness perturbed, and summing the
results with Gaussian weights. Every fuzzed layer multiplies the iteration
count, so two layers with 5 steps each means 25 full simulations -- it is the
most expensive option in the program.

The quantiles come from the inverse normal CDF, arranged so the perturbations
are symmetric about zero:

.. math::
    \Delta t = \pm\, \text{fuzz} \cdot \Phi^{-1}\!\left(
        \frac{j(j-1) + \tfrac12}{w}\right) \cdot \tfrac{1}{\sqrt 2}

The :math:`1/\sqrt2` (the C's literal ``.7071067``) means the resulting
distribution has standard deviation ``fuzz/sqrt(2)``, not ``fuzz`` -- so the
parameter is not the roughness sigma directly.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import ndtri

#: The C's literal for 1/sqrt(2) (creatr.c:682).
_INV_SQRT2 = 0.7071067


@dataclass(frozen=True, slots=True)
class FuzzStep:
    """One perturbation of a fuzzed layer."""

    delta: float
    """Thickness offset, in the layer's own units."""

    weight: float
    """Relative weight of this replica; weights over a layer sum to 1."""


def fuzz_steps(amount: float, steps: int) -> list[FuzzStep]:
    """Perturbations and weights for one fuzzed layer (creatr.c:668-684).

    ``steps`` of 0 or 1 means no fuzzing: a single unperturbed replica.
    """
    if steps <= 1 or amount == 0.0:
        return [FuzzStep(delta=0.0, weight=1.0)]

    half = steps // 2
    # Normalisation differs for even and odd step counts (creatr.c:672-676).
    size = 0.5 / (half * half) if 2 * half == steps else 1.0 / (2 * half * half + steps)

    out: list[FuzzStep] = []
    for k in range(steps):
        j = min(k + 1, steps - k)
        y = (j * (j - 1) + 0.5) * size
        x = float(ndtri(y))
        if k >= half:
            x = -x
        out.append(FuzzStep(delta=amount * x * _INV_SQRT2, weight=(2 * j - 1) * size))
    return out


def iteration_count(steps_per_layer: list[int]) -> int:
    """Total replicas: the product over fuzzed layers (creatr.c:433-436)."""
    total = 1
    for steps in steps_per_layer:
        if steps > 0:
            total *= steps
    return total


def replica_thicknesses(
    thicknesses: np.ndarray,
    amounts: list[float],
    steps_per_layer: list[int],
) -> list[tuple[np.ndarray, float]]:
    """Enumerate every replica as ``(thicknesses, amplitude)``.

    The C walks a mixed-radix counter (``runner``/``runner/fuzzs``,
    creatr.c:666-685) so each iteration selects one step from each fuzzed layer;
    this reproduces that enumeration explicitly.

    Amplitudes over all replicas sum to 1, so a fuzzed simulation preserves the
    total yield of an unfuzzed one.
    """
    thicknesses = np.asarray(thicknesses, dtype=np.float64)
    per_layer = [
        fuzz_steps(amounts[i], steps_per_layer[i]) for i in range(thicknesses.size)
    ]

    out: list[tuple[np.ndarray, float]] = []
    for iteration in range(iteration_count(steps_per_layer)):
        adjusted = thicknesses.copy()
        amplitude = 1.0
        runner = iteration
        for index, steps in enumerate(steps_per_layer):
            if steps <= 0:
                continue
            step = per_layer[index][runner % steps]
            adjusted[index] += step.delta
            amplitude *= step.weight
            runner //= steps
        out.append((adjusted, amplitude))
    return out
