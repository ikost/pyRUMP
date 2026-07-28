"""Bragg-rule mixing of stopping powers across a slab's composition.

Port of ``CalcAverageStop`` (creatr.c:1096-1122).

RUMP applies plain linear additivity of elemental stopping cross-sections,
weighted by areal density -- **no compound (CAB) correction**. Because the weights
are areal densities in 1e15 at/cm^2 rather than fractions, the summed coefficients
have units of *eV through the slab*, not eV/(1e15 at/cm^2). The depth stepping
then needs no further thickness factor.

Since the mixing is linear and the stopping model is a polynomial in a fixed
variable, the sum can be taken over the coefficients once per slab rather than
over energies at every step. That precomputed ``(n_slab, 6)`` array is the hottest
data structure in the simulation.
"""

from __future__ import annotations

import numpy as np

from .table import NDEG, StoppingTable


def bragg_coefficients(
    table: StoppingTable,
    areal_density: np.ndarray,
    element_z: list[int],
) -> np.ndarray:
    """Combine per-element polynomials into per-slab polynomials.

    Parameters
    ----------
    table:
        Fitted stopping polynomials, one per target element.
    areal_density:
        ``(n_slab, n_element)`` areal densities in 1e15 atoms/cm^2.
    element_z:
        Atomic numbers, one per column of ``areal_density``.

    Returns
    -------
    ``(n_slab, NDEG)`` coefficients whose evaluation gives energy loss in eV
    across the whole slab.
    """
    areal_density = np.atleast_2d(np.asarray(areal_density, dtype=np.float64))
    if areal_density.shape[1] != len(element_z):
        raise ValueError(
            f"areal_density has {areal_density.shape[1]} columns but "
            f"{len(element_z)} elements were given"
        )

    # (n_element, NDEG) stack of elemental coefficients.
    per_element = np.array(
        [table.coefficients[z] for z in element_z], dtype=np.float64
    )
    # The C skips elements with zero areal density (creatr.c:1098); with a plain
    # matrix product those contribute nothing anyway.
    return areal_density @ per_element


def evaluate_slab_loss(
    table: StoppingTable, coefficients: np.ndarray, energy_keV
) -> np.ndarray:
    """Evaluate slab-summed coefficients at one energy per slab.

    Returns energy loss in eV for a normal-incidence traversal; the caller
    applies the path-length secant.
    """
    x = table.transform(energy_keV)
    x = np.atleast_1d(x)
    powers = np.vander(x, NDEG, increasing=True)
    return np.einsum("ij,ij->i", powers, np.atleast_2d(coefficients))
