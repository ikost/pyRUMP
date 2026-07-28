"""Slab discretization: turning user layers into the simulation's sublayers.

Port of the geometry/thickness half of ``FillSimStructure`` (creatr.c:548-710).

The step size is **path-length based, not depth based**:

.. code-block:: c

    maxpath = sample->maxpth / max(|secin|, |secout|)   /* creatr.c:587 */
    num_sublayers = (int)(1 + cm2_thick / maxpath)      /* creatr.c:702 */

so a tilted sample is automatically cut into more, thinner slabs. ``maxpth``
defaults to 200 (1e15 at/cm^2), which the 1985 paper justifies as the point where
the third-order energy-loss expansion reaches ~1e-5 fractional error.

Note the ``int()`` truncation: a layer exactly ``maxpath`` thick gets 2 sublayers,
not 1. Reproduced deliberately.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..profiles.equations import (
    mix_composition,
    recommended_sublayers as _recommended_sublayers,
    species_fraction,
)

#: Default ``maxpth`` in 1e15 atoms/cm^2 (sample.h:16).
DEFAULT_MAXPATH = 200.0


@dataclass(slots=True)
class SlabGrid:
    """Flattened sublayers, the array the simulation actually marches over."""

    areal_density: np.ndarray
    """``(n_slab,)`` total areal density per slab, 1e15 atoms/cm^2."""

    composition: np.ndarray
    """``(n_slab, n_element)`` per-element areal density, 1e15 atoms/cm^2."""

    layer_index: np.ndarray
    """``(n_slab,)`` index of the originating user layer."""

    element_z: list[int]
    """Atomic numbers, one per column of :attr:`composition`."""

    @property
    def n_slab(self) -> int:
        return int(self.areal_density.size)

    @property
    def depth(self) -> np.ndarray:
        """Cumulative areal density at each slab's *back* face."""
        return np.cumsum(self.areal_density)

    @property
    def total_thickness(self) -> float:
        return float(self.areal_density.sum())


def sublayer_count(
    areal_thickness: float,
    maxpath: float,
    *,
    explicit: int = 0,
    sub_thickness: float = 0.0,
    equation_sublayers: int = 0,
) -> int:
    """How many sublayers a layer is cut into (creatr.c:693-703).

    Precedence: an explicit count, then an explicit sub-thickness, then the
    equation's *recommended* count, then the ``maxpath`` default.

    That third case is easy to miss: once a layer carries an equation, RUMP
    ignores ``maxpath`` entirely and uses the per-equation value from
    ``eqlist`` (sim2.c:101) -- 5 for CONSTANT, 30 for THINFILM, and so on.
    """
    if explicit:
        return explicit
    if sub_thickness:
        return int(1 + areal_thickness / sub_thickness)
    if equation_sublayers > 0:
        return equation_sublayers
    return int(1 + areal_thickness / maxpath)


def path_limited_maxpath(maxpth: float, sec_in: float, sec_out: float) -> float:
    """``maxpth`` reduced by the longer of the two path secants (creatr.c:587)."""
    return maxpth / max(abs(sec_in), abs(sec_out))


def build_uniform_grid(
    layer_thicknesses: np.ndarray,
    layer_compositions: np.ndarray,
    element_z: list[int],
    *,
    maxpth: float = DEFAULT_MAXPATH,
    sec_in: float = 1.0,
    sec_out: float = 1.0,
    explicit_sublayers: list[int] | None = None,
    profiles: "list | None" = None,
    layer_species: np.ndarray | None = None,
    layer_densities: np.ndarray | None = None,
    species_densities: np.ndarray | None = None,
) -> SlabGrid:
    """Build a slab grid from layers of uniform composition.

    Parameters
    ----------
    layer_thicknesses:
        ``(n_layer,)`` areal thickness per layer, 1e15 atoms/cm^2.
    layer_compositions:
        ``(n_layer, n_element)`` *fractional* composition per layer. Rows are
        normalised, matching RUMP's treatment of matrix and species.
    element_z:
        Atomic numbers, one per composition column.

    A layer may carry a :class:`~pyrump.profiles.equations.Profile`, in which
    case its composition varies with depth and its sublayer count comes from the
    equation rather than from ``maxpth``.
    """
    thicknesses = np.atleast_1d(np.asarray(layer_thicknesses, dtype=np.float64))
    compositions = np.atleast_2d(np.asarray(layer_compositions, dtype=np.float64))
    if compositions.shape[0] != thicknesses.size:
        raise ValueError(
            f"{thicknesses.size} layers but {compositions.shape[0]} composition rows"
        )
    if compositions.shape[1] != len(element_z):
        raise ValueError(
            f"{compositions.shape[1]} composition columns but {len(element_z)} elements"
        )

    maxpath = path_limited_maxpath(maxpth, sec_in, sec_out)
    explicit = explicit_sublayers or [0] * thicknesses.size
    profiles = profiles or [None] * thicknesses.size

    areal: list[float] = []
    per_element: list[np.ndarray] = []
    origin: list[int] = []

    for index, thickness in enumerate(thicknesses):
        if thickness <= 0:
            continue  # creatr.c:689 skips empty layers outright
        profile = profiles[index]
        count = sublayer_count(
            float(thickness),
            maxpath,
            explicit=explicit[index],
            equation_sublayers=(
                _recommended_sublayers(profile.type) if profile is not None else 0
            ),
        )
        slab_thickness = float(thickness) / count

        row = compositions[index]
        total = row.sum()
        fractions = row / total if total > 0 else row

        if profile is None:
            per_element.extend([fractions * slab_thickness] * count)
        else:
            density = (
                float(layer_densities[index]) if layer_densities is not None else 1.0
            )
            # Layer thickness in cm, needed by the depth-dependent forms
            # (creatr.c:692): areal density / atomic density / 1e8.
            thickness_cm = float(thickness) / density / 1e8
            fraction = species_fraction(
                profile,
                count,
                thickness_cm,
                areal_thickness=float(thickness),
                matrix_density=density,
                species_density=(
                    float(species_densities[index])
                    if species_densities is not None
                    else 1.0
                ),
            )
            species = (
                layer_species[index]
                if layer_species is not None
                else np.zeros_like(row)
            )
            blended = mix_composition(fraction, row, species)
            per_element.extend(list(blended * slab_thickness))

        areal.extend([slab_thickness] * count)
        origin.extend([index] * count)

    if not areal:
        raise ValueError("sample has no layers of positive thickness")

    return SlabGrid(
        areal_density=np.array(areal, dtype=np.float64),
        composition=np.array(per_element, dtype=np.float64),
        layer_index=np.array(origin, dtype=np.int64),
        element_z=list(element_z),
    )
