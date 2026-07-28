"""Absorber layers: dead layers, windows and stopping foils.

Port of the ``fsurf`` handling in ``FillSimStructure`` (creatr.c:815-838) and
``SimFlyout`` (creatr.c:1971).

The first ``absorber_layers`` layers of a sample are not part of the sample at
all -- they represent material between it and the detector, such as a silicon
dead layer or a Mylar window. Three consequences, all of which the simulation
must honour:

1. **Nothing scatters from them.** The depth loop starts at ``samm->fsurf``.
2. **They are traversed on the exit path only**, so they attenuate the outgoing
   particle without contributing yield.
3. **They are not tilted with the sample.** ``SimFlyout`` forces normal
   incidence through them (``if (layer < samm->fsurf) secout = 1E-3f``), because
   a detector window does not rotate when the sample does. This is easy to miss
   and matters as soon as the sample is tilted.

They also do not advance the depth counters (creatr.c:838), so a profile in a
following layer measures depth from the true sample surface.

``fres_only_absorber`` restricts the absorber to forward-recoil spectra
(creatr.c:816); not implemented, since ERD is not yet supported.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, slots=True)
class AbsorberSpec:
    """How many leading layers are absorber rather than sample."""

    layers: int = 0
    fres_only: bool = False

    def __post_init__(self) -> None:
        if self.fres_only:
            raise NotImplementedError(
                "fres_only_absorber applies to forward-recoil spectra, which "
                "pyRUMP does not yet simulate"
            )


def first_sample_slab(layer_index: np.ndarray, absorber_layers: int) -> int:
    """Index of the first slab that belongs to the sample proper (``fsurf``).

    Slabs before this are absorber: traversed on the way out, never scattered
    from.
    """
    if absorber_layers <= 0:
        return 0
    beyond = np.flatnonzero(layer_index >= absorber_layers)
    return int(beyond[0]) if beyond.size else int(layer_index.size)


def sample_depth(areal_density: np.ndarray, first_slab: int) -> np.ndarray:
    """Cumulative depth measured from the sample surface, not the absorber's.

    Absorber slabs do not advance the depth counters (creatr.c:838).
    """
    depth = np.zeros_like(areal_density)
    depth[first_slab:] = np.cumsum(areal_density[first_slab:])
    return depth
