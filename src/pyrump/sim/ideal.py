"""The per-isotope depth loop: slabs in, bricks out.

Port of ``SimCideal`` (creatr.c:1594-1930), restricted to the backscattering
(RBS) case with Rutherford or screened-Rutherford cross-sections. Recoil and
tabulated-resonance paths are later milestones.

Yield, per creatr.c:1809-1823:

.. code-block:: text

    r      = N_slab * ratde * (sec_in / [eps]) * isotope_fraction
    height = BSCALE * sigma(E) * r

``BSCALE = 6.241507`` converts barns/sr into RUMP's internal
counts/uC/msr/(1e15 at/cm^2) -- essentially 1/e in picocoulombs.

The loop is sequential per slab because ``flyout`` walks the whole overlying
stack each time. It is the O(n_slab^2) term the 1985 paper identifies as
dominating, and it is why RUMP works hard to keep the slab count low.
"""

from __future__ import annotations

import numpy as np

from ..physics.kinematics import kinematic_factor
from ..physics.xsec.rutherford import CrossSection
from ..stopping.table import StoppingTable
from .bricks import Bricks
from .outbound import epsilon_factor, flyout
from .precal import InboundPath

#: barns/sr -> counts/uC/msr/(1e15 at/cm^2) (creatr.c:1613).
BSCALE = 6.241507


def simulate_isotope(
    table: StoppingTable,
    coefficients_in: np.ndarray,
    coefficients_out: np.ndarray,
    slab_element_density: np.ndarray,
    inbound: InboundPath,
    cross_section: CrossSection,
    *,
    m_beam: float,
    m_target: float,
    scattering_angle_deg: float,
    isotope_fraction: float,
    sec_in: float,
    sec_out: float,
    cutoff_keV: float,
    straggle_geometry: float | None = None,
    first_slab: int = 0,
) -> Bricks:
    """Produce one brick per slab for a single target isotope.

    Parameters
    ----------
    slab_element_density:
        ``(n_slab,)`` areal density of *this element* per slab, 1e15 at/cm^2.
    inbound:
        Result of :func:`pyrump.sim.precal.march_inbound`.
    straggle_geometry:
        The ``stragc = sec_in*K + sec_out`` factor. ``None`` disables straggling,
        which is RUMP's default.
    """
    kinematic = kinematic_factor(m_beam, m_target, scattering_angle_deg)
    n_slab = slab_element_density.size

    rows: list[tuple] = []
    front_valid = False
    e_front = h_front = 0.0
    ratde = 1.0

    # Absorber slabs sit in front of the sample: traversed on the way out by
    # flyout, but never scattered from (creatr.c:1033).
    for slab in range(first_slab, n_slab):
        if slab >= inbound.reached:
            break
        if slab_element_density[slab] <= 0:
            # Nothing of this element here; the cached front edge goes stale.
            front_valid = False
            continue

        energy_in = inbound.energy[slab]

        # Only the *geometry* of the front edge is reused between adjacent
        # slabs: e_front and ratde carry over because the outward path through
        # the overlying stack is unchanged (creatr.c:1791-1795, guarded by `ok`).
        if not front_valid:
            e_front, ratde = flyout(
                table,
                coefficients_out,
                slab - 1,
                kinematic * energy_in,
                sec_out=sec_out,
                cutoff_keV=cutoff_keV,
                first_surface=first_slab,
            )
            if e_front <= cutoff_keV:
                break

        # The front *height*, by contrast, is recomputed every slab from that
        # slab's own areal density (creatr.c:1804-1823, outside the `ok` guard).
        # With uniform composition this is indistinguishable from carrying it
        # over; with a depth profile it is not, and reusing it shifts the whole
        # spectrum by one slab.
        eps = epsilon_factor(
            table,
            coefficients_in[slab],
            coefficients_out[slab],
            energy_in,
            kinematic,
            sec_in=sec_in,
            sec_out=sec_out,
        )
        r_front = (
            slab_element_density[slab] * ratde * (sec_in / eps) * isotope_fraction
        )
        h_front = BSCALE * float(cross_section(energy_in)[0]) * r_front

        # Back edge of this slab.
        energy_out = inbound.energy[slab + 1]
        e_back, ratde = flyout(
            table,
            coefficients_out,
            slab,
            kinematic * energy_out,
            sec_out=sec_out,
            cutoff_keV=cutoff_keV,
            first_surface=first_slab,
        )
        if e_back < cutoff_keV:
            break

        eps_back = epsilon_factor(
            table,
            coefficients_in[slab],
            coefficients_out[slab],
            energy_out,
            kinematic,
            sec_in=sec_in,
            sec_out=sec_out,
        )
        r_back = (
            slab_element_density[slab] * ratde * (sec_in / eps_back) * isotope_fraction
        )
        h_back = BSCALE * float(cross_section(energy_out)[0]) * r_back

        # The paper's exact area, from the pre-computed Rutherford integral.
        # The shipped fill stage ignores it; kept for the parabolic mode.
        area = (
            BSCALE
            * cross_section.csigma
            * slab_element_density[slab]
            * isotope_fraction
            * sec_in
        )

        if straggle_geometry is None:
            sig_front = sig_back = 0.0
        else:
            sig_front = straggle_geometry * inbound.straggle[slab]
            sig_back = straggle_geometry * inbound.straggle[slab + 1]

        rows.append((e_front, e_back, h_front, h_back, area, sig_front, sig_back))

        # The next slab's front edge is this slab's back edge -- bricks tile in
        # energy. Only the geometry carries; the height is recomputed above.
        e_front = e_back
        front_valid = True

    return Bricks.from_list(rows) if rows else Bricks.empty(0)


def straggle_geometry_factor(kinematic: float, sec_in: float, sec_out: float) -> float:
    """RUMP's in/out straggling combination (creatr.c:1661).

    .. warning::
       This is ``sec_in*K + sec_out`` applied *linearly* to the inbound Bohr
       variance. The physically correct combination is
       ``K^2 * sigma_in^2 + sigma_out^2``. Reproduced as-is; see
       ``docs/rump-quirks.md`` entry 9.
    """
    return sec_in * kinematic + sec_out
