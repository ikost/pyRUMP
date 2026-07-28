"""ZBL85 (Ziegler-Biersack-Littmark) stopping powers.

Port of ``ziegler.c`` (``zstop``/``pstop``/``hestop``/``histop``), which is itself a
transliteration of Ziegler's FORTRAN. All coefficients come from ``pscoef.dat``
(the TRIM/SRIM ``SCOEF`` table); nothing is hard-coded here.

Three projectile regimes, selected by Z1:

* ``Z1 == 1`` - Andersen-Ziegler proton fit, the base for everything else
* ``Z1 == 2`` - helium, via an He/H stopping-power *ratio*
* ``Z1 <= 92`` - heavy ions, via Brandt-Kitagawa effective charge

Plus universal ZBL nuclear stopping, which is projectile-independent.

Energies are in keV; the natural result unit is eV/(1e15 atoms/cm^2), which is what
the simulation kernel works in.

.. note::
   Accuracy is Ziegler's own: "mean accuracy 9%", valid below 100 000 keV/amu,
   and only for Z = 1..92. RUMP's Mylar pseudo-element (Z=93) has no Ziegler
   path -- see :mod:`pyrump.stopping.mylar`.
"""

from __future__ import annotations

from enum import IntEnum

import numpy as np

from ..model.element import Element

#: Ziegler's tables stop at uranium.
MAX_Z = 92

#: Reduced energy above which the parameterisation is not valid (keV/amu).
MAX_REDUCED_ENERGY = 100_000.0


class StoppingUnits(IntEnum):
    """``zstop``'s ``units`` argument (ziegler.c:259-278)."""

    EV_1E15_ATOMS = 1
    """eV/(1e15 atoms/cm^2) -- the simulation's internal unit."""

    MEV_CM2_MG = 2
    EV_PER_ANGSTROM = 3
    LSS_REDUCED = 4


def _as_array(energy_keV) -> np.ndarray:
    return np.atleast_1d(np.asarray(energy_keV, dtype=np.float64))


def proton_stopping(e_per_amu: np.ndarray, z2: int, pcoef: np.ndarray) -> np.ndarray:
    """Andersen-Ziegler proton electronic stopping (``pstop``, ziegler.c:394).

    ``e_per_amu`` is keV/amu. Below 25 keV/amu the fit is replaced by a
    velocity-proportional extrapolation.
    """
    peo = 25.0
    pe = np.maximum(peo, e_per_amu)

    # Low- and high-energy branches joined harmonically.
    sl = pcoef[0] * pe ** pcoef[1] + pcoef[2] * pe ** pcoef[3]
    sh = pcoef[4] / pe ** pcoef[5] * np.log(pcoef[6] / pe + pcoef[7] * pe)
    se = sl * sh / (sl + sh)

    # Velocity-proportional below peo; the exponent switches at Z2 = 6.
    velocity_power = 0.45 if z2 > 6 else 0.25
    return np.where(
        e_per_amu <= peo,
        se * (e_per_amu / peo) ** velocity_power,
        se,
    )


def helium_stopping(
    e_per_amu: np.ndarray, z1: int, z2: int, pcoef: np.ndarray
) -> np.ndarray:
    """Helium electronic stopping via the He/H ratio (``hestop``, ziegler.c:291)."""
    heo = 1.0
    he = np.maximum(heo, e_per_amu)
    b = np.log(he)

    a = (
        0.744647
        + 0.142913 * b
        + 0.0156235 * b * b
        - 0.0026665 * b**3
        + 1.32512e-6 * b**8
    )
    heh = 1.0 - np.exp(-np.minimum(30.0, a))

    # Z1^3 (Barkas) correction to the ratio.
    a = (7.6 - np.maximum(0.0, np.log(he))) ** 2
    heh = heh * (1.0 + (0.007 + 0.00005 * z2) * np.exp(-a))

    se = proton_stopping(he, z2, pcoef) * (z1 * heh) ** 2
    return np.where(e_per_amu <= heo, se * np.sqrt(e_per_amu / heo), se)


def heavy_ion_stopping(
    e_per_amu: np.ndarray,
    z1: int,
    z2: int,
    pcoef: np.ndarray,
    fermi_velocity: float,
    lambda_screening: float,
) -> np.ndarray:
    """Heavy-ion stopping via Brandt-Kitagawa effective charge (``histop``).

    ``fermi_velocity`` is the target's; ``lambda_screening`` is the *projectile's*
    (ziegler.c:215-221) -- an easy detail to get backwards.
    """
    yrmin, vrmin = 0.13, 1.0

    v = np.sqrt(e_per_amu / 25.0) / fermi_velocity
    vr = np.where(
        v < 1.0,
        (3 * fermi_velocity / 4) * (1 + 2 * v * v / 3 - v**4 / 15),
        v * fermi_velocity * (1 + 1 / (5 * np.maximum(v, 1e-30) ** 2)),
    )

    z1_23 = z1**0.6667
    yr = np.maximum(yrmin, vr / z1_23)
    yr = np.maximum(yr, vrmin / z1_23)

    a = -0.003845 / yr - 0.09876 + 1.0406 * yr - 0.08483 * yr * yr + 0.01294 * yr**3
    q = np.minimum(1.0, np.maximum(0.0, 1.0 - np.exp(-np.minimum(a, 50.0))))

    # Ionization level -> effective charge (screening length l).
    b = 0.26 - 0.004 * z1
    l0 = 2 * 0.24 * (1.0 - q) ** 0.6667 / (z1**0.3333 * (1.0 - 0.143 * (1.0 - q)))

    q1 = max(0.0, 0.6 - 0.015 * z1)
    q2 = max(0.0, 0.8 - 0.02 * z1)
    q3 = max(0.0, 1.0 - 0.02 * z1)

    l1_mid = 2.0 * 0.24 * (1.0 - q2) ** 0.6667 / (z1**0.3333 * (1.0 - 0.143 * (1.0 - q2)))
    denom = (q2 - q1) if (q2 - q1) != 0 else 1.0
    l1 = np.select(
        [q < q1, q < q2, q < q3],
        [l0, l1_mid + (b - l1_mid) * (q - q1) / denom, np.full_like(q, b)],
        default=b * (1.0 - q) / (0.02 * z1),
    )
    length = np.maximum(l1, l0 * lambda_screening)

    zeta = q + (1.0 / (2.0 * fermi_velocity**2)) * (1.0 - q) * np.log(
        1 + (4 * length * fermi_velocity / 1.919) ** 2
    )
    # Z1^3 effect (ref. 779).
    zeta = zeta * (
        1.0
        + (1.0 / (z1 * z1))
        * (0.18 + 0.0015 * z2)
        * np.exp(-((7.6 - np.maximum(0.0, np.log(e_per_amu))) ** 2))
    )

    high = proton_stopping(e_per_amu, z2, pcoef) * (zeta * z1) ** 2

    # Below yrmin, fall back to velocity-proportional stopping anchored at eee.
    vrmin_low = max(vrmin, yrmin * z1_23)
    vmin = 0.5 * (
        vrmin_low + np.sqrt(max(0.0, vrmin_low**2 - 0.8 * fermi_velocity**2))
    )
    eee = 25 * vmin * vmin
    sp_at_eee = proton_stopping(np.array([eee]), z2, pcoef)[0]
    ratio = e_per_amu / eee
    low = sp_at_eee * (zeta * z1) ** 2 * np.sqrt(ratio)
    if z2 == 6:
        # Special correction for low-energy ions in carbon.
        low = low * ratio ** (0.75 * 0.5) / np.sqrt(ratio)

    return np.where(yr > max(yrmin, vrmin / z1_23), high, low)


def nuclear_stopping(
    energy_keV: np.ndarray, z1: int, m1: float, z2: int, m2: float
) -> np.ndarray:
    """Universal ZBL nuclear stopping (ziegler.c:246-256), eV/(1e15 at/cm^2)."""
    screen = z1**0.23 + z2**0.23
    epsilon = 32.53 * m2 * energy_keV / (z1 * z2 * (m1 + m2) * screen)

    a = 0.01321 * epsilon**0.21226 + 0.19593 * np.sqrt(epsilon)
    sn = np.where(
        epsilon < 30.0,
        0.5 * np.log(1 + 1.1383 * epsilon) / (epsilon + a),
        np.log(np.maximum(epsilon, 1.0 + 1e-30)) / (2 * epsilon),
    )
    return sn * z1 * z2 * m1 * 8.462 / ((m1 + m2) * screen)


class ZieglerStopping:
    """ZBL85 stopping for a given periodic table."""

    def __init__(self, elements: list[Element]):
        self._elements = elements

    def _params(self, z: int):
        element = self._elements[z - 1]
        if element.ziegler is None:
            raise ValueError(f"no Ziegler parameters for Z={z} (table covers 1..{MAX_Z})")
        return element.ziegler

    def __call__(
        self,
        z1: int,
        m1: float,
        z2: int,
        energy_keV,
        units: StoppingUnits = StoppingUnits.EV_1E15_ATOMS,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(electronic, nuclear)`` stopping.

        Out-of-range Z reproduces RUMP's behaviour: rather than raising, it
        returns the finite dummies ``(40, 5)`` so callers do not blow up
        (ziegler.c:199-203).
        """
        energy = _as_array(energy_keV)

        if not (1 <= z1 <= MAX_Z and 1 <= z2 <= MAX_Z):
            return np.full_like(energy, 40.0), np.full_like(energy, 5.0)

        projectile = self._params(z1)
        target = self._params(z2)

        if m1 == 0.0:
            m1 = projectile.mass_most_abundant or projectile.mass_average
        m2 = target.mass_average
        pcoef = np.asarray(target.proton_coefficients, dtype=np.float64)

        e_per_amu = energy / m1
        if np.any(e_per_amu > MAX_REDUCED_ENERGY):
            raise ValueError(
                f"reduced energy above {MAX_REDUCED_ENERGY} keV/amu is out of range"
            )

        if z1 == 1:
            se = proton_stopping(e_per_amu, z2, pcoef)
        elif z1 == 2:
            se = helium_stopping(e_per_amu, z1, z2, pcoef)
        else:
            se = heavy_ion_stopping(
                e_per_amu,
                z1,
                z2,
                pcoef,
                target.fermi_velocity,
                projectile.lambda_screening,
            )

        sn = nuclear_stopping(energy, z1, m1, z2, m2)
        return self._convert(se, sn, units, m2, target.atomic_density_e22)

    @staticmethod
    def _convert(se, sn, units, m2, atomic_density_e22):
        """Apply zstop's unit conversion (ziegler.c:259-278).

        .. warning::
           Case 2 reproduces a bug in the C: it computes ``sn`` from the
           *already converted* ``se`` rather than from ``sn``
           (``sn=se*.60222/m2;`` at ziegler.c:263). Preserved for fidelity --
           only ``EV_1E15_ATOMS`` (a no-op) is used by the simulation.
        """
        if units == StoppingUnits.EV_1E15_ATOMS:
            return se, sn
        if units == StoppingUnits.MEV_CM2_MG:
            se_converted = se * 0.60222 / m2
            return se_converted, se_converted * 0.60222 / m2
        if units == StoppingUnits.EV_PER_ANGSTROM:
            # zatrho is scaled to absolute atoms/cm^3 on load, hence the 1e22.
            scale = atomic_density_e22 * 1e22 * 1e-23
            return se * scale, sn * scale
        if units == StoppingUnits.LSS_REDUCED:
            raise NotImplementedError(
                "LSS reduced units need z1/m1, which this helper does not receive"
            )
        raise ValueError(f"unknown stopping units: {units}")
