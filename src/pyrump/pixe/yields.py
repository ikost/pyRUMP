r"""Line yields from a layered thin-film sample.

The SIM sample is cut into sublayers (slabs) the way the RBS simulation cuts
it, only finer (:data:`PIXE_MAXPTH`), and the beam's energy in each slab
comes from the same stopping powers and the same inward march. The yield of
line :math:`\ell` of element Z is then summed over the film's slabs k:

.. math::
    Y_\ell = H\,N_\text{ion}\,\frac{\Omega}{4\pi}\,\varepsilon(E_\ell)
             \sum_k \frac{n_{Z,k}}{\cos\alpha}\,\sigma_\ell(\bar E_k)\,A_{k,\ell}

.. math::
    A_{k,\ell} = \exp\!\Big(-\sum_{j<k}\frac{\mu_j\,\rho t_j}{\cos\theta_d}\Big)
                 \,\frac{1-e^{-x_k}}{x_k},\qquad
    x_k = \frac{\mu_k\,\rho t_k}{\cos\theta_d}

with :math:`n_{Z,k}` the element's areal density in the slab,
:math:`\bar E_k` the beam's mean energy across it, :math:`\mu` the slab's
mass attenuation coefficient at the line energy (mixture rule),
:math:`\rho t` its mass thickness, :math:`\alpha` the beam's and
:math:`\theta_d` the detector's angle to the sample normal, and
:math:`\varepsilon` the detector efficiency including any filters. The
second factor of :math:`A` is the exact self-absorption of a uniform slab.
For a truly thin film :math:`A \to 1` and :math:`\bar E_k \to E_0`.

What is left out, on purpose:

* **the substrate** -- the last SIM layer. Its peaks come from the measured
  bare-substrate spectrum, not from this sum;
* **RBS absorber layers** -- foils in front of the RBS detector, which the
  beam never crosses;
* **roughness** (``FUZZ``): the slabs are flat.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from ..model.geometry import Geometry, GeometryKind
from ..sim.absorber import first_sample_slab
from ..sim.engine import build_sample_grid
from ..sim.precal import march_inbound
from ..stopping.bragg import bragg_coefficients
from ..stopping.table import StoppingTable
from .atomic import AtomicData, XrayLine, atomic_data
from .detector import Absorber, PixeDetector
from .production import line_production
from .xsect import ion_name, ionisation_table

AVOGADRO = 6.02214076e23
#: The elementary charge in µC.
ELEMENTARY_CHARGE_UC = 1.602176634e-13
BARN_CM2 = 1e-24

#: Largest sublayer, in 1e15 atoms/cm^2 -- finer than RBS's default (1000),
#: since sigma changes steeply with energy and the X-rays are absorbed on
#: their way out. About 15 nm of Pt, 20 nm of Si.
PIXE_MAXPTH = 100.0

#: Line energies the attenuation table covers for every element, keV
#: (xraylib's data start between 0.1007 and 0.1037 keV).
LINE_ENERGY_RANGE = (0.11, 99.0)


@dataclass(frozen=True, slots=True)
class Exposure:
    """What sets the absolute number of X-rays: the ions that arrived and
    the instrumental constants."""

    charge_uC: float
    charge_state: int = 1
    correction: float = 1.0
    """RUMP's CORR, applied as for RBS: the charge is divided by it."""

    live_fraction: float = 1.0
    """The PIXE detector's live time over real time."""

    h: tuple[float, float, float] = (1.0, 1.0, 1.0)
    """The instrumental constant H for K, L and M lines."""

    @property
    def ions(self) -> float:
        return (
            self.charge_uC
            / (ELEMENTARY_CHARGE_UC * self.charge_state * self.correction)
            * self.live_fraction
        )

    def h_for(self, family: str) -> float:
        return self.h["KLM".index(family)]


@dataclass(slots=True)
class LineYield:
    """One line's simulated yield: counts in the detector."""

    z: int
    symbol: str
    line: XrayLine
    counts: float
    sigma_barn: float
    """Production cross section at the beam's full energy."""

    efficiency: float
    by_layer: dict[int, float]
    """Counts per SIM layer (0 = layer 1, the surface)."""

    @property
    def family(self) -> str:
        return self.line.shell[0]

    @property
    def energy_keV(self) -> float:
        return self.line.energy_keV


def element_density_g_cm3(element) -> float:
    return element.atomic_density * element.mass / AVOGADRO


def absorber_mu(absorber: Absorber, energy_keV, periodic_table, atomic: AtomicData):
    """``(mu in cm^2/g, density in g/cm^3)`` of an absorber's material: an
    element, or a compound by the mixture rule (mass fractions from pyRUMP's
    atomic masses)."""
    compound = absorber.compound
    if compound is None:
        element = periodic_table.by_symbol(absorber.material)
        return atomic.mu(element.z, energy_keV), element_density_g_cm3(element)
    elements = {symbol: periodic_table.by_symbol(symbol) for symbol in compound.atoms}
    grams = {s: n * elements[s].mass for s, n in compound.atoms.items()}
    total = sum(grams.values())
    mu = sum(grams[s] / total * atomic.mu(e.z, energy_keV) for s, e in elements.items())
    return mu, compound.density_g_cm3


def transmission(absorber: Absorber, energy_keV, periodic_table, atomic: AtomicData) -> np.ndarray:
    """Fraction of X-rays getting through ``absorber`` at normal incidence."""
    mu, density = absorber_mu(absorber, energy_keV, periodic_table, atomic)
    solid = np.exp(-mu * density * absorber.thickness_um * 1e-4)
    hole = absorber.hole_percent / 100.0
    return hole + (1.0 - hole) * solid


def efficiency(detector: PixeDetector, energy_keV, periodic_table, atomic: AtomicData) -> np.ndarray:
    """Intrinsic efficiency (window transmission times absorption in the
    crystal), times the filters' transmission."""
    energy = np.asarray(energy_keV, dtype=np.float64)
    mu, density = absorber_mu(detector.crystal, energy, periodic_table, atomic)
    crystal_mass = density * detector.crystal.thickness_um * 1e-4
    result = transmission(detector.window, energy, periodic_table, atomic)
    result = result * -np.expm1(-mu * crystal_mass)
    for absorber in detector.filters:
        result = result * transmission(absorber, energy, periodic_table, atomic)
    return result


def pixe_grid(sample, alpha_deg: float, exit_deg: float, periodic_table):
    """The sample's slabs for PIXE: RBS's sublayering, only finer."""
    fine = replace(sample, maxpth=min(sample.maxpth, PIXE_MAXPTH))
    geometry = Geometry(theta=alpha_deg, phi=10.0, psi=exit_deg, kind=GeometryKind.GENERAL)
    return fine, geometry, build_sample_grid(fine, geometry, periodic_table)


def simulate_lines(
    sample,
    beam,
    theta_deg: float,
    detector: PixeDetector,
    exposure: Exposure,
    registry,
    periodic_table,
    *,
    atomic: AtomicData | None = None,
    faithful: bool = True,
) -> list[LineYield]:
    """Yield of every line of every film element, strongest first.

    ``sample`` is the RBS simulation's :class:`~pyrump.sim.engine.UniformSample`
    and ``beam`` its :class:`~pyrump.sim.engine.Beam`; ``theta_deg`` is the
    sample's tilt, the RBS geometry's signed THETA. The beam comes in at
    ``|theta_deg|`` to the normal, the X-rays leave at
    :meth:`PixeDetector.exit_angle`.
    """
    atomic = atomic or atomic_data()
    ionisation = ionisation_table(ion_name(beam.z, beam.mass))
    if not ionisation.e_min <= beam.e0_MeV <= ionisation.e_max:
        raise ValueError(
            f"beam energy {beam.e0_MeV:g} MeV outside the cross-section tables "
            f"({ionisation.e_min:g}-{ionisation.e_max:g} MeV)"
        )
    alpha_deg, exit_deg = abs(theta_deg), detector.exit_angle(theta_deg)
    fine, geometry, grid = pixe_grid(sample, alpha_deg, exit_deg, periodic_table)
    cos_in = np.cos(np.radians(alpha_deg))
    cos_out = np.cos(np.radians(exit_deg))
    if cos_in <= 0 or cos_out <= 0:
        raise ValueError("the beam and the detector must face the sample's front")

    # The beam's energy through the slabs, exactly as the RBS march has it.
    table = StoppingTable.build(registry, beam.z, beam.mass, beam.e0_MeV, list(grid.element_z))
    coefficients = bragg_coefficients(table, grid.composition, grid.element_z)
    first = first_sample_slab(grid.layer_index, fine.absorber_layers)
    inbound = march_inbound(
        table, coefficients, grid.composition, grid.element_z,
        e0_keV=beam.e0_MeV * 1000.0, sec_in=geometry.sec_in,
        cutoff_keV=table.cutoff * 1000.0, z_beam=beam.z, first_slab=first,
        faithful=faithful,
    )

    # Film slabs: past the absorber, before the substrate, reached by the beam.
    substrate = int(grid.layer_index.max())
    film = np.arange(first, min(inbound.reached, grid.n_slab))
    film = film[grid.layer_index[film] != substrate]
    if film.size == 0:
        return []
    mean_energy_MeV = 0.5e-3 * (inbound.energy[film] + inbound.energy[film + 1])
    if mean_energy_MeV.min() < ionisation.e_min:
        raise ValueError(
            f"the beam slows below {ionisation.e_min:g} MeV inside the film, "
            "under the cross-section tables"
        )

    # Mass thickness and mass fractions of each film slab.
    masses = np.array([periodic_table.by_z(z).mass for z in grid.element_z])
    atoms = grid.composition[film] * 1e15  # atoms/cm^2 per element
    grams = atoms * masses / AVOGADRO
    mass_thickness = grams.sum(axis=1)  # g/cm^2
    fractions = grams / mass_thickness[:, None]

    def attenuation(energy_keV: float) -> np.ndarray:
        """A_k for one X-ray energy: absorption above slab k, times the
        slab's own self-absorption."""
        mu = sum(
            fractions[:, i] * atomic.mu(z, energy_keV)
            for i, z in enumerate(grid.element_z)
            if fractions[:, i].any()
        )
        x = mu * mass_thickness / cos_out
        above = np.concatenate(([0.0], np.cumsum(x)[:-1]))
        self_absorption = np.where(x > 1e-12, -np.expm1(-x) / np.maximum(x, 1e-300), 1.0)
        return np.exp(-above) * self_absorption

    solid = detector.solid_angle_msr * 1e-3 / (4 * np.pi)
    e_min, e_max = LINE_ENERGY_RANGE
    # Cross sections at each slab's mean energy, and at the full beam energy
    # last (reported in LINES, not used in the sum).
    energies = np.append(mean_energy_MeV, beam.e0_MeV)
    yields: list[LineYield] = []
    for column, z in enumerate(grid.element_z):
        per_slab_atoms = atoms[:, column] / cos_in
        if not per_slab_atoms.any():
            continue
        element = periodic_table.by_z(z)
        for produced in line_production(z, energies, ionisation, atomic):
            energy = produced.line.energy_keV
            if not e_min <= energy <= e_max:
                continue
            sigma = produced.sigma_barn[:-1] * BARN_CM2
            eff = float(efficiency(detector, energy, periodic_table, atomic))
            per_slab = per_slab_atoms * sigma * attenuation(energy)
            scale = exposure.h_for(produced.family) * exposure.ions * solid * eff
            by_layer: dict[int, float] = {}
            for layer, value in zip(grid.layer_index[film], per_slab * scale):
                by_layer[int(layer)] = by_layer.get(int(layer), 0.0) + float(value)
            yields.append(LineYield(
                z=z, symbol=element.symbol, line=produced.line,
                counts=float(per_slab.sum() * scale),
                sigma_barn=float(produced.sigma_barn[-1]), efficiency=eff,
                by_layer=by_layer,
            ))
    return sorted(yields, key=lambda y: -y.counts)
