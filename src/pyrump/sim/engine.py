"""Top-level simulation entry point.

Corresponds to ``SimCreateDetails`` (creatr.c:267), assembling the stages:

``StoppingTable`` -> ``SlabGrid`` -> Bragg coefficients -> inbound march ->
per-isotope depth loop -> bricks.

The brick -> channel fill, detector convolution and pileup stages are later
milestones; :func:`simulate_bricks` stops at the point the C calls
``SimFillSpectrum``.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

from ..atomic.density import layer_atomic_density
from ..atomic.tables import PeriodicTable
from ..model.detector import Measurement, yield_normalisation
from ..model.geometry import Geometry
from ..model.spectrum import Calibration, Spectrum
from ..physics.kinematics import kinematic_factor
from ..physics.xsec.rutherford import setup_scatter
from ..stopping.bragg import bragg_coefficients
from ..stopping.registry import StoppingRegistry
from ..stopping.table import StoppingTable
from .absorber import first_sample_slab
from .bricks import Bricks
from .convolve import convolve_detector
from .fill.straggled import fill_straggled
from .fill.trapezoid import fill_trapezoid
from .fuzz import replica_thicknesses
from .multiscatter import add_multiple_scattering
from .pileup import new_pileup
from .ideal import simulate_isotope, straggle_geometry_factor
from .precal import march_inbound
from .slabs import DEFAULT_MAXPATH, build_uniform_grid


@dataclass(slots=True)
class Beam:
    """Incident beam parameters."""

    e0_MeV: float = 2.0
    z: int = 2
    mass: float = 4.0026
    charge_state: int = 1


@dataclass(slots=True)
class UniformSample:
    """A stack of layers, each of uniform or depth-varying composition."""

    thicknesses: list[float]
    """Areal thickness per layer, 1e15 atoms/cm^2."""

    element_z: list[int]
    compositions: list[list[float]]
    """``(n_layer, n_element)``; rows are normalised internally."""

    sublayers: list[int] | None = None
    maxpth: float = DEFAULT_MAXPATH
    straggle: float = 0.0
    """Multiplier on Bohr straggling. RUMP's default of 0 disables it."""

    profiles: list | None = None
    """Optional per-layer :class:`~pyrump.profiles.equations.Profile`."""

    species: list[list[float]] | None = None
    """``(n_layer, n_element)`` composition each profile blends toward."""

    densities: list[float] | None = None
    """Per-layer matrix density in 1e23 at/cm^3. Computed if omitted."""

    absorber_layers: int = 0
    """Leading layers that are absorber (dead layer / window), not sample."""

    multiple: float = 0.0
    """Ad-hoc multiple-scattering strength. No physical basis; 0 disables."""

    fuzz_amounts: list[float] | None = None
    """Per-layer thickness roughness, in the layer's own units."""

    fuzz_steps: list[int] | None = None
    """Per-layer replica count. Every fuzzed layer multiplies the total."""

    tags: dict = field(default_factory=dict)


def simulate_bricks(
    sample: UniformSample,
    beam: Beam,
    geometry: Geometry,
    registry: StoppingRegistry,
    periodic_table: PeriodicTable,
    *,
    screening: bool = True,
) -> Bricks:
    """Run the forward model up to the brick stage.

    Bricks are emitted one block per target isotope, heaviest first, matching
    the order the C produces so captures can be compared element-wise.
    """
    geometry.validate()

    # Layer densities drive the depth scale of every depth-dependent profile, so
    # they are derived the same way the C does (inverse-density averaging over
    # the layer's own composition) rather than being left to the caller.
    atomic_densities = np.array(
        [periodic_table.by_z(z).atomic_density for z in sample.element_z],
        dtype=np.float64,
    )
    compositions = np.atleast_2d(np.asarray(sample.compositions, dtype=np.float64))
    if sample.densities is not None:
        matrix_densities = np.asarray(sample.densities, dtype=np.float64)
    else:
        matrix_densities = np.array(
            [layer_atomic_density(row, atomic_densities) for row in compositions],
            dtype=np.float64,
        )
    if sample.species is not None:
        species_densities = np.array(
            [
                layer_atomic_density(row, atomic_densities)
                for row in np.atleast_2d(np.asarray(sample.species, dtype=np.float64))
            ],
            dtype=np.float64,
        )
    else:
        species_densities = None

    table = StoppingTable.build(
        registry, beam.z, beam.mass, beam.e0_MeV, list(sample.element_z)
    )
    grid = build_uniform_grid(
        np.asarray(sample.thicknesses, dtype=np.float64),
        np.asarray(sample.compositions, dtype=np.float64),
        list(sample.element_z),
        maxpth=sample.maxpth,
        sec_in=geometry.sec_in,
        sec_out=geometry.sec_out,
        explicit_sublayers=sample.sublayers,
        profiles=sample.profiles,
        layer_species=(
            np.asarray(sample.species, dtype=np.float64)
            if sample.species is not None
            else None
        ),
        layer_densities=matrix_densities,
        species_densities=species_densities,
    )
    coefficients = bragg_coefficients(table, grid.composition, grid.element_z)
    cutoff_keV = table.cutoff * 1000.0
    first_slab = first_sample_slab(grid.layer_index, sample.absorber_layers)

    inbound = march_inbound(
        table,
        coefficients,
        grid.composition,
        grid.element_z,
        e0_keV=beam.e0_MeV * 1000.0,
        sec_in=geometry.sec_in,
        cutoff_keV=cutoff_keV,
        straggle_scale=sample.straggle,
        z_beam=beam.z,
        first_slab=first_slab,
    )

    blocks: list[np.ndarray] = []
    for column, z_target in enumerate(sample.element_z):
        element = periodic_table.by_z(z_target)
        # Monoisotopic elements still get one pass, at the average mass.
        isotopes = sorted(element.isotopes, key=lambda i: -i.mass) or [
            type("Iso", (), {"mass": element.mass, "fraction": 1.0})()
        ]
        for isotope in isotopes:
            cross_section = setup_scatter(
                beam.z,
                beam.mass,
                z_target,
                isotope.mass,
                geometry.scattering_angle,
                screening=screening,
            )
            geometry_factor = None
            if sample.straggle:
                kinematic = kinematic_factor(
                    beam.mass, isotope.mass, geometry.scattering_angle
                )
                geometry_factor = straggle_geometry_factor(
                    kinematic, geometry.sec_in, geometry.sec_out
                )

            block = simulate_isotope(
                table,
                coefficients,
                coefficients,  # backscatter: in and out share the projectile table
                grid.composition[:, column],
                inbound,
                cross_section,
                m_beam=beam.mass,
                m_target=isotope.mass,
                scattering_angle_deg=geometry.scattering_angle,
                isotope_fraction=isotope.fraction,
                sec_in=geometry.sec_in,
                sec_out=geometry.sec_out,
                cutoff_keV=cutoff_keV,
                straggle_geometry=geometry_factor,
                first_slab=first_slab,
            )
            if len(block):
                blocks.append(block.data)

    if not blocks:
        return Bricks.empty(0)
    return Bricks(data=np.concatenate(blocks))


def simulate(
    sample: UniformSample,
    beam: Beam,
    geometry: Geometry,
    registry: StoppingRegistry,
    periodic_table: PeriodicTable,
    calibration: Calibration,
    measurement: Measurement | None = None,
    *,
    screening: bool = True,
    convolve_edge: str = "rump",
) -> Spectrum:
    """Full forward model: sample in, channel spectrum out.

    ``convolve_edge`` selects RUMP's non-count-conserving edge handling
    (default) or a renormalising variant; see :mod:`pyrump.sim.convolve`.
    """
    measurement = measurement or Measurement()

    # Stage order follows SimCreateDetails (creatr.c:307-345) exactly, and it
    # matters: resolution is convolved before normalisation, pile-up needs real
    # counts so it comes after, and the multiple-scattering tail is last.
    counts = np.zeros(calibration.npt, dtype=np.float64)

    for thicknesses, amplitude in _fuzz_replicas(sample):
        replica = sample if thicknesses is None else replace(
            sample, thicknesses=list(thicknesses)
        )
        bricks = simulate_bricks(
            replica, beam, geometry, registry, periodic_table, screening=screening
        )
        if len(bricks) == 0:
            continue
        # SimAnlyz (anlyz.c:207) picks the fill routine per brick on whether
        # either straggling width is non-zero; with straggling the trapezoid is
        # abandoned for two Gaussian-broadened triangles.
        fill = fill_straggled if bricks.has_straggling else fill_trapezoid
        counts += amplitude * fill(bricks, calibration)

    counts = convolve_detector(
        counts, calibration, measurement.fwhm_keV, mode=convolve_edge
    )
    counts *= yield_normalisation(measurement)

    if measurement.tau_us > 0 and measurement.current_nA > 0:
        counts = new_pileup(
            counts,
            tau_us=measurement.tau_us,
            current_nA=measurement.current_nA,
            charge_uC=measurement.charge_uC,
        )[: calibration.npt]

    if sample.multiple:
        counts = add_multiple_scattering(
            counts,
            strength=sample.multiple,
            charge_uC=measurement.charge_uC,
            omega_msr=measurement.omega_msr,
        )

    return Spectrum(counts=counts, calibration=calibration)


def _fuzz_replicas(sample: UniformSample):
    """Yield ``(thicknesses, amplitude)`` for each fuzz replica.

    Without fuzzing this is a single unperturbed pass, so the common path costs
    nothing.
    """
    if not sample.fuzz_steps or not any(s > 1 for s in sample.fuzz_steps):
        yield None, 1.0
        return
    amounts = sample.fuzz_amounts or [0.0] * len(sample.thicknesses)
    yield from replica_thicknesses(
        np.asarray(sample.thicknesses, dtype=np.float64), amounts, sample.fuzz_steps
    )
