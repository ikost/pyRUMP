#!/usr/bin/env python3
"""Runnable demonstration of pyRUMP — ``python main.py``.

Exercises the whole package end to end and prints what it found, so you can
check the install works and see the API in use without reading the test suite.

    python main.py              # run everything
    python main.py --list       # show the individual demos
    python main.py simulate     # run just one

Needs the RUMP data tables. They are found automatically at ``C-code/rump/data``
or via ``PYRUMP_DATA`` / ``PYRUMP_C_REFERENCE``.
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings
from pathlib import Path

import numpy as np


# --------------------------------------------------------------------- setup


def find_data() -> Path:
    """Locate the data tables, or explain how to point at them."""
    candidates = []
    for key in ("PYRUMP_DATA",):
        if os.environ.get(key):
            candidates.append(Path(os.environ[key]))
    if os.environ.get("PYRUMP_C_REFERENCE"):
        candidates.append(Path(os.environ["PYRUMP_C_REFERENCE"]) / "rump" / "data")
    candidates.append(Path(__file__).parent / "C-code" / "rump" / "data")

    for path in candidates:
        if (path / "atom4.dat").is_file():
            return path
    raise SystemExit(
        "Could not find the data tables (atom4.dat, pscoef.dat, newstop.kal,\n"
        "density.tab). Set PYRUMP_DATA to the directory holding them."
    )


def load(data: Path):
    """Build the periodic table, stopping registry and compound densities.

    Building the registry refits stopping polynomials, so do it once and pass
    it around rather than rebuilding inside a loop.
    """
    from pyrump.atomic.density import DensityTable
    from pyrump.atomic.tables import PeriodicTable
    from pyrump.io.kalbitzer import parse_kalbitzer
    from pyrump.stopping.kalbitzer import KalbitzerStopping
    from pyrump.stopping.registry import StoppingRegistry
    from pyrump.stopping.ziegler import ZieglerStopping

    table = PeriodicTable.load(data / "atom4.dat", data / "pscoef.dat")
    registry = StoppingRegistry(
        table.elements,
        kalbitzer=KalbitzerStopping(parse_kalbitzer(data / "newstop.kal"), table.elements),
        ziegler=ZieglerStopping(table.elements),
    )
    return table, registry, DensityTable.load(data / "density.tab")


def heading(text: str) -> None:
    print(f"\n\033[1m{text}\033[0m\n" + "-" * len(text))


# --------------------------------------------------------------------- demos


def demo_tables(data, table, registry, densities) -> None:
    """Atomic data and stopping powers."""
    heading("Atomic data and stopping powers")

    silicon = table.by_z(14)
    print(f"  silicon: mass {silicon.mass}, density {silicon.atomic_density:.4e} at/cm^3")
    print(f"           isotopes " + ", ".join(
        f"{i.mass_number}Si {i.fraction * 100:.2f}%" for i in
        sorted(silicon.isotopes, key=lambda i: i.mass)
    ))

    print("\n  stopping power of 2 MeV ions in silicon:")
    for symbol, z, mass in (("H", 1, 1.00797), ("He", 2, 4.0026), ("Li", 3, 6.939)):
        result = registry(z, mass, 14, [2000.0])
        print(f"    {symbol:2s} {result.values[0]:7.2f} eV/(1e15 at/cm^2)"
              f"   [{result.source.name.lower()}]")
    print("\n  note He and H use Konac/Kalbitzer, not Ziegler -- it takes priority")


def demo_kinematics(data, table, registry, densities) -> None:
    """Where each element's surface edge lands."""
    from pyrump.physics.kinematics import kinematic_factor

    heading("Kinematics: surface edges for 2 MeV He at 170 degrees")
    print(f"  {'element':8s} {'K':>8s} {'edge (keV)':>12s}")
    for symbol in ("C", "O", "Si", "Ti", "Cu", "Ag", "Au"):
        element = table.by_symbol(symbol)
        mass = max(element.isotopes, key=lambda i: i.fraction).mass
        k = kinematic_factor(4.0026, mass, 170.0)
        print(f"  {symbol:8s} {k:8.4f} {k * 2000:12.1f}")
    print("\n  heavier target -> less energy lost -> higher edge")


def demo_simulate(data, table, registry, densities) -> None:
    """Simulate a spectrum from a layered sample."""
    from pyrump.model.detector import Measurement
    from pyrump.model.geometry import Geometry
    from pyrump.model.spectrum import Calibration
    from pyrump.sim.engine import Beam, UniformSample, simulate

    heading("Simulating a gold marker buried in silicon")

    sample = UniformSample(
        thicknesses=[500.0, 100.0, 2000.0],      # 1e15 atoms/cm^2
        element_z=[14, 79],                      # silicon, gold
        compositions=[[1.0, 0.0], [0.0, 1.0], [1.0, 0.0]],
    )
    calibration = Calibration(kevch=5.0, kev0=0.0, npt=1024)
    spectrum = simulate(
        sample,
        Beam(e0_MeV=2.0, z=2, mass=4.0026),
        Geometry(theta=0.0, phi=10.0),           # phi = 180 - scattering angle
        registry, table, calibration,
        Measurement(omega_msr=1.0, charge_uC=10.0, fwhm_keV=15.0),
    )

    counts = spectrum.counts
    occupied = np.flatnonzero(counts)
    print(f"  {spectrum.total():.0f} counts across channels "
          f"{occupied.min()}..{occupied.max()}")

    peak = int(np.argmax(counts))
    print(f"  highest channel {peak} at {calibration.edge_energy(peak):.0f} keV "
          f"({counts[peak]:.0f} counts) -- the gold marker")
    _sparkline(counts[occupied.min(): occupied.max() + 1])


def demo_identify(data, table, registry, densities) -> None:
    """Identify elements in a measured spectrum from its edges."""
    from pyrump.io.rbs import read_rbs
    from pyrump.physics.kinematics import kinematic_factor

    path = data / "Fixed" / "2A.rbs"
    if not path.is_file():
        print("\n  (skipping: no example .rbs files present)")
        return

    heading("Identifying elements in a measured spectrum")
    measured = read_rbs(path)
    print(f"  {measured.identifier}")
    print(f"  {measured.e0_MeV} MeV Z={measured.zbeam}, scattering angle "
          f"{measured.geometry.scattering_angle}, {measured.calibration.kevch} keV/ch")

    counts = measured.counts
    smoothed = np.convolve(counts, np.ones(5) / 5, mode="same")
    slope = np.diff(smoothed)
    edges = [
        i for i in range(30, 700)
        if slope[i] < -0.02 * smoothed.max() and smoothed[i] > 0.02 * smoothed.max()
    ]
    grouped: list[list[int]] = []
    for index in edges:
        if grouped and index - grouped[-1][-1] <= 6:
            grouped[-1].append(index)
        else:
            grouped.append([index])

    predictions = []
    for z in range(3, 84):
        element = table.by_z(z)
        if not element.isotopes:
            continue
        mass = max(element.isotopes, key=lambda i: i.fraction).mass
        channel = measured.calibration.channel_of(
            kinematic_factor(measured.mbeam, mass, measured.geometry.scattering_angle)
            * measured.e0_MeV * 1000
        )
        predictions.append((element.symbol, float(channel)))

    print("\n  falling edges found, and the nearest predicted surface edges:")
    for group in grouped:
        channel = int(np.mean(group))
        nearest = sorted(predictions, key=lambda p: abs(p[1] - channel))[:3]
        print(f"    channel {channel:4d}  ->  " +
              ", ".join(f"{sym} ({ch:.0f})" for sym, ch in nearest))
    print("\n  C, O and In/Sn: indium tin oxide on a polymer")
    print("  (In and Sn are 0.6 channels apart -- not separable)")


def demo_fit(data, table, registry, densities) -> None:
    """Recover a known thickness from synthetic data."""
    from pyrump.fit.lm import fit
    from pyrump.fit.parameters import FitInputs, parameter, thickness
    from pyrump.fit.windows import Window, WindowSet
    from pyrump.model.detector import Measurement
    from pyrump.model.geometry import Geometry
    from pyrump.model.spectrum import Calibration
    from pyrump.sim.engine import Beam, UniformSample, simulate

    heading("Fitting: recovering a thickness from noisy data")

    calibration = Calibration(kevch=5.0, kev0=0.0, npt=1024)
    geometry = Geometry(theta=0.0, phi=10.0)
    measurement = Measurement(omega_msr=1.0, charge_uC=10.0, fwhm_keV=15.0)

    def run(inputs: FitInputs) -> np.ndarray:
        return simulate(
            inputs.sample, inputs.beam, inputs.geometry, registry, table,
            inputs.calibration, inputs.measurement,
        ).counts

    truth = 1200.0
    reference = FitInputs(
        sample=UniformSample([truth], [14], [[1.0]]), beam=Beam(),
        geometry=geometry, calibration=calibration, measurement=measurement,
    )
    noisy = np.random.default_rng(3).poisson(np.maximum(run(reference), 0)).astype(float)

    guess = 900.0
    inputs = FitInputs(
        sample=UniformSample([guess], [14], [[1.0]]), beam=Beam(),
        geometry=geometry, calibration=calibration, measurement=measurement,
    )
    print(f"  truth {truth:.0f}, starting guess {guess:.0f} (1e15 atoms/cm^2)")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        result = fit(
            run, noisy, inputs, [thickness(0), parameter("fwhm")],
            windows=WindowSet(error=[Window(200, 260)]),
        )

    for name, value in result.parameters.items():
        sigma = result.uncertainties.get(name)
        print(f"    {name:16s} {value:10.2f}" + (f"  +/- {sigma:.2f}" if sigma else ""))
    print(f"  reduced chi-square {result.reduced_chi_square:.3f} on {result.dof} dof, "
          f"{result.n_evaluations} evaluations")
    if result.correlation is not None:
        print(f"  thickness/fwhm correlation {result.correlation[0, 1]:+.2f}"
              "  -- check this, not just the error bars")


def demo_lcm(data, table, registry, densities) -> None:
    """Read a RUMP sample description."""
    from pyrump.script.lcm import read_lcm, to_sample, write_lcm

    path = data / "Fixed" / "ITO.lcm"
    if not path.is_file():
        print("\n  (skipping: no example .lcm files present)")
        return

    heading("Reading a RUMP sample description")
    script = read_lcm(path)
    print(f"  {len(script.layers)} layers, elements {', '.join(script.elements)}")
    for index, layer in enumerate(script.layers):
        composition = " ".join(f"{k}{v:g}" for k, v in layer.composition.items())
        print(f"    layer {index}: {layer.thickness:g} {layer.unit:5s} {composition}")

    sample = to_sample(script, table, densities)
    print("\n  as areal density (1e15 atoms/cm^2): " +
          ", ".join(f"{t:.1f}" for t in sample.thicknesses))
    identical = write_lcm(script) == path.read_text()
    print(f"  round-trips byte-identically: {identical}")


DEMOS = {
    "tables": demo_tables,
    "kinematics": demo_kinematics,
    "simulate": demo_simulate,
    "identify": demo_identify,
    "fit": demo_fit,
    "lcm": demo_lcm,
}


def _sparkline(values: np.ndarray, width: int = 64) -> None:
    """A rough shape of the spectrum, for a terminal."""
    blocks = " ▁▂▃▄▅▆▇█"
    if values.size == 0 or values.max() <= 0:
        return
    binned = np.array_split(values, min(width, values.size))
    heights = np.array([b.mean() for b in binned])
    scaled = (heights / heights.max() * (len(blocks) - 1)).round().astype(int)
    print("  " + "".join(blocks[i] for i in scaled))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("demo", nargs="*", choices=list(DEMOS) + [], default=None,
                        help="which demo(s) to run; default is all")
    parser.add_argument("--list", action="store_true", help="list the demos")
    parser.add_argument("--data", help="directory holding atom4.dat and friends")
    args = parser.parse_args(argv)

    if args.list:
        for name, function in DEMOS.items():
            print(f"  {name:12s} {function.__doc__.splitlines()[0]}")
        return 0

    data = Path(args.data) if args.data else find_data()
    print(f"pyRUMP demo — data tables from {data}")
    table, registry, densities = load(data)

    for name in (args.demo or list(DEMOS)):
        DEMOS[name](data, table, registry, densities)

    print("\nAll demos completed.")
    print("Run the test suite with:  pytest -q")
    return 0


if __name__ == "__main__":
    sys.exit(main())
