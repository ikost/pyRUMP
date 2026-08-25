"""Command-line interface.

``pyrump [simulate | fit | convert | plot | shell]``

A thin wrapper over the library. Running ``pyrump`` with no subcommand enters
the interactive RUMP shell (:mod:`pyrump.shell`); the subcommands here exist for
one-off jobs and for scripting.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from .. import __version__
from ._common import data_dir as _data_dir
from ._common import load_tables as _load
from ._common import read_spectrum as _read_spectrum
from ._common import resolve_beam


def _build(args, table, registry, densities):
    """Assemble the simulation inputs from CLI arguments."""
    from pyrump.model.detector import Measurement
    from pyrump.model.geometry import Geometry, GeometryKind
    from pyrump.model.spectrum import Calibration
    from pyrump.script.lcm import read_lcm, to_sample
    from pyrump.sim.engine import Beam

    script = read_lcm(args.sample)
    sample = to_sample(script, table, densities)

    z, mass = resolve_beam(table, args.beam)
    beam = Beam(e0_MeV=args.energy, z=z, mass=mass)
    geometry = Geometry(
        theta=args.theta, phi=args.phi, psi=args.psi,
        kind=GeometryKind[args.geometry.upper()],
    )
    calibration = Calibration(
        kevch=args.kevch, kev0=args.kev0, first=0.0, npt=args.channels
    )
    measurement = Measurement(
        omega_msr=args.omega, charge_uC=args.charge, fwhm_keV=args.fwhm,
        current_nA=args.current, tau_us=args.tau,
    )
    return sample, beam, geometry, calibration, measurement


def _simulation_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("sample", type=Path, help="sample description (.lcm)")
    parser.add_argument("--beam", default="4He", help="beam species (default: 4He)")
    parser.add_argument("--energy", type=float, default=2.0, help="MeV")
    parser.add_argument("--theta", type=float, default=0.0, help="incidence angle, deg")
    parser.add_argument(
        "--phi", type=float, default=10.0,
        help="RUMP's phi: 180 minus the scattering angle (default 10 = 170 deg)",
    )
    parser.add_argument("--psi", type=float, default=0.0, help="exit angle, deg")
    parser.add_argument(
        "--geometry", default="cornell", choices=["cornell", "ibm", "general"]
    )
    parser.add_argument("--kevch", type=float, default=5.0, help="keV per channel")
    parser.add_argument("--kev0", type=float, default=0.0, help="keV at channel 0")
    parser.add_argument("--channels", type=int, default=1024)
    parser.add_argument("--fwhm", type=float, default=15.0, help="detector FWHM, keV")
    parser.add_argument("--omega", type=float, default=1.0, help="solid angle, msr")
    parser.add_argument("--charge", type=float, default=10.0, help="uC")
    parser.add_argument("--current", type=float, default=0.0, help="nA (pile-up)")
    parser.add_argument("--tau", type=float, default=0.0, help="shaping time, us")


def command_simulate(args) -> int:
    from pyrump.io.ascii import write_ascii
    from pyrump.io.rbs import RbsSpectrum, write_rbs
    from pyrump.sim.engine import simulate

    table, registry, densities = _load(_data_dir(args.data))
    sample, beam, geometry, calibration, measurement = _build(
        args, table, registry, densities
    )
    spectrum = simulate(
        sample, beam, geometry, registry, table, calibration, measurement
    )

    if args.output is None:
        for energy, counts in zip(spectrum.energies, spectrum.counts):
            print(f"{energy:.3f} {counts:.6f}")
    elif Path(args.output).suffix.lower() == ".rbs":
        write_rbs(
            args.output,
            RbsSpectrum(
                counts=spectrum.counts, calibration=calibration, geometry=geometry,
                measurement=measurement, e0_MeV=beam.e0_MeV, zbeam=beam.z,
                mbeam=beam.mass, identifier=f"pyRUMP simulation of {args.sample.name}",
            ),
        )
        print(f"wrote {args.output}", file=sys.stderr)
    else:
        write_ascii(args.output, spectrum.counts, two_column=args.two_column)
        print(f"wrote {args.output}", file=sys.stderr)

    print(
        f"total {spectrum.total():.1f} counts in "
        f"{np.count_nonzero(spectrum.counts)} channels",
        file=sys.stderr,
    )
    return 0


def command_convert(args) -> int:
    from pyrump.io.ascii import write_ascii
    from pyrump.io.rbs import RbsSpectrum, write_rbs

    source = _read_spectrum(args.input)
    counts = source.counts
    target = Path(args.output)

    if target.suffix.lower() == ".rbs":
        if not isinstance(source, RbsSpectrum):
            raise SystemExit(
                "writing .rbs needs the beam and geometry metadata that ASCII "
                "files do not carry; convert from a .rbs source instead"
            )
        write_rbs(target, source)
    else:
        write_ascii(
            target, counts,
            identifier=getattr(source, "identifier", ""),
            two_column=args.two_column,
        )
    print(f"{args.input} -> {target} ({counts.size} channels)", file=sys.stderr)
    return 0


def command_plot(args) -> int:
    import matplotlib

    if args.output:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from pyrump.model.spectrum import Calibration, Spectrum
    from pyrump.plot.spectra import plot_comparison, plot_spectrum

    first = _read_spectrum(args.input)
    calibration = getattr(first, "calibration", Calibration(npt=first.counts.size))
    data = Spectrum(counts=np.asarray(first.counts), calibration=calibration)

    if args.compare:
        other = _read_spectrum(args.compare)
        model = Spectrum(
            counts=np.asarray(other.counts),
            calibration=getattr(other, "calibration", calibration),
        )
        figure = plot_comparison(data, model)
    else:
        figure = plot_spectrum(data, label=args.input.name).figure

    if args.output:
        figure.savefig(args.output, dpi=150, bbox_inches="tight")
        print(f"wrote {args.output}", file=sys.stderr)
    else:
        plt.show()
    return 0


def command_fit(args) -> int:
    from pyrump.fit.lm import fit
    from pyrump.fit.parameters import FitInputs, parameter, thickness
    from pyrump.fit.windows import Window, WindowSet
    from pyrump.sim.engine import simulate

    table, registry, densities = _load(_data_dir(args.data))
    sample, beam, geometry, calibration, measurement = _build(
        args, table, registry, densities
    )
    observed = _read_spectrum(args.data_file).counts

    parameters = []
    for name in args.vary:
        if name.lower().startswith("thickness"):
            layer = int(name.split(":")[1]) if ":" in name else 0
            parameters.append(thickness(layer))
        else:
            parameters.append(parameter(name))
    if not parameters:
        raise SystemExit("nothing to fit: pass --vary at least once")

    windows = WindowSet(
        error=[Window(low, high) for low, high in (args.window or [])]
    )
    inputs = FitInputs(
        sample=sample, beam=beam, geometry=geometry,
        calibration=calibration, measurement=measurement,
    )

    def run(current: FitInputs) -> np.ndarray:
        return simulate(
            current.sample, current.beam, current.geometry, registry, table,
            current.calibration, current.measurement,
        ).counts

    result = fit(run, observed, inputs, parameters, windows=windows)

    print(f"reduced chi-square {result.reduced_chi_square:.4f} on {result.dof} dof")
    print(f"{result.n_evaluations} evaluations, {result.message}")
    if result.n_invalid:
        print(f"warning: {result.n_invalid} windowed channels had zero predicted counts")
    for name, value in result.parameters.items():
        sigma = result.uncertainties.get(name)
        print(f"  {name:24s} {value:14.5g}" + (f"  +/- {sigma:.5g}" if sigma else ""))
    return 0


def _shell_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "macro", nargs="?", type=Path,
        help="command file to XEQ at startup, after the rc file",
    )
    parser.add_argument(
        "--norc", action="store_true", help="skip ~/.pyrumprc"
    )
    parser.add_argument(
        "--batch", action="store_true",
        help="run the macro and exit instead of dropping to the prompt",
    )


def _shell_defaults(args) -> None:
    """Fill in shell options for the bare ``pyrump`` invocation."""
    for name, value in (("macro", None), ("norc", False), ("batch", False)):
        if not hasattr(args, name):
            setattr(args, name, value)


def command_shell(args) -> int:
    from pyrump.shell.repl import run_shell

    return run_shell(
        data=args.data, macro=args.macro, norc=args.norc, batch=args.batch
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="pyrump",
        description="Rutherford backscattering simulation and analysis.",
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"pyrump {__version__}"
    )
    parser.add_argument("--data", help="directory holding atom4.dat and friends")
    sub = parser.add_subparsers(dest="command")

    shell_parser = sub.add_parser("shell", help="interactive RUMP shell (the default)")
    _shell_arguments(shell_parser)
    shell_parser.set_defaults(func=command_shell)

    simulate_parser = sub.add_parser("simulate", help="simulate a spectrum")
    _simulation_arguments(simulate_parser)
    simulate_parser.add_argument("-o", "--output", type=Path)
    simulate_parser.add_argument("--two-column", action="store_true")
    simulate_parser.set_defaults(func=command_simulate)

    fit_parser = sub.add_parser("fit", help="fit a sample to measured data")
    _simulation_arguments(fit_parser)
    fit_parser.add_argument("data_file", type=Path, help="measured spectrum")
    fit_parser.add_argument(
        "--vary", action="append", default=[],
        help="parameter to fit, e.g. thickness:0 or fwhm (repeatable)",
    )
    fit_parser.add_argument(
        "--window", action="append", nargs=2, type=int, metavar=("LOW", "HIGH"),
        help="error window in channels (repeatable, up to 10)",
    )
    fit_parser.set_defaults(func=command_fit)

    convert_parser = sub.add_parser("convert", help="convert between spectrum formats")
    convert_parser.add_argument("input", type=Path)
    convert_parser.add_argument("output", type=Path)
    convert_parser.add_argument("--two-column", action="store_true")
    convert_parser.set_defaults(func=command_convert)

    plot_parser = sub.add_parser("plot", help="plot a spectrum")
    plot_parser.add_argument("input", type=Path)
    plot_parser.add_argument("--compare", type=Path, help="overlay a second spectrum")
    plot_parser.add_argument("-o", "--output", type=Path, help="save instead of showing")
    plot_parser.set_defaults(func=command_plot)

    args = parser.parse_args(argv)
    if args.command is None:
        # Bare "pyrump" is the interactive shell -- the way the original was used.
        _shell_defaults(args)
        return command_shell(args)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
