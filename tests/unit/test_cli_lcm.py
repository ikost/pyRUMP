"""M13 acceptance: the ``.lcm`` subset, plotting helpers and the CLI.

The `.lcm` format is a subset of RUMP's interactive command language, so the
sharpest test available is a **byte-identical round-trip** of the file RUMP
itself wrote (``data/Fixed/ITO.lcm``).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from pyrump.atomic.density import DensityTable  # noqa: E402
from pyrump.atomic.tables import PeriodicTable  # noqa: E402
from pyrump.cli.__main__ import main  # noqa: E402
from pyrump.model.spectrum import Calibration, Spectrum  # noqa: E402
from pyrump.plot.spectra import (  # noqa: E402
    plot_comparison,
    plot_depth_profile,
    plot_spectrum,
)
from pyrump.profiles.equations import EquationType  # noqa: E402
from pyrump.script.lcm import parse_lcm, read_lcm, to_sample, write_lcm  # noqa: E402


def _data_dir() -> Path | None:
    env = os.environ.get("PYRUMP_C_REFERENCE")
    roots = [Path(env)] if env else []
    roots.append(Path(__file__).resolve().parents[2] / "C-code")
    for root in roots:
        if (root / "rump" / "data" / "atom4.dat").is_file():
            return root / "rump" / "data"
    return None


DATA = _data_dir()
ITO = DATA / "Fixed" / "ITO.lcm" if DATA else None
needs_data = pytest.mark.skipif(DATA is None, reason="legacy data tables unavailable")

SIMPLE = """\
Sim Reset
Layer 1
 Thick 2000 A
 Composition Si 1 /
Maxpth 200
Foil disable
"""


# ------------------------------------------------------------------- .lcm


def test_parses_a_simple_sample():
    script = parse_lcm(SIMPLE)
    assert len(script.layers) == 1
    assert script.layers[0].thickness == 2000.0
    assert script.layers[0].unit == "A"
    assert script.layers[0].composition == {"Si": 1.0}
    assert script.maxpth == 200.0
    assert script.elements == ["Si"]


def test_next_starts_a_new_layer():
    script = parse_lcm(
        "Sim Reset\nLayer 1\n Thick 100 A\n Composition Si 1 /\n"
        "Next\n Thick 200 A\n Composition Au 1 /\n"
    )
    assert len(script.layers) == 2
    assert script.layers[1].composition == {"Au": 1.0}
    assert script.elements == ["Si", "Au"]


def test_composition_list_terminated_by_slash():
    script = parse_lcm(
        "Sim Reset\nLayer 1\n Thick 1 A\n Composition In 2 O 3 Sn 0.1 /\n"
    )
    assert script.layers[0].composition == {"In": 2.0, "O": 3.0, "Sn": 0.1}


def test_element_without_amount_is_rejected():
    with pytest.raises(ValueError, match="no amount"):
        parse_lcm("Sim Reset\nLayer 1\n Thick 1 A\n Composition Si /\n")


def test_equation_and_species():
    script = parse_lcm(
        "Sim Reset\nLayer 1\n Thick 1000 A\n Composition Si 1 /\n"
        " Equation Linear 0 0.2\n Species Au 1 /\n"
    )
    layer = script.layers[0]
    assert layer.profile is not None
    assert layer.profile.type is EquationType.LINEAR
    assert layer.profile.parameters == (0.0, 0.2)
    assert layer.species == {"Au": 1.0}


def test_equation_aliases():
    """RUMP accepts Error for ERFC, Implant for Gaussian, Thicfilm for ThickFilm."""
    for name, expected in (
        ("Error", EquationType.ERFC),
        ("Implant", EquationType.GAUSSIAN),
        ("Thicfilm", EquationType.THICKFILM),
    ):
        script = parse_lcm(
            f"Sim Reset\nLayer 1\n Thick 1 A\n Composition Si 1 /\n Equation {name} 1 2 3\n"
        )
        assert script.layers[0].profile.type is expected


def test_unknown_equation_rejected():
    with pytest.raises(ValueError, match="unknown equation"):
        parse_lcm("Sim Reset\nLayer 1\n Thick 1 A\n Composition Si 1 /\n Equation Wobble 1\n")


def test_unrecognised_commands_are_collected_not_fatal():
    """Real files carry plotting commands we do not implement."""
    script = parse_lcm(SIMPLE + "Plot\nOverlay theory\n")
    assert script.ignored == ["Plot", "Overlay theory"]
    assert len(script.layers) == 1


def test_fuzzy_and_absorber():
    script = parse_lcm(
        "Sim Reset\nLayer 1\n Thick 100 A\n Composition Si 1 /\n Fuzzy 25 5\n"
        "Absorber 1\nStraggle 1.5\n"
    )
    assert script.layers[0].fuzz_amount == 25.0
    assert script.layers[0].fuzz_steps == 5
    assert script.absorber_layers == 1
    assert script.straggle == 1.5


@needs_data
def test_ito_round_trips_byte_identically():
    """The strongest check: reproduce RUMP's own output exactly."""
    original = ITO.read_text()
    assert write_lcm(parse_lcm(original)) == original


@needs_data
def test_ito_converts_to_a_sample():
    table = PeriodicTable.load(DATA / "atom4.dat", DATA / "pscoef.dat")
    densities = DensityTable.load(DATA / "density.tab")
    sample = to_sample(read_lcm(ITO), table, densities)

    assert sample.element_z == [49, 8, 50, 6, 1]  # In O Sn C H
    assert len(sample.thicknesses) == 3
    # 151 A of ITO; ITO is absent from density.tab so RUMP falls back to
    # silicon's 0.4997, giving 151 * 0.4997.
    assert sample.thicknesses[0] == pytest.approx(151 * 0.4997, rel=1e-3)
    # 10 um of a light polymer is a very large areal density.
    assert sample.thicknesses[1] > 10_000


@needs_data
@pytest.mark.parametrize(
    "unit, magnitude, expected_ratio",
    [("A", 1000.0, 1.0), ("nm", 100.0, 1.0), ("um", 0.1, 1.0)],
)
def test_length_units_agree(unit, magnitude, expected_ratio):
    """1000 A, 100 nm and 0.1 um are the same thickness."""
    table = PeriodicTable.load(DATA / "atom4.dat", DATA / "pscoef.dat")
    densities = DensityTable.load(DATA / "density.tab")
    script = parse_lcm(
        f"Sim Reset\nLayer 1\n Thick {magnitude} {unit}\n Composition Si 1 /\n"
    )
    sample = to_sample(script, table, densities)
    reference = to_sample(parse_lcm(SIMPLE.replace("2000 A", "1000 A")), table, densities)
    assert sample.thicknesses[0] == pytest.approx(
        reference.thicknesses[0] * expected_ratio, rel=1e-9
    )


@needs_data
def test_atomic_unit_is_already_areal():
    table = PeriodicTable.load(DATA / "atom4.dat", DATA / "pscoef.dat")
    densities = DensityTable.load(DATA / "density.tab")
    script = parse_lcm("Sim Reset\nLayer 1\n Thick 500 /CM2\n Composition Si 1 /\n")
    assert to_sample(script, table, densities).thicknesses[0] == pytest.approx(500.0)


# --------------------------------------------------------------- plotting


def _spectrum(n=200, peak=100):
    counts = np.zeros(n)
    counts[50:150] = np.linspace(10, peak, 100)
    return Spectrum(counts=counts, calibration=Calibration(kevch=5.0, npt=n))


def test_plot_spectrum_returns_axes():
    ax = plot_spectrum(_spectrum(), label="test")
    assert ax.get_xlabel() == "Energy (keV)"
    assert len(ax.lines) == 1


def test_plot_spectrum_channel_axis():
    ax = plot_spectrum(_spectrum(), energy_axis=False)
    assert ax.get_xlabel() == "Channel"


def test_plot_comparison_has_residual_panel():
    data = _spectrum()
    model = _spectrum(peak=110)
    figure = plot_comparison(data, model)
    assert len(figure.axes) == 2
    assert r"residual ($\sigma$)" in figure.axes[1].get_ylabel()


def test_plot_comparison_without_residuals():
    figure = plot_comparison(_spectrum(), _spectrum(peak=110), residuals=False)
    assert len(figure.axes) == 1


def test_plot_comparison_shades_windows():
    window = np.zeros(200, dtype=bool)
    window[60:100] = True
    figure = plot_comparison(_spectrum(), _spectrum(peak=110), window=window)
    assert figure.axes[0].patches  # the axvspan


@needs_data
def test_plot_depth_profile():
    from pyrump.sim.slabs import build_uniform_grid

    grid = build_uniform_grid([1000.0], [[1.0, 2.0]], [14, 8])
    ax = plot_depth_profile(grid, ["Si", "O"])
    assert len(ax.lines) == 2
    assert "Depth" in ax.get_xlabel()


# -------------------------------------------------------------------- CLI


@needs_data
def test_cli_simulate_writes_ascii(tmp_path, capsys):
    sample = tmp_path / "si.lcm"
    sample.write_text(SIMPLE)
    out = tmp_path / "out.dat"
    assert main(["--data", str(DATA), "simulate", str(sample), "-o", str(out)]) == 0
    assert out.exists()
    values = [float(line) for line in out.read_text().split() if line]
    assert sum(values) > 0


@needs_data
def test_cli_simulate_writes_rbs_that_reads_back(tmp_path):
    from pyrump.io.rbs import read_rbs

    sample = tmp_path / "si.lcm"
    sample.write_text(SIMPLE)
    out = tmp_path / "out.rbs"
    assert main(["--data", str(DATA), "simulate", str(sample), "-o", str(out)]) == 0
    spectrum = read_rbs(out)
    assert spectrum.counts.sum() > 0
    assert spectrum.e0_MeV == pytest.approx(2.0)
    assert spectrum.zbeam == 2


@needs_data
def test_cli_convert_rbs_to_ascii(tmp_path):
    out = tmp_path / "converted.dat"
    assert main(["convert", str(DATA / "Fixed" / "2A.rbs"), str(out)]) == 0
    assert len(out.read_text().splitlines()) >= 2048


@needs_data
def test_cli_plot_saves_a_file(tmp_path):
    out = tmp_path / "plot.png"
    assert main(["plot", str(DATA / "Fixed" / "2A.rbs"), "-o", str(out)]) == 0
    assert out.stat().st_size > 1000


@needs_data
def test_cli_fit_recovers_a_thickness(tmp_path, capsys):
    """End to end: simulate a thicker layer, then fit a thinner start to it."""
    truth = tmp_path / "truth.lcm"
    truth.write_text(SIMPLE.replace("2000 A", "2400 A"))
    start = tmp_path / "start.lcm"
    start.write_text(SIMPLE)
    data = tmp_path / "data.rbs"

    main(["--data", str(DATA), "simulate", str(truth), "-o", str(data)])
    capsys.readouterr()

    assert main([
        "--data", str(DATA), "fit", str(start), str(data),
        "--vary", "thickness:0", "--window", "190", "226",
    ]) == 0
    output = capsys.readouterr().out
    assert "thickness[0]" in output
    fitted = float(output.split("thickness[0]")[1].split()[0])
    # 2400 A of silicon is 2400 * 0.4977 = 1194 in 1e15 at/cm^2.
    assert fitted == pytest.approx(1194.0, rel=0.02)


def test_cli_fit_requires_something_to_vary(tmp_path):
    sample = tmp_path / "s.lcm"
    sample.write_text(SIMPLE)
    data = tmp_path / "d.dat"
    data.write_text("1\n2\n3\n")
    with pytest.raises(SystemExit, match="nothing to fit"):
        main(["fit", str(sample), str(data)])


def test_cli_reports_missing_data_directory(tmp_path, monkeypatch):
    monkeypatch.delenv("PYRUMP_DATA", raising=False)
    monkeypatch.delenv("PYRUMP_C_REFERENCE", raising=False)
    monkeypatch.chdir(tmp_path)
    sample = tmp_path / "s.lcm"
    sample.write_text(SIMPLE)
    with pytest.raises(SystemExit, match="Could not find the data tables"):
        main(["simulate", str(sample)])


def test_cli_requires_a_subcommand():
    with pytest.raises(SystemExit):
        main([])
