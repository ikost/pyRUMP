"""PERT: parameter selection, windows, GO, and write-back into the sample.

The acceptance test is a round trip -- simulate a known Au marker, add Poisson
noise, then drive the shell exactly as a user would and check both that the
thickness comes back and that it lands in the SIM sample description, where
``SIM SAVE`` would pick it up.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from pyrump.model.detector import Measurement  # noqa: E402
from pyrump.model.geometry import Geometry, GeometryKind  # noqa: E402
from pyrump.model.spectrum import Calibration  # noqa: E402
from pyrump.shell.dispatch import CommandError  # noqa: E402
from pyrump.shell.repl import execute_line  # noqa: E402
from pyrump.shell.session import Buffer, Session  # noqa: E402
from pyrump.sim.engine import Beam, UniformSample, simulate  # noqa: E402


from conftest import data_dir

DATA = data_dir()
needs_data = pytest.mark.skipif(DATA is None, reason="legacy data tables unavailable")

#: Truth for the synthetic sample, 1e15 at/cm^2 of Au on Si.
TRUTH = 300.0
GUESS = 200.0

SAMPLE = """Sim Reset
Layer 1
 Thick {guess} /cm2
 Composition Au 1 /
Next
 Thick 5000 /cm2
 Composition Si 1 /
Maxpth 200
"""


def run(session: Session, *lines: str) -> None:
    stack = ["rump"]
    for line in lines:
        execute_line(session, line, stack)


@pytest.fixture(scope="module")
def synthetic():
    """A noisy Au-on-Si spectrum, and the parameters that produced it."""
    if DATA is None:
        pytest.skip("legacy data tables unavailable")
    session = Session.create(str(DATA))
    calibration = Calibration(kevch=5.0, kev0=0.0, npt=512)
    geometry = Geometry(theta=0.0, phi=10.0, kind=GeometryKind.CORNELL)
    measurement = Measurement(omega_msr=1.0, charge_uC=10.0, fwhm_keV=15.0)
    beam = Beam(e0_MeV=2.0, z=2, mass=4.0026)
    truth = UniformSample(
        thicknesses=[TRUTH, 5000.0],
        element_z=[79, 14],
        compositions=[[1.0, 0.0], [0.0, 1.0]],
    )
    clean = simulate(
        truth, beam, geometry, session.registry, session.table,
        calibration, measurement,
    )
    counts = np.random.default_rng(7).poisson(
        np.clip(clean.counts, 0, None)
    ).astype(float)
    return counts, calibration, geometry, measurement, beam


@pytest.fixture
def session(synthetic, tmp_path):
    counts, calibration, geometry, measurement, beam = synthetic
    from pyrump.model.spectrum import Spectrum

    built = Session.create(str(DATA))
    built.buffers.load(
        Buffer(
            spectrum=Spectrum(counts=counts.copy(), calibration=calibration),
            beam=beam, geometry=geometry, measurement=measurement, name="au",
        ),
        1,
    )
    built.buffers.active = 1
    sample = tmp_path / "au.lcm"
    sample.write_text(SAMPLE.format(guess=GUESS))
    run(built, f"sim get {sample}")
    return built


# -- selection and windows -------------------------------------------------


@needs_data
def test_window_and_normalization_are_recorded(session):
    run(session, "pert", "window 355 375", "norm 140 200")
    state = session.pert
    assert [(w.low, w.high) for w in state.windows.error] == [(355, 375)]
    assert (state.windows.normalisation.low, state.windows.normalisation.high) == (
        140, 200
    )


@needs_data
def test_thickness_selects_a_layer_by_one_based_number(session):
    run(session, "pert", "thick 1")
    assert [v.name for v in session.pert.varying] == ["thickness[0]"]


@needs_data
def test_selecting_a_layer_outside_the_sample_is_rejected(session):
    with pytest.raises(CommandError, match="outside 1-2"):
        run(session, "pert", "thick 9")


@needs_data
def test_composition_needs_an_element_in_the_sample(session):
    with pytest.raises(CommandError, match="not in the sample"):
        run(session, "pert", "composition 1 Xe")


@needs_data
def test_the_same_parameter_cannot_be_selected_twice(session):
    with pytest.raises(CommandError, match="already being varied"):
        run(session, "pert", "thick 1", "thick 1")


@needs_data
def test_go_without_a_selection_is_rejected(session):
    with pytest.raises(CommandError, match="nothing selected"):
        run(session, "pert", "window 355 375", "go")


@needs_data
def test_single_and_multi_modes(session):
    # One run() call, because each starts a fresh mode stack at the RUMP level
    # -- SINGLE and MULTI only exist inside PERT.
    run(session, "pert", "single")
    assert session.pert.multi is False
    run(session, "pert", "multi")
    assert session.pert.multi is True


@needs_data
def test_a_normalisation_window_with_a_free_correction_is_rejected(session):
    """Degenerate: both absorb the same scale factor (pert.c:1163)."""
    with pytest.raises(CommandError):
        run(
            session, "pert", "window 355 375", "norm 140 200",
            "correction", "go",
        )


# -- the fit itself --------------------------------------------------------


@needs_data
def test_go_recovers_the_thickness_and_writes_it_back(session, capsys):
    run(session, "pert", "window 355 375", "norm 140 200", "thick 1", "go")

    # Recovered to within a few percent of truth, from a 33% low guess.
    fitted = session.script.layers[0].thickness
    assert fitted == pytest.approx(TRUTH, rel=0.05)
    assert fitted != GUESS

    output = capsys.readouterr().out
    assert "reduced chi-square" in output
    assert "thickness[0]" in output


@needs_data
def test_offset_recovers_a_calibration_shift_and_writes_it_back(tmp_path, capsys):
    """A sample-charging-style energy shift: OFFSET alone should recover kev0.

    Calibration feeds the channel binning inside ``simulate()`` (not just axis
    labels), so a genuine kev0 mismatch between the buffer's guess and the
    data that produced it is exactly what a charging shift looks like.
    """
    from pyrump.model.spectrum import Spectrum

    true_kev0 = 8.0
    true_calibration = Calibration(kevch=5.0, kev0=true_kev0, npt=512)
    geometry = Geometry(theta=0.0, phi=10.0, kind=GeometryKind.CORNELL)
    measurement = Measurement(omega_msr=1.0, charge_uC=10.0, fwhm_keV=15.0)
    beam = Beam(e0_MeV=2.0, z=2, mass=4.0026)
    truth = UniformSample(
        thicknesses=[TRUTH, 5000.0],
        element_z=[79, 14],
        compositions=[[1.0, 0.0], [0.0, 1.0]],
    )

    session = Session.create(str(DATA))
    clean = simulate(
        truth, beam, geometry, session.registry, session.table,
        true_calibration, measurement,
    )
    counts = np.random.default_rng(11).poisson(
        np.clip(clean.counts, 0, None)
    ).astype(float)

    # The buffer starts out unshifted -- the wrong calibration for this data.
    guess_calibration = Calibration(kevch=5.0, kev0=0.0, npt=512)
    session.buffers.load(
        Buffer(
            spectrum=Spectrum(counts=counts, calibration=guess_calibration),
            beam=beam, geometry=geometry, measurement=measurement, name="au",
        ),
        1,
    )
    session.buffers.active = 1
    sample = tmp_path / "au.lcm"
    sample.write_text(SAMPLE.format(guess=TRUTH))
    run(session, f"sim get {sample}")

    run(session, "pert", "window 355 375", "norm 140 200", "offset", "go")

    fitted = session.buffers[1].calibration.kev0
    assert fitted == pytest.approx(true_kev0, abs=1.0)
    assert fitted != 0.0

    output = capsys.readouterr().out
    assert "kev(0)" in output


@needs_data
def test_the_fit_leaves_the_simulation_stale_so_compare_redraws(session):
    run(session, "pert", "window 355 375", "norm 140 200", "thick 1", "go")
    assert session.dirty is True


@needs_data
def test_the_fitted_sample_survives_sim_save(session, tmp_path):
    run(session, "pert", "window 355 375", "norm 140 200", "thick 1", "go")
    out = tmp_path / "fitted.lcm"
    run(session, f"sim save {out}")

    from pyrump.script.lcm import read_lcm

    assert read_lcm(out).layers[0].thickness == pytest.approx(TRUTH, rel=0.05)


@needs_data
def test_thickness_write_back_preserves_the_unit(session):
    run(session, "pert", "window 355 375", "norm 140 200", "thick 1", "go")
    # The script keeps magnitude + unit; only the magnitude should move.
    assert session.script.layers[0].unit == "/cm2"
