"""PERT: parameter selection, windows, GO, and write-back into the sample.

The acceptance test is a round trip -- simulate a known Au marker, add Poisson
noise, then drive the shell exactly as a user would and check both that the
thickness comes back and that it lands in the SIM sample description, where
``SIM SAVE`` would pick it up.
"""

from __future__ import annotations

import re
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
    """The displayed name echoes back the 1-based layer number the user
    typed, not the fit engine's internal 0-based index."""
    run(session, "pert", "thick 1")
    assert [v.name for v in session.pert.varying] == ["layer 1 thickness"]


@needs_data
def test_sim_delete_warns_about_a_pert_selection_on_a_later_layer(session, capsys):
    """Deleting layer 1 shifts layer 2 down to index 0, so the earlier
    ``thickness 2`` selection (layer index 1) now silently points at
    whatever ended up there -- SIM must say so."""
    run(session, "pert", "thickness 2", "return", "sim", "layer 1", "delete")
    output = capsys.readouterr().out
    assert "WARNING" in output
    assert "layer 2 thickness" in output
    assert [layer.thickness for layer in session.script.layers] == [5000.0]


@needs_data
def test_sim_delete_takes_a_layer_number_directly(session, capsys):
    """``DELETE n`` selects that layer and removes it in one step, without
    a preceding ``LAYER n`` -- same effect, same shift warning."""
    run(session, "pert", "thickness 2", "return", "sim", "delete 1")
    output = capsys.readouterr().out
    assert "WARNING" in output
    assert "layer 2 thickness" in output
    assert [layer.thickness for layer in session.script.layers] == [5000.0]


@needs_data
def test_sim_delete_of_a_later_layer_does_not_warn(session, capsys):
    """Deleting layer 2 never shifts layer 1, so a selection on layer 1
    (the earlier, unaffected layer) needs no warning."""
    run(session, "pert", "thickness 1", "return", "sim", "layer 2", "delete")
    output = capsys.readouterr().out
    assert "WARNING" not in output


@needs_data
def test_sim_open_warns_about_a_pert_selection_on_the_current_layer(session, capsys):
    """Inserting a blank layer above layer 1 shifts it (and the
    already-selected ``thickness 1``) down to index 1."""
    run(session, "pert", "thickness 1", "return", "sim", "layer 1", "open")
    output = capsys.readouterr().out
    assert "WARNING" in output
    assert "layer 1 thickness" in output


@needs_data
def test_selecting_a_layer_outside_the_sample_is_rejected(session):
    with pytest.raises(CommandError, match="outside 1-2"):
        run(session, "pert", "thick 9")


@needs_data
def test_sim_delete_of_an_undefined_layer_number_is_rejected(session):
    with pytest.raises(CommandError, match="layer 9 is not defined"):
        run(session, "sim", "delete 9")
    assert [layer.thickness for layer in session.script.layers] == [GUESS, 5000.0]


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


# -- numbering, WINDOW REMOVE, CLEAR <n> ------------------------------------


@needs_data
def test_window_remove_deletes_the_numbered_window(session):
    run(session, "pert", "window 100 200", "window 300 400", "window remove 1")
    assert [(w.low, w.high) for w in session.pert.windows.error] == [(300, 400)]


@needs_data
def test_window_remove_out_of_range_is_rejected(session):
    run(session, "pert", "window 100 200")
    with pytest.raises(CommandError, match="outside 1-1"):
        run(session, "pert", "window remove 5")


@needs_data
def test_window_remove_with_no_windows_is_rejected(session):
    with pytest.raises(CommandError, match="no error windows are set"):
        run(session, "pert", "window remove 1")


@needs_data
def test_window_remove_prints_the_resulting_numbered_list(session, capsys):
    run(session, "pert", "window 100 200", "window 300 400")
    capsys.readouterr()
    run(session, "pert", "window remove 1")
    output = capsys.readouterr().out
    assert "[1] 300-400" in output
    assert "100-200" not in output


@needs_data
def test_clear_n_removes_one_varying_parameter(session):
    run(session, "pert", "thick 1", "composition 1 Au", "clear 1")
    assert [v.name for v in session.pert.varying] == ["layer 1 composition Au"]


@needs_data
def test_clear_n_out_of_range_is_rejected(session):
    run(session, "pert", "thick 1")
    with pytest.raises(CommandError, match="outside 1-1"):
        run(session, "pert", "clear 5")


@needs_data
def test_clear_n_with_nothing_selected_is_rejected(session):
    with pytest.raises(CommandError, match="nothing selected to clear"):
        run(session, "pert", "clear 1")


@needs_data
def test_clear_with_no_args_still_wipes_everything(session):
    run(session, "pert", "window 100 200", "norm 140 200", "thick 1", "clear")
    state = session.pert
    assert state.varying == []
    assert state.windows.error == []
    assert state.windows.normalisation is None


@needs_data
def test_clear_n_prints_the_resulting_numbered_varying_list(session, capsys):
    run(session, "pert", "thick 1", "composition 1 Au")
    capsys.readouterr()
    run(session, "pert", "clear 1")
    output = capsys.readouterr().out
    assert "varying:" in output
    assert "[1] layer 1 composition Au" in output
    assert "layer 1 thickness" not in output


@needs_data
def test_parms_numbers_windows_and_varying_parameters(session, capsys):
    run(session, "pert", "window 100 200", "window 300 400", "thick 1", "parms")
    output = capsys.readouterr().out
    assert "[1] 100-200" in output
    assert "[2] 300-400" in output
    assert "[1] layer 1 thickness" in output


# -- clearer names and per-layer element validation -------------------------


@needs_data
def test_composition_display_name_uses_the_element_symbol(session):
    run(session, "pert", "composition 1 Au")
    assert session.pert.varying[0].name == "layer 1 composition Au"


@needs_data
def test_composition_rejects_an_element_not_declared_in_this_layer(session):
    """Si only appears in layer 2 of the fixture's sample."""
    with pytest.raises(CommandError, match="not part of layer 1"):
        run(session, "pert", "composition 1 Si")


@needs_data
def test_species_rejects_an_element_not_declared_as_a_species_in_this_layer(
    session, tmp_path
):
    sample = tmp_path / "species.lcm"
    sample.write_text(
        "Sim Reset\n"
        "Layer 1\n"
        " Thick 200 /cm2\n"
        " Composition Au 1 /\n"
        " Species Mn 1 /\n"
        "Next\n"
        " Thick 5000 /cm2\n"
        " Composition Si 1 /\n"
        "Maxpth 200\n"
    )
    run(session, f"sim get {sample}")
    with pytest.raises(CommandError, match="not part of layer 2"):
        run(session, "pert", "species 2 Mn")


@needs_data
def test_species_accepts_an_element_declared_in_this_layer(session, tmp_path):
    sample = tmp_path / "species.lcm"
    sample.write_text(
        "Sim Reset\n"
        "Layer 1\n"
        " Thick 200 /cm2\n"
        " Composition Au 1 /\n"
        " Species Mn 1 /\n"
        "Next\n"
        " Thick 5000 /cm2\n"
        " Composition Si 1 /\n"
        "Maxpth 200\n"
    )
    run(session, f"sim get {sample}")
    run(session, "pert", "species 1 Mn")
    assert session.pert.varying[0].name == "layer 1 species Mn"


@needs_data
def test_show_prints_the_sample_description(session, capsys):
    run(session, "pert", "show")
    output = capsys.readouterr().out
    assert "Au" in output
    assert "Si" in output


# -- optional search bounds --------------------------------------------------


@needs_data
def test_thickness_accepts_an_optional_bound(session):
    run(session, "pert", "thickness 1 100 500")
    entry = session.pert.varying[0]
    assert entry.bounds == (100.0, 500.0)
    assert entry.parameter.lower == 100.0
    assert entry.parameter.upper == 500.0


@needs_data
def test_thickness_with_no_bound_keeps_the_default_range(session):
    run(session, "pert", "thickness 1")
    entry = session.pert.varying[0]
    assert entry.bounds is None
    assert entry.parameter.lower == 0.0  # thickness's own physical floor


@needs_data
def test_a_reversed_bound_is_rejected(session):
    with pytest.raises(CommandError, match="empty bound"):
        run(session, "pert", "thickness 1 500 100")


@needs_data
def test_a_bound_with_only_one_number_is_rejected(session):
    with pytest.raises(CommandError, match="expected a maximum bound"):
        run(session, "pert", "thickness 1 100")


@needs_data
def test_a_simple_parameter_accepts_a_bound(session):
    run(session, "pert", "mev 1.5 2.5")
    entry = session.pert.varying[0]
    assert entry.bounds == (1.5, 2.5)
    assert entry.parameter.lower == 1.5
    assert entry.parameter.upper == 2.5


@needs_data
def test_composition_accepts_a_bound(session):
    run(session, "pert", "composition 1 Au 0.5 1.0")
    entry = session.pert.varying[0]
    assert entry.bounds == (0.5, 1.0)


@needs_data
def test_parms_shows_a_bound_only_for_a_parameter_that_has_one(session, capsys):
    run(session, "pert", "thickness 1 100 500", "mev", "parms")
    output = capsys.readouterr().out
    assert "layer 1 thickness  bounds 100-500" in output
    assert "mev" in output
    assert "mev  bounds" not in output


@needs_data
def test_bound_overrides_the_parameter_own_default_floor(session):
    """Thickness normally floors at 0 by default; an explicit negative
    lower bound should win outright, matching the original -- pyRUMP's
    built-in floor is a convenience default, not a hard physical rule the
    user cannot override."""
    run(session, "pert", "thickness 1 -50 500")
    assert session.pert.varying[0].parameter.lower == -50.0


@needs_data
def test_bound_is_enforced_by_the_fit(session, capsys):
    """Bounding thickness 1 to [190, 210] -- below the true 300 -- should
    stop the solver at the boundary rather than reaching the true value,
    proving the bound reaches the solver and isn't just cosmetic."""
    run(session, "pert", "window 355 375", "norm 140 200", "thickness 1 190 210", "go")
    fitted = session.script.layers[0].thickness
    assert fitted <= 210.5
    assert fitted != pytest.approx(TRUTH, rel=0.05)


@needs_data
def test_save_and_get_round_trip_a_bound(session, tmp_path):
    pert_file = tmp_path / "bounded.pert"
    run(
        session, "pert", "thickness 1 100 500",
        f"save {pert_file}", "clear", f"get {pert_file}",
    )
    entry = session.pert.varying[0]
    assert entry.bounds == (100.0, 500.0)
    assert entry.parameter.lower == 100.0
    assert entry.parameter.upper == 500.0


# -- GET/SAVE ----------------------------------------------------------------


@needs_data
def test_save_and_get_round_trip(session, tmp_path):
    pert_file = tmp_path / "usual.pert"
    run(
        session, "pert",
        "window 100 200", "window 300 400", "norm 140 200",
        "thick 1", "composition 1 Au", "mev",
        f"save {pert_file}",
        "clear",
        f"get {pert_file}",
    )
    assert pert_file.exists()
    state = session.pert
    assert [(w.low, w.high) for w in state.windows.error] == [(100, 200), (300, 400)]
    assert (state.windows.normalisation.low, state.windows.normalisation.high) == (
        140, 200,
    )
    assert [v.name for v in state.varying] == [
        "layer 1 thickness", "layer 1 composition Au", "mev",
    ]


@needs_data
def test_save_and_get_preserve_single_mode(session, tmp_path):
    pert_file = tmp_path / "single.pert"
    run(
        session, "pert", "single", "thick 1",
        f"save {pert_file}", "clear", f"get {pert_file}",
    )
    assert session.pert.multi is False


@needs_data
def test_get_replaces_the_current_selection_rather_than_merging(session, tmp_path):
    pert_file = tmp_path / "usual.pert"
    run(
        session, "pert", "thick 1", "mev",
        f"save {pert_file}",
        "clear",
        "fwhm",
        f"get {pert_file}",
    )
    assert [v.name for v in session.pert.varying] == ["layer 1 thickness", "mev"]


@needs_data
def test_get_missing_file_is_rejected(session, tmp_path):
    missing = tmp_path / "nope.pert"
    with pytest.raises(CommandError, match="no such file"):
        run(session, "pert", f"get {missing}")


@needs_data
def test_save_appends_the_pert_extension(session, tmp_path):
    base = tmp_path / "usual"
    run(session, "pert", "thick 1", f"save {base}")
    assert (tmp_path / "usual.pert").exists()


@needs_data
def test_save_and_get_round_trip_a_species_selection(session, tmp_path):
    sample = tmp_path / "species.lcm"
    sample.write_text(
        "Sim Reset\n"
        "Layer 1\n"
        " Thick 200 /cm2\n"
        " Composition Au 1 /\n"
        " Species Mn 1 /\n"
        "Next\n"
        " Thick 5000 /cm2\n"
        " Composition Si 1 /\n"
        "Maxpth 200\n"
    )
    pert_file = tmp_path / "species.pert"
    run(
        session, f"sim get {sample}",
        "pert", "species 1 Mn",
        f"save {pert_file}",
        "clear",
        f"get {pert_file}",
    )
    assert [v.name for v in session.pert.varying] == ["layer 1 species Mn"]


# -- PERT as a one-shot from RUMP level --------------------------------------


@needs_data
def test_pert_command_works_as_a_one_shot_from_rump_level(session, capsys):
    run(session, "pert thick 1")
    assert [v.name for v in session.pert.varying] == ["layer 1 thickness"]

    run(session, "pert parms")
    output = capsys.readouterr().out
    assert "layer 1 thickness" in output


@needs_data
def test_pert_one_shot_does_not_stay_in_pert_mode(session):
    run(session, "pert thick 1")
    with pytest.raises(CommandError, match="unrecognized command: go"):
        run(session, "go")


@needs_data
def test_pert_one_shot_rejects_an_unrecognized_command(session):
    with pytest.raises(CommandError, match="unrecognized PERT command"):
        run(session, "pert bogus")


@needs_data
def test_pert_get_go_runs_as_a_single_line(session, tmp_path, capsys):
    pert_file = tmp_path / "usual.pert"
    run(
        session, "pert", "window 355 375", "norm 140 200", "thick 1",
        f"save {pert_file}",
    )

    run(session, f"pert get {pert_file} go")

    fitted = session.script.layers[0].thickness
    assert fitted == pytest.approx(TRUTH, rel=0.05)
    output = capsys.readouterr().out
    assert "reduced chi-square" in output
    assert "fitted stack" in output


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
    assert "layer 1 thickness" in output
    assert "fitted stack" in output
    assert "Au" in output
    assert "Si" in output


@needs_data
def test_volume_on_prints_per_evaluation_progress_during_go(session, capsys):
    """A fit with several varying parameters can take seconds with nothing
    printed in between otherwise, easy to mistake for a hung prompt."""
    run(session, "pert", "window 355 375", "norm 140 200", "thick 1", "volume", "go")
    output = capsys.readouterr().out
    assert "eval   1" in output
    assert "chi2/dof" in output


@needs_data
def test_volume_off_by_default_prints_no_per_evaluation_progress(session, capsys):
    run(session, "pert", "window 355 375", "norm 140 200", "thick 1", "go")
    output = capsys.readouterr().out
    assert "eval " not in output


@needs_data
def test_go_prints_the_fitted_thickness_rounded_to_a_whole_unit(session, capsys):
    """Sub-angstrom precision is meaningless for RBS, so both the fitted
    value and the "(was ...)" comparison should be whole numbers -- unlike
    other fit parameters (e.g. an uncertainty), which keep full precision."""
    run(session, "pert", "window 355 375", "norm 140 200", "thick 1", "go")
    output = capsys.readouterr().out
    # The GO result line is the last "layer 1 thickness" mention -- earlier
    # ones are just the "varying ..." selection echo, with nothing after.
    value_field = re.findall(r"layer 1 thickness[ \t]+(\S+)", output)[-1]
    assert "." not in value_field
    assert "(was 200)" in output


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
