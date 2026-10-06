"""SNAPSHOT: saving the session as it stands, and picking it up again.

The acceptance test is a round trip: drive a session as a user would --
load, fit or tune by hand, SNAPSHOT -- then run the restore macro in a fresh
session started somewhere else, and check that it is the same session
(:func:`pyrump.shell.snapshot.fingerprint`, which covers the data, every
buffer parameter, the settings, the sample, and the PERT and PIXE setup).
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from pyrump.cli.__main__ import _implicit_shell  # noqa: E402
from pyrump.io.rbs import read_rbs, write_rbs  # noqa: E402
from pyrump.model.detector import Measurement  # noqa: E402
from pyrump.model.geometry import Geometry, GeometryKind  # noqa: E402
from pyrump.model.spectrum import Calibration, Spectrum  # noqa: E402
from pyrump.shell import snapshot  # noqa: E402
from pyrump.shell.commands import rump  # noqa: E402
from pyrump.shell.dispatch import CommandError  # noqa: E402
from pyrump.shell.repl import execute_file, execute_line  # noqa: E402
from pyrump.shell.session import Buffer, Session  # noqa: E402
from pyrump.sim.engine import Beam, UniformSample, simulate  # noqa: E402

from conftest import data_dir

DATA = data_dir()
EXAMPLES = Path(__file__).resolve().parents[2] / "examples"
needs_data = pytest.mark.skipif(DATA is None, reason="legacy data tables unavailable")

SAMPLE = """Sim Reset
Layer 1
 Thick 200 /cm2
 Composition Au 1 /
Next
 Thick 5000 /cm2
 Composition Si 1 /
Maxpth 200
"""


def run(session: Session, *lines: str, stack: list[str] | None = None) -> list[str]:
    stack = stack if stack is not None else ["rump"]
    for line in lines:
        execute_line(session, line, stack)
    return stack


def restored(macro: Path) -> Session:
    """A fresh session that has run ``macro``, as ``pyrump <macro>`` would."""
    fresh = Session.create(str(DATA))
    execute_file(fresh, macro, ["rump"])
    return fresh


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    import matplotlib.pyplot as plt

    plt.close("all")


@pytest.fixture(scope="module")
def au_spectrum(tmp_path_factory):
    """300e15 Au/cm^2 on Si, with Poisson noise, as a binary .rbs file."""
    if DATA is None:
        pytest.skip("legacy data tables unavailable")
    tables = Session.create(str(DATA))
    calibration = Calibration(kevch=5.0, kev0=0.0, npt=512)
    geometry = Geometry(theta=0.0, phi=10.0, kind=GeometryKind.CORNELL)
    measurement = Measurement(omega_msr=1.0, charge_uC=10.0, fwhm_keV=15.0)
    beam = Beam(e0_MeV=2.0, z=2, mass=4.0026)
    truth = UniformSample(
        thicknesses=[300.0, 5000.0], element_z=[79, 14],
        compositions=[[1.0, 0.0], [0.0, 1.0]],
    )
    clean = simulate(
        truth, beam, geometry, tables.registry, tables.table, calibration, measurement
    )
    counts = np.random.default_rng(7).poisson(np.clip(clean.counts, 0, None)).astype(float)
    buffer = Buffer(
        spectrum=Spectrum(counts=counts, calibration=calibration),
        beam=beam, geometry=geometry, measurement=measurement, identifier="au",
    )
    path = tmp_path_factory.mktemp("data") / "au.rbs"
    write_rbs(path, buffer.to_rbs())
    return path


@pytest.fixture
def work(tmp_path, monkeypatch):
    """The folder the user works in -- where snapshots land."""
    folder = tmp_path / "work"
    folder.mkdir()
    monkeypatch.chdir(folder)
    return folder


@pytest.fixture
def session(au_spectrum, work) -> Session:
    """The Au spectrum read with GET, and a starting sample read with SIM GET."""
    built = Session.create(str(DATA))
    (work / "start.lcm").write_text(SAMPLE)
    run(built, f"get {au_spectrum}", "sim get start.lcm")
    return built


def elsewhere(tmp_path, monkeypatch) -> Path:
    folder = tmp_path / "elsewhere"
    folder.mkdir()
    monkeypatch.chdir(folder)
    return folder


# -- what a snapshot needs ----------------------------------------------------


@needs_data
def test_snapshot_without_data_or_sample_says_what_is_missing(work):
    empty = Session.create(str(DATA))
    with pytest.raises(CommandError) as error:
        run(empty, "snapshot")
    message = str(error.value)
    assert "no measured spectrum" in message and "XEQ" in message
    assert "no SIM sample" in message and "SIM GET" in message
    assert not list(work.iterdir())


@needs_data
def test_snapshot_without_a_sample_names_only_that(au_spectrum, work):
    loaded = Session.create(str(DATA))
    run(loaded, f"get {au_spectrum}")
    with pytest.raises(CommandError) as error:
        run(loaded, "snap")
    assert "no SIM sample" in str(error.value)
    assert "no measured spectrum" not in str(error.value)


@needs_data
def test_snapshot_of_a_buffer_not_read_from_a_file_refuses(session, work):
    run(session, "copy 1 2", "pointat 2")
    session.buffers[2].path = None
    with pytest.raises(CommandError, match="not read from a file.*WRITE"):
        run(session, "snapshot")
    assert not (work / "au_fit.xeq").exists()


# -- the round trip -----------------------------------------------------------


@needs_data
def test_hand_tuned_session_restores_exactly_from_another_folder(
    session, work, tmp_path, monkeypatch
):
    run(
        session, "correction 0.93", "conversion 5.01 1.5", "sim thick 250",
        "pert window 355 375", "pert thick 1", "region 100 450", "log", "snap",
    )
    for name in ("au_fit.xeq", "au.lcm", "au.pert", "au.png", "au.report"):
        assert (work / name).exists(), name
    saved = snapshot.fingerprint(session)
    assert session.saved_state == saved and not snapshot.unsaved(session)

    elsewhere(tmp_path, monkeypatch)
    again = restored(work / "au_fit.xeq")
    assert snapshot.fingerprint(again) == saved
    assert again.buffers.active_buffer.measurement.correction == 0.93
    assert again.script.layers[0].thickness == 250
    assert (again.plot.low, again.plot.high, again.plot.yscale) == (100, 450, "log")
    assert not snapshot.unsaved(again)  # its last line marks it saved


@needs_data
def test_restore_macro_paths_are_relative_to_its_own_folder(session, work, au_spectrum):
    run(session, "snap")
    text = (work / "au_fit.xeq").read_text()
    relative = Path(os.path.relpath(au_spectrum, work.resolve())).as_posix()
    assert relative.startswith("../")
    assert f'GET "{relative}"' in text
    assert 'SIM GET "au.lcm"' in text
    assert str(work) not in text


@needs_data
def test_xeq_finds_the_restore_macro_by_its_bare_name(session, work):
    run(session, "snap")
    again = Session.create(str(DATA))
    run(again, "xeq au_fit")
    assert snapshot.fingerprint(again) == snapshot.fingerprint(session)


@needs_data
def test_snapshot_into_a_named_folder(session, work):
    (work / "results").mkdir()
    run(session, "snapshot results/gold")
    assert (work / "results" / "gold_fit.xeq").exists()
    assert (work / "results" / "gold.report").exists()
    assert snapshot.fingerprint(restored(work / "results" / "gold_fit.xeq")) == (
        snapshot.fingerprint(session)
    )


@needs_data
def test_snapshot_into_a_missing_folder_is_an_error(session, work):
    with pytest.raises(CommandError, match="no such folder"):
        run(session, "snapshot nowhere/gold")


@needs_data
def test_changed_counts_are_saved_with_the_snapshot(session, work, tmp_path, monkeypatch):
    """BACKGROUND, SMOOTH and the like change the counts after reading:
    reloading the original file would undo that, so the spectrum is saved
    as it is now, next to the macro."""
    run(session, "smooth")
    run(session, "snap")
    assert (work / "au_fit.rbs").exists()
    assert 'GET "au_fit.rbs"' in (work / "au_fit.xeq").read_text()

    elsewhere(tmp_path, monkeypatch)
    again = restored(work / "au_fit.xeq")
    np.testing.assert_allclose(
        again.buffers.active_buffer.spectrum.counts,
        session.buffers.active_buffer.spectrum.counts, rtol=1e-6,
    )
    assert snapshot.stem(again.buffers.active_buffer) == "au"


@needs_data
def test_acquisition_macro_and_paired_pixe_restore(work, tmp_path, monkeypatch):
    """The common case: an RC43 acquisition .RBS read with XEQ, its .PIX
    paired, a fit setup -- restored from another folder."""
    data = tmp_path / "data"
    data.mkdir()
    for name in ("MnPt.RBS", "MnPt.PIX", "MnPt.lcm", "MnPt.pert"):
        shutil.copy(EXAMPLES / name, data / name)
    first = Session.create(str(DATA))
    run(
        first, "pixe pair on", f"xeq {data / 'MnPt.RBS'}", f"sim get {data / 'MnPt.lcm'}",
        f"pert get {data / 'MnPt.pert'}", "pixe fwhm 140", "correction 0.97", "snap",
    )
    text = (work / "MA8408_fit.xeq").read_text()
    assert 'XEQ "../data/MnPt.RBS"' in text
    assert 'PIXE GET "../data/MnPt.PIX"' in text
    assert "PIXE fwhm 140" in text

    elsewhere(tmp_path, monkeypatch)
    again = restored(work / "MA8408_fit.xeq")
    assert snapshot.fingerprint(again) == snapshot.fingerprint(first)
    buffer = again.buffers.active_buffer
    assert buffer.macro == (data / "MnPt.RBS").resolve()
    assert buffer.pixe is not None and again.pixe.pair


@needs_data
def test_a_corrected_identifier_and_date_are_restored(session, work, tmp_path, monkeypatch):
    """Reloading the data file alone would bring back the typo."""
    run(session, "identifier 'au  marker, corrected'", "date '12:00  10-06-2026'", "snap")
    text = (work / "au_fit.xeq").read_text()
    assert "IDENTIFIER 'au  marker, corrected'" in text
    elsewhere(tmp_path, monkeypatch)
    again = restored(work / "au_fit.xeq").buffers.active_buffer
    assert again.identifier == "au  marker, corrected"
    assert again.date == "12:00  10-06-2026"


@needs_data
def test_an_identifier_with_a_quote_uses_the_other_kind(session, work, tmp_path, monkeypatch):
    run(session, 'identifier "Ta 2\' annealed"', "snap")
    assert 'IDENTIFIER "Ta 2\' annealed"' in (work / "Ta_fit.xeq").read_text()
    elsewhere(tmp_path, monkeypatch)
    assert restored(work / "Ta_fit.xeq").buffers.active_buffer.identifier == "Ta 2' annealed"


@needs_data
def test_correcting_the_identifier_is_unsaved_but_keeps_the_fit_a_fit(session, work):
    run(session, "pert window 355 375", "pert thick 1", "pert go", "snap")
    run(session, "identifier 'au fixed'")
    assert snapshot.unsaved(session)
    run(session, "snap")
    assert "fit (GO)" in (work / "au.report").read_text().split("---")[-2]


# -- a name another dataset already has -----------------------------------------


def acquisition(folder: Path, name: str, identifier: str, pix: bool = True) -> Path:
    """The MnPt example as ``<name>.RBS`` (and ``.PIX``), identified as given."""
    folder.mkdir(parents=True, exist_ok=True)
    lines = (EXAMPLES / "MnPt.RBS").read_text().splitlines()
    lines = [f"Identifier '{identifier}" if l.startswith("Identifier") else l for l in lines]
    (folder / f"{name}.RBS").write_text("\n".join(lines) + "\n")
    if pix:
        shutil.copy(EXAMPLES / "MnPt.PIX", folder / f"{name}.PIX")
    return folder / f"{name}.RBS"


def load(path: Path, *lines: str) -> Session:
    loaded = Session.create(str(DATA))
    run(loaded, "pixe pair on", f"xeq {path}", f"sim get {EXAMPLES / 'MnPt.lcm'}", *lines)
    return loaded


@needs_data
def test_an_identifier_another_dataset_has_falls_back_to_the_file_name(
    tmp_path, monkeypatch, capsys
):
    """A123 is the bad spot; the fresh spot A123b, renamed A123 by IDENTIFIER,
    is snapshotted in the same folder."""
    data = tmp_path / "data"
    acquisition(data, "A123", "A123  bad spot")
    good = acquisition(data, "A123b", "A123b  fresh spot")
    monkeypatch.chdir(data)
    session = load(good, "identifier 'A123  good spot'", "snap")
    output = capsys.readouterr().out
    assert "A123 is taken by another dataset (A123.RBS)" in output
    assert "named after the data file instead: A123b" in output
    assert (data / "A123b_fit.xeq").exists() and (data / "A123b.report").exists()
    assert not (data / "A123_fit.xeq").exists() and not (data / "A123.report").exists()

    elsewhere(tmp_path, monkeypatch)
    again = restored(data / "A123b_fit.xeq")
    assert snapshot.fingerprint(again) == snapshot.fingerprint(session)
    assert again.buffers.active_buffer.identifier == "A123  good spot"
    assert again.buffers.active_buffer.macro == good.resolve()


@needs_data
def test_a_namesake_next_to_the_data_counts_with_snapshots_elsewhere(tmp_path, work, capsys):
    data = tmp_path / "data"
    acquisition(data, "A123", "A123  bad spot")
    good = acquisition(data, "A123b", "A123b  fresh spot")
    load(good, "identifier 'A123'", "snap")
    assert "named after the data file instead: A123b" in capsys.readouterr().out
    assert (work / "A123b_fit.xeq").exists()


@needs_data
def test_an_earlier_snapshot_of_another_file_takes_the_name(tmp_path, work, capsys):
    """The two spots in folders of their own: only the bad spot's snapshot
    in the work folder holds the name A123."""
    bad = acquisition(tmp_path / "bad", "A123", "A123  bad spot")
    good = acquisition(tmp_path / "good", "A123b", "A123b  fresh spot")
    load(bad, "snap")
    assert (work / "A123_fit.xeq").exists()
    before = (work / "A123_fit.xeq").read_text()
    load(good, "identifier 'A123'", "snap")
    output = capsys.readouterr().out
    assert "A123_fit.xeq, the snapshot of A123.RBS" in output
    assert (work / "A123b_fit.xeq").exists()
    assert (work / "A123_fit.xeq").read_text() == before


@needs_data
def test_snapshotting_the_same_dataset_again_updates_it(tmp_path, work, capsys):
    """Its own earlier snapshot is no clash -- after SMOOTH (saved as
    _fit.rbs) and after restoring from that, too."""
    good = acquisition(tmp_path / "data", "A123", "A123")
    session = load(good, "snap", "correction 0.9", "snap", "smooth", "snap")
    assert "taken" not in capsys.readouterr().out
    assert (work / "A123_fit.rbs").exists()
    again = restored(work / "A123_fit.xeq")
    run(again, "correction 0.8", "snap")
    assert "taken" not in capsys.readouterr().out
    assert snapshot.fingerprint(restored(work / "A123_fit.xeq")) == snapshot.fingerprint(again)
    assert snapshot.stem(session.buffers.active_buffer) == "A123"


@needs_data
def test_when_both_names_are_taken_it_asks_for_one(tmp_path, monkeypatch):
    data = tmp_path / "data"
    acquisition(data, "A123", "A123  bad spot")
    good = acquisition(data, "A123b", "A123b")
    acquisition(tmp_path / "other", "A123b", "A123b  elsewhere")
    # an A123b snapshot here, of another A123b.RBS
    (data / "A123b_fit.xeq").write_text('! source: "../other/A123b.RBS"\n')
    monkeypatch.chdir(data)
    with pytest.raises(CommandError, match="choose one: SNAPSHOT <name>"):
        load(good, "identifier 'A123'", "snap")


@needs_data
def test_an_explicit_name_is_the_users_choice_but_not_over_another_snapshot(
    tmp_path, work
):
    data = tmp_path / "data"
    acquisition(data, "A123", "A123  bad spot")
    good = acquisition(data, "A123b", "A123b")
    load(good, "snap A123")  # a raw A123.RBS elsewhere: allowed, as asked
    assert (work / "A123_fit.xeq").exists()
    bad = acquisition(tmp_path / "bad", "B7", "B7")
    with pytest.raises(CommandError, match="A123 is taken by another dataset"):
        load(bad, "snap A123")


@needs_data
def test_a_users_own_lcm_of_that_name_is_kept(session, work, capsys):
    (work / "au.lcm").write_text(SAMPLE)
    (work / "au.pert").write_text("window 355 375\n")
    run(session, "pert thick 1", "snap")
    assert "kept your au.lcm as au.lcm.bak" in capsys.readouterr().out
    assert (work / "au.lcm.bak").read_text() == SAMPLE
    assert (work / "au.pert.bak").exists()
    run(session, "snap")  # now its own snapshot's files: replaced, no more backups
    assert not (work / "au.lcm.2.bak").exists()


# -- what the report says -----------------------------------------------------


@needs_data
def test_snapshot_right_after_go_records_the_fit(session, work):
    run(session, "pert window 355 375", "pert thick 1", "pert go", "snap")
    report = (work / "au.report").read_text()
    assert "fit (GO)" in report
    assert "+/-" in report
    assert "set by hand" not in report


@needs_data
def test_snapshot_after_tuning_a_fit_by_hand_says_so(session, work):
    run(session, "pert window 355 375", "pert thick 1", "pert go", "correction 0.9", "snap")
    report = (work / "au.report").read_text()
    assert "set by hand" in report
    assert "uncertainties no longer apply" in report
    assert "+/-" not in report
    assert "reduced chi-square" in report and "PERT error windows 355-375" in report


@needs_data
def test_snapshot_without_any_fit(session, work):
    run(session, "snap")
    report = (work / "au.report").read_text()
    assert "not by GO, so there are no uncertainties" in report
    assert "REGION channels" in report  # no PERT windows: scored over the plot region


@needs_data
def test_report_on_writes_the_restore_macro_too(session, work, capsys):
    run(session, "pert window 355 375", "pert thick 1", "pert report", "pert go")
    output = capsys.readouterr().out
    for line in ("updated au.report", "wrote au.pert", "wrote au.lcm", "wrote au.png",
                 "wrote au_fit.xeq"):
        assert line in output
    assert "fit (GO)" in (work / "au.report").read_text()
    assert not snapshot.unsaved(session)


@needs_data
def test_report_on_for_a_buffer_not_read_from_a_file_still_records_the_fit(
    session, work, capsys
):
    session.buffers.active_buffer.path = None
    run(session, "pert window 355 375", "pert thick 1", "pert report", "pert go")
    assert (work / "au.report").exists()
    assert not (work / "au_fit.xeq").exists()
    assert "no restore macro: buffer 1 was not read from a file" in capsys.readouterr().out


# -- at the prompt ------------------------------------------------------------


@needs_data
def test_snap_at_the_sim_prompt_stays_there(session, work):
    stack = run(session, "sim")
    run(session, "snap", stack=stack)
    assert stack == ["rump", "sim"]
    assert (work / "au_fit.xeq").exists()


@needs_data
def test_quit_warns_about_an_unsaved_session(session, work, monkeypatch):
    questions = []
    monkeypatch.setattr(rump, "_interactive", lambda session: True)
    monkeypatch.setattr("builtins.input", lambda text: questions.append(text) or "n")
    run(session, "quit")
    run(session, "snap", "quit")
    run(session, "correction 0.5", "quit")
    assert "no SNAPSHOT" in questions[0]
    assert questions[1] == "Really quit pyRUMP? [y/N] "
    assert "not saved" in questions[2]


@needs_data
def test_buffers_listing_names_the_acquisition_file(work, tmp_path):
    shutil.copy(EXAMPLES / "MnPt.RBS", tmp_path / "MnPt.RBS")
    loaded = Session.create(str(DATA))
    run(loaded, f"xeq {tmp_path / 'MnPt.RBS'}")
    assert "MnPt.RBS (XEQ)" in loaded.buffers.active_buffer.describe()


# -- finding files from inside a macro ----------------------------------------


@needs_data
def test_a_macro_finds_files_in_its_own_folder_first(tmp_path, monkeypatch):
    macros = tmp_path / "macros"
    macros.mkdir()
    (macros / "here.lcm").write_text(SAMPLE)
    (macros / "load.cmd").write_text("sim get here.lcm\n")
    elsewhere(tmp_path, monkeypatch)
    loaded = Session.create(str(DATA))
    run(loaded, f"xeq {macros / 'load'}")
    assert len(loaded.script.layers) == 2


@needs_data
def test_a_macro_still_finds_files_in_the_working_directory(tmp_path, monkeypatch):
    macros = tmp_path / "macros"
    macros.mkdir()
    (macros / "load.cmd").write_text("sim get there.lcm\n")
    folder = elsewhere(tmp_path, monkeypatch)
    (folder / "there.lcm").write_text(SAMPLE)
    loaded = Session.create(str(DATA))
    run(loaded, f"xeq {macros / 'load'}")
    assert len(loaded.script.layers) == 2


@needs_data
def test_a_bare_name_tries_rump_macro_extensions_in_rump_order(tmp_path, monkeypatch):
    """lexp.c's LexMacroExtList: .mac, .xeq, .cmd -- then pyRUMP's .rbs."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "go.cmd").write_text("mev 1.0\n")
    (tmp_path / "go.xeq").write_text("mev 2.0\n")
    loaded = Session.create(str(DATA))
    run(loaded, "xeq go")
    assert loaded.settings.experiment_defaults.beam.e0_MeV == 2.0
    (tmp_path / "go.mac").write_text("mev 3.0\n")
    run(loaded, "xeq go")
    assert loaded.settings.experiment_defaults.beam.e0_MeV == 3.0


@needs_data
def test_a_cmd_macro_still_runs_by_its_bare_name(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "old.cmd").write_text("mev 1.5\n")
    loaded = Session.create(str(DATA))
    run(loaded, "xeq old")
    assert loaded.settings.experiment_defaults.beam.e0_MeV == 1.5


# -- the command line ---------------------------------------------------------


@pytest.mark.parametrize(
    "argv, expected",
    [
        (["au_fit.xeq"], ["shell", "au_fit.xeq"]),
        (["au_fit"], ["shell", "au_fit"]),  # a bare name, as XEQ takes it
        (["--batch", "au_fit.xeq"], ["shell", "--batch", "au_fit.xeq"]),
        (["--data", "tables", "au_fit.xeq"], ["--data", "tables", "shell", "au_fit.xeq"]),
        (["--norc"], ["shell", "--norc"]),
        (["shell", "au_fit.xeq"], ["shell", "au_fit.xeq"]),
        (["fit", "x.lcm", "y.rbs"], ["fit", "x.lcm", "y.rbs"]),
        (["nonsense"], ["nonsense"]),  # no such macro: left for argparse to reject
        (["--version"], ["--version"]),
        ([], []),
    ],
)
def test_pyrump_with_a_macro_means_the_shell(argv, expected, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "au_fit.xeq").write_text("")
    assert _implicit_shell(argv) == expected


def test_binary_write_keeps_the_identifier_the_stem_comes_from(au_spectrum):
    assert read_rbs(au_spectrum).identifier == "au"
    assert os.path.basename(au_spectrum) == "au.rbs"
