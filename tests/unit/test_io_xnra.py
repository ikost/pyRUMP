"""SIMNRA ``.xnra`` reading and writing, and the GETNRA/WRITENRA commands.

Expected values come from the files themselves: ``examples/MnPt.xnra`` was
saved by SIMNRA 7.04 from the same measurement as ``examples/MnPt.RBS``, so
the RUMP header of the one is the known answer for the other.
"""

from __future__ import annotations

import math
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from pyrump.io.rbs import read_rbs  # noqa: E402
from pyrump.io.xnra import (  # noqa: E402
    ELEMENTARY_CHARGE,
    IDF_NS,
    SIMNRA_NS,
    XnraFormatError,
    XnraLayer,
    exit_angle,
    fluence,
    format_number,
    read_xnra,
    write_xnra,
)
from pyrump.model.geometry import Geometry, GeometryKind  # noqa: E402
from pyrump.physics.xsec.rutherford import ScreeningModel  # noqa: E402
from pyrump.shell.commands import rump  # noqa: E402
from pyrump.shell.dispatch import CommandError  # noqa: E402
from pyrump.shell.repl import execute_line  # noqa: E402
from pyrump.shell.session import Buffer, Session  # noqa: E402

from conftest import data_dir

EXAMPLES = Path(__file__).resolve().parents[2] / "examples"
MNPT = EXAMPLES / "MnPt.xnra"
NS = {"i": IDF_NS, "s": SIMNRA_NS}

DATA = data_dir()
needs_data = pytest.mark.skipif(DATA is None, reason="legacy data tables unavailable")


def _raw(path: Path, xpath: str) -> ET.Element:
    return ET.parse(path).getroot().find(xpath, NS)


def _rbs_header(path: Path) -> dict[str, str]:
    """The keyword lines of an RC43 EMPTY/SWALLOW macro."""
    header = {}
    for line in path.read_text().splitlines():
        if line.strip().lower() == "swallow":
            break
        key, _, rest = line.strip().partition(" ")
        header[key.lower()] = rest.strip()
    return header


def _edited(tmp_path: Path, edit) -> Path:
    """A copy of MnPt.xnra with ``edit(root)`` applied."""
    tree = ET.parse(MNPT)
    edit(tree.getroot())
    target = tmp_path / "edited.xnra"
    tree.write(target, encoding="utf-8", xml_declaration=True)
    return target


# -- reading ----------------------------------------------------------------


def test_mnpt_spectrum_metadata():
    spectrum = read_xnra(MNPT).spectrum
    assert spectrum.e0_MeV == pytest.approx(1.902)
    assert spectrum.zbeam == 2
    assert spectrum.mbeam == pytest.approx(4.00260325413)
    assert spectrum.calibration.kev0 == pytest.approx(50.0)
    assert spectrum.calibration.kevch == pytest.approx(1.59)
    assert spectrum.calibration.first == 0.0
    assert spectrum.measurement.fwhm_keV == pytest.approx(24.0)
    assert spectrum.identifier == "MnPt.RBS"  # the data's legend in SIMNRA
    y = _raw(MNPT, "i:sample/i:spectra/i:spectrum/i:data/i:simpledata/i:y").text.split()
    assert spectrum.counts.size == len(y)
    assert spectrum.counts.sum() == pytest.approx(sum(float(v) for v in y))


def test_mnpt_sample_layers():
    script = read_xnra(MNPT).script
    assert [layer.unit for layer in script.layers] == ["/CM2"] * 4
    assert [list(layer.composition) for layer in script.layers] == [
        ["Ru"], ["Mn", "Pt"], ["Ru"], ["Si"]
    ]
    assert script.layers[1].thickness == pytest.approx(241.491562542493)
    assert script.layers[1].composition["Mn"] == pytest.approx(0.733127677322798)
    assert script.layers[1].composition["Pt"] == pytest.approx(0.266872322677202)
    assert script.layers[3].thickness == pytest.approx(5000.0)


def test_dose_folds_the_solid_angle_into_the_charge():
    """SIMNRA's fluence is particles x solid angle: MnPt.RBS's 16 uC at
    2.7 msr must come back as the same product at OMEGA 1."""
    header = _rbs_header(EXAMPLES / "MnPt.RBS")
    rump_product = float(header["charge"]) * float(header["omega"])
    m = read_xnra(MNPT).spectrum.measurement
    assert m.omega_msr == 1.0
    assert m.charge_state == 1 and m.correction == 1.0
    assert m.omega_msr * m.charge_uC == pytest.approx(rump_product)
    assert fluence(m) == pytest.approx(rump_product * 1e-6 / ELEMENTARY_CHARGE)


def test_ibm_geometry_takes_theta_sign_from_the_exit_angle():
    g = read_xnra(MNPT).geometry
    assert g.kind is GeometryKind.IBM
    assert (g.theta, g.phi, g.psi) == pytest.approx((9.0, 11.0, 20.0))
    assert g.sec_out == pytest.approx(1 / math.cos(math.radians(20.0)))


def test_unsupported_features_are_reported(tmp_path):
    def edit(root):
        spectrum = root.find("i:sample/i:spectra/i:spectrum", NS)
        quad = spectrum.findall(".//i:calibrationparameter", NS)[2]
        quad.text = format_number(1e-5)
        layer = root.find(".//i:layers/i:layer", NS)
        layer.find("s:hasroughness", NS).text = "true"
        spectrum.find("i:reactions/i:technique", NS).text = "ERDA"

    notices = read_xnra(_edited(tmp_path, edit)).notices
    assert any("quadratic calibration term" in n for n in notices)
    assert any("layer 1: roughness ignored" in n for n in notices)
    assert any("technique ERDA" in n for n in notices)


def test_isotopes_fold_into_their_element(tmp_path):
    def edit(root):
        names = root.findall(".//i:layers/i:layer", NS)[3].findall(".//i:name", NS)
        names[0].text = "28Si"

    source = read_xnra(_edited(tmp_path, edit))
    assert source.script.layers[3].composition == {"Si": pytest.approx(1.0)}
    assert any("isotope 28Si read as natural Si" in n for n in source.notices)


def test_generic_idf_is_refused(tmp_path):
    def edit(root):
        attributes = root.find("i:attributes", NS)
        attributes.remove(attributes.find("s:filetype", NS))

    with pytest.raises(XnraFormatError, match="generic IDF"):
        read_xnra(_edited(tmp_path, edit))


# -- writing ----------------------------------------------------------------


def _buffer_from(path: Path) -> Buffer:
    source = read_xnra(path)
    buffer = Buffer.from_rbs(source.spectrum, path)
    buffer.source_xml = source.xml
    return buffer


def test_round_trip_keeps_data_metadata_and_simnra_settings(tmp_path):
    original = read_xnra(MNPT)
    layers = [
        XnraLayer(layer.thickness, dict(layer.composition)) for layer in original.script.layers
    ]
    target = tmp_path / "out.xnra"
    write_xnra(target, _buffer_from(MNPT), layers=layers, base=original.xml)

    again = read_xnra(target)
    np.testing.assert_array_equal(again.spectrum.counts, original.spectrum.counts)
    assert again.spectrum.calibration == original.spectrum.calibration
    assert again.spectrum.geometry == original.spectrum.geometry
    assert again.spectrum.measurement == original.spectrum.measurement
    assert [(l.thickness, l.composition) for l in again.script.layers] == pytest.approx(
        [(l.thickness, l.composition) for l in original.script.layers]
    )
    # Settings pyRUMP does not model survive.
    assert _raw(target, ".//s:scaling/s:xaxis/s:max").text == _raw(
        MNPT, ".//s:scaling/s:xaxis/s:max"
    ).text
    assert _raw(target, ".//i:simulation/i:physics/i:crosssections") is not None
    # Counts agree with the lists they count.
    root = ET.parse(target).getroot()
    assert int(root.find(".//i:nlayers", NS).text) == len(root.findall(".//i:layers/i:layer", NS))
    assert int(root.find(".//i:nelements", NS).text) == len(
        root.findall(".//i:elements/i:element", NS)
    )


@pytest.mark.parametrize(
    "geometry",
    [
        Geometry(theta=9.0, phi=11.0, kind=GeometryKind.CORNELL),
        Geometry(theta=9.0, phi=11.0, kind=GeometryKind.IBM),
        Geometry(theta=-9.0, phi=11.0, kind=GeometryKind.IBM),
        Geometry(theta=30.0, phi=15.0, psi=35.0, kind=GeometryKind.GENERAL),
    ],
    ids=["cornell", "ibm+", "ibm-", "general"],
)
def test_geometry_round_trip(tmp_path, geometry):
    buffer = _buffer_from(MNPT)
    buffer.geometry = geometry
    target = tmp_path / "geometry.xnra"
    write_xnra(target, buffer)

    written = float(_raw(target, ".//i:geometry/i:exitangle").text)
    assert written == pytest.approx(exit_angle(geometry))
    back = read_xnra(target).geometry
    assert back.kind is geometry.kind
    assert (back.theta, back.phi) == pytest.approx((geometry.theta, geometry.phi))
    assert back.sec_out == pytest.approx(geometry.sec_out)


def test_fresh_file_is_a_valid_xnra(tmp_path):
    """No source document: the skeleton alone must satisfy SIMNRA's rules."""
    source = read_rbs(EXAMPLES / "2A.rbs")
    buffer = Buffer.from_rbs(source)
    target = tmp_path / "fresh.xnra"
    write_xnra(target, buffer, layers=[XnraLayer(1000.0, {"Si": 1.0})])

    root = ET.parse(target).getroot()
    assert root.tag == f"{{{IDF_NS}}}idf"
    assert root.find("i:attributes/s:filetype", NS).text == "xnra"
    assert root.find("i:attributes/s:xnraversionnr", NS).text == "1.1"
    assert len(root.findall("i:sample", NS)) == 1
    assert len(root.findall(".//i:spectra/i:spectrum", NS)) == 1

    back = read_xnra(target)
    np.testing.assert_allclose(back.spectrum.counts, source.counts)
    assert back.spectrum.calibration.kevch == pytest.approx(source.calibration.kevch)
    assert back.spectrum.e0_MeV == pytest.approx(source.e0_MeV)
    assert back.script.layers[0].composition == {"Si": 1.0}


def test_empty_sample_reads_as_none(tmp_path):
    target = tmp_path / "empty.xnra"
    write_xnra(target, Buffer.from_rbs(read_rbs(EXAMPLES / "2A.rbs")))
    assert read_xnra(target).script is None
    # SIMNRA's own form: one layer of zero thickness.
    assert _raw(target, ".//i:nlayers").text == "1"


def test_screening_is_written_in_simnra_terms(tmp_path):
    buffer = _buffer_from(MNPT)
    target = tmp_path / "screening.xnra"
    write_xnra(target, buffer, screening=ScreeningModel.NONE, base=buffer.source_xml)
    assert _raw(target, ".//i:crosssectiondefault/i:Rutherford").text == "true"
    notices = write_xnra(target, buffer, screening=ScreeningModel.LECUYER, base=buffer.source_xml)
    assert any("L'Ecuyer" in n for n in notices)


# -- the shell --------------------------------------------------------------


@pytest.fixture
def session() -> Session:
    if DATA is None:
        pytest.skip("legacy data tables unavailable")
    return Session.create(str(DATA))


def run(session: Session, *lines: str) -> None:
    stack = ["rump"]
    for line in lines:
        execute_line(session, line, stack)


def test_command_names_resolve():
    table = rump.TABLE
    assert table.match("GN").handler is rump.cmd_getnra
    assert table.match("getn").handler is rump.cmd_getnra
    assert table.match("WN").handler is rump.cmd_writenra
    assert table.match("writen").handler is rump.cmd_writenra
    # Existing abbreviations are untouched.
    assert table.match("GET").handler is rump.cmd_get
    assert table.match("READ").handler is rump.cmd_get
    assert table.match("WRITE").handler is rump.cmd_write
    assert table.match("GE") is None or table.match("GE").handler is not rump.cmd_getnra


@needs_data
def test_getnra_loads_spectrum_and_sample(session, capsys):
    run(session, f"GN {MNPT}")
    out = capsys.readouterr().out
    assert "OMEGA 1 msr, CHARGE 43.2 uC" in out
    assert "geometry IBM: theta 9, phi 11" in out
    buffer = session.buffers.active_buffer
    assert session.buffers.active == 1
    assert buffer.path == MNPT.resolve()
    assert buffer.source_xml
    assert [list(layer.composition) for layer in session.script.layers][1] == ["Mn", "Pt"]


@needs_data
def test_getnra_data_only_leaves_sim_alone(session):
    run(session, "SIM", "LAYER 1", "THICK 100 A", "COMP Au 1 /", "RETURN")
    run(session, f"GN {MNPT} -data")
    assert session.script.layers[0].composition == {"Au": 1.0}


@needs_data
@pytest.mark.parametrize("answer, kept", [("n", True), ("y", False)])
def test_getnra_asks_before_replacing_a_sample(session, monkeypatch, answer, kept):
    run(session, "SIM", "LAYER 1", "THICK 100 A", "COMP Au 1 /", "RETURN")
    monkeypatch.setattr(rump.sys.stdin, "isatty", lambda: True, raising=False)
    questions = []
    monkeypatch.setattr("builtins.input", lambda q: questions.append(q) or answer)
    run(session, f"GN {MNPT}")
    assert len(questions) == 1 and "Replace SIM sample" in questions[0]
    assert (session.script.layers[0].composition == {"Au": 1.0}) is kept


@needs_data
def test_getnra_in_a_macro_replaces_without_asking(session, tmp_path, monkeypatch):
    run(session, "SIM", "LAYER 1", "THICK 100 A", "COMP Au 1 /", "RETURN")
    monkeypatch.setattr("builtins.input", lambda q: pytest.fail("asked inside a macro"))
    macro = tmp_path / "load.cmd"
    macro.write_text(f"GN {MNPT}\n")
    run(session, f"XEQ {macro}")
    assert list(session.script.layers[1].composition) == ["Mn", "Pt"]


@needs_data
def test_getnra_simulation_goes_into_its_own_buffer(session, capsys):
    run(session, f"GN {MNPT} -sim")
    assert "SIMNRA simulation -> buffer 2" in capsys.readouterr().out
    assert session.buffers.active == 1
    assert session.buffers[1].path == MNPT.resolve()
    simulated = session.buffers[2]
    np.testing.assert_array_equal(
        simulated.spectrum.counts, read_xnra(MNPT).simulation.counts
    )


@needs_data
def test_simulation_only_file_needs_the_flag(session, tmp_path):
    def edit(root):
        data = root.find(".//i:spectrum/i:data/i:simpledata", NS)
        data.find("i:x", NS).text = ""
        data.find("i:y", NS).text = ""

    path = _edited(tmp_path, edit)
    with pytest.raises(CommandError, match="-simulation"):
        run(session, f"GN {path}")
    run(session, f"GN {path} -simulation")
    assert session.buffers.active_buffer.name.startswith("SIMNRA simulation")


@needs_data
def test_get_reads_xnra_as_a_spectrum(session):
    run(session, f"GET {MNPT}")
    assert session.buffers.active_buffer.source_xml
    assert not session.script.layers  # GET never touches SIM


@needs_data
def test_writenra_from_the_rbs_macro_matches_simnra(session, tmp_path):
    """XEQ MnPt.RBS (RUMP's own header), then WN: the fluence must be the one
    SIMNRA itself computed from the same header."""
    target = tmp_path / "from-rbs.xnra"
    run(session, f"XEQ {EXAMPLES / 'MnPt.RBS'}", f"WN {target}")
    written = float(_raw(target, ".//i:beam/i:beamfluence").text)
    expected = float(_raw(MNPT, ".//i:beam/i:beamfluence").text)
    assert written == pytest.approx(expected)


@needs_data
def test_writenra_then_getnra_round_trips_through_the_session(session, tmp_path):
    target = tmp_path / "again"
    run(session, f"GN {MNPT}", f"WN {target}", "NEWALL", f"GN {target}.xnra")
    assert session.buffers.active_buffer.spectrum.counts.sum() == pytest.approx(
        read_xnra(MNPT).spectrum.counts.sum()
    )
    root = ET.parse(f"{target}.xnra").getroot()
    total = root.find(".//i:simulation/i:simpledata/i:y", NS).text.split()
    assert len(total) == session.buffers.active_buffer.n_channels  # pyRUMP's own simulation


@needs_data
def test_pyrump_agrees_with_simnra_on_peak_positions(session):
    """Simulate SIMNRA's own sample with SIMNRA's own settings: the peak
    centroids must agree within half a channel. This catches errors in the
    energy calibration, the beam energy or the angle conversion, which move
    peaks by channels. Heights are not compared -- the stopping and
    straggling models differ."""
    source = read_xnra(MNPT)
    run(session, f"GN {MNPT}", "SCREENING ANDERSEN")
    ours = session.simulation().spectrum.counts
    theirs = source.simulation.counts

    def centroid(y, low, high):
        channels = np.arange(low, high)
        return float((channels * y[low:high]).sum() / y[low:high].sum())

    # Windows around the Pt, surface Ru and Mn signals, read off SIMNRA's
    # simulated spectrum.
    for low, high in ((1030, 1090), (990, 1030), (850, 900)):
        assert centroid(ours, low, high) == pytest.approx(centroid(theirs, low, high), abs=0.5)
