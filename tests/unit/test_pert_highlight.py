"""HIGHLIGHT: PERT's windows shaded on the RBS and PIXE plots.

Driven on the MnPt example (RBS acquisition macro with its .PIX), through
the shell as a user would; the shading is read back off the axes.
"""

from __future__ import annotations

import contextlib
import io
from pathlib import Path

import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from matplotlib.colors import to_rgba  # noqa: E402

from pyrump.shell import snapshot  # noqa: E402
from pyrump.shell.plotting import FIT_TONE, NORM_TONE  # noqa: E402
from pyrump.shell.repl import execute_line  # noqa: E402
from pyrump.shell.session import Session  # noqa: E402

from conftest import data_dir

DATA = data_dir()
EXAMPLES = Path(__file__).resolve().parents[2] / "examples"
needs_data = pytest.mark.skipif(DATA is None, reason="legacy data tables unavailable")


def run(session: Session, *lines: str) -> None:
    stack = ["rump"]
    with contextlib.redirect_stdout(io.StringIO()):
        for line in lines:
            execute_line(session, line, stack)


def spans(ax, tone) -> list[tuple[float, float]]:
    """The x ranges shaded in ``tone`` on ``ax``."""
    color = to_rgba(tone["color"], tone["alpha"])
    found = []
    for patch in ax.patches:
        if patch.get_facecolor() == pytest.approx(color):
            x = [vertex[0] for vertex in patch.get_path().vertices]
            low, high = patch.get_transform().transform([(min(x), 0), (max(x), 0)])[:, 0]
            found.append(ax.transData.inverted().transform([(low, 0), (high, 0)])[:, 0])
    return [tuple(round(v, 3) for v in span) for span in found]


@pytest.fixture
def session():
    if DATA is None:
        pytest.skip("legacy data tables unavailable")
    built = Session.create(str(DATA))
    run(built, "pixe pair on", f"xeq {EXAMPLES / 'MnPt.RBS'}", f"sim get {EXAMPLES / 'MnPt.lcm'}")
    yield built
    import matplotlib.pyplot as plt

    plt.close("all")


@needs_data
def test_compare_shades_the_windows_in_both_panels(session):
    run(session, "pert window 450 700", "pert window 900 1150", "pert normalize 1180 1250",
        "compare")
    top, bottom = session.figure.axes[:2]
    for ax in (top, bottom):
        assert spans(ax, FIT_TONE) == [(449.5, 700.5), (899.5, 1150.5)]
        assert spans(ax, NORM_TONE) == [(1179.5, 1250.5)]


@needs_data
def test_highlight_off_and_on_redraws_at_once(session):
    run(session, "compare", "pert window 450 700")
    assert spans(session.figure.axes[0], FIT_TONE)  # adding the window redrew
    run(session, "pert highlight off")
    assert spans(session.figure.axes[0], FIT_TONE) == []
    run(session, "pert highlight")
    assert spans(session.figure.axes[0], FIT_TONE) == [(449.5, 700.5)]
    run(session, "pert window clear")
    assert spans(session.figure.axes[0], FIT_TONE) == []


@needs_data
def test_clear_and_get_redraw_the_shading(session, tmp_path):
    run(session, "pert window 450 700", "pert pixwin 575 600", "compare", "pixe", "return")
    run(session, "pert clear")
    assert spans(session.figure.axes[0], FIT_TONE) == []
    assert spans(session.pixe.figure.axes[0], FIT_TONE) == []
    (tmp_path / "setup.pert").write_text("window 900 1150\nnormalize 1180 1250\n")
    run(session, f"pert get {tmp_path / 'setup.pert'}")
    assert spans(session.figure.axes[0], FIT_TONE) == [(899.5, 1150.5)]
    assert spans(session.figure.axes[0], NORM_TONE) == [(1179.5, 1250.5)]
    (tmp_path / "empty.pert").write_text("")
    run(session, f"pert get {tmp_path / 'empty.pert'}")
    assert spans(session.figure.axes[0], FIT_TONE) == []


@needs_data
def test_plot_shades_too_and_follows_the_energy_axis(session):
    run(session, "pert window 450 700", "plot")
    ax = session.figure.axes[0]
    assert spans(ax, FIT_TONE) == [(449.5, 700.5)]
    run(session, "energy")
    calibration = session.buffers.active_buffer.calibration
    expected = tuple(round(float(calibration.edge_energy(c)), 3) for c in (449.5, 700.5))
    assert spans(session.figure.axes[0], FIT_TONE) == [expected]


@needs_data
def test_a_window_outside_region_does_not_widen_the_view(session):
    run(session, "region 300 800", "compare")
    before = session.figure.axes[0].get_xlim()
    run(session, "pert window 1000 1100")
    assert session.figure.axes[0].get_xlim() == pytest.approx(before)


@needs_data
def test_pixe_windows_are_shaded_in_the_pixe_window(session):
    run(session, "compare", "pixe", "return", "pert pixwin 575 600")
    calibration = session.buffers.active_buffer.pixe.calibration
    first = round(calibration.first)
    expected = tuple(
        round(float(calibration.edge_energy(c - first)), 3) for c in (574.5, 600.5)
    )
    for ax in session.pixe.figure.axes[:2]:
        assert spans(ax, FIT_TONE) == [expected]
    run(session, "pert highlight off")
    assert spans(session.pixe.figure.axes[0], FIT_TONE) == []


@needs_data
def test_highlight_is_a_standing_preference(session, tmp_path):
    run(session, "pert highlight off", "pert clear")
    assert session.pert.highlight is False
    (tmp_path / "setup.pert").write_text("window 450 700\n")
    run(session, f"pert get {tmp_path / 'setup.pert'}")
    assert session.pert.highlight is False
    assert "highlight" not in (tmp_path / "setup.pert").read_text()


@needs_data
def test_a_snapshot_keeps_highlight_off(session, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    run(session, "pert window 450 700", "pert highlight off", "snap")
    macro = next(tmp_path.glob("*_fit.xeq")).read_text()
    assert "PERT HIGHLIGHT OFF" in macro
    assert snapshot.fingerprint(session) is not None
