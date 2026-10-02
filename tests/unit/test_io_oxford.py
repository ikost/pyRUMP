"""Oxford Instruments ASCII spectra (RC43's ``.PIX``).

Expected values are read off the files themselves, so the test can't drift
from what the header actually says.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

from pyrump.io.oxford import OxfordFormatError, is_oxford, read_oxford
from pyrump.model.spectrum import Calibration
from pyrump.pixe.data import PixeData

EXAMPLES = Path(__file__).resolve().parents[2] / "examples"
MNPT_PIX = EXAMPLES / "MnPt.PIX"


def _header_value(path: Path, pattern: str) -> float:
    match = re.search(pattern, path.read_text(encoding="latin-1"))
    return float(match.group(1))


def _rows(path: Path) -> list[tuple[int, float]]:
    lines = path.read_text(encoding="latin-1").splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip().lower().startswith("channel"))
    return [(int(a), float(b)) for a, b in (line.split() for line in lines[start + 1 :] if line.strip())]


def test_rc43_file():
    spectrum = read_oxford(MNPT_PIX)
    rows = _rows(MNPT_PIX)
    assert spectrum.counts.size == len(rows) == 2048
    assert spectrum.first_channel == rows[0][0] == 1
    assert spectrum.counts.sum() == pytest.approx(sum(count for _, count in rows))
    assert spectrum.live_time_s == _header_value(MNPT_PIX, r"Live Time:\s*([\d.]+)")
    assert spectrum.real_time_s == _header_value(MNPT_PIX, r"Real Time:\s*([\d.]+)")
    gain = _header_value(MNPT_PIX, r"A =\s*([-\d.E+]+)")
    offset = _header_value(MNPT_PIX, r"B =\s*([-\d.E+]+)")
    assert spectrum.calibration == Calibration(kevch=gain, kev0=offset, first=1.0, npt=2048)
    # Channel N starts at N*A + B: the file's own channel numbers, not indices.
    assert spectrum.calibration.edge_energy(0) == pytest.approx(gain + offset)
    assert "PIXE X-Ray" in spectrum.identifier
    assert spectrum.notices == []


def _write_old_export(path: Path, *, first: int = 0, channels=None, heading: bool = True) -> Path:
    """The layout of GUPIX's 1998 sample: tabs, "Course", 0-based channels,
    and A/B that are not a keV/channel calibration."""
    channels = channels if channels is not None else range(first, first + 8)
    body = "\n".join(f"{c}\t\t{10 * (c - first)}" for c in channels)
    path.write_text(
        "Oxford Instruments, Inc\n\n\nID: \n\n"
        "Acquisition Date: 16-Dec-1998 16:45\n"
        "Elapsed Real Time:  \t600.00\n"
        "Elapsed Live Time:\t599.00\n"
        "Conversion Gain:\t2048\t\t\tCalibration\n"
        "High Voltage:\t\t0\t\t\tA = -2.560E+001\n"
        "Course Gain:\t\t0\t\t\tB = 1.903E-004\n"
        "Fine Gain:\t\t0.00\t\t\tC = -1.263E-004\n\n\n"
        + ("Channel\t\tContents\n" if heading else "")
        + body + "\n",
        encoding="latin-1",
    )
    return path


def test_older_export_without_a_usable_calibration(tmp_path):
    spectrum = read_oxford(_write_old_export(tmp_path / "old.txt"))
    assert spectrum.first_channel == 0
    np.testing.assert_array_equal(spectrum.counts, 10.0 * np.arange(8))
    assert (spectrum.live_time_s, spectrum.real_time_s) == (599.0, 600.0)
    assert spectrum.date == "16-Dec-1998 16:45"
    assert spectrum.calibration is None
    assert any("not a keV/channel" in notice for notice in spectrum.notices)


def test_default_calibration_goes_on_the_files_own_channel_numbers(tmp_path):
    default = Calibration(kevch=0.01, kev0=-0.05, first=1.0, npt=2048)
    data, notices = PixeData.read(_write_old_export(tmp_path / "old.txt"), default)
    assert data.calibration == Calibration(kevch=0.01, kev0=-0.05, first=0.0, npt=8)
    assert data.live_fraction == pytest.approx(599.0 / 600.0)
    assert notices


def test_rejects_other_files(tmp_path):
    assert not is_oxford(EXAMPLES / "MnPt.RBS")
    with pytest.raises(OxfordFormatError, match="not an Oxford"):
        read_oxford(EXAMPLES / "MnPt.RBS")
    with pytest.raises(OxfordFormatError, match="Channel Contents"):
        read_oxford(_write_old_export(tmp_path / "a.txt", heading=False))
    with pytest.raises(OxfordFormatError, match="not consecutive"):
        read_oxford(_write_old_export(tmp_path / "b.txt", channels=[0, 1, 3]))
