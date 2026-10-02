"""Oxford Instruments ASCII spectra, as NEC's RC43 writes them for PIXE.

The file is a free-form header, a ``Channel  Contents`` heading, then one
``channel count`` pair per line::

    Oxford Instruments, Inc
    ID:           MnPt.PIX  PIXE X-Ray  LT =  951.132 RT  953.434 Gain  61.6
    Acquisition Date:           08-19-2026 16:43:34
    Elapsed Real Time:           953.434
    Elapsed Live Time:           951.132
    Conversion Gain:             2048 Calibration
    High Voltage:   -120        A =            1.009699E-02
    Coarse Gain:   44           B =           -3.864563E-02
    Fine Gain:   1.10           C = -1.263E-004

    Channel   Contents
     1             0
     2             0

Exports differ in the details: RC43 numbers channels from 1 and uses spaces,
an older export (GUPIX's sample ``Oxford.txt``) numbers from 0, uses tabs
and spells "Course Gain". The header is therefore read by label, never by
position, and the channel numbers are taken from the file.

**Calibration.** In RC43's files ``A`` is the gain in keV per channel and
``B`` the offset in keV. The older export writes ``A = -25.6``,
``B = 1.9E-4`` -- not a keV/channel gain at all -- so the header calibration
is accepted only when it looks like one (:data:`GAIN_RANGE`), and is never
more than a starting value: PIXE fits refine gain and offset. ``C`` is the
same ``-1.263E-004`` in every file seen, the 1998 one included, so it is a
template constant rather than a quadratic term, and is ignored.

pyRUMP's calibration gives a channel's *lower edge*
(:mod:`pyrump.model.spectrum`), with ``first`` set to the file's first
channel number, so channel *N* starts at ``N*A + B``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..model.spectrum import Calibration

#: The first line that identifies the format.
SIGNATURE = "oxford instruments"

#: Header gains (keV per channel) accepted as a calibration: 1 to 50 eV/ch.
GAIN_RANGE = (0.001, 0.05)

_NUMBER = r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"
_HEADER = {
    "real": re.compile(r"elapsed\s+real\s+time\s*:\s*" + _NUMBER, re.I),
    "live": re.compile(r"elapsed\s+live\s+time\s*:\s*" + _NUMBER, re.I),
    "gain": re.compile(r"\bA\s*=\s*" + _NUMBER),
    "offset": re.compile(r"\bB\s*=\s*" + _NUMBER),
}
_ID = re.compile(r"^\s*ID\s*:\s*(.*?)\s*$", re.I)
_DATE = re.compile(r"acquisition\s+date\s*:\s*(.*?)\s*$", re.I)
_HEADING = re.compile(r"^\s*channel\s+contents\s*$", re.I)


class OxfordFormatError(ValueError):
    """The file is not an Oxford Instruments spectrum pyRUMP can read."""


@dataclass(slots=True)
class OxfordSpectrum:
    """What :func:`read_oxford` found."""

    counts: np.ndarray
    first_channel: int
    """Number the file gives its first channel (1 for RC43, 0 for others)."""

    calibration: Calibration | None
    """From the header's ``A``/``B``, or ``None`` when they don't look like a
    keV/channel calibration (see :data:`GAIN_RANGE`)."""

    live_time_s: float | None = None
    real_time_s: float | None = None
    identifier: str = ""
    date: str = ""
    notices: list[str] = field(default_factory=list)


def is_oxford(path: str | Path) -> bool:
    """Whether the file's first non-blank line is the Oxford signature."""
    try:
        with open(path, encoding="latin-1") as handle:
            for line in handle:
                if line.strip():
                    return line.strip().lower().startswith(SIGNATURE)
    except OSError:
        return False
    return False


def read_oxford(path: str | Path) -> OxfordSpectrum:
    """Read an Oxford Instruments ASCII spectrum."""
    path = Path(path)
    lines = path.read_text(encoding="latin-1").splitlines()
    if not any(line.strip() for line in lines) or not is_oxford(path):
        raise OxfordFormatError(
            f"{path.name} is not an Oxford Instruments spectrum "
            f"(its first line is not {SIGNATURE!r})"
        )

    heading = next((i for i, line in enumerate(lines) if _HEADING.match(line)), None)
    if heading is None:
        raise OxfordFormatError(f"{path.name} has no 'Channel Contents' heading")
    header = lines[:heading]

    values: dict[str, float] = {}
    identifier = date = ""
    for line in header:
        for key, pattern in _HEADER.items():
            match = pattern.search(line)
            if match and key not in values:
                values[key] = float(match.group(1))
        if not identifier and (match := _ID.match(line)) and match.group(1):
            identifier = match.group(1)
        if not date and (match := _DATE.search(line)):
            date = match.group(1)

    channels, counts = _channel_data(lines[heading + 1 :], path)
    notices: list[str] = []
    calibration = None
    gain, offset = values.get("gain"), values.get("offset")
    if gain is not None and offset is not None and GAIN_RANGE[0] <= gain <= GAIN_RANGE[1]:
        calibration = Calibration(
            kevch=gain, kev0=offset, first=float(channels[0]), npt=counts.size
        )
    else:
        notices.append(
            f"header calibration A = {gain}, B = {offset} is not a keV/channel "
            "gain and offset: the PIXE calibration setting is used instead"
        )

    return OxfordSpectrum(
        counts=counts,
        first_channel=int(channels[0]),
        calibration=calibration,
        live_time_s=values.get("live"),
        real_time_s=values.get("real"),
        identifier=identifier,
        date=date,
        notices=notices,
    )


def _channel_data(lines: list[str], path: Path) -> tuple[np.ndarray, np.ndarray]:
    """``(channel numbers, counts)`` from the ``channel count`` lines."""
    channels: list[int] = []
    counts: list[float] = []
    for number, line in enumerate(lines, start=1):
        fields = line.split()
        if not fields:
            continue
        if len(fields) != 2:
            raise OxfordFormatError(
                f"{path.name}: data line {number} is {line.strip()!r}, "
                "expected 'channel count'"
            )
        try:
            channels.append(int(float(fields[0])))
            counts.append(float(fields[1]))
        except ValueError:
            raise OxfordFormatError(
                f"{path.name}: data line {number} is {line.strip()!r}, "
                "expected two numbers"
            ) from None
    if not counts:
        raise OxfordFormatError(f"{path.name} has no channel data")
    channel_array = np.array(channels)
    if np.any(np.diff(channel_array) != 1):
        raise OxfordFormatError(f"{path.name}: channel numbers are not consecutive")
    return channel_array, np.array(counts, dtype=np.float64)
