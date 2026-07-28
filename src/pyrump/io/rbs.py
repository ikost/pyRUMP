"""RUMP's native ``.RBS`` binary spectrum format.

Implements the "RUMP Binary [v 1.1] Data Format Specification"
(``html/RUMP/rbs_inf.htm``, Aug 1994), cross-checked against ``rump/rbs_rdwr.c``.

Unlike the rest of the legacy distribution, this format is **explicitly open**:
the specification ships a reference C API and states it "may be freely
incorporated into user written applications", with record numbers above 1000h
reserved for third parties.

Structure
---------

A TIFF-inspired record stream of 32-bit **big-endian** words::

    [length_in_words] [record_type] [data...] [checksum]

``length`` counts every word including itself and the checksum, so a record
carries ``length - 3`` data words and may hold at most 1024 of them.

**Every word in a record, checksum included, sums to zero modulo 2^32.** That is
the format's integrity check and it is verified on read.

Data compression
----------------

Counts are stored in one of four encodings, selected by the record type:

===== ==========================================================
mode  encoding
===== ==========================================================
0     IEEE-754 single-precision floats, unpacked
1     32-bit integers, unpacked
2     differential integers
3     differential integers plus zero-run compression
===== ==========================================================

Differential encoding stores the first value absolutely, then one signed byte
per delta, escaping to a 16-bit delta via ``0x80`` and to a fresh 32-bit
absolute via ``0x8000``.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..model.detector import Measurement
from ..model.geometry import Geometry, GeometryKind
from ..model.spectrum import Calibration, Spectrum

#: Program identifier written by RUMP in the first record.
RUMP_PROGRAM_ID = 0x10211210

#: Record types (rbs_rdwr.c:99-118).
PROGRAM_ID_REC = 0x0000
COMMENT_REC = 0x0001
SILENT_REC = 0x0002
DATA_INIT_REC = 0x0010
DATA_REC = 0x0011
DATA_REC_0 = 0x0012
DATA_REC_1 = 0x0013
DATA_REC_2 = 0x0014
DATA_REC_3 = 0x0015
ARRAY_INIT_REC = 0x0020
ID_REC = 0x0101
LTCT_REC = 0x0102
DATE_REC = 0x0103
CORR_REC = 0x0110
ACCEL_REC = 0x0111
MCA_REC = 0x0112
RBS_REC = 0x0120
FRES_REC = 0x0121
PIXE_REC = 0x0122
NREAC_REC = 0x0123

#: Maximum data words in one record, and items per data record.
MAX_DATA_WORDS = 1024

#: Marks a zero-run-compressed byte stream (rbs_rdwr.c).
ZERO_COMPRESS_FLAG = 0x80

_UINT32_MODULUS = 1 << 32


class RbsFormatError(ValueError):
    """The file is not a well-formed RUMP binary spectrum."""


@dataclass(slots=True)
class RbsSpectrum:
    """A spectrum plus the metadata the format carries."""

    counts: np.ndarray
    calibration: Calibration
    geometry: Geometry
    measurement: Measurement
    e0_MeV: float = 0.0
    zbeam: int = 0
    mbeam: float = 0.0
    identifier: str = ""
    date: str = ""
    livetime: str = ""
    spectrum_type: str = "RBS"
    nspectra: int = 1
    comments: list[str] = field(default_factory=list)

    def to_spectrum(self) -> Spectrum:
        return Spectrum(counts=self.counts.copy(), calibration=self.calibration)


# ------------------------------------------------------------------ reading


class _Reader:
    """Cursor over one record's data words, with byte-level access."""

    def __init__(self, payload: bytes):
        self.data = payload
        self.pos = 0

    def uint32(self) -> int:
        value = struct.unpack_from(">I", self.data, self.pos)[0]
        self.pos += 4
        return value

    def int32(self) -> int:
        value = struct.unpack_from(">i", self.data, self.pos)[0]
        self.pos += 4
        return value

    def real(self) -> float:
        value = struct.unpack_from(">f", self.data, self.pos)[0]
        self.pos += 4
        return float(value)

    def int16(self) -> int:
        value = struct.unpack_from(">h", self.data, self.pos)[0]
        self.pos += 2
        return value

    def int8(self) -> int:
        value = struct.unpack_from(">b", self.data, self.pos)[0]
        self.pos += 1
        return value

    def string(self) -> str:
        """Length in bytes, then characters packed four per word."""
        length = self.int32()
        raw = self.data[self.pos : self.pos + length]
        # Advance over the whole words the characters occupy.
        self.pos += ((length + 3) // 4) * 4
        return raw.decode("latin-1").rstrip("\x00")

    @property
    def remaining(self) -> int:
        return len(self.data) - self.pos


def _records(blob: bytes):
    """Walk the record stream, verifying each checksum."""
    offset = 0
    total = len(blob)
    while offset + 12 <= total:
        length = struct.unpack_from(">I", blob, offset)[0]
        if not 3 <= length <= MAX_DATA_WORDS + 3:
            raise RbsFormatError(
                f"record at byte {offset} declares {length} words, outside 3..1027"
            )
        end = offset + 4 * length
        if end > total:
            raise RbsFormatError(f"record at byte {offset} runs past end of file")

        words = struct.unpack_from(f">{length}I", blob, offset)
        if sum(words) % _UINT32_MODULUS != 0:
            raise RbsFormatError(
                f"checksum failure at byte {offset}: record words sum to "
                f"{sum(words) % _UINT32_MODULUS:#010x}, expected 0"
            )

        record_type = words[1]
        payload = blob[offset + 8 : end - 4]
        yield record_type, payload
        offset = end

    if offset != total:
        raise RbsFormatError(f"{total - offset} trailing bytes after last record")


def _unzero(data: bytes) -> bytes:
    """Expand zero-run compression (``UnZeroCompress``, rbs_rdwr.c:889-916).

    A stream is marked by a leading ``0x80``; the next byte is the repeat
    marker. ``<marker><n>`` expands to *n* zeros, and ``<marker><0>`` is a
    literal marker byte.
    """
    if not data or data[0] != ZERO_COMPRESS_FLAG:
        return data

    marker = data[1]
    out = bytearray()
    i = 2
    while i < len(data):
        byte = data[i]
        i += 1
        if byte != marker or i >= len(data):
            out.append(byte)
            continue
        count = data[i]
        i += 1
        if count == 0:
            out.append(marker)
        else:
            out.extend(b"\x00" * count)
    return bytes(out)


def _read_differential(reader: _Reader, wanted: int) -> list[float]:
    """Decode differential integers (``read_compress``, rbs_rdwr.c:818).

    One record holds at most 1024 values: an absolute first value followed by
    signed deltas, escaping to wider forms as needed.
    """
    out: list[float] = []
    last = reader.int32()
    out.append(float(last))

    count = 1
    while len(out) < wanted and count != MAX_DATA_WORDS:
        count += 1
        if reader.remaining < 1:
            break
        byte = reader.int8()
        if (byte & 0xFF) != 0x80:
            last += byte
        else:
            word = reader.int16()
            if (word & 0xFFFF) != 0x8000:
                last += word
            else:
                last = reader.int32()
        out.append(float(last))
    return out


def read_rbs(path: str | Path) -> RbsSpectrum:
    """Read a ``.RBS`` binary spectrum."""
    blob = Path(path).read_bytes()

    counts: list[float] = []
    npt = nspectra = 0
    compression = 0
    calibration = Calibration()
    geometry = Geometry(theta=0.0, phi=0.0)
    kevch = kev0 = first = fwhm = 0.0
    corr = 1.0
    e0 = mbeam = q = current = 0.0
    zbeam = cbeam = 0
    geom = GeometryKind.CORNELL
    theta = phi = psi = omega = 0.0
    identifier = date = livetime = ""
    spectrum_type = "RBS"
    comments: list[str] = []
    saw_program_id = False

    for record_type, payload in _records(blob):
        reader = _Reader(payload)

        if record_type == PROGRAM_ID_REC:
            program = reader.uint32()
            if program != RUMP_PROGRAM_ID:
                raise RbsFormatError(
                    f"not a RUMP file: program id {program:#010x}, "
                    f"expected {RUMP_PROGRAM_ID:#010x}"
                )
            saw_program_id = True
        elif record_type in (COMMENT_REC, SILENT_REC):
            comments.append(reader.string())
        elif record_type == ID_REC:
            identifier = reader.string()
        elif record_type == LTCT_REC:
            livetime = reader.string()
        elif record_type == DATE_REC:
            date = reader.string()
        elif record_type == CORR_REC:
            corr = reader.real()
        elif record_type == ACCEL_REC:
            e0 = reader.real()
            zbeam = reader.int32()
            mbeam = reader.real()
            cbeam = reader.int32()
            q = reader.real()
            current = reader.real()
        elif record_type == MCA_REC:
            kevch = reader.real()
            kev0 = reader.real()
            first = reader.real()
            fwhm = reader.real()
        elif record_type in (RBS_REC, FRES_REC):
            spectrum_type = "RBS" if record_type == RBS_REC else "FRES"
            geom = GeometryKind(reader.int32())
            theta = reader.real()
            phi = reader.real()
            psi = reader.real()
            omega = reader.real()
        elif record_type in (DATA_INIT_REC, ARRAY_INIT_REC):
            compression = reader.int32()
            if not 0 <= compression <= 3:
                raise RbsFormatError(f"unknown compression mode {compression}")
            npt = reader.int32()
            nspectra = reader.int32() if record_type == ARRAY_INIT_REC else 1
        elif record_type in (DATA_REC, DATA_REC_0, DATA_REC_1, DATA_REC_2, DATA_REC_3):
            # 0x0011 is the generic data record: its encoding comes from the
            # preceding init record. 0x0012-0x0015 name an encoding explicitly
            # and override it (rbs_rdwr.c:774-779). The shipped fixtures all
            # use the generic form, so a reader that only handles 0x0012-0x0015
            # parses every header correctly and silently returns zeros.
            mode = compression if record_type == DATA_REC else record_type - DATA_REC_0
            wanted = npt * nspectra - len(counts)
            if wanted <= 0:
                continue
            if mode == 0:
                counts.extend(
                    reader.real() for _ in range(min(MAX_DATA_WORDS, wanted))
                )
            elif mode == 1:
                counts.extend(
                    float(reader.int32()) for _ in range(min(MAX_DATA_WORDS, wanted))
                )
            else:
                if mode == 3:
                    reader = _Reader(_unzero(payload))
                counts.extend(_read_differential(reader, wanted))
        elif record_type in (PIXE_REC, NREAC_REC):
            spectrum_type = "PIXE" if record_type == PIXE_REC else "NUCLEAR"
        # Unknown record types are skipped: the spec reserves >=1000h for
        # third-party use, so an unrecognised record is not an error.

    if not saw_program_id:
        raise RbsFormatError("missing program-identifier record")
    if npt == 0:
        raise RbsFormatError("file contains no data records")

    values = np.array(counts[: npt * nspectra], dtype=np.float64)
    return RbsSpectrum(
        counts=values,
        calibration=Calibration(kevch=kevch, kev0=kev0, first=first, npt=npt),
        geometry=Geometry(theta=theta, phi=phi, psi=psi, kind=geom),
        measurement=Measurement(
            omega_msr=omega,
            charge_uC=q,
            correction=corr,
            charge_state=cbeam,
            current_nA=current,
            fwhm_keV=fwhm,
        ),
        e0_MeV=e0,
        zbeam=zbeam,
        mbeam=mbeam,
        identifier=identifier,
        date=date,
        livetime=livetime,
        spectrum_type=spectrum_type,
        nspectra=nspectra,
        comments=comments,
    )


# ------------------------------------------------------------------ writing


def _record(record_type: int, payload: bytes) -> bytes:
    """Frame a payload as a record, appending the zeroing checksum."""
    if len(payload) % 4:
        raise ValueError("record payload must be a whole number of words")
    length = len(payload) // 4 + 3
    header = struct.pack(">II", length, record_type)
    words = struct.unpack(f">{length - 1}I", header + payload)
    checksum = (-sum(words)) % _UINT32_MODULUS
    return header + payload + struct.pack(">I", checksum)


def _string_payload(text: str) -> bytes:
    raw = text.encode("latin-1")
    padded = raw + b"\x00" * ((-len(raw)) % 4)
    return struct.pack(">i", len(raw)) + padded


def write_rbs(path: str | Path, spectrum: RbsSpectrum, *, version: str = "1.0") -> None:
    """Write a ``.RBS`` binary spectrum.

    Writes compression mode 0 (unpacked floats), which is what RUMP itself
    defaults to: it writes 1.0-level files unless ``CONFIG WRITE_LEVEL`` raises
    it, because 1.1 features break 1.0 readers.
    """
    if version not in ("1.0", "1.1"):
        raise ValueError(f"unknown format version {version!r}")

    major, minor = 1, (0 if version == "1.0" else 1)
    out = bytearray()
    out += _record(PROGRAM_ID_REC, struct.pack(">IHH", RUMP_PROGRAM_ID, major, minor))

    for comment in spectrum.comments:
        out += _record(COMMENT_REC, _string_payload(comment))
    if spectrum.identifier:
        out += _record(ID_REC, _string_payload(spectrum.identifier))
    if spectrum.livetime:
        out += _record(LTCT_REC, _string_payload(spectrum.livetime))
    if spectrum.date:
        out += _record(DATE_REC, _string_payload(spectrum.date))

    measurement = spectrum.measurement
    out += _record(CORR_REC, struct.pack(">f", measurement.correction))
    out += _record(
        ACCEL_REC,
        struct.pack(
            ">fififf",
            spectrum.e0_MeV,
            spectrum.zbeam,
            spectrum.mbeam,
            measurement.charge_state,
            measurement.charge_uC,
            measurement.current_nA,
        ),
    )
    calibration = spectrum.calibration
    out += _record(
        MCA_REC,
        struct.pack(
            ">ffff",
            calibration.kevch,
            calibration.kev0,
            calibration.first,
            measurement.fwhm_keV,
        ),
    )
    geometry = spectrum.geometry
    out += _record(
        RBS_REC,
        struct.pack(
            ">iffff",
            int(geometry.kind),
            geometry.theta,
            geometry.phi,
            geometry.psi,
            measurement.omega_msr,
        ),
    )

    npt = int(calibration.npt)
    out += _record(DATA_INIT_REC, struct.pack(">ii", 0, npt))

    values = np.asarray(spectrum.counts, dtype=np.float32)
    for start in range(0, npt, MAX_DATA_WORDS):
        chunk = values[start : start + MAX_DATA_WORDS]
        out += _record(DATA_REC_0, chunk.astype(">f4").tobytes())

    Path(path).write_bytes(bytes(out))
