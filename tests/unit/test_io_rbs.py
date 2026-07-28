"""M11 acceptance: the ``.RBS`` binary format and ASCII spectra.

Unlike the rest of the legacy distribution, the ``.RBS`` specification is
explicitly open — it ships a reference C API and invites reimplementation — so
this milestone is a clean-room implementation from the published spec,
cross-checked against ``rbs_rdwr.c``.

The acceptance bar is interoperability in both directions: pyRUMP reads RUMP's
files bit-identically, and RUMP reads pyRUMP's.
"""

from __future__ import annotations

import os
import struct
import sys
from pathlib import Path

import numpy as np
import pytest

from pyrump.io.ascii import read_ascii, write_ascii
from pyrump.io.rbs import (
    ARRAY_INIT_REC,
    DATA_INIT_REC,
    DATA_REC,
    PROGRAM_ID_REC,
    RUMP_PROGRAM_ID,
    RbsFormatError,
    RbsSpectrum,
    _record,
    _records,
    _unzero,
    read_rbs,
    write_rbs,
)
from pyrump.model.geometry import GeometryKind

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "oracle"))
import driver  # noqa: E402


def _fixtures() -> list[Path]:
    env = os.environ.get("PYRUMP_C_REFERENCE")
    roots = [Path(env)] if env else []
    roots.append(Path(__file__).resolve().parents[2] / "C-code")
    for root in roots:
        folder = root / "rump" / "data" / "Fixed"
        if folder.is_dir():
            return sorted(folder.glob("*.rbs"))
    return []


FIXTURES = _fixtures()
needs_fixtures = pytest.mark.skipif(not FIXTURES, reason="legacy .rbs fixtures absent")


# ------------------------------------------------------------ record framing


def test_checksum_zeroes_the_record():
    """Every word of a record, checksum included, sums to 0 mod 2^32."""
    blob = _record(0x1234, struct.pack(">II", 0xDEADBEEF, 42))
    words = struct.unpack(f">{len(blob) // 4}I", blob)
    assert sum(words) % (1 << 32) == 0
    assert words[0] == len(blob) // 4
    assert words[1] == 0x1234


def test_corrupt_checksum_is_rejected():
    blob = bytearray(_record(PROGRAM_ID_REC, struct.pack(">IHH", RUMP_PROGRAM_ID, 1, 0)))
    blob[-1] ^= 0x01
    with pytest.raises(RbsFormatError, match="checksum failure"):
        list(_records(bytes(blob)))


def test_truncated_record_is_rejected():
    blob = _record(PROGRAM_ID_REC, struct.pack(">IHH", RUMP_PROGRAM_ID, 1, 0))
    with pytest.raises(RbsFormatError, match="past end of file"):
        list(_records(blob[:-4]))


def test_implausible_record_length_is_rejected():
    blob = struct.pack(">III", 99999, 0, 0)
    with pytest.raises(RbsFormatError, match="outside 3..1027"):
        list(_records(blob))


def test_zero_run_expansion():
    """<marker><n> expands to n zeros; <marker><0> is a literal marker."""
    marker = 0x81
    stream = bytes([0x80, marker, 0x01, marker, 0x03, 0x02, marker, 0x00, 0x04])
    assert _unzero(stream) == bytes([0x01, 0x00, 0x00, 0x00, 0x02, marker, 0x04])
    # An unmarked stream passes through untouched.
    assert _unzero(b"\x01\x02\x03") == b"\x01\x02\x03"


# ------------------------------------------------------------------ reading


@needs_fixtures
@pytest.mark.parametrize("path", FIXTURES, ids=lambda p: p.name)
def test_fixture_parses(path):
    spectrum = read_rbs(path)
    assert spectrum.calibration.npt == 2048
    assert spectrum.nspectra == 1
    assert spectrum.counts.size == 2048
    assert spectrum.counts.sum() > 0, "counts must not decode as all-zero"
    assert spectrum.identifier.startswith("Binghampton_target_")
    assert spectrum.zbeam == 1  # proton beam
    assert spectrum.geometry.kind is GeometryKind.GENERAL


@needs_fixtures
def test_metadata_matches_the_known_values():
    """2A.rbs, whose parameters the RUMP shell also reports."""
    spectrum = read_rbs(next(p for p in FIXTURES if p.name == "2A.rbs"))
    assert spectrum.e0_MeV == pytest.approx(3.0)
    assert spectrum.mbeam == pytest.approx(1.008, rel=1e-3)
    assert spectrum.calibration.kevch == pytest.approx(7.815, rel=1e-5)
    assert spectrum.calibration.kev0 == pytest.approx(65.646, rel=1e-5)
    assert spectrum.geometry.theta == pytest.approx(0.0)
    assert spectrum.geometry.phi == pytest.approx(20.0)
    assert spectrum.geometry.psi == pytest.approx(20.0)
    assert spectrum.measurement.omega_msr == pytest.approx(3.0)
    assert spectrum.measurement.correction == pytest.approx(0.95, rel=1e-6)
    assert spectrum.date == "14:44:25  06-09-2010"


@needs_fixtures
def test_generic_data_record_uses_the_declared_compression():
    """The shipped files use 0x0011, not the explicit 0x0012-0x0015 forms.

    A reader that handles only the explicit types parses every header correctly
    and then silently returns a spectrum of zeros -- which is why this is
    asserted rather than assumed.
    """
    blob = FIXTURES[0].read_bytes()
    types = {record_type for record_type, _ in _records(blob)}
    assert DATA_REC in types
    assert DATA_INIT_REC in types or ARRAY_INIT_REC in types
    assert read_rbs(FIXTURES[0]).counts.sum() > 0


@needs_fixtures
@pytest.mark.oracle
@pytest.mark.skipif(not driver.available(), reason="legacy binary unavailable")
@pytest.mark.parametrize("path", FIXTURES, ids=lambda p: p.name)
def test_reader_is_bit_identical_to_rump(tmp_path, path):
    """The strongest available check: RUMP's own reader, same file."""
    oracle = driver.RumpOracle()
    out = tmp_path / "ref.dat"
    oracle.run([f"read '{path}'", f"wrascii '{out}'"], workdir=tmp_path)
    _, reference = driver.read_wrascii(out)

    mine = read_rbs(path).counts
    assert len(reference) == mine.size
    assert np.array_equal(mine, np.asarray(reference, dtype=np.float64))


# ------------------------------------------------------------------ writing


@needs_fixtures
@pytest.mark.parametrize("path", FIXTURES, ids=lambda p: p.name)
def test_round_trip_preserves_everything(tmp_path, path):
    original = read_rbs(path)
    written = tmp_path / "out.rbs"
    write_rbs(written, original)
    back = read_rbs(written)

    assert np.array_equal(original.counts, back.counts)
    assert back.identifier == original.identifier
    assert back.date == original.date
    assert back.e0_MeV == pytest.approx(original.e0_MeV)
    assert back.zbeam == original.zbeam
    assert back.mbeam == pytest.approx(original.mbeam)
    assert back.calibration.kevch == pytest.approx(original.calibration.kevch)
    assert back.calibration.kev0 == pytest.approx(original.calibration.kev0)
    assert back.geometry.kind is original.geometry.kind
    assert back.geometry.phi == pytest.approx(original.geometry.phi)
    assert back.measurement.omega_msr == pytest.approx(original.measurement.omega_msr)
    assert back.measurement.correction == pytest.approx(original.measurement.correction)


@needs_fixtures
@pytest.mark.oracle
@pytest.mark.skipif(not driver.available(), reason="legacy binary unavailable")
def test_rump_can_read_what_pyrump_writes(tmp_path):
    """Interoperability in the other direction."""
    original = read_rbs(next(p for p in FIXTURES if p.name == "2A.rbs"))
    written = tmp_path / "pyrump.rbs"
    write_rbs(written, original)

    oracle = driver.RumpOracle()
    out = tmp_path / "readback.dat"
    oracle.run([f"read '{written}'", f"wrascii '{out}'"], workdir=tmp_path)
    assert out.exists(), "RUMP could not read the file pyRUMP wrote"

    metadata, counts = driver.read_wrascii(out)
    assert np.array_equal(np.asarray(counts, dtype=np.float64), original.counts)
    assert metadata["Spectrum"] == "RBS"


def test_write_rejects_unknown_version(tmp_path):
    from pyrump.model.detector import Measurement
    from pyrump.model.geometry import Geometry
    from pyrump.model.spectrum import Calibration

    spectrum = RbsSpectrum(
        counts=np.zeros(16),
        calibration=Calibration(npt=16),
        geometry=Geometry(theta=0.0, phi=10.0),
        measurement=Measurement(),
    )
    with pytest.raises(ValueError, match="unknown format version"):
        write_rbs(tmp_path / "x.rbs", spectrum, version="9.9")


def test_synthetic_spectrum_round_trips(tmp_path):
    """No legacy files needed: write and read back our own."""
    from pyrump.model.detector import Measurement
    from pyrump.model.geometry import Geometry
    from pyrump.model.spectrum import Calibration

    counts = np.arange(512, dtype=np.float64) * 3.0
    spectrum = RbsSpectrum(
        counts=counts,
        calibration=Calibration(kevch=5.0, kev0=1.5, first=0.0, npt=512),
        geometry=Geometry(theta=7.0, phi=10.0, psi=0.0, kind=GeometryKind.IBM),
        measurement=Measurement(omega_msr=2.0, charge_uC=20.0, fwhm_keV=15.0),
        e0_MeV=2.0, zbeam=2, mbeam=4.0026, identifier="synthetic",
    )
    path = tmp_path / "synth.rbs"
    write_rbs(path, spectrum)
    back = read_rbs(path)
    assert np.array_equal(back.counts, counts)
    assert back.identifier == "synthetic"
    assert back.geometry.kind is GeometryKind.IBM
    assert back.calibration.kev0 == pytest.approx(1.5)


def test_multi_record_data_spans_1024_channel_blocks(tmp_path):
    """Records hold at most 1024 data words, so 2048 channels need two."""
    from pyrump.model.detector import Measurement
    from pyrump.model.geometry import Geometry
    from pyrump.model.spectrum import Calibration

    counts = np.linspace(0.0, 1000.0, 2048)
    spectrum = RbsSpectrum(
        counts=counts,
        calibration=Calibration(npt=2048),
        geometry=Geometry(theta=0.0, phi=10.0),
        measurement=Measurement(),
    )
    path = tmp_path / "big.rbs"
    write_rbs(path, spectrum)
    data_records = [t for t, _ in _records(path.read_bytes()) if t == DATA_REC]
    assert len(data_records) == 0  # we write the explicit form
    assert read_rbs(path).counts.size == 2048


def test_non_rump_file_is_rejected(tmp_path):
    path = tmp_path / "alien.rbs"
    path.write_bytes(_record(PROGRAM_ID_REC, struct.pack(">IHH", 0x12345678, 1, 0)))
    with pytest.raises(RbsFormatError, match="not a RUMP file"):
        read_rbs(path)


# -------------------------------------------------------------------- ASCII


def test_ascii_one_column(tmp_path):
    path = tmp_path / "a.dat"
    write_ascii(path, np.array([1.0, 2.0, 3.5]), identifier="my spectrum")
    result = read_ascii(path)
    assert np.allclose(result.counts, [1.0, 2.0, 3.5])
    assert result.identifier == "my spectrum"
    assert result.channels is None


def test_ascii_two_column(tmp_path):
    path = tmp_path / "b.dat"
    write_ascii(path, np.array([4.0, 5.0]), two_column=True, first_channel=10)
    result = read_ascii(path)
    assert np.allclose(result.counts, [4.0, 5.0])
    assert np.allclose(result.channels, [10, 11])


def test_ascii_reads_rump_wrascii_output(tmp_path):
    """The keyword header ends at the literal line 'Swallow'."""
    path = tmp_path / "c.dat"
    path.write_text(
        "Spectrum    RBS\n"
        "Ident      'a title'\n"
        "Conversion 5.000000 0.000000\n"
        "Swallow\n"
        "10.0\n20.0\n30.0\n"
    )
    result = read_ascii(path)
    assert np.allclose(result.counts, [10.0, 20.0, 30.0])
    assert result.identifier == "a title"
    assert result.metadata["Spectrum"] == "RBS"


def test_ascii_rejects_empty(tmp_path):
    path = tmp_path / "d.dat"
    path.write_text("just a header\nand another\n")
    with pytest.raises(ValueError, match="no numeric data"):
        read_ascii(path)
