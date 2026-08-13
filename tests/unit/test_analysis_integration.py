"""M15 acceptance: INTEGRAL/THICKNESS/INTSET's shared region-integration algorithm."""

from __future__ import annotations

import numpy as np
import pytest

from pyrump.analysis.integration import integrate_region
from pyrump.model.detector import Measurement
from pyrump.model.geometry import Geometry
from pyrump.model.spectrum import Calibration, Spectrum
from pyrump.shell.session import Buffer
from pyrump.sim.engine import Beam

from conftest import data_dir

DATA = data_dir()
needs_data = pytest.mark.skipif(DATA is None, reason="legacy data tables unavailable")


def _buffer(counts, kevch=5.0, e0_MeV=2.0):
    return Buffer(
        spectrum=Spectrum(counts=np.asarray(counts, dtype=np.float64),
                           calibration=Calibration(kevch=kevch, npt=len(counts))),
        beam=Beam(e0_MeV=e0_MeV, z=2, mass=4.0026),
        geometry=Geometry(theta=0.0, phi=10.0),
        measurement=Measurement(),  # omega=1, charge_uC=10, correction=1, charge_state=1
    )


# With the Measurement defaults above, NormK*kevch == charge_state*correction /
# (omega*charge_uC) == 1*1/(1*10) == 0.1, independent of kevch (it cancels).
_NORM = 0.1


def test_discrete_mode_matches_the_hand_computed_bugfix_formula():
    counts = [0, 0, 10, 20, 30, 20, 10, 0, 0, 0]
    buffer = _buffer(counts)
    result = integrate_region(buffer, 2, 6, interp=False)
    # gross = sum(counts[2:7]) = 90; net = gross - (counts[2]+counts[6])*(6+1-2)/2
    #       = 90 - (10+10)*5/2 = 40 -- the "+1" is anlytc.c's 1984 bugfix.
    assert result.gross == pytest.approx(90 * _NORM)
    assert result.net == pytest.approx(40 * _NORM)
    assert result.ist == 2
    assert result.iend == 6


def test_interpolated_mode_matches_the_hand_computed_triangular_correction():
    counts = [0, 0, 10, 20, 30, 20, 10, 0, 0, 0]
    buffer = _buffer(counts)
    result = integrate_region(buffer, 1.7, 6.3, interp=True)
    # ist=round(1.7+.5)=2, iend=round(6.3-.5)=6, dxst=dxend=0.3; boundary
    # corrections -2.45 each side -> gross=85.1; net=85.1-0.5*14*4.6=52.9.
    assert result.gross == pytest.approx(85.1 * _NORM)
    assert result.net == pytest.approx(52.9 * _NORM)


def test_reversed_bounds_are_reordered_first():
    counts = [0, 0, 10, 20, 30, 20, 10, 0, 0, 0]
    buffer = _buffer(counts)
    forward = integrate_region(buffer, 2, 6, interp=False)
    backward = integrate_region(buffer, 6, 2, interp=False)
    assert backward.gross == pytest.approx(forward.gross)
    assert backward.net == pytest.approx(forward.net)


def test_target_token_none_stops_after_gross_net():
    counts = [0, 0, 10, 20, 30, 20, 10, 0, 0, 0]
    buffer = _buffer(counts)
    result = integrate_region(buffer, 2, 6, interp=False)
    assert result.thickness is None


@needs_data
def test_thickness_surface_approximation():
    from pyrump.cli._common import load_tables

    table, registry, _ = load_tables(DATA)
    counts = [0.0] * 10
    counts[3:7] = [100.0, 200.0, 200.0, 100.0]
    buffer = _buffer(counts)
    result = integrate_region(
        buffer, 3, 6, interp=False, registry=registry, table=table,
        target_token="Si", qmode=0,
    )
    assert result.thickness is not None
    assert result.thickness.compensated is None
    assert result.thickness.gross_atoms_cm2 > 0
    assert result.thickness.density_g_cc == pytest.approx(2.32, abs=0.05)


@needs_data
def test_thickness_compensated_pass_needs_alpha_in_query_mode():
    from pyrump.cli._common import load_tables

    table, registry, _ = load_tables(DATA)
    counts = [0.0] * 10
    counts[3:7] = [100.0, 200.0, 200.0, 100.0]
    buffer = _buffer(counts)
    with pytest.raises(ValueError, match="alpha"):
        integrate_region(
            buffer, 3, 6, interp=False, registry=registry, table=table,
            target_token="Si", qmode=2,
        )
    result = integrate_region(
        buffer, 3, 6, interp=False, registry=registry, table=table,
        target_token="Si", qmode=2, alpha_override=0.85,
    )
    assert result.thickness.compensated is not None


@needs_data
def test_thickness_rejects_an_empty_region():
    from pyrump.cli._common import load_tables

    table, registry, _ = load_tables(DATA)
    counts = [100.0] * 10
    buffer = _buffer(counts)
    with pytest.raises(ValueError, match="no region"):
        integrate_region(
            buffer, 5, 5, interp=False, registry=registry, table=table, target_token="Si",
        )
