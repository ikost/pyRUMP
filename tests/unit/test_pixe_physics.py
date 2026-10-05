"""PIXE physics: cross sections, the vacancy cascade, line yields and the
detector response.

Expected values come from closed forms built out of the same tables the
code uses (the thin-target formula, Beer-Lambert attenuation, the cascade
written out by hand), not from numbers remembered from the literature.
"""

from __future__ import annotations

import numpy as np
import pytest

from pyrump.model.spectrum import Calibration
from pyrump.pixe import yields as yields_module
from pyrump.pixe.atomic import atomic_data
from pyrump.pixe.detector import PixeDetector
from pyrump.pixe.production import line_production
from pyrump.pixe.spectrum import resolution_sigma_keV, synthesize
from pyrump.pixe.xsect import ion_name, ionisation_table
from pyrump.pixe.yields import (
    AVOGADRO,
    BARN_CM2,
    Exposure,
    efficiency,
    element_density_g_cm3,
    simulate_lines,
)
from pyrump.sim.engine import Beam, UniformSample

from conftest import data_dir

DATA = data_dir()
needs_data = pytest.mark.skipif(DATA is None, reason="legacy data tables unavailable")
PROTON = Beam(e0_MeV=2.0, z=1, mass=1.00728)
HELIUM = Beam(e0_MeV=1.9, z=2, mass=4.0026)


@pytest.fixture(scope="module")
def tables():
    if DATA is None:
        pytest.skip("legacy data tables unavailable")
    from pyrump.shell.session import Session

    session = Session.create(str(DATA))
    return session.table, session.registry


def _film(table, layers: list[tuple[str, float]], substrate: str = "Si") -> UniformSample:
    """Single-element layers, thicknesses in 1e15 at/cm^2, on a substrate."""
    symbols = list(dict.fromkeys([s for s, _ in layers] + [substrate]))
    rows = [[1.0 if s == symbol else 0.0 for s in symbols] for symbol, _ in layers]
    rows.append([1.0 if s == substrate else 0.0 for s in symbols])
    return UniformSample(
        thicknesses=[t for _, t in layers] + [1e5],
        element_z=[table.by_symbol(s).z for s in symbols],
        compositions=rows,
    )


def _counts(lines, symbol: str, iupac: str) -> float:
    return next(line.counts for line in lines if line.symbol == symbol and line.line.line == iupac)


# -- cross sections -----------------------------------------------------------


def test_cross_sections_hit_the_table_and_refuse_to_extrapolate():
    table = ionisation_table("He4")
    log_e, log_s, _ = table._curves[26, "K"]
    np.testing.assert_allclose(table.sigma(26, "K", np.exp(log_e)), np.exp(log_s), rtol=1e-12)
    with pytest.raises(ValueError, match="outside"):
        table.sigma(26, "K", 20.0)
    assert table.sigma(26, "K", [1.0, 2.0])[1] > table.sigma(26, "K", 1.0)


def test_only_protons_and_helium_have_tables():
    assert ion_name(1, 1.00728) == "H1" and ion_name(2, 4.0026) == "He4"
    with pytest.raises(ValueError, match="only protons and 4He"):
        ion_name(2, 3.016)


def test_l_shell_cascade_written_out():
    """sigma(L3 line) = [s3 + s2 f23 + s1 (f13 + f12 f23)] omega3 F --
    Maxwell et al. 1989, eq. 7."""
    z, energy = 78, np.array([2.0])
    ion, atomic = ionisation_table("H1"), atomic_data()
    s1, s2, s3 = (ion.sigma(z, s, energy)[0] for s in ("L1", "L2", "L3"))
    f12, f13, f23 = (atomic.coster_kronig(z, t) for t in ("L12", "L13", "L23"))
    alpha1 = next(l for l in atomic.lines(z) if l.line == "L3M5")
    expected = (s3 + s2 * f23 + s1 * (f13 + f12 * f23)) * atomic.fluorescence_yield(z, "L3") * alpha1.rate
    produced = next(p for p in line_production(z, energy, ion, atomic) if p.line.line == "L3M5")
    assert produced.sigma_barn[0] == pytest.approx(expected, rel=1e-12)


# -- yields -------------------------------------------------------------------


@needs_data
def test_thin_film_limit(tables):
    """A very thin layer: no energy loss, no absorption -- the textbook
    thin-target yield."""
    table, registry = tables
    thickness = 0.5  # 1e15 at/cm^2
    detector = PixeDetector()
    exposure = Exposure(charge_uC=10.0)
    lines = simulate_lines(_film(table, [("Ti", thickness)]), PROTON, 0.0, detector,
                           exposure, registry, table)
    produced = next(p for p in line_production(22, np.array([PROTON.e0_MeV]),
                                               ionisation_table("H1"), atomic_data())
                    if p.line.line == "KL3")
    expected = (
        exposure.ions * detector.solid_angle_msr * 1e-3 / (4 * np.pi)
        * float(efficiency(detector, produced.line.energy_keV, table, atomic_data()))
        * thickness * 1e15 * produced.sigma_barn[0] * BARN_CM2
    )
    assert _counts(lines, "Ti", "KL3") == pytest.approx(expected, rel=2e-3)


@needs_data
def test_overlayer_attenuates_by_beer_lambert(tables):
    """Tilting the detector lengthens only the X-rays' way out through a Pt
    cap: the Ti Kα ratio between two angles is exp(-mu rho t (1/cos2 - 1/cos1)),
    the beam's path being unchanged."""
    table, registry = tables
    cap = 400.0
    sample = _film(table, [("Pt", cap), ("Ti", 0.5)])
    exposure = Exposure(charge_uC=10.0)
    counts = {
        angle: _counts(simulate_lines(sample, PROTON, 0.0, PixeDetector(angle_deg=angle),
                                      exposure, registry, table), "Ti", "KL3")
        for angle in (0.0, 60.0)
    }
    pt = table.by_symbol("Pt")
    mass_thickness = cap * 1e15 * pt.mass / AVOGADRO
    ti_kalpha = next(l.energy_keV for l in atomic_data().lines(22) if l.line == "KL3")
    mu = float(atomic_data().mu(78, ti_kalpha))
    expected = np.exp(-mu * mass_thickness * (1 / np.cos(np.radians(60.0)) - 1))
    assert counts[60.0] / counts[0.0] == pytest.approx(expected, rel=2e-3)
    assert expected < 0.9  # the test means something


@needs_data
def test_yield_is_linear_in_charge_solid_angle_h_and_live_time(tables):
    table, registry = tables
    sample = _film(table, [("Mn", 50.0)])

    def mn_kalpha(exposure, detector=PixeDetector()):
        return _counts(simulate_lines(sample, HELIUM, 9.0, detector, exposure, registry, table),
                       "Mn", "KL3")

    base = mn_kalpha(Exposure(charge_uC=10.0))
    assert mn_kalpha(Exposure(charge_uC=20.0)) == pytest.approx(2 * base, rel=1e-12)
    assert mn_kalpha(Exposure(charge_uC=10.0, live_fraction=0.5)) == pytest.approx(base / 2, rel=1e-12)
    assert mn_kalpha(Exposure(charge_uC=10.0, h=(3.0, 1.0, 1.0))) == pytest.approx(3 * base, rel=1e-12)
    assert mn_kalpha(Exposure(charge_uC=10.0), PixeDetector(solid_angle_msr=2.5)) == pytest.approx(
        2.5 * base, rel=1e-12)
    assert mn_kalpha(Exposure(charge_uC=10.0, correction=2.0)) == pytest.approx(base / 2, rel=1e-12)


@needs_data
def test_finer_sublayers_converge(tables, monkeypatch):
    table, registry = tables
    sample = _film(table, [("Ru", 40.0), ("Pt", 300.0), ("Mn", 300.0)])

    def run():
        lines = simulate_lines(sample, HELIUM, 9.0, PixeDetector(), Exposure(charge_uC=10.0),
                               registry, table)
        return {(l.symbol, l.line.line): l.counts for l in lines}

    coarse = run()
    monkeypatch.setattr(yields_module, "PIXE_MAXPTH", 10.0)
    fine = run()
    for key in (("Pt", "L3M5"), ("Mn", "KL3"), ("Ru", "L3M5")):
        assert fine[key] == pytest.approx(coarse[key], rel=5e-3), key


@needs_data
def test_substrate_and_absorber_layers_make_no_lines(tables):
    table, registry = tables
    sample = _film(table, [("Al", 1000.0), ("Ti", 50.0)], substrate="Cu")
    sample.absorber_layers = 1
    lines = simulate_lines(sample, PROTON, 0.0, PixeDetector(), Exposure(charge_uC=1.0),
                           registry, table)
    assert {line.symbol for line in lines} == {"Ti"}
    assert all(set(line.by_layer) == {1} for line in lines)


# -- the detector response ----------------------------------------------------


def test_resolution_reproduces_the_fwhm_at_mn_kalpha():
    detector = PixeDetector(fwhm_eV=131.0, fano=0.13)
    fwhm = 2 * np.sqrt(2 * np.log(2)) * resolution_sigma_keV(detector, 5.8988)
    assert fwhm * 1000 == pytest.approx(131.0, rel=1e-3)
    assert resolution_sigma_keV(detector, 1.74) < resolution_sigma_keV(detector, 9.44)


def _one_line(energy_keV: float, counts: float):
    from pyrump.pixe.atomic import XrayLine
    from pyrump.pixe.yields import LineYield

    line = XrayLine(25, "KL3", "K", energy_keV, 1.0)
    return LineYield(25, "Mn", line, counts, 1.0, 1.0, {0: counts})


def test_peaks_keep_their_area_and_centroid():
    calibration = Calibration(kevch=0.005, kev0=0.0, first=0.0, npt=4000)
    spectrum = synthesize([_one_line(5.8988, 1000.0)], calibration, PixeDetector(), escape=False)
    assert spectrum.counts.sum() == pytest.approx(1000.0, rel=1e-9)
    centres = calibration.edge_energy(np.arange(4000)) + calibration.kevch / 2
    assert (centres * spectrum.counts).sum() / 1000.0 == pytest.approx(5.8988, abs=1e-4)


def test_escape_peak_moves_area_1_74_kev_down():
    calibration = Calibration(kevch=0.005, kev0=0.0, first=0.0, npt=4000)
    detector = PixeDetector()
    with_escape = synthesize([_one_line(5.8988, 1000.0)], calibration, detector)
    assert with_escape.counts.sum() == pytest.approx(1000.0, rel=1e-9)
    escape_region = (calibration.edge_energy(np.arange(4000)) > 3.9) & (
        calibration.edge_energy(np.arange(4000)) < 4.4)
    fraction = with_escape.counts[escape_region].sum() / 1000.0
    assert 1e-3 < fraction < 0.05
    below_edge = synthesize([_one_line(1.74, 1000.0)], calibration, detector)
    assert below_edge.counts[calibration.edge_energy(np.arange(4000)) < 0.5].sum() < 1e-9


@needs_data
def test_efficiency_is_window_times_crystal(tables):
    table, _ = tables
    atomic = atomic_data()
    detector = PixeDetector()
    energy = 2.0
    be, si = table.by_symbol("Be"), table.by_symbol("Si")
    window = np.exp(-atomic.mu(4, energy) * element_density_g_cm3(be) * 8.0e-4)
    crystal = 1 - np.exp(-atomic.mu(14, energy) * element_density_g_cm3(si) * 450e-4)
    assert float(efficiency(detector, energy, table, atomic)) == pytest.approx(
        float(window * crystal), rel=1e-12)
