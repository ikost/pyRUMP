"""M7 acceptance: straggling and the Gaussian-broadened triangle fill.

pyRUMP does **not** port ``SimStragf``'s seven-regime rational approximation; it
uses a closed ``erf`` form instead. That is a deliberate departure from
faithfulness, so it is justified by measurement rather than assertion:

* the closed form is checked against numerical quadrature (machine precision);
* ``SimStragf`` is checked against the same quadrature, exposing its ~1.7e-6
  approximation error;
* the composed spectrum is checked against the C.

The one trap: ``SimStragf`` rescales its argument as ``x*(1+3*sig)`` before doing
anything (anlyz.c:387), so its ``x`` is in units of the broadened width, not the
triangle base. Overlooking that makes the function appear wrong by 0.23 out of a
total range of 0.5.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import integrate
from scipy.special import erf

from pyrump.atomic.tables import PeriodicTable
from pyrump.io.kalbitzer import parse_kalbitzer
from pyrump.model.detector import Measurement
from pyrump.model.geometry import Geometry
from pyrump.model.spectrum import Calibration
from pyrump.sim.engine import Beam, UniformSample, simulate
from pyrump.sim.fill.straggled import (
    add_triangle,
    stragf,
    triangle_gaussian_integral,
)
from pyrump.stopping.kalbitzer import KalbitzerStopping
from pyrump.stopping.registry import StoppingRegistry
from pyrump.stopping.ziegler import ZieglerStopping

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "oracle"))
import oracle as ora  # noqa: E402


from conftest import data_dir

DATA = data_dir()
CAL = Calibration(kevch=5.0, kev0=0.0, first=0.0, npt=1024)
SIGMAS = [0.05, 0.2, 0.5, 1.0, 3.0]


def _quadrature(x: float, s: float) -> float:
    """Independent reference: integral of (triangle * Gaussian) up to x."""
    integrand = lambda u: (1 - u) * 0.5 * (1 + erf((x - u) / s))  # noqa: E731
    value, _ = integrate.quad(integrand, 0.0, 1.0, limit=200)
    return value - 0.25


# ------------------------------------------------- the convolution integral


@pytest.mark.parametrize("s", SIGMAS)
def test_closed_form_matches_quadrature(s):
    """Machine precision -- this validates the derivation itself."""
    for x in (-2.0, -0.5, 0.0, 0.3, 0.7, 1.0, 1.5, 3.0):
        assert triangle_gaussian_integral(x, s)[0] == pytest.approx(
            _quadrature(x, s), abs=1e-12
        )


def test_zero_width_degenerates_to_the_bare_triangle():
    """anlyz.c:378-382's special case: no Gaussian, just the triangle integral."""
    assert triangle_gaussian_integral(-1.0, 0.0)[0] == -0.25
    assert triangle_gaussian_integral(0.0, 0.0)[0] == -0.25
    # int_0^x (1-t) dt - 1/4
    assert triangle_gaussian_integral(0.5, 0.0)[0] == pytest.approx(0.5 - 0.125 - 0.25)
    assert triangle_gaussian_integral(1.0, 0.0)[0] == pytest.approx(0.25)


@pytest.mark.parametrize("s", SIGMAS)
def test_limits_are_plus_and_minus_a_quarter(s):
    assert triangle_gaussian_integral(-50.0, s)[0] == pytest.approx(-0.25, abs=1e-12)
    assert triangle_gaussian_integral(50.0, s)[0] == pytest.approx(0.25, abs=1e-12)


@pytest.mark.parametrize("s", SIGMAS)
def test_integral_is_monotonic(s):
    """It is the integral of a non-negative function."""
    values = triangle_gaussian_integral(np.linspace(-5.0, 6.0, 400), s)
    # The flat tails sit at exactly +/-0.25, where float64 round-off in the
    # erf/exp combination is ~1e-15 -- larger than eps(0.25) alone.
    assert np.all(np.diff(values) >= -5e-15)


def test_broadening_smooths_the_transition():
    """Larger sigma spreads the rise over a wider range."""
    sharp = triangle_gaussian_integral(np.linspace(-2, 3, 200), 0.05)
    broad = triangle_gaussian_integral(np.linspace(-2, 3, 200), 2.0)
    assert np.max(np.diff(sharp)) > np.max(np.diff(broad))


# --------------------------------------------------- comparison against C

oracle_only = pytest.mark.skipif(
    not ora.available() or ora.data_dir() is None, reason="oracle unavailable"
)


@pytest.fixture
def oracle() -> ora.Oracle:
    handle = ora.Oracle.load()
    handle.reset_stopping_tables()
    return handle


@oracle_only
@pytest.mark.parametrize("s", SIGMAS)
def test_stragf_contract_matches_the_c(oracle, s):
    """Our ``stragf`` reproduces RUMP's, including the (1+3*sig) rescaling."""
    x = np.linspace(-1.5, 1.7, 25)
    ours = stragf(x, s)
    theirs = oracle.stragf(x, s)
    assert np.allclose(ours, theirs, atol=5e-6)


@oracle_only
def test_rump_approximation_error_is_quantified(oracle):
    """Record what replacing SimStragf actually costs -- or saves.

    The closed form is exact to machine precision; RUMP's rational fit carries a
    small but real error. This test documents the gap rather than hiding it.
    """
    worst_ours = worst_theirs = 0.0
    for s in SIGMAS:
        for x in np.linspace(-1.5, 1.7, 21):
            scaled = x * (1.0 + 3.0 * s)
            reference = _quadrature(scaled, s)
            worst_ours = max(worst_ours, abs(stragf(x, s)[0] - reference))
            worst_theirs = max(worst_theirs, abs(oracle.stragf(x, s)[0] - reference))

    assert worst_ours < 1e-12, "closed form should be exact"
    # RUMP's own claim is "rapidly and accurately"; measured, that is ~2e-6.
    assert 1e-8 < worst_theirs < 1e-4
    assert worst_theirs > worst_ours * 1e5


# ------------------------------------------------------ triangle deposition


def test_triangle_conserves_area():
    """A broadened triangle deposits height*width/2 -- it is a triangle.

    Broadening redistributes counts but must never create or destroy them, so
    the total is independent of sigma. (The two triangles of a split brick then
    sum to de*(h_front+h_back)/2, exactly the trapezoid area, which is why the
    straggled and unstraggled spectra have the same integral.)
    """
    expected = 0.5 * 3.0 * 40.0
    for sigma in (0.0, 0.5, 2.0, 10.0, 40.0):
        counts = np.zeros(CAL.npt)
        add_triangle(counts, CAL, height=3.0, energy=500.0, de=40.0, sigma=sigma)
        assert counts.sum() == pytest.approx(expected, rel=1e-9), f"sigma={sigma}"


def test_triangle_direction():
    """Negative de mirrors the triangle about its peak."""
    rising = np.zeros(CAL.npt)
    add_triangle(rising, CAL, 2.0, 500.0, 50.0, 1.0)
    falling = np.zeros(CAL.npt)
    add_triangle(falling, CAL, 2.0, 500.0, -50.0, 1.0)

    centroid = lambda c: float((np.arange(c.size) * c).sum() / c.sum())  # noqa: E731
    assert centroid(rising) > centroid(falling)
    assert rising.sum() == pytest.approx(falling.sum(), rel=1e-9)


def test_broadening_widens_the_deposit():
    narrow, wide = np.zeros(CAL.npt), np.zeros(CAL.npt)
    add_triangle(narrow, CAL, 2.0, 500.0, 40.0, 0.5)
    add_triangle(wide, CAL, 2.0, 500.0, 40.0, 20.0)
    assert np.count_nonzero(wide > 1e-9) > np.count_nonzero(narrow > 1e-9)
    assert wide.sum() == pytest.approx(narrow.sum(), rel=1e-9)


# ------------------------------------------------------- full spectrum


STRAGGLE_CASES = [
    ("Si straggle=1", UniformSample([1000.0], [14], [[1.0]], straggle=1.0)),
    ("Si thick", UniformSample([5000.0], [14], [[1.0]], straggle=1.0)),
    ("Au film", UniformSample([500.0], [79], [[1.0]], straggle=1.0)),
    ("Si straggle=3", UniformSample([1000.0], [14], [[1.0]], straggle=3.0)),
    ("SiO2", UniformSample([1000.0], [14, 8], [[1.0, 2.0]], straggle=1.0)),
]


@pytest.fixture(scope="module")
def table() -> PeriodicTable:
    assert DATA is not None
    return PeriodicTable.load(DATA / "atom4.dat", DATA / "pscoef.dat")


@pytest.fixture(scope="module")
def registry(table) -> StoppingRegistry:
    assert DATA is not None
    return StoppingRegistry(
        table.elements,
        kalbitzer=KalbitzerStopping(parse_kalbitzer(DATA / "newstop.kal"), table.elements),
        ziegler=ZieglerStopping(table.elements),
    )


@oracle_only
@pytest.mark.skipif(DATA is None, reason="legacy tables unavailable")
@pytest.mark.parametrize("label, sample", STRAGGLE_CASES, ids=[c[0] for c in STRAGGLE_CASES])
def test_straggled_spectrum_matches_oracle(oracle, registry, table, label, sample):
    geometry = Geometry(theta=0.0, phi=10.0)
    measurement = Measurement(omega_msr=1.0, charge_uC=10.0)
    mine = simulate(sample, Beam(), geometry, registry, table, CAL, measurement)

    oracle.set_beam(
        e0_MeV=2.0, phi=10.0, theta=0.0, kevch=CAL.kevch, kev0=CAL.kev0,
        first=CAL.first, npt=CAL.npt, fwhm=0.0, q_uC=10.0, omega=1.0, corr=1.0,
    )
    oracle.set_sample(
        sample.thicknesses, sample.element_z, sample.compositions,
        straggle=sample.straggle,
    )
    theirs = oracle.simulate_spectrum()

    n = min(len(theirs), len(mine.counts))
    ours, ref = mine.counts[:n], theirs[:n]
    assert ours.sum() == pytest.approx(ref.sum(), rel=1e-5)
    assert np.max(np.abs(ours - ref)) < 1e-4 * ref.max()


@pytest.mark.skipif(DATA is None, reason="legacy tables unavailable")
def test_straggling_broadens_the_edge(registry, table):
    """Physical check independent of the oracle: edges soften, counts survive."""
    geometry = Geometry(theta=0.0, phi=10.0)
    sharp = simulate(
        UniformSample([1000.0], [14], [[1.0]]), Beam(), geometry, registry, table, CAL
    )
    broad = simulate(
        UniformSample([1000.0], [14], [[1.0]], straggle=1.0),
        Beam(), geometry, registry, table, CAL,
    )
    assert np.count_nonzero(broad.counts) > np.count_nonzero(sharp.counts)
    assert broad.total() == pytest.approx(sharp.total(), rel=1e-3)
