"""M8 acceptance: detector-resolution convolution.

The last stage of the forward model. Two things here are easy to get wrong and
both change results:

* the kernel is **channel-integrated**, not point-sampled -- weight *j* is the
  Gaussian's integral across channel *j*;
* the convolution **loses counts at both edges**, because contributions falling
  outside the channel range are simply discarded.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.special import ndtr

from pyrump.atomic.tables import PeriodicTable
from pyrump.io.kalbitzer import parse_kalbitzer
from pyrump.model.detector import Measurement
from pyrump.model.geometry import Geometry
from pyrump.model.spectrum import Calibration
from pyrump.sim.convolve import (
    FWHM_TO_SIGMA,
    LMAX,
    convolve_detector,
    full_kernel,
    gaussian_kernel,
)
from pyrump.sim.engine import Beam, UniformSample, simulate
from pyrump.stopping.kalbitzer import KalbitzerStopping
from pyrump.stopping.registry import StoppingRegistry
from pyrump.stopping.ziegler import ZieglerStopping

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "oracle"))
import oracle as ora  # noqa: E402


from conftest import data_dir

DATA = data_dir()
CAL = Calibration(kevch=5.0, kev0=0.0, first=0.0, npt=1024)


# ------------------------------------------------------------------- kernel


def test_kernel_is_channel_integrated():
    """Weights are CDF differences, not point samples of the Gaussian."""
    fwhm, kevch = 30.0, 5.0
    sigma = fwhm / FWHM_TO_SIGMA
    half = gaussian_kernel(fwhm, kevch)

    # Central channel spans +/- half a channel about zero.
    expected0 = ndtr(0.5 * kevch / sigma) - ndtr(-0.5 * kevch / sigma)
    assert half[0] == pytest.approx(expected0, rel=1e-12)

    # Channel j spans (j-0.5, j+0.5) channel widths out.
    for j in (1, 2, 3):
        expected = ndtr((j + 0.5) * kevch / sigma) - ndtr((j - 0.5) * kevch / sigma)
        assert half[j] == pytest.approx(expected, rel=1e-12)


def test_kernel_truncates_at_three_sigma():
    """lgauss = int(3*sigma/kevch + 2), so the tail beyond 3 sigma is dropped."""
    fwhm, kevch = 30.0, 5.0
    sigma = fwhm / FWHM_TO_SIGMA
    half = gaussian_kernel(fwhm, kevch)
    assert half.size == int(3.0 * sigma / kevch + 2)

    total = full_kernel(half).sum()
    # A 3-sigma Gaussian holds ~99.7%; the shortfall is the truncation.
    assert 0.99 < total < 1.0
    assert 1.0 - total == pytest.approx(8.49e-4, rel=0.1)


def test_kernel_is_clamped_to_lmax():
    """A fine calibration would want a huge kernel; the C caps it at 100."""
    half = gaussian_kernel(fwhm_keV=500.0, kevch=0.1)
    assert half.size == LMAX


def test_zero_resolution_is_a_no_op():
    counts = np.arange(100.0)
    assert np.array_equal(convolve_detector(counts, CAL, 0.0), counts)
    assert np.array_equal(convolve_detector(counts, CAL, -1.0), counts)


def test_kernel_is_symmetric_and_peaked():
    kernel = full_kernel(gaussian_kernel(30.0, 5.0))
    assert np.allclose(kernel, kernel[::-1])
    assert kernel.argmax() == kernel.size // 2


# ------------------------------------------------------------ edge behaviour


def test_counts_are_lost_at_both_edges():
    """RUMP discards contributions that fall outside the channel range.

    The head loop (creatr.c:1233-1236) drops the ``i-k`` half once ``k > i``,
    and the tail loop (creatr.c:1284-1291) drops the ``i+k`` half once
    ``k >= npt-i``. The behaviour is symmetric -- both ends leak.
    """
    calibration = Calibration(kevch=5.0, kev0=0.0, first=0.0, npt=200)
    truncation = 1.0 - full_kernel(gaussian_kernel(30.0, 5.0)).sum()

    middle = np.zeros(200)
    middle[100] = 1000.0
    interior = convolve_detector(middle, calibration, 30.0).sum()
    # Away from the edges only the 3-sigma truncation is lost.
    assert interior == pytest.approx(1000.0 * (1.0 - truncation), rel=1e-9)

    low, high = np.zeros(200), np.zeros(200)
    low[2] = 1000.0
    high[197] = 1000.0
    low_total = convolve_detector(low, calibration, 30.0).sum()
    high_total = convolve_detector(high, calibration, 30.0).sum()

    assert low_total < 0.9 * interior
    assert high_total < 0.9 * interior
    # Symmetric positions leak identically.
    assert low_total == pytest.approx(high_total, rel=1e-9)


def test_renormalize_mode_conserves_counts():
    calibration = Calibration(kevch=5.0, kev0=0.0, first=0.0, npt=200)
    for position in (2, 100, 197):
        counts = np.zeros(200)
        counts[position] = 1000.0
        renormalised = convolve_detector(
            counts, calibration, 30.0, mode="renormalize"
        )
        assert renormalised.sum() == pytest.approx(1000.0, rel=1e-9)


def test_unknown_mode_rejected():
    with pytest.raises(ValueError, match="unknown convolution mode"):
        convolve_detector(np.zeros(10), CAL, 10.0, mode="bogus")


# ------------------------------------------------------------ physical shape


def test_convolution_broadens_and_preserves_centroid():
    counts = np.zeros(400)
    counts[200] = 1000.0
    calibration = Calibration(kevch=5.0, kev0=0.0, first=0.0, npt=400)
    out = convolve_detector(counts, calibration, 40.0)

    assert np.count_nonzero(out) > 1
    centroid = float((np.arange(out.size) * out).sum() / out.sum())
    assert centroid == pytest.approx(200.0, abs=1e-9)


def test_wider_resolution_broadens_more():
    counts = np.zeros(400)
    counts[200] = 1000.0
    calibration = Calibration(kevch=5.0, kev0=0.0, first=0.0, npt=400)
    narrow = convolve_detector(counts, calibration, 15.0)
    wide = convolve_detector(counts, calibration, 60.0)
    assert np.count_nonzero(wide) > np.count_nonzero(narrow)
    assert wide.max() < narrow.max()


# --------------------------------------------------------- against the C

pytestmark = pytest.mark.skipif(
    DATA is None or not ora.available() or ora.data_dir() is None, reason="legacy tables or oracle unavailable"
)


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


@pytest.fixture
def oracle() -> ora.Oracle:
    handle = ora.Oracle.load()
    handle.reset_stopping_tables()
    return handle


CASES = [
    ("Si fwhm=15", UniformSample([1000.0], [14], [[1.0]]), 15.0),
    ("Si fwhm=30", UniformSample([1000.0], [14], [[1.0]]), 30.0),
    ("Au fwhm=15", UniformSample([500.0], [79], [[1.0]]), 15.0),
    ("SiO2 fwhm=20", UniformSample([1000.0], [14, 8], [[1.0, 2.0]]), 20.0),
    ("straggle + fwhm", UniformSample([1000.0], [14], [[1.0]], straggle=1.0), 15.0),
    ("thick Si fwhm=25", UniformSample([5000.0], [14], [[1.0]]), 25.0),
]


@pytest.mark.parametrize("label, sample, fwhm", CASES, ids=[c[0] for c in CASES])
def test_convolved_spectrum_matches_oracle(oracle, registry, table, label, sample, fwhm):
    geometry = Geometry(theta=0.0, phi=10.0)
    measurement = Measurement(omega_msr=1.0, charge_uC=10.0, fwhm_keV=fwhm)
    mine = simulate(sample, Beam(), geometry, registry, table, CAL, measurement)

    oracle.set_beam(
        e0_MeV=2.0, phi=10.0, theta=0.0, kevch=CAL.kevch, kev0=CAL.kev0,
        first=CAL.first, npt=CAL.npt, fwhm=fwhm, q_uC=10.0, omega=1.0, corr=1.0,
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


def test_resolution_softens_the_surface_edge(registry, table):
    """A physical check independent of the oracle."""
    geometry = Geometry(theta=0.0, phi=10.0)
    sample = UniformSample([1000.0], [14], [[1.0]])

    sharp = simulate(
        sample, Beam(), geometry, registry, table, CAL,
        Measurement(omega_msr=1.0, charge_uC=10.0, fwhm_keV=0.0),
    )
    blurred = simulate(
        sample, Beam(), geometry, registry, table, CAL,
        Measurement(omega_msr=1.0, charge_uC=10.0, fwhm_keV=30.0),
    )

    assert np.count_nonzero(blurred.counts) > np.count_nonzero(sharp.counts)
    # The edge is well inside the range, so almost nothing is lost.
    assert blurred.total() == pytest.approx(sharp.total(), rel=2e-3)
    # The steepest gradient must soften.
    assert np.max(np.abs(np.diff(blurred.counts))) < np.max(
        np.abs(np.diff(sharp.counts))
    )
