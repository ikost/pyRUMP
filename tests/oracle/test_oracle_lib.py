"""M0.3 acceptance: the unit-level physics oracle is callable and sane.

Like test_oracle_smoke.py these exercise the C, not pyRUMP. They establish that
the numbers later milestones will be compared against are trustworthy.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
import oracle as ora  # noqa: E402

pytestmark = [
    pytest.mark.oracle,
    pytest.mark.skipif(
        not ora.available() or ora.data_dir() is None,
        reason="oracle library not built (run tests/oracle/build_oracle.py)",
    ),
]


@pytest.fixture(scope="module")
def oracle() -> ora.Oracle:
    return ora.Oracle.load()


def test_tables_loaded(oracle):
    assert oracle.num_elements == 93
    assert oracle.element_symbol(14) == "Si"
    assert oracle.element_mass(14) == pytest.approx(28.086, rel=1e-5)
    # STOP_SQRT is what ships (stopping.c:73); the whole M3 fit depends on it.
    assert oracle.stop_type == ora.STOP_SQRT


def test_ziegler_block_matches_the_data_file(oracle):
    params = oracle.ziegler_params(14)
    assert params[0] == 28  # most abundant mass number
    assert params[1] == pytest.approx(27.977, rel=1e-4)
    assert params[2] == pytest.approx(28.086, rel=1e-4)
    assert params[3] == pytest.approx(2.3212, rel=1e-4)  # g/cm^3
    assert params[4] == pytest.approx(4.977e22, rel=1e-3)  # scaled on load


@pytest.mark.parametrize(
    "z1, m1, z2, energy, expected",
    [
        # Literature values for electronic stopping in eV/(1e15 at/cm^2).
        (2, 4.0026, 14, 2000.0, 48.9),  # 2 MeV He in Si
        (1, 1.008, 14, 2000.0, 5.33),  # 2 MeV H in Si
    ],
)
def test_zstop_is_physically_correct(oracle, z1, m1, z2, energy, expected):
    se, sn = oracle.zstop(z1, m1, z2, energy)
    assert se == pytest.approx(expected, rel=0.02)
    assert 0 < sn < se, "nuclear stopping is small but positive at MeV energies"


def test_zstop_rejects_out_of_range_z(oracle):
    """Ziegler covers Z=1..92; RUMP returns finite dummies rather than failing."""
    oracle.clear_log()
    se, sn = oracle.zstop(2, 4.0026, 93, 2000.0)
    assert (se, sn) == (40.0, 5.0)
    assert "not in range" in oracle.log()


def test_stopping_polynomial_is_the_simulation_path(oracle):
    """M3's crux: RUMP simulates with fitted coefficients, not raw Ziegler."""
    coef, scale = oracle.stopping_coefficients(2, 4.0026, 2.0, 14)
    assert len(coef) == ora.NDEG == 6
    assert scale == pytest.approx(1.0), "no energy rescaling for 4He"

    energies = [500.0, 1000.0, 2000.0, 2200.0]
    fitted = oracle.stopping(2, 4.0026, 2.0, 14, energies)
    assert all(s > 0 for s in fitted)
    # He in Si peaks near 500 keV and falls thereafter.
    assert fitted[0] > fitted[-1]


def test_fit_window_follows_the_beam_energy(oracle):
    """emin/emax scale with E_beam (stopping.c:316-319), so they are not constants.

    The port must derive the window the same way: it changes with the beam and
    therefore shifts the whole depth scale.
    """
    window = oracle.stopping_range(2, 4.0026, 2.0)
    assert window["type"] == ora.STOP_SQRT
    assert window["emin"] == pytest.approx(0.04 * 2.0, rel=1e-5)
    assert window["emax"] == pytest.approx(1.15 * 2.0, rel=1e-5)
    assert window["cutoff"] == pytest.approx(0.03 * 2.0, rel=1e-5)


def test_polynomial_tracks_the_model_inside_the_window_only(oracle):
    """Inside [emin, emax] the fit is close; outside it diverges fast.

    This is a trap for the port: evaluating the polynomial beyond emax silently
    returns a plausible-looking but wrong number, as RUMP does.
    """
    beam_MeV = 2.0
    window = oracle.stopping_range(2, 4.0026, beam_MeV)
    emax_keV = window["emax"] * 1000.0

    # Note: He on Si resolves to the Konac/Kalbitzer table, which takes priority
    # over Ziegler (stopping.c:481 before :483), so agreement with zstop is only
    # approximate even inside the window.
    for energy in (500.0, 1000.0, 2000.0):
        assert energy < emax_keV
        fitted = oracle.stopping(2, 4.0026, beam_MeV, 14, [energy])[0]
        raw, _ = oracle.zstop(2, 4.0026, 14, energy)
        assert fitted == pytest.approx(raw, rel=0.15)

    far_outside = 3000.0
    assert far_outside > emax_keV
    fitted = oracle.stopping(2, 4.0026, beam_MeV, 14, [far_outside])[0]
    raw, _ = oracle.zstop(2, 4.0026, 14, far_outside)
    assert abs(fitted - raw) / raw > 0.2, "extrapolation beyond emax should diverge"


def test_stopping_derivatives_are_consistent(oracle):
    """dS/dE from the analytic macros must match a numerical derivative."""
    energy = 2000.0
    step = 1.0
    _, ds, _ = oracle.stopping(2, 4.0026, 2.0, 14, [energy], derivatives=True)
    lo = oracle.stopping(2, 4.0026, 2.0, 14, [energy - step])[0]
    hi = oracle.stopping(2, 4.0026, 2.0, 14, [energy + step])[0]
    assert ds[0] == pytest.approx((hi - lo) / (2 * step), rel=1e-3)


def test_double_build_is_rejected_not_silently_wrong():
    """The -DREAL_IS_DOUBLE build loads corrupt tables; it must not be usable.

    RUMP's readers use scanf "%f" against REAL fields (ziegler.c:100,114;
    atomio.c:162-173), so doubling REAL writes 4 bytes into 8-byte fields.
    """
    if not ora.available("double"):
        pytest.skip("double build not present")
    with pytest.raises(RuntimeError, match="corrupt tables"):
        ora.Oracle("double")
