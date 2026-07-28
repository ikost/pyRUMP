"""cffi binding to the RUMP physics oracle.

Loads ``libpyrump_oracle_float`` (see :mod:`build_oracle`) and exposes the legacy
routines as ordinary Python functions.

.. warning::
   **The ``double`` build does not work, and this is a defect in RUMP itself.**

   The original plan was to build the C twice -- once as shipped (``REAL`` is
   ``float``) and once with ``-DREAL_IS_DOUBLE`` -- and diff the two to measure
   RUMP's own single-precision rounding. That is not possible: every table reader
   hard-codes ``%f`` in its ``scanf`` format while writing into a ``REAL`` field
   (``ziegler.c:100,114``; ``atomio.c:162,163,173``). With ``REAL`` as ``double``
   those calls write four bytes into an eight-byte field, so the Ziegler and
   atomic tables load as denormal garbage and ``zstop`` silently returns 0.

   ``REAL_IS_DOUBLE`` is never set by any shipped makefile, so the bug was latent.
   Fixing it means patching the readers, which would make the oracle a *modified*
   RUMP rather than a reference.

   Consequence for pyRUMP: the ``float`` build is the authoritative oracle, and
   agreement tolerances must be set from single-precision reasoning (~1e-6
   relative on a single value, looser where cancellation occurs) rather than
   measured against a double-precision twin.
"""

from __future__ import annotations

import os
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BUILD = HERE / "build"

_CDEF = """
int    OracleInit(const char *data_dir);
int    OracleGetStopType(void);
void   OracleSetStopType(int type);

void   OracleZStop(int z1, double m1, int z2, double energy_keV,
                   int units, double *se, double *sn);

int    OracleStoppingCoefficients(int z_beam, double m_beam, double e_beam_MeV,
                                  int z_target, double *coef, double *e_scale);
int    OracleStoppingEvaluate(int z_beam, double m_beam, double e_beam_MeV,
                              int z_target, const double *energies, int n,
                              double *out_s, double *out_ds, double *out_dds);

int    OracleStoppingRange(int z_beam, double m_beam, double e_beam_MeV,
                           double *emin, double *emax, double *cutoff, int *type);

int    OracleCrossSection(int recoil, int z1, double m1, int z2, double m2,
                          double phi_deg, double kev_max,
                          const double *energies, int n, double *out);
int    OracleSigmaConstants(int recoil, int z1, double m1, int z2, double m2,
                            double phi_deg, double kev_max, double *out);

int    OracleSetSample(int n_layer, const double *thickness,
                       int n_element, const int *element_z,
                       const double *composition, const int *sublayers,
                       double maxpth, double straggle, double multiple);
int    OracleSetBeam(double e0_MeV, int zbeam, double mbeam, int cbeam,
                     double q_uC, double current_nA,
                     double kevch, double kev0, double first, int npt,
                     double fwhm, double tau, int geom,
                     double phi, double theta, double psi,
                     double omega, double corr);
int    OracleSetLayerEquation(int layer_index, int eqn_index,
                              const double *params, int n_param,
                              const double *species, int n_element);
int    OracleEquationSublayers(int eqn_index);
int    OracleSetAbsorber(int n_layers);
int    OracleSetLayerFuzz(int layer_index, double amount, int steps);
int    OracleSimulate(int capture_only);
int    OracleBrickCount(void);
int    OracleBrickOverflow(void);
void   OracleBricks(double *out, int n);
int    OracleSpectrum(double *out, int n);

double OracleStragf(double x, double sig);
double OracleNdtri(double p);

void   OracleResetStoppingTables(void);

int    OracleNumElements(void);
double OracleElementMass(int z);
double OracleElementDensity(int z);
double OracleRealMass(int z, int iso);
void   OracleElementSymbol(int z, char *out);
void   OracleZieglerParams(int z, double *out);

const char *OracleLog(void);
void   OracleLogClear(void);
"""

#: ``zstop`` unit codes (ziegler.c).
UNITS_EV_1E15_ATOMS = 1
UNITS_KEV_PER_MICRON = 2
UNITS_EV_PER_ANGSTROM = 3
UNITS_MEV_MG_CM2 = 4

#: ``STOPPING_TYPE`` values; RUMP ships with STOP_SQRT active (stopping.c:73).
STOP_LINEAR = 0
STOP_SQRT = 1

#: Degree-5 polynomial => 6 coefficients (``NDEG``).
NDEG = 6


def library_path(precision: str) -> Path:
    ext = ".dylib" if sys.platform == "darwin" else ".so"
    return BUILD / f"libpyrump_oracle_{precision}{ext}"


#: The only usable build; see the module docstring for why "double" is not.
DEFAULT_PRECISION = "float"


def available(precision: str = DEFAULT_PRECISION) -> bool:
    return library_path(precision).is_file()


def data_dir() -> Path | None:
    env = os.environ.get("PYRUMP_C_REFERENCE")
    roots = [Path(env)] if env else []
    roots.append(HERE.parents[1] / "C-code")
    for root in roots:
        if (root / "rump" / "data" / "atom4.dat").is_file():
            return root / "rump" / "data"
    return None


class Oracle:
    """Handle on one build of the physics library."""

    def __init__(self, precision: str = DEFAULT_PRECISION):
        import cffi

        lib_path = library_path(precision)
        if not lib_path.is_file():
            raise RuntimeError(
                f"{lib_path} not built. Run: python tests/oracle/build_oracle.py"
            )
        data = data_dir()
        if data is None:
            raise RuntimeError("legacy data tables not found")

        self.precision = precision
        self._ffi = cffi.FFI()
        self._ffi.cdef(_CDEF)
        self._lib = self._ffi.dlopen(str(lib_path))

        status = self._lib.OracleInit(str(data).encode())
        if status != 0:
            raise RuntimeError(
                f"OracleInit failed with status {status}: {self.log()}"
            )
        self._verify_tables()

    def _verify_tables(self) -> None:
        """Fail loudly if the tables loaded as garbage.

        The scanf/REAL width mismatch corrupts data silently rather than
        erroring, and a corrupt oracle would invalidate every test that trusts
        it. Silicon's Ziegler block is a cheap canary: 28.086 amu, 2.32 g/cc.
        """
        params = self.ziegler_params(14)
        mass_average, density = params[2], params[3]
        if not (28.0 < mass_average < 28.2 and 2.2 < density < 2.5):
            raise RuntimeError(
                f"{self.precision} oracle loaded corrupt tables "
                f"(Si mass={mass_average!r}, density={density!r}). "
                "This is expected for the 'double' build: RUMP's readers use "
                "scanf %f against REAL fields, so -DREAL_IS_DOUBLE writes 4 bytes "
                "into 8-byte doubles. Use precision='float'."
            )

    # ---------------------------------------------------------------- logging
    def log(self) -> str:
        return self._ffi.string(self._lib.OracleLog()).decode("utf-8", "replace")

    def clear_log(self) -> None:
        self._lib.OracleLogClear()

    # ------------------------------------------------------------- stop type
    @property
    def stop_type(self) -> int:
        return self._lib.OracleGetStopType()

    @stop_type.setter
    def stop_type(self, value: int) -> None:
        self._lib.OracleSetStopType(value)

    # ------------------------------------------------------- raw ZBL stopping
    def zstop(
        self,
        z1: int,
        m1: float,
        z2: int,
        energy_keV: float,
        units: int = UNITS_EV_1E15_ATOMS,
    ) -> tuple[float, float]:
        """Raw ZBL85 stopping, bypassing RUMP's polynomial refit.

        Returns ``(electronic, nuclear)``. Default units are
        eV/(1e15 atoms/cm^2), which is what the simulation kernel works in.
        """
        se = self._ffi.new("double *")
        sn = self._ffi.new("double *")
        self._lib.OracleZStop(z1, m1, z2, energy_keV, units, se, sn)
        return se[0], sn[0]

    # ------------------------------------------------ fitted polynomial layer
    def stopping_coefficients(
        self, z_beam: int, m_beam: float, e_beam_MeV: float, z_target: int
    ) -> tuple[np.ndarray, float]:
        """The degree-5 coefficients RUMP actually uses, plus the energy scale."""
        coef = self._ffi.new(f"double[{NDEG}]")
        scale = self._ffi.new("double *")
        ok = self._lib.OracleStoppingCoefficients(
            z_beam, m_beam, e_beam_MeV, z_target, coef, scale
        )
        if not ok:
            raise RuntimeError(
                f"no stopping table for Z1={z_beam} m={m_beam} Z2={z_target}: {self.log()}"
            )
        return np.frombuffer(self._ffi.buffer(coef), dtype=np.float64).copy(), scale[0]

    def stopping(
        self,
        z_beam: int,
        m_beam: float,
        e_beam_MeV: float,
        z_target: int,
        energies_keV,
        *,
        derivatives: bool = False,
    ):
        """Evaluate the fitted stopping polynomial as the simulation kernel does.

        With ``derivatives=True`` also returns dS/dE and d2S/dE2, which the
        3rd-order energy-loss expansion needs.
        """
        energies = np.ascontiguousarray(energies_keV, dtype=np.float64)
        n = energies.size
        out_s = np.empty(n, dtype=np.float64)
        out_ds = np.empty(n, dtype=np.float64) if derivatives else None
        out_dds = np.empty(n, dtype=np.float64) if derivatives else None

        def ptr(array):
            if array is None:
                return self._ffi.NULL
            return self._ffi.cast("double *", array.ctypes.data)

        ok = self._lib.OracleStoppingEvaluate(
            z_beam,
            m_beam,
            e_beam_MeV,
            z_target,
            self._ffi.cast("const double *", energies.ctypes.data),
            n,
            ptr(out_s),
            ptr(out_ds),
            ptr(out_dds),
        )
        if not ok:
            raise RuntimeError(
                f"no stopping table for Z1={z_beam} m={m_beam} Z2={z_target}: {self.log()}"
            )
        return (out_s, out_ds, out_dds) if derivatives else out_s

    def stopping_range(
        self, z_beam: int, m_beam: float, e_beam_MeV: float
    ) -> dict[str, float | int]:
        """The fit window RUMP derives from the beam energy, in MeV.

        The fitted polynomial is only valid inside ``[emin, emax]``; beyond
        ``emax`` it diverges from the underlying model within a few hundred keV.
        """
        emin = self._ffi.new("double *")
        emax = self._ffi.new("double *")
        cutoff = self._ffi.new("double *")
        kind = self._ffi.new("int *")
        ok = self._lib.OracleStoppingRange(
            z_beam, m_beam, e_beam_MeV, emin, emax, cutoff, kind
        )
        if not ok:
            raise RuntimeError(f"no stopping table for Z1={z_beam}: {self.log()}")
        return {
            "emin": emin[0],
            "emax": emax[0],
            "cutoff": cutoff[0],
            "type": kind[0],
        }

    # ---------------------------------------------------------- cross sections
    def cross_section(
        self,
        z1: int,
        m1: float,
        z2: int,
        m2: float,
        phi_deg: float,
        energies_keV,
        *,
        recoil: bool = False,
        kev_max: float | None = None,
    ) -> np.ndarray:
        """Cross-section in barns/sr.

        ``phi_deg`` is the **true** scattering angle, not RUMP's supplement:
        ``creatr.c:1629`` converts with ``sp.phi = 180 - samm->phi`` before
        calling in.
        """
        energies = np.ascontiguousarray(energies_keV, dtype=np.float64)
        out = np.empty(energies.size, dtype=np.float64)
        ok = self._lib.OracleCrossSection(
            1 if recoil else 0,
            z1,
            m1,
            z2,
            m2,
            phi_deg,
            float(energies.max()) if kev_max is None else kev_max,
            self._ffi.cast("const double *", energies.ctypes.data),
            energies.size,
            self._ffi.cast("double *", out.ctypes.data),
        )
        if not ok:
            raise RuntimeError(f"cross-section setup failed: {self.log()}")
        return out

    def sigma_constants(
        self,
        z1: int,
        m1: float,
        z2: int,
        m2: float,
        phi_deg: float,
        *,
        recoil: bool = False,
        kev_max: float = 3000.0,
    ) -> dict[str, float]:
        """The setup constants ``(csigma, csig_0, csig_f)``, for term-by-term checks."""
        out = self._ffi.new("double[3]")
        ok = self._lib.OracleSigmaConstants(
            1 if recoil else 0, z1, m1, z2, m2, phi_deg, kev_max, out
        )
        if not ok:
            raise RuntimeError(f"cross-section setup failed: {self.log()}")
        return {"csigma": out[0], "csig_0": out[1], "csig_f": out[2]}

    # --------------------------------------------------------- simulation
    def set_sample(
        self,
        thicknesses,
        element_z: list[int],
        compositions,
        *,
        sublayers: list[int] | None = None,
        maxpth: float = 200.0,
        straggle: float = 0.0,
        multiple: float = 0.0,
    ) -> None:
        """Define a sample of uniform layers, thicknesses in 1e15 at/cm^2."""
        thickness = np.ascontiguousarray(thicknesses, dtype=np.float64)
        zs = np.ascontiguousarray(element_z, dtype=np.int32)
        comp = np.ascontiguousarray(
            np.atleast_2d(compositions), dtype=np.float64
        )
        subs = np.ascontiguousarray(
            sublayers if sublayers is not None else [0] * thickness.size,
            dtype=np.int32,
        )
        ok = self._lib.OracleSetSample(
            thickness.size,
            self._ffi.cast("const double *", thickness.ctypes.data),
            zs.size,
            self._ffi.cast("const int *", zs.ctypes.data),
            self._ffi.cast("const double *", comp.ctypes.data),
            self._ffi.cast("const int *", subs.ctypes.data),
            maxpth,
            straggle,
            multiple,
        )
        if not ok:
            raise RuntimeError(f"OracleSetSample rejected the sample: {self.log()}")

    def set_beam(
        self,
        *,
        e0_MeV: float = 2.0,
        zbeam: int = 2,
        mbeam: float = 4.0026,
        cbeam: int = 1,
        q_uC: float = 10.0,
        current_nA: float = 0.0,
        kevch: float = 5.0,
        kev0: float = 0.0,
        first: float = 0.0,
        npt: int = 1024,
        fwhm: float = 0.0,
        tau: float = 5.0,
        geom: int = 0,
        phi: float = 10.0,
        theta: float = 0.0,
        psi: float = 0.0,
        omega: float = 1.0,
        corr: float = 1.0,
    ) -> None:
        """Configure the beam and detector.

        Note ``phi`` follows RUMP's convention: it is the *supplement* of the
        scattering angle, so 170 degrees is entered as 10.
        """
        self._lib.OracleSetBeam(
            e0_MeV, zbeam, mbeam, cbeam, q_uC, current_nA,
            kevch, kev0, first, npt, fwhm, tau, geom,
            phi, theta, psi, omega, corr,
        )

    #: probe_eqns order in sim_probe.c.
    EQUATIONS = [
        "none", "constant", "linear", "erfc", "exponential", "semi-infinite",
        "thinfilm", "buriedthinfilm", "thickfilm", "timedependent",
        "gaussian", "edgeworth",
    ]

    def set_layer_equation(
        self, layer_index: int, equation: str | None, params, species
    ) -> None:
        """Attach a depth-profile equation to a layer (or None to clear)."""
        index = -1 if equation is None else self.EQUATIONS.index(equation.lower())
        p = np.ascontiguousarray(params or [0.0], dtype=np.float64)
        sp = np.ascontiguousarray(species, dtype=np.float64)
        ok = self._lib.OracleSetLayerEquation(
            layer_index, index,
            self._ffi.cast("const double *", p.ctypes.data), p.size,
            self._ffi.cast("const double *", sp.ctypes.data), sp.size,
        )
        if not ok:
            raise RuntimeError(f"OracleSetLayerEquation rejected {equation!r}")

    def equation_sublayers(self, equation: str) -> int:
        """The equation's recommended sublayer count, which overrides maxpath."""
        return self._lib.OracleEquationSublayers(self.EQUATIONS.index(equation.lower()))

    def set_absorber(self, n_layers: int) -> None:
        """Mark the first n layers as absorber rather than sample."""
        if not self._lib.OracleSetAbsorber(n_layers):
            raise RuntimeError("OracleSetAbsorber rejected the layer count")

    def set_layer_fuzz(self, layer_index: int, amount: float, steps: int) -> None:
        """Attach lateral thickness roughness to a layer."""
        if not self._lib.OracleSetLayerFuzz(layer_index, amount, steps):
            raise RuntimeError("OracleSetLayerFuzz rejected the layer")

    def simulate_bricks(self) -> np.ndarray:
        """Run the engine and capture every brick it emits.

        Returns an ``(n, 9)`` array of
        ``(z, mass, efront, eback, hfront, hback, qqq, sigf, sigb)``.

        Works by swapping ``SimFillSpectrum`` -- a function pointer in the C --
        for a recorder, so `creatr.c` runs completely unmodified.
        """
        if not self._lib.OracleSimulate(1):
            raise RuntimeError(f"simulation failed: {self.log()}")
        if self._lib.OracleBrickOverflow():
            raise RuntimeError("brick capture buffer overflowed")
        count = self._lib.OracleBrickCount()
        out = np.empty((count, 9), dtype=np.float64)
        self._lib.OracleBricks(self._ffi.cast("double *", out.ctypes.data), count)
        return out

    def simulate_spectrum(self) -> np.ndarray:
        """Run the engine normally and return the finished channel spectrum."""
        if not self._lib.OracleSimulate(0):
            raise RuntimeError(f"simulation failed: {self.log()}")
        size = self._lib.OracleSpectrum(self._ffi.NULL, 0)
        out = np.empty(max(size, 1), dtype=np.float64)
        self._lib.OracleSpectrum(self._ffi.cast("double *", out.ctypes.data), out.size)
        return out[:size]

    def stragf(self, x, sig: float) -> np.ndarray:
        """RUMP's ``SimStragf``: integral of (unit triangle * Gaussian).

        The triangle is 1 at x=0 falling to 0 at x=1; ``sig`` uses the
        exp(-t^2/sig^2) convention, i.e. sqrt(2) times the true sigma. The
        result runs from -0.25 (far left) to +0.25 (far right).
        """
        values = np.atleast_1d(np.asarray(x, dtype=np.float64))
        return np.array([self._lib.OracleStragf(float(v), sig) for v in values])

    def ndtri(self, p) -> np.ndarray:
        """RUMP's inverse normal CDF, used to place FUZZ replicas."""
        values = np.atleast_1d(np.asarray(p, dtype=np.float64))
        return np.array([self._lib.OracleNdtri(float(v)) for v in values])

    def reset_stopping_tables(self) -> None:
        """Flush the C's session cache of fitted stopping tables.

        ``RbsStpfind`` reuses a cached table whenever the new beam energy fits
        inside an existing window (stopping.c:274-279), which makes results
        order-dependent. Call this before any comparison that assumes a fit was
        made at a specific beam energy.
        """
        self._lib.OracleResetStoppingTables()

    # ---------------------------------------------------------------- element
    @property
    def num_elements(self) -> int:
        return self._lib.OracleNumElements()

    def element_mass(self, z: int) -> float:
        return self._lib.OracleElementMass(z)

    def element_density(self, z: int) -> float:
        return self._lib.OracleElementDensity(z)

    def real_mass(self, z: int, mass_number: int = 0) -> float:
        return self._lib.OracleRealMass(z, mass_number)

    def element_symbol(self, z: int) -> str:
        buf = self._ffi.new("char[8]")
        self._lib.OracleElementSymbol(z, buf)
        return self._ffi.string(buf).decode()

    def ziegler_params(self, z: int) -> np.ndarray:
        out = self._ffi.new("double[7]")
        self._lib.OracleZieglerParams(z, out)
        return np.frombuffer(self._ffi.buffer(out), dtype=np.float64).copy()

    @staticmethod
    @lru_cache(maxsize=4)
    def load(precision: str = DEFAULT_PRECISION) -> "Oracle":
        """Cached loader -- the C keeps global table state, so reuse one handle."""
        return Oracle(precision)
