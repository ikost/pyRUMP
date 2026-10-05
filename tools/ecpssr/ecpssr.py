"""
ECPSSR ionization cross sections (Brandt & Lapicki) for K and L1-L3 subshells.

sigma_ECPSSR = C_s(x_C) * f_s(z) * sigma_PWBA(m_R * eta, zeta * theta)

Implementation notes (atomic units unless stated):
  Z_s   = Z2 - 0.3 (K), Z2 - 4.15 (L)         screened charge (Slater)
  theta = n^2 U_s / (Z_s^2 Ry)                  U_s = observed binding energy
  eta   = (v1 / Z_s)^2
  xi    = 2 n sqrt(eta) / theta
  zeta  = 1 + 2 Z1 / (Z_s theta) * (g_s(xi) - h_s(xi))      binding + polarization
  h_s   = 2 n I(c_s / xi) / (theta xi^3),  c = 1.5 n (K, L1), 1.25 n (L2, L3)
  m_R   = sqrt(1 + 1.1 y^2) + y,  y = 0.40 (Z_s/c)^2 / (n xi/zeta) (K, L1),
                                   y = 0.15 (Z_s/c)^2 / (xi/zeta)   (L2, L3)
  energy loss   z = sqrt(1 - 4 zeta / (M theta xi^2)),
                f(z) = 2^-p (p-1)^-1 [(pz - 1)(1+z)^p + (pz + 1)(1-z)^p],  p = 9 (K, L1), 11 (L2, L3)
  Coulomb       x_C = 2 pi d q0 zeta / (z (1+z)),  d q0 = Z1 Z2 U / (M v1^3)
                C = p E_{p+1}(x_C)
  M = reduced mass of projectile + target atom (electron masses)
  Binding energies U_s: xraylib EdgeEnergy.  zeta*theta must stay >= 0.3, the
  lower limit of the PWBA form-factor grid; below it pwba raises rather than
  clamping (the shipped tables never reach it).
"""
import numpy as np
from scipy.special import expn
from pwba import ReducedPWBA

HARTREE_EV = 27.211386245988
RY_EV = HARTREE_EV / 2
AMU_ME = 1822.888486209
A0SQ_BARN = 0.529177210903e-8 ** 2 * 1e24          # a0^2 in barn
C_AU = 137.035999084

SHELL_PARAMS = {
    #       n  screen  c_pol   p   y-coef  y uses n
    "K":  (1, 0.30, 1.5,  9, 0.40, True),
    "L1": (2, 4.15, 1.5,  9, 0.40, True),
    "L2": (2, 4.15, 1.25, 11, 0.15, False),
    "L3": (2, 4.15, 1.25, 11, 0.15, False),
}


def g_func(shell, xi):
    if shell == "K":
        num = 1 + 9*xi + 31*xi**2 + 98*xi**3 + 12*xi**4 + 25*xi**5 + 4.2*xi**6 + 0.515*xi**7
        return num / (1 + xi)**9
    if shell == "L1":
        num = 1 + 9*xi + 31*xi**2 + 49*xi**3 + 162*xi**4 + 63*xi**5 + 18*xi**6 + 1.97*xi**7
        return num / (1 + xi)**9
    num = (1 + 10*xi + 45*xi**2 + 102*xi**3 + 331*xi**4 + 6.7*xi**5 + 58*xi**6
           + 7.8*xi**7 + 0.888*xi**8)
    return num / (1 + xi)**10


def I_func(x):
    if x <= 0.035:
        return 0.75 * np.pi * (np.log(1.0 / x**2) - 1.0)
    if x <= 3.1:
        return np.exp(-2*x) / (0.031 + 0.213*x**0.5 + 0.005*x - 0.069*x**1.5 + 0.324*x**2)
    if x <= 11.0:
        return 2.0 * np.exp(-2*x) / x**1.6
    return 0.0


_PWBA = {}


def pwba(shell):
    if shell not in _PWBA:
        _PWBA[shell] = ReducedPWBA(shell)
    return _PWBA[shell]


def ecpssr(shell, Z1, A1, Z2, A2, E_MeV, U_keV, parts=False):
    """Ionization cross section (barn) of subshell `shell` of element Z2 by an
    ion (Z1, mass A1 in u) with lab kinetic energy E_MeV.  U_keV = binding energy."""
    n, scr, cpol, p, ycoef, yn = SHELL_PARAMS[shell]
    Zs = Z2 - scr
    U = U_keV * 1e3 / HARTREE_EV                 # Hartree
    M1 = A1 * AMU_ME
    M = M1 * (A2 * AMU_ME) / (M1 + A2 * AMU_ME)  # reduced mass (m_e)
    E = E_MeV * 1e6 / HARTREE_EV
    v1 = np.sqrt(2.0 * E / M1)
    eta = v1**2 / Zs**2
    theta = n**2 * U_keV * 1e3 / (Zs**2 * RY_EV)
    xi = 2.0 * n * np.sqrt(eta) / theta
    h = 2.0 * n * I_func(cpol * n / xi) / (theta * xi**3)
    g = g_func(shell, xi)
    zeta = 1.0 + 2.0 * Z1 / (Zs * theta) * (g - h)
    y = ycoef * (Zs / C_AU)**2 / ((n if yn else 1.0) * xi / zeta)
    mR = np.sqrt(1.0 + 1.1 * y**2) + y
    sigma0 = 8.0 * np.pi * Z1**2 / Zs**4 * A0SQ_BARN
    eta_R = mR * eta
    tz = zeta * theta                    # pwba raises below its theta >= 0.3 grid
    s_pssr = sigma0 / eta_R * pwba(shell).I(eta_R, tz)
    dl = 4.0 * zeta / (M * theta * xi**2)
    if dl >= 1.0:
        return 0.0
    z = np.sqrt(1.0 - dl)
    fz = 2.0**(-p) / (p - 1) * ((p*z - 1) * (1 + z)**p + (p*z + 1) * (1 - z)**p)
    dq0 = Z1 * Z2 * U / (M * v1**3)
    xC = 2.0 * np.pi * dq0 * zeta / (z * (1 + z))
    C = p * expn(p + 1, xC)
    sig = C * fz * s_pssr
    if parts:
        return dict(sigma=sig, pwba=sigma0 / eta * pwba(shell).I(eta, theta),
                    pssr=s_pssr, zeta=zeta, mR=mR, C=C, fz=fz, xi=xi, theta=theta, eta=eta)
    return sig
