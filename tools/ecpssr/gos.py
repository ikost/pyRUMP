"""
Hydrogenic generalized form factors for bound -> continuum transitions,
computed numerically (atomic units, nuclear charge Z = 1; for a screened
charge Z_s lengths scale as a0/Z_s, momenta as Z_s/a0, energies as Z_s^2 Hartree).

  G_nl(k^2, q) = (1/(2l+1)) sum_{m,m',l'} |<k l' m'| exp(i q.r) |n l m>|^2
               = sum_{l', lam} (2 lam+1)(2 l'+1) (l' lam l; 0 0 0)^2 N_{l'}^2 O_{lam}(k^2, q)^2

  O_lam = int u_{k l'}(r) j_lam(q r) R_nl(r) r dr, with u the regular Coulomb solution
  normalised as u = r^(l'+1)(1 + ...) at the origin, and N^2 the energy-normalisation
  (per Hartree) factor, which is analytic in k^2:

  N_l^2 = 2^(2l+2) prod_{s=1..l} (1 + s^2 k^2) / ((2l+1)!)^2 * 1/(1 - exp(-2 pi / k)).

For k^2 < 0 (energy transfer below the hydrogenic threshold, needed for the
Merzbacher-Lewis "outer screening" with theta < 1) the expression is continued
analytically in k^2 and the factor 1/(1 - exp(-2 pi/k)) is replaced by 1.
This reproduces exactly the tabulated Khandelwal-Choi-Merzbacher / Brandt-Lapicki
PWBA (checked against the Bethe coefficient C1(theta) for the K shell).
G is per electron.
"""
import numpy as np
from math import lgamma, factorial
from scipy.special import spherical_jn


def bound_radial(n, l, r):
    """Hydrogenic bound radial function R_nl(r) for Z = 1 (atomic units)."""
    if (n, l) == (1, 0):
        return 2.0 * np.exp(-r)
    if (n, l) == (2, 0):
        return (1.0 / np.sqrt(2.0)) * (1.0 - r / 2.0) * np.exp(-r / 2.0)
    if (n, l) == (2, 1):
        return (1.0 / (2.0 * np.sqrt(6.0))) * r * np.exp(-r / 2.0)
    if (n, l) == (3, 0):
        return 2.0 / (3.0 * np.sqrt(3.0)) * (1 - 2 * r / 3 + 2 * r**2 / 27) * np.exp(-r / 3)
    if (n, l) == (3, 1):
        return 8.0 / (27.0 * np.sqrt(6.0)) * r * (1 - r / 6) * np.exp(-r / 3)
    if (n, l) == (3, 2):
        return 4.0 / (81.0 * np.sqrt(30.0)) * r**2 * np.exp(-r / 3)
    raise ValueError((n, l))


def norm2(ls, k2s):
    """Energy-normalisation factor N_l^2 (analytic continuation for k^2 < 0)."""
    k2s = np.asarray(k2s, float)
    res = np.empty((len(ls), len(k2s)))
    k = np.sqrt(np.abs(k2s))
    fac = np.ones_like(k2s)
    pos = k2s > 0
    fac[pos] = 1.0 / (-np.expm1(-2 * np.pi / k[pos]))
    for i, l in enumerate(ls):
        p = np.ones_like(k2s)
        for s_ in range(1, l + 1):
            p = p * (1.0 + s_ * s_ * k2s)
        res[i] = np.exp((2 * l + 2) * np.log(2.0) - 2 * lgamma(2 * l + 2)) * p * fac
    return res


def coulomb_series_solution(ls, k2s, r):
    """Regular solution of u'' = [l(l+1)/r^2 - 2/r - k^2] u with u = r^(l+1)(1 + ...)
    at the origin, on the uniform grid r (r[0] = h).  Returns [len(ls), len(k2s), len(r)].
    Power series near the origin (where Numerov is inaccurate for large l), Numerov beyond."""
    h = r[1] - r[0]
    L = np.asarray(ls, float)[:, None]
    K2 = np.asarray(k2s, float)[None, :]
    nL, nK, nR = len(ls), K2.shape[1], len(r)
    out = np.empty((nL, nK, nR))

    def f(ri):
        return L * (L + 1) / ri**2 - 2.0 / ri - K2

    rs = np.maximum(2 * h, 12.0 * h * np.sqrt(L * (L + 1)))
    nser = int(np.ceil(rs.max() / h)) + 2
    rser = r[:nser]
    J = 30
    a = [np.ones((nL, nK)), (-2.0 / (2 * L + 2)) * np.ones((nL, nK))]
    for j in range(2, J):
        a.append((-2.0 * a[j - 1] - K2 * a[j - 2]) / (j * (j + 2 * L + 1)))
    ser = np.zeros((nL, nK, nser))
    pw = np.ones(nser)
    for j in range(J):
        ser += a[j][:, :, None] * pw[None, None, :]
        pw = pw * rser
    ser *= np.exp((L + 1) * np.log(rser)[None, :])[:, None, :]
    use_ser = (rser[None, :] <= rs)

    u0 = ser[:, :, 0].copy()
    u1 = ser[:, :, 1].copy()
    out[:, :, 0] = u0
    out[:, :, 1] = u1
    h12 = h * h / 12.0
    w0 = 1 - h12 * f(r[0])
    w1 = 1 - h12 * f(r[1])
    for i in range(2, nR):
        w2 = 1 - h12 * f(r[i])
        u2 = ((12 - 10 * w1) * u1 - w0 * u0) / w2
        if i < nser:
            m = use_ser[:, i][:, None]
            u2 = np.where(m, ser[:, :, i], u2)
        out[:, :, i] = u2
        u0, u1, w0, w1 = u1, u2, w1, w2
    return out


def threej_sq(lp, lam, l):
    """(l' lam l; 0 0 0)^2."""
    J = lp + lam + l
    if J % 2 == 1 or lam < abs(lp - l) or lam > lp + l:
        return 0.0
    g = J // 2
    num = factorial(J - 2 * lp) * factorial(J - 2 * lam) * factorial(J - 2 * l)
    den = factorial(J + 1)
    return num / den * (factorial(g) / (factorial(g - lp) * factorial(g - lam) * factorial(g - l)))**2


def gos_table(n, l, k2s, qs, lmax=45, rmax=None, h=0.002, chunk=4, verbose=False):
    """Numerical G_nl(k^2, q) (per electron) on grids k2s (signed) and qs."""
    k2s = np.asarray(k2s, float)
    if rmax is None:
        rmax = {1: 40.0, 2: 80.0, 3: 120.0}[n]
    r = np.arange(1, int(round(rmax / h)) + 1) * h
    Rb = bound_radial(n, l, r)
    w = np.ones_like(r)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    w *= h / 3.0
    base = Rb * r * w
    G = np.zeros((len(k2s), len(qs)))
    qr = np.outer(qs, r)
    jcache = {}
    for c0 in range(0, lmax + 1, chunk):
        lchunk = list(range(c0, min(c0 + chunk, lmax + 1)))
        U = coulomb_series_solution(lchunk, k2s, r)
        N2 = norm2(lchunk, k2s)
        for lam_old in [x for x in jcache if x < lchunk[0] - l]:
            del jcache[lam_old]
        for i, lp in enumerate(lchunk):
            for lam in range(abs(lp - l), lp + l + 1):
                t = threej_sq(lp, lam, l)
                if t == 0.0:
                    continue
                if lam not in jcache:
                    jcache[lam] = spherical_jn(lam, qr) * base[None, :]
                O = U[i] @ jcache[lam].T
                G += (2 * lam + 1) * (2 * lp + 1) * t * N2[i][:, None] * O**2
        if verbose:
            print("  l' <=", lchunk[-1], flush=True)
    return G


def gos_1s_analytic(k2, q):
    """Closed form (Walske / Merzbacher-Lewis) for 1s, per electron, same
    normalisation as gos_table:  G = 2 Q A(k^2, Q),  Q = q^2,
      A = 2^7 [Q + (k^2+1)/3] exp(-(2/k) atan(2k/(Q-k^2+1)))
          / ([(Q-k^2+1)^2 + 4k^2]^3 (1 - exp(-2 pi/k)))
    continued analytically for k^2 = -kappa^2 < 0:
      A = 2^7 [Q + (1-kappa^2)/3] exp(-(2/kappa) artanh(2 kappa/(Q+1+kappa^2)))
          / ([(Q+1+kappa^2)^2 - 4 kappa^2]^3)."""
    k2 = np.asarray(k2, float)[:, None]
    Q = np.asarray(q, float)[None, :] ** 2
    A = np.empty(np.broadcast(k2, Q).shape)
    pos = (k2 > 0)[:, 0]
    if pos.any():
        kk2 = k2[pos]
        k = np.sqrt(kk2)
        a = Q - kk2 + 1.0
        ang = np.arctan2(2 * k, a)
        A[pos] = (2**7 * (Q + (kk2 + 1) / 3.0) * np.exp(-(2.0 / k) * ang)
                  / ((a**2 + 4 * kk2) ** 3 * (-np.expm1(-2 * np.pi / k))))
    neg = ~pos
    if neg.any():
        ka2 = -k2[neg]
        ka = np.sqrt(ka2)
        b = Q + 1.0 + ka2
        A[neg] = (2**7 * (Q + (1 - ka2) / 3.0) * np.exp(-(2.0 / ka) * np.arctanh(2 * ka / b))
                  / ((b**2 - 4 * ka2) ** 3))
    return 2.0 * Q * A
