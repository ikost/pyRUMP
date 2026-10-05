"""
Hydrogenic PWBA reduced cross sections I_s(eta, theta) for the K (1s), L1 (2s)
and L2/L3 (2p) subshells in the Khandelwal-Choi-Merzbacher / Brandt-Lapicki
formulation (Merzbacher-Lewis outer screening).

  sigma_PWBA = 8 pi a0^2 Z1^2 / (eta Z_s^4) * I_s(eta, theta)

  I_s = N_e * Int_{eps_min}^{inf} d(eps) Int_{q0}^{inf} dq/q^3 G_nl(k^2, q)

Scaled hydrogenic units (charge Z_s, energies in Z_s^2 Hartree, momenta in Z_s/a0):
  eps     = energy transfer,  k^2/2 = eps - 1/(2 n^2)   (hydrogenic continuum energy)
  eps_min = theta / (2 n^2)   (observed binding energy; k^2 < 0 below the hydrogenic
                               threshold -> analytic continuation, see gos.py)
  q0      = eps / sqrt(eta),  eta = (v1/(Z_s v0))^2,  theta = n^2 U / (Z_s^2 Ry)

Universal function of Brandt-Lapicki / ISICS:  F_s(theta, eta/theta^2) = theta I_s / eta.
"""
import os
import numpy as np
from gos import gos_table, gos_1s_analytic

HERE = os.path.dirname(os.path.abspath(__file__))

SHELLS = {            # name: (n, l, electrons)
    "K": (1, 0, 2),
    "L1": (2, 0, 2),
    "L2": (2, 1, 2),
    "L3": (2, 1, 4),
    "M1": (3, 0, 2),
    "M2": (3, 1, 2),     # 3p1/2: 1/3 of the 3p shell
    "M3": (3, 1, 4),     # 3p3/2: 2/3
    "M4": (3, 2, 4),     # 3d3/2: 0.4 of the 3d shell
    "M5": (3, 2, 6),     # 3d5/2: 0.6
}
THETA_MIN = {1: 0.30, 2: 0.30, 3: 0.30}   # lowest reduced binding energy supported
Q_GRID = np.geomspace(1e-3, 80.0, 230)


def e_grid(n):
    """Continuum kinetic energies e = k^2/2 (scaled Hartree), signed.
    Linear spacing where the integrand is steep (near threshold), geometric above."""
    e_lo = (THETA_MIN[n] - 1.0) / (2 * n * n)
    if n == 1:
        neg = np.linspace(e_lo, 0.0, 141)[:-1]
        pos = np.concatenate([np.linspace(0.0, 2.0, 201), np.geomspace(2.0, 80.0, 70)[1:]])
    elif n == 2:
        neg = np.linspace(e_lo, 0.0, 89)[:-1]
        pos = np.concatenate([np.linspace(0.0, 1.0, 161), np.geomspace(1.0, 20.0, 50)[1:]])
    else:
        # n = 3 binds 1/18 Hartree (scaled): the steep part sits at lower
        # continuum energies than for n = 2.
        neg = np.linspace(e_lo, 0.0, 71)[:-1]
        pos = np.concatenate([np.linspace(0.0, 0.6, 151), np.geomspace(0.6, 20.0, 55)[1:]])
    return np.concatenate([neg, pos])


def gos_grid(n, l, cache=True, verbose=False):
    """G_nl (per electron) on (e_grid(n), Q_GRID); numeric results cached as .npz."""
    E = e_grid(n)
    k2 = 2.0 * E
    k2 = np.where(k2 == 0.0, 1e-12, k2)
    if (n, l) == (1, 0):
        return gos_1s_analytic(k2, Q_GRID)
    fn = os.path.join(HERE, f"gos_{n}{l}.npz")
    if cache and os.path.exists(fn):
        d = np.load(fn)
        if d["G"].shape == (len(E), len(Q_GRID)) and np.allclose(d["E"], E):
            return d["G"]
    neg = k2 < 0
    G = np.zeros((len(E), len(Q_GRID)))
    if verbose:
        print(f"GOS {n}{l}: {neg.sum()} sub-threshold energies (analytic continuation)", flush=True)
    G[neg] = gos_table(n, l, k2[neg], Q_GRID, rmax=400.0, h=0.004, chunk=2, verbose=verbose)
    if verbose:
        print(f"GOS {n}{l}: {(~neg).sum()} continuum energies", flush=True)
    G[~neg] = gos_table(n, l, k2[~neg], Q_GRID, rmax=80.0, h=0.002, chunk=2, verbose=verbose)
    if cache:
        np.savez_compressed(fn, G=G, E=E, Q=Q_GRID)
    return G


def _seg_int(ga, gb, qa, qb, lo, hi):
    """Integral over [lo, hi] (inside [qa, qb]) of the power law through (qa, ga), (qb, gb);
    linear interpolation where the two values are not both positive."""
    ga = np.asarray(ga, float)
    gb = np.asarray(gb, float)
    ok = (ga > 0) & (gb > 0)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        s = np.where(ok, -np.log(np.where(ok, gb, 1.0) / np.where(ok, ga, 1.0)) / np.log(qb / qa), 0.0)
        one = np.abs(1.0 - s) < 1e-9
        pl = np.where(one, ga * qa * np.log(hi / lo),
                      ga * qa ** s * (hi ** (1 - s) - lo ** (1 - s)) / np.where(one, 1.0, 1 - s))
    glo = ga + (gb - ga) * (lo - qa) / (qb - qa)
    ghi = ga + (gb - ga) * (hi - qa) / (qb - qa)
    lin = 0.5 * (glo + ghi) * (hi - lo)
    return np.where(ok, pl, lin)


class ReducedPWBA:
    """I_s(eta, theta) evaluator for one subshell."""

    def __init__(self, shell, **kw):
        n, l, ne = SHELLS[shell]
        self.n, self.l, self.ne = n, l, ne
        self.E = e_grid(n)
        G = gos_grid(n, l, **kw) * ne
        q = Q_GRID
        g = G / q[None, :] ** 3                              # integrand G/q^3
        seg = np.empty((G.shape[0], len(q) - 1))
        for j in range(len(q) - 1):
            seg[:, j] = _seg_int(g[:, j], g[:, j + 1], q[j], q[j + 1], q[j], q[j + 1])
        Gm, Gm1 = np.maximum(G[:, -1], 1e-300), np.maximum(G[:, -2], 1e-300)
        self.p = np.maximum(-np.log(Gm / Gm1) / np.log(q[-1] / q[-2]), 0.5)   # tail G ~ q^-p
        tail = G[:, -1] / q[-1] ** 2 / (self.p + 2.0)
        H = np.zeros_like(G)
        H[:, -1] = tail
        H[:, :-1] = tail[:, None] + np.cumsum(seg[:, ::-1], axis=1)[:, ::-1]
        self.G, self.g, self.H = G, g, H

    def _Hq(self, i, q0):
        """int_{q0}^{inf} G_i(q)/q^3 dq for energy row i (piecewise power law in q)."""
        q = Q_GRID
        if q0 <= q[0]:
            return self.H[i, 0] + self.G[i, 0] / q[0] ** 2 * np.log(q[0] / q0)   # G ~ q^2
        if q0 >= q[-1]:
            return self.H[i, -1] * (q[-1] / q0) ** (self.p[i] + 2.0)
        j = np.searchsorted(q, q0) - 1
        return float(self.H[i, j + 1] + _seg_int(self.g[i, j], self.g[i, j + 1],
                                                  q[j], q[j + 1], q0, q[j + 1]))

    def I(self, eta, theta):
        n = self.n
        if theta < THETA_MIN[n]:
            raise ValueError(f"theta={theta:.3f} below supported minimum {THETA_MIN[n]}")
        e0 = (theta - 1.0) / (2 * n * n)             # lowest continuum energy (signed)
        E = self.E
        i0 = np.searchsorted(E, e0)                   # first grid point >= e0
        eb = 1.0 / (2 * n * n)
        se = np.sqrt(eta)
        vals = np.array([self._Hq(i, (E[i] + eb) / se) for i in range(i0, len(E))])
        total = np.sum(0.5 * (vals[1:] + vals[:-1]) * np.diff(E[i0:]))
        if i0 > 0 and E[i0] > e0:                     # partial first interval [e0, E[i0]]
            q0 = (e0 + eb) / se
            t = (e0 - E[i0 - 1]) / (E[i0] - E[i0 - 1])
            v0 = self._Hq(i0 - 1, q0) * (1 - t) + self._Hq(i0, q0) * t
            total += 0.5 * (v0 + vals[0]) * (E[i0] - e0)
        return float(total)

    def F(self, theta, eta_over_theta2):
        eta = eta_over_theta2 * theta ** 2
        return theta * self.I(eta, theta) / eta
