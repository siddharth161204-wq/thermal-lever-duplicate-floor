"""Verification of the finite-volume cell solver in two dimensions against an analytic series solution.

Problem: finite anisotropic cylinder (k_r, k_z), uniform source q''', linear Robin h on the lateral surface and
both ends, no shell. Series: theta(r, z) = sum_n a_n / (k_r lam_n^2) J0(lam_n r) [1 - C_n cosh(m_n (z - L/2))],
k_r lam_n J1(lam_n R) = h J0(lam_n R), m_n = lam_n sqrt(k_r / k_z),
C_n = h / (k_z m_n sinh(m_n L/2) + h cosh(m_n L/2)),
a_n = q''' [R J1(lam_n R) / lam_n] / [(R^2/2)(J0(lam_n R)^2 + J1(lam_n R)^2)].
Also reports the effect of the can shell on the registered quantities (shell on versus off, same grid).
"""
import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import brentq
from scipy.special import j0, j1

from common import K_R, K_Z, L, RAD, q_cell, save
import resolved_cell as RC


def roots(h, n_terms):
    Bi = h * RAD / K_R
    f = lambda x: x * j1(x) - Bi * j0(x)          # x = lam R
    out = []
    x = 1e-9
    step = 0.05
    fa = f(x)
    while len(out) < n_terms:
        xb = x + step
        fb = f(xb)
        if fa * fb < 0:
            out.append(brentq(f, x, xb, xtol=1e-15, rtol=1e-15, maxiter=200))
        x, fa = xb, fb
    return np.array(out) / RAD


def series(h, q3, r, z, n_terms=4000):
    lam = roots(h, n_terms)
    x = lam * RAD
    a = q3 * (RAD * j1(x) / lam) / (0.5 * RAD ** 2 * (j0(x) ** 2 + j1(x) ** 2))
    m = lam * math.sqrt(K_R / K_Z)
    # C_n written with exp to avoid overflow: cosh(m(z-L/2))/[...] with large m
    th = 0.0
    for an, ln, mn in zip(a, lam, m):
        u = mn * L / 2
        denom = K_Z * mn * math.tanh(u) + h            # divided by cosh(u)
        ratio = math.exp(mn * (abs(z - L / 2) - L / 2)) * (1 + math.exp(-2 * mn * abs(z - L / 2))) / (1 + math.exp(-2 * u))
        th += an / (K_R * ln ** 2) * j0(ln * r) * (1 - h / denom * ratio)
    return th


def fv_linear(nr, h, q3, shell):
    nz = int(round(L / RAD * nr))
    nz += nz % 2
    g = RC.Grid(nr, nz, shell=shell)
    A_air = g.A_lat + g.A_top + g.A_bot
    J = (g.K + sp.diags(h * A_air)).tocsc()
    th = spla.spsolve(J, q3 * g.V)
    i_c = 0 + (nz // 2) * (nr + 1)        # r = 0, z = L/2
    i_rim = nr                             # r = R, z = 0
    return g, th, th[i_c], th[i_rim]


def main():
    h = 15.3
    Q = q_cell(0.5)
    q3 = Q / (math.pi * RAD ** 2 * L)
    centre = series(h, q3, 0.0, L / 2)
    rim = series(h, q3, RAD, 0.0)
    out = {"h_linear": h, "series_centre_rise_K": centre, "series_rim_rise_K": rim,
           "series_centre_minus_rim_K": centre - rim, "fv": []}
    for nr in (20, 40, 80, 160):
        g, th, c, rmn = fv_linear(nr, h, q3, shell=False)
        out["fv"].append({"nr": nr, "centre_rise_K": c, "rim_rise_K": rmn, "centre_minus_rim_K": c - rmn,
                          "max_minus_min_K": float(th.max() - th.min()),
                          "rel_err_centre_minus_rim": (c - rmn - (centre - rim)) / (centre - rim)})
    # effect of the shell on the registered quantity (C/2, h = 10, radiation on), nr = 160
    g = RC.Grid(160, 1064, shell=True)
    g0 = RC.Grid(160, 1064, shell=False)
    th1, _ = RC.solve(g, Q, 10.0)
    th0, _ = RC.solve(g0, Q, 10.0)
    m1 = RC.metrics(g, th1, Q, 10.0)
    m0 = RC.metrics(g0, th0, Q, 10.0)
    out["shell_effect_C2_h10"] = {"with_shell": {k: m1[k] for k in ("dT_cell_K", "dT_hot_minus_mean_K", "cold_r_mm", "cold_z_mm")},
                                  "without_shell": {k: m0[k] for k in ("dT_cell_K", "dT_hot_minus_mean_K", "cold_r_mm", "cold_z_mm")}}
    # D6: "the within-cell spread would be about 5 % larger with 1.24 W/mK" (C/2, h = 10)
    gk = RC.Grid(160, 1064, kr=1.24)
    thk, _ = RC.solve(gk, Q, 10.0)
    mk = RC.metrics(gk, thk, Q, 10.0)
    out["kr_1.24_vs_1.30_C2_h10"] = {"dT_cell_kr_1.30": m1["dT_cell_K"], "dT_cell_kr_1.24": mk["dT_cell_K"],
                                     "ratio": mk["dT_cell_K"] / m1["dT_cell_K"]}
    print(out["kr_1.24_vs_1.30_C2_h10"])
    save("verify_fv_series.json", out)
    print("series centre", centre, "rim", rim, "diff", centre - rim)
    for f in out["fv"]:
        print(f)
    print(out["shell_effect_C2_h10"])


if __name__ == "__main__":
    main()
