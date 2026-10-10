"""Independent steady axisymmetric finite-volume model of the resolved cell (D6, D12, D13, D14, D15 V3).

Geometry (learned from the two-cell study's mesh files, read for geometry only): a solid cylinder R = 10.55 mm,
L = 70.1 mm, all of it jelly roll (no core, no separate end caps), k_r = 1.3 W/mK, k_z = 25 W/mK; a 0.25 mm steel
can at 15 W/mK represented as a zero-thickness shell with in-plane conductance k t on the lateral surface and on
both end discs (the shells meet at the rims; no through-thickness resistance). Uniform heat per unit volume of the
cylinder, q''' = Q / (pi R^2 L). Losses from the whole outer surface: h (T - T_amb) + eps sigma (T^4 - T_amb^4),
35 C, eps 0.8. Base-cooled case (D12): the base (z = 0) exchanges only with a 35 C plate through 600 W/(m2 K);
the lateral surface and the top keep convection and radiation.

Method: vertex-centred finite volumes on a structured (r, z) grid with nodes on the axis and on every boundary;
face conductances k A / dx; shell links between neighbouring boundary nodes; Newton iteration on the radiation;
each Newton step solved directly (SuperLU) or by conjugate gradients with an algebraic-multigrid preconditioner for
the largest grids. Nothing here comes from the study's scikit-fem code.
"""
import math
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from common import (A_END, A_OUTER, EPS, H_BASE, H_LIST, K_CAN, K_R, K_Z, L, RAD, RATES, RATE_NAMES, SIGMA,
                    T_AMB, T_CAN, q_cell, save)

KT = K_CAN * T_CAN


class Grid:
    def __init__(self, nr, nz, shell=True, kr=K_R, kz=K_Z, ends_adiabatic=False):
        self.nr, self.nz = nr, nz
        r = np.linspace(0.0, RAD, nr + 1)
        z = np.linspace(0.0, L, nz + 1)
        self.r, self.z = r, z
        rw = np.concatenate([[0.0], 0.5 * (r[:-1] + r[1:])])
        re = np.concatenate([0.5 * (r[:-1] + r[1:]), [RAD]])
        zs = np.concatenate([[0.0], 0.5 * (z[:-1] + z[1:])])
        zn = np.concatenate([0.5 * (z[:-1] + z[1:]), [L]])
        ring = math.pi * (re ** 2 - rw ** 2)          # (nr+1,)
        dzc = zn - zs                                 # (nz+1,)
        self.V = np.outer(dzc, ring).ravel()          # node volumes, index j*(nr+1)+i
        n = (nr + 1) * (nz + 1)
        self.n = n
        idx = np.arange(n).reshape(nz + 1, nr + 1)
        rows, cols, vals = [], [], []

        def link(a, b, g):
            rows.extend([a, b, a, b])
            cols.extend([a, b, b, a])
            vals.extend([g, g, -g, -g])

        # radial links in the roll
        dr = np.diff(r)
        g_r = kr * 2 * math.pi * np.outer(dzc, re[:-1] / dr)          # (nz+1, nr)
        link(idx[:, :-1].ravel(), idx[:, 1:].ravel(), g_r.ravel())
        # axial links in the roll
        dz = np.diff(z)
        g_z = kz * np.outer(1.0 / dz, ring)                           # (nz, nr+1)
        link(idx[:-1, :].ravel(), idx[1:, :].ravel(), g_z.ravel())
        if shell:
            # lateral shell along z at r = R
            link(idx[:-1, -1], idx[1:, -1], KT * 2 * math.pi * RAD / dz)
            # end shells along r at z = 0 and z = L
            g_e = KT * 2 * math.pi * re[:-1] / dr
            link(idx[0, :-1], idx[0, 1:], g_e)
            link(idx[-1, :-1], idx[-1, 1:], g_e)
        rows = np.concatenate([np.atleast_1d(np.asarray(x)) for x in rows])
        cols = np.concatenate([np.atleast_1d(np.asarray(x)) for x in cols])
        vals = np.concatenate([np.atleast_1d(np.asarray(x, dtype=float)) for x in vals])
        self.K = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
        # boundary areas per node
        self.A_lat = np.zeros(n)
        self.A_lat[idx[:, -1]] = 2 * math.pi * RAD * dzc
        self.A_bot = np.zeros(n)
        self.A_top = np.zeros(n)
        if not ends_adiabatic:
            self.A_bot[idx[0, :]] = ring
            self.A_top[idx[-1, :]] = ring
        self.rr = np.tile(r, nz + 1)
        self.zz = np.repeat(z, nr + 1)


def solve(grid, Q, h, base_h=None, tol=1e-12, use_amg=None, verbose=False):
    """Steady solve. Returns the rise field theta = T - T_amb and Newton/linear-solver diagnostics."""
    q3 = Q / (math.pi * RAD ** 2 * L)
    S = q3 * grid.V
    A_air = grid.A_lat + grid.A_top + (grid.A_bot if base_h is None else 0.0)
    g_plate = (base_h * grid.A_bot) if base_h is not None else np.zeros(grid.n)
    # initial guess: lumped rise with the same losses, linearised
    A_tot = A_air.sum()
    th = np.full(grid.n, Q / (h * A_tot + 4 * EPS * SIGMA * A_tot * T_AMB ** 3 + g_plate.sum()))
    if use_amg is None:
        use_amg = grid.n > 250000
    info = {"newton": 0, "linear": "amg-cg" if use_amg else "superlu", "cg_converged": True}
    for it in range(1, 40):
        T = T_AMB + th
        rad = EPS * SIGMA * A_air * (T ** 4 - T_AMB ** 4)
        F = grid.K @ th + h * A_air * th + rad + g_plate * th - S
        dJ = h * A_air + 4 * EPS * SIGMA * A_air * T ** 3 + g_plate
        J = (grid.K + sp.diags(dJ)).tocsr()
        if use_amg:
            import pyamg
            ml = pyamg.smoothed_aggregation_solver(J, symmetry="symmetric")
            M = ml.aspreconditioner(cycle="V")
            x, cginfo = spla.cg(J, -F, rtol=1e-13, atol=0.0, M=M, maxiter=2000)
            if cginfo != 0:
                info["cg_converged"] = False
            d = x
        else:
            d = spla.spsolve(J.tocsc(), -F)
        th = th + d
        info["newton"] = it
        if verbose:
            print("   newton", it, "max|d| =", np.abs(d).max())
        if np.abs(d).max() < tol * max(1.0, np.abs(th).max()):
            break
    T = T_AMB + th
    losses = h * A_air * th + EPS * SIGMA * A_air * (T ** 4 - T_AMB ** 4) + g_plate * th
    info["energy_residual_rel"] = float((S.sum() - losses.sum()) / S.sum())
    return th, info


def lumped_rise(Q, h, base_h=None):
    """Lumped cell with the same losses: one temperature for the whole outer surface."""
    a_air = A_OUTER - (A_END if base_h is not None else 0.0)
    g_pl = base_h * A_END if base_h is not None else 0.0
    th = Q / (h * a_air + g_pl + 4 * EPS * SIGMA * a_air * T_AMB ** 3)
    for _ in range(100):
        T = T_AMB + th
        f = h * a_air * th + g_pl * th + EPS * SIGMA * a_air * (T ** 4 - T_AMB ** 4) - Q
        df = h * a_air + g_pl + 4 * EPS * SIGMA * a_air * T ** 3
        d = -f / df
        th += d
        if abs(d) < 1e-15:
            break
    return th


def metrics(grid, th, Q, h, base_h=None):
    i_hot, i_cold = int(np.argmax(th)), int(np.argmin(th))
    A_surf = grid.A_lat + grid.A_top + grid.A_bot
    vol_mean = float(grid.V @ th / grid.V.sum())
    surf_mean = float(A_surf @ th / A_surf.sum())
    th_l = lumped_rise(Q, h, base_h)
    zc = grid.zz[i_cold]
    return {
        "T_max_C": T_AMB + th[i_hot] - 273.15, "T_min_C": T_AMB + th[i_cold] - 273.15,
        "T_vol_mean_C": T_AMB + vol_mean - 273.15, "T_surface_mean_C": T_AMB + surf_mean - 273.15,
        "T_lumped_C": T_AMB + th_l - 273.15,
        "dT_cell_K": float(th[i_hot] - th[i_cold]), "dT_hot_minus_mean_K": float(th[i_hot] - vol_mean),
        "hot_r_mm": 1e3 * grid.rr[i_hot], "hot_z_mm": 1e3 * grid.zz[i_hot],
        "cold_r_mm": 1e3 * grid.rr[i_cold], "cold_z_mm": 1e3 * zc,
        "cold_dist_from_nearest_end_mm": 1e3 * min(zc, L - zc),
        "V3b_surface_mean_minus_lumped_over_rise": (surf_mean - th_l) / th_l,
        "V3c_volume_mean_minus_lumped_over_rise": (vol_mean - th_l) / th_l,
    }


def grid_sizes(levels):
    """Nr doubles per level; Nz chosen so that dz = dr (aspect L/R = 6.645), and kept even (node at mid-height)."""
    out = []
    for nr in levels:
        nz = int(round(L / RAD * nr))
        nz += nz % 2
        out.append((nr, nz))
    return out


def analytic_check():
    """Infinite-cylinder check: ends adiabatic, no shell, linear Robin with h_eff; exact T(r) = Ts + q R^2/(4k)(1-r^2/R^2)."""
    out = []
    h_eff = 15.3
    Q = q_cell(0.5)
    q3 = Q / (math.pi * RAD ** 2 * L)
    Ts = Q / (h_eff * 2 * math.pi * RAD * L)
    for nr in (10, 20, 40):
        g = Grid(nr, 40, shell=False, ends_adiabatic=True)
        # linear Robin: emulate by solving with radiation removed
        S = q3 * g.V
        J = (g.K + sp.diags(h_eff * g.A_lat)).tocsc()
        th = spla.spsolve(J, S)
        exact = Ts + q3 * RAD ** 2 / (4 * K_R) * (1 - (g.rr / RAD) ** 2)
        out.append({"nr": nr, "max_abs_err_K": float(np.abs(th - exact).max()),
                    "centre_minus_surface_numeric": float(th[0] - th[nr]),
                    "centre_minus_surface_exact": q3 * RAD ** 2 / (4 * K_R)})
    return out


def richardson(vals, ratio=2.0):
    """Observed order and extrapolate from the last three values of a sequence on grids refined by `ratio`."""
    f1, f2, f3 = vals[-3], vals[-2], vals[-1]
    d1, d2 = f2 - f1, f3 - f2
    if d1 == 0 or d2 == 0 or (d1 / d2) <= 0:
        return {"order": float("nan"), "extrapolated": f3, "rel_change_last": (f3 - f2) / f3}
    p = math.log(d1 / d2) / math.log(ratio)
    ext = f3 + d2 / (ratio ** p - 1)
    return {"order": p, "extrapolated": ext, "rel_change_last": (f3 - f2) / f3}


def main(levels=(20, 40, 80, 160, 320), all_rates_level=160):
    t0 = time.time()
    out = {"analytic_infinite_cylinder": analytic_check(), "convergence": {}, "all_cases": [], "grids": {}}
    print("analytic check", out["analytic_infinite_cylinder"])
    sizes = grid_sizes(levels)
    grids = {}
    for nr, nz in sizes:
        grids[nr] = Grid(nr, nz)
        out["grids"][str(nr)] = {"nr": nr, "nz": nz, "nodes": grids[nr].n, "dr_mm": 1e3 * RAD / nr, "dz_mm": 1e3 * L / nz}
    # grid convergence: C/2 and 2C at h = 5, 10, 25, plus the base-cooled case at h = 10 for every rate
    conv_cases = [(c, h, None) for c in (0.5, 2.0) for h in H_LIST] + [(c, 10.0, H_BASE) for c in RATES]
    for c, h, bh in conv_cases:
        key = f"c={c}_h={h}_base={bh}"
        seq = []
        for nr, nz in sizes:
            th, info = solve(grids[nr], q_cell(c), h, bh)
            m = metrics(grids[nr], th, q_cell(c), h, bh)
            m.update(nr=nr, nz=nz, nodes=grids[nr].n, **info)
            seq.append(m)
            print(f"{key:28s} nr={nr:4d} nz={nz:5d} dT={m['dT_cell_K']:.8f} hot-mean={m['dT_hot_minus_mean_K']:.8f} "
                  f"Tvol={m['T_vol_mean_C']:.8f} Tsurf={m['T_surface_mean_C']:.8f} Tl={m['T_lumped_C']:.8f} "
                  f"cold=({m['cold_r_mm']:.3f},{m['cold_z_mm']:.3f}) hot=({m['hot_r_mm']:.3f},{m['hot_z_mm']:.3f}) "
                  f"newton={info['newton']} res={info['energy_residual_rel']:.1e} [{time.time() - t0:.0f}s]", flush=True)
        rich = {q: richardson([s[q] for s in seq]) for q in ("dT_cell_K", "dT_hot_minus_mean_K", "T_vol_mean_C",
                                                             "T_surface_mean_C", "T_max_C", "T_min_C")}
        out["convergence"][key] = {"crate": c, "h": h, "base_h": bh, "q_W": q_cell(c), "sequence": seq, "richardson": rich}
    # every registered rate and h on one fine grid
    g = grids[all_rates_level]
    for c in RATES:
        for h in H_LIST:
            th, info = solve(g, q_cell(c), h, None)
            m = metrics(g, th, q_cell(c), h, None)
            m.update(crate=c, h=h, base_h=None, q_W=q_cell(c), nr=g.nr, nz=g.nz, **info)
            out["all_cases"].append(m)
    save("resolved_cell.json", out)
    print("done in", round(time.time() - t0), "s")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "quick":
        main(levels=(10, 20, 40), all_rates_level=40)
    else:
        main()
