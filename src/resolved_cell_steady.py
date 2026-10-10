#!/usr/bin/env python3
"""Step 2(b): the resolved cell at steady state, with a uniform internal heat and no neighbour (DECISIONS.md D13, D14).

Reuses the two-cell study's scikit-fem cell: src/twocell/resolved_cell2.py, mesh_cell2.py and lumped_two_cell.py,
copied unchanged from github.com/siddharth161204-wq/m50-two-cell-pilot at tag v1.0.0 (commit 2b35ac7; doi
10.5281/zenodo.23253738; SHA-256 in DECISIONS.md D13). From that model come the half-cylinder mesh (diameter 21.1 mm,
length 70.1 mm), the anisotropic jelly-roll conduction operator (1.3 W/mK through the layers, 25 W/mK along them,
the two-cell study's D11), the 0.25 mm steel can as a conducting shell (15 W/mK) and the boundary groups.

What this driver changes, without editing the copied files:
  1. No neighbour: no strips, no air gap and no radiation from a neighbour. Convection acts on the WHOLE outer surface
     (the two-cell study removed the gap-facing quarter because a neighbour sat there).
  2. Ambient 35 C; h swept over 5, 10 and 25 W/(m2 K); grey radiation at emissivity 0.8 to 35 C from the whole outer
     surface, the two-cell model's room-radiation term.
  3. The uniform internal heat is spread per unit volume of the cylinder.
  4. The steady state is solved directly by Newton iteration on the same assembled operators instead of by time
     stepping, because only the steady field is wanted.
  5. Secondary and exploratory (D12): the bottom face sits on a cold plate at 35 C through 600 W/(m2 K); the other
     surfaces keep convection and radiation.

Usage:
  python3 src/resolved_cell_steady.py --out results            (the sweep, after the predictions are public)
  python3 src/resolved_cell_steady.py --verify --out verification   (V3)
"""
import argparse, json, math, os, sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "twocell"))
sys.path.insert(0, HERE)
import lumped_two_cell as LP                       # noqa: E402  (constants of the two-cell cell, unchanged)
from resolved_cell2 import Cell2Model              # noqa: E402
import thermal_lever as TL                         # noqa: E402  (rates, resistance, heat per cell)

SIGMA = TL.SIGMA
EPS = TL.EPS
T_AMB = TL.T_AMB
R_CELL = LP.D_CELL / 2
A_END = math.pi * R_CELL ** 2


class SteadyCell:
    def __init__(self, refine=1.0):
        t0 = time.time()
        self.m = Cell2Model(refine)
        self.refine = refine
        m = self.m
        self.full = lambda k: 2.0 * m.w[k].sum()
        self.n = m.C.shape[0]
        self.build_s = time.time() - t0

    def operators(self, h, base_h=None):
        m = self.m
        if base_h is None:
            w = m.w["outer"]
            area = LP.A_CELL                              # exact outer area; weights normalised as in the two-cell V0
            d_conv = h * area * w / (2.0 * w.sum())
            w_rad = EPS * SIGMA * area * w / (2.0 * w.sum())
            d_base = np.zeros(self.n)
        else:
            w = m.w["top"] + m.w["lateral"]
            area = LP.A_CELL - A_END
            d_conv = h * area * w / (2.0 * w.sum())
            w_rad = EPS * SIGMA * area * w / (2.0 * w.sum())
            wb = m.w["bottom"]
            d_base = base_h * A_END * wb / (2.0 * wb.sum())
        return d_conv, w_rad, d_base

    def solve(self, q_W, h, T_amb=T_AMB, base_h=None, tol=1e-11, maxit=60):
        """Newton iteration on the rise u = T - T_amb (the conduction operator annihilates a uniform field), each step
        solved by conjugate gradients with a smoothed-aggregation algebraic multigrid preconditioner (pyamg) built once
        from the linearised operator at u = 0, to a relative residual of 1e-13."""
        import pyamg
        m = self.m
        d_conv, w_rad, d_base = self.operators(h, base_h)
        src = 0.5 * q_W * m.w_vol / m.w_vol.sum()                 # half model, uniform per unit volume
        K = m.K.tocsr()
        J0 = (K + sp.diags(d_conv + d_base + 4.0 * w_rad * T_amb ** 3)).tocsr()
        ml = pyamg.smoothed_aggregation_solver(J0, symmetry="symmetric")
        M = ml.aspreconditioner(cycle="V")
        u = np.zeros(self.n)
        self.cg_info = []
        for it in range(maxit):
            T = T_amb + u
            F = K @ u + (d_conv + d_base) * u + w_rad * (T ** 4 - T_amb ** 4) - src
            J = (K + sp.diags(d_conv + d_base + 4.0 * w_rad * T ** 3)).tocsr()
            du, info = spla.cg(J, -F, rtol=1e-13, atol=0.0, maxiter=5000, M=M)
            self.cg_info.append(info)
            u = u + du
            if np.abs(du).max() < tol:
                break
        T = T_amb + u
        loss = 2.0 * float((d_conv + d_base) @ (T - T_amb) + w_rad @ (T ** 4 - T_amb ** 4))
        vol_mean = float(m.w_vol @ T / m.w_vol.sum())
        w_out = m.w["outer"]
        surf_mean = float(w_out @ T / w_out.sum())
        i_hot, i_cold = int(np.argmax(T)), int(np.argmin(T))
        loc = lambda i: {"r_mm": 1e3 * float(m.r_node[i]), "z_mm": 1e3 * float(m.z_node[i]),
                         "theta_deg": math.degrees(float(m.th_node[i]))}
        return {"q_W": q_W, "h": h, "base_h": base_h, "refine": self.refine, "nodes": self.n, "newton_iterations": it + 1,
                "T_max_C": float(T.max() - 273.15), "T_min_C": float(T.min() - 273.15),
                "T_vol_mean_C": vol_mean - 273.15, "T_surface_mean_C": surf_mean - 273.15,
                "dT_cell_K": float(T.max() - T.min()), "dT_hot_minus_mean_K": float(T.max() - vol_mean),
                "hot": loc(i_hot), "cold": loc(i_cold),
                "energy_residual_rel": (q_W - loss) / q_W if q_W else 0.0}


def lumped_cell(q_W, h, T_amb=T_AMB, base_h=None):
    """A single lumped cell with the same losses: h A (T - Ta) + eps sigma A (T^4 - Ta^4) (+ plate) = q."""
    area = LP.A_CELL if base_h is None else LP.A_CELL - A_END
    g_base = 0.0 if base_h is None else base_h * A_END
    T = T_amb + 1.0
    for _ in range(100):
        F = q_W - h * area * (T - T_amb) - EPS * SIGMA * area * (T ** 4 - T_amb ** 4) - g_base * (T - T_amb)
        dF = -h * area - 4 * EPS * SIGMA * area * T ** 3 - g_base
        dT = -F / dF
        T += dT
        if abs(dT) < 1e-13:
            break
    return T


def sweep(out, refine=1.0):
    cell = SteadyCell(refine)
    rows = []
    for c in TL.CRATES:
        q = TL.heat_per_cell(c)
        for h in TL.H_LIST:
            r = cell.solve(q, h)
            r["crate"] = c
            r["T_lumped_C"] = lumped_cell(q, h) - 273.15
            rows.append(r)
        r = cell.solve(q, 10.0, base_h=TL.H_BASE)            # secondary: cold plate under the base, h = 10 elsewhere
        r["crate"] = c
        r["T_lumped_C"] = lumped_cell(q, 10.0, base_h=TL.H_BASE) - 273.15
        rows.append(r)
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "resolved_cell.json"), "w") as f:
        json.dump({"mesh_nodes_half_model": cell.n, "refine": refine, "rows": rows}, f, indent=1)
    import csv
    keys = ["crate", "q_W", "h", "base_h", "T_max_C", "T_min_C", "T_vol_mean_C", "T_surface_mean_C", "T_lumped_C",
            "dT_cell_K", "dT_hot_minus_mean_K", "energy_residual_rel", "newton_iterations"]
    with open(os.path.join(out, "resolved_cell.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys + ["hot_r_mm", "hot_z_mm", "cold_r_mm", "cold_z_mm"]); w.writeheader()
        for r in rows:
            d = {k: r[k] for k in keys}
            d.update({"hot_r_mm": r["hot"]["r_mm"], "hot_z_mm": r["hot"]["z_mm"], "cold_r_mm": r["cold"]["r_mm"], "cold_z_mm": r["cold"]["z_mm"]})
            w.writerow(d)
    return rows


def verify_v3(out):
    """V3 as fixed in DECISIONS.md D15 before any resolved run."""
    res = {}
    base = SteadyCell(1.0)
    fine = SteadyCell(2.0)
    res["mesh"] = {"refine1_nodes": base.n, "refine2_nodes": fine.n}
    cases = []
    q = TL.heat_per_cell(0.5)
    for h in TL.H_LIST:
        r1 = base.solve(q, h)
        r2 = fine.solve(q, h)
        TL_ = lumped_cell(q, h)
        rise = TL_ - T_AMB
        offset_hand = q / (8 * math.pi * LP.K_RADIAL * LP.L_CELL)           # radial conduction, infinite cylinder
        cases.append({
            "h": h, "q_W": q,
            "a_energy_residual_rel": r1["energy_residual_rel"], "a_pass": abs(r1["energy_residual_rel"]) < 1e-8,
            "b_surface_mean_minus_lumped_over_rise": (r1["T_surface_mean_C"] + 273.15 - TL_) / rise,
            "b_pass": abs((r1["T_surface_mean_C"] + 273.15 - TL_) / rise) < 0.01,
            "c_volume_mean_minus_lumped_over_rise": (r1["T_vol_mean_C"] + 273.15 - TL_) / rise,
            "c_brief_literal_pass_1pct": abs((r1["T_vol_mean_C"] + 273.15 - TL_) / rise) < 0.01,
            "c_hand_offset_over_rise_stated_before_the_run": offset_hand / rise,
            "d_dT_cell_refine1_K": r1["dT_cell_K"], "d_dT_cell_refine2_K": r2["dT_cell_K"],
            "d_rel_change": (r2["dT_cell_K"] - r1["dT_cell_K"]) / r1["dT_cell_K"],
            "d_pass": abs((r2["dT_cell_K"] - r1["dT_cell_K"]) / r1["dT_cell_K"]) < 0.05,
            "d_dT_hot_minus_mean_rel_change": (r2["dT_hot_minus_mean_K"] - r1["dT_hot_minus_mean_K"]) / r1["dT_hot_minus_mean_K"],
        })
    res["cases"] = cases
    # e: analytic radial check at these conditions (roll operator alone, lateral Robin, ends adiabatic)
    m = base.m
    h_eff = 10.0 + 4 * EPS * SIGMA * T_AMB ** 3
    q3 = q / (math.pi * R_CELL ** 2 * LP.L_CELL)
    w_lat = m.w["lateral"]
    import pyamg
    A = (m.K_roll + sp.diags(h_eff * w_lat)).tocsr()
    ml = pyamg.smoothed_aggregation_solver(A, symmetry="symmetric")
    x = T_AMB + ml.solve(q3 * m.w_vol, tol=1e-13, accel="cg", maxiter=2000)
    ana = T_AMB + q3 * R_CELL / (2 * h_eff) + q3 * (R_CELL ** 2 - m.r_node ** 2) / (4 * LP.K_RADIAL)
    res["e_analytic_radial"] = {"h_eff": h_eff, "q3_W_per_m3": q3, "max_rise_analytic_K": float(ana.max() - T_AMB),
                                "max_abs_err_K": float(np.abs(x - ana).max()),
                                "rel_err": float(np.abs(x - ana).max() / (ana.max() - T_AMB)),
                                "center_to_surface_analytic_K": float(q3 * R_CELL ** 2 / (4 * LP.K_RADIAL)),
                                "center_to_surface_numeric_K": float(x.max() - x[np.argmax(m.r_node)])}
    res["e_pass"] = res["e_analytic_radial"]["rel_err"] < 0.01
    res["V3_numerical_pass"] = all(c["a_pass"] and c["b_pass"] and c["d_pass"] for c in cases) and res["e_pass"]
    res["V3_brief_literal_volume_mean_pass"] = all(c["c_brief_literal_pass_1pct"] for c in cases)
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "v3_resolved_cell.json"), "w") as f:
        json.dump(res, f, indent=1)
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--refine", type=float, default=1.0)
    a = ap.parse_args()
    if a.verify:
        print(json.dumps(verify_v3(a.out), indent=1))
    else:
        rows = sweep(a.out, a.refine)
        print("Step 2(b) written:", len(rows), "cases")
