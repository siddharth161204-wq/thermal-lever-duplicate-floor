#!/usr/bin/env python3
"""Independent check of Stage 4, part (c): the resolved cell at 10 C ambient, 1C and 2C, h = 10 (D59).

As D59 states: resolved_cell_steady.SteadyCell(1.0).solve(q, 10.0, T_amb=283.15) with q = (c x 5 A)^2 x R(10 C);
the within-cell spread is its hottest minus coldest node. Stage 1's two cases at 35 C (1C and 2C, h = 10) are rerun
beside it against results/resolved_cell.csv. As the checker's own bracket, the analytic centre-to-surface rise of an
infinite cylinder with the same uniform heat and the through-layer conductivity, q / (4 pi k_r L), is reported with it.

Usage: python3 verification/stage4/independent/resolved.py
Output: out/resolved_cell.csv, out/resolved_cell.json
"""
import csv, json, math, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = os.path.join(HERE, "out")
sys.path.insert(0, os.path.join(ROOT, "src"))

R35 = 0.01871279570954531
K_R, L_CELL = 1.3, 0.0701


def main():
    import resolved_cell_steady as RC
    with open(os.path.join(OUT, "resistance.json")) as f:
        R10 = json.load(f)["ambients"]["10"]["R_registered_ohm"]
    with open(os.path.join(OUT, "module.csv")) as f:
        mod = list(csv.DictReader(f))
    with open(os.path.join(ROOT, "results", "resolved_cell.csv")) as f:
        arch = [r for r in csv.DictReader(f) if r["base_h"] == "" and float(r["h"]) == 10.0]
    t0 = time.time()
    cell = RC.SteadyCell(1.0)
    rows = []
    for Ta, R in ((10.0, R10), (35.0, R35)):
        for c in (1.0, 2.0):
            q = (c * 5.0) ** 2 * R
            t1 = time.time()
            r = cell.solve(q, 10.0, T_amb=Ta + 273.15)
            m = [x for x in mod if float(x["T_amb_C"]) == Ta and float(x["crate"]) == c and float(x["h"]) == 10.0][0]
            ms = float(m["spread_K"])
            rec = {"T_amb_C": Ta, "crate": c, "h": 10.0, "R_ohm": R, "q_W": q, "nodes": r["nodes"],
                   "newton_iterations": r["newton_iterations"], "cg_all_converged": r["cg_all_converged"],
                   "energy_residual_rel": r["energy_residual_rel"], "T_max_C": r["T_max_C"], "T_min_C": r["T_min_C"],
                   "T_vol_mean_C": r["T_vol_mean_C"], "T_surface_mean_C": r["T_surface_mean_C"],
                   "dT_cell_K": r["dT_cell_K"], "dT_hot_minus_mean_K": r["dT_hot_minus_mean_K"],
                   "hot_r_mm": r["hot"]["r_mm"], "hot_z_mm": r["hot"]["z_mm"], "cold_r_mm": r["cold"]["r_mm"],
                   "cold_z_mm": r["cold"]["z_mm"], "module_spread_K": ms, "ratio_cell_over_module": r["dT_cell_K"] / ms,
                   "below_one_third": r["dT_cell_K"] < ms / 3.0,
                   "analytic_infinite_cylinder_K": q / (4 * math.pi * K_R * L_CELL),
                   "analytic_ratio_over_module": q / (4 * math.pi * K_R * L_CELL) / ms, "seconds": time.time() - t1}
            if Ta == 35.0:
                a = [x for x in arch if float(x["crate"]) == c][0]
                rec["archived_dT_cell_K"] = float(a["dT_cell_K"])
                rec["rerun_minus_archived_K"] = r["dT_cell_K"] - float(a["dT_cell_K"])
                rec["rerun_bitwise_equal"] = r["dT_cell_K"] == float(a["dT_cell_K"])
                rec["archived_q_equal"] = float(a["q_W"]) == q
            rows.append(rec)
            print(rec, flush=True)
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(OUT, "resolved_cell.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: repr(v) if isinstance(v, float) else v for k, v in r.items()})
    with open(os.path.join(OUT, "resolved_cell.json"), "w") as f:
        json.dump({"mesh_build_s": cell.build_s, "rows": rows, "total_s": time.time() - t0}, f, indent=1)


if __name__ == "__main__":
    main()
