#!/usr/bin/env python3
"""Stage 4: the eight-cell module at 35, 25 and 10 C ambient and the resolved cell at 10 C (DECISIONS.md D58, D59),
with V1 (D62): the 35 C cases against Stage 1's archived values.

Stage 1's network (thermal_lever.module_steady) and resolved cell (resolved_cell_steady.SteadyCell) are imported
unchanged. Each cell's heat is (c x 5 A)^2 x R(ambient), R from results/stage4_resistance.json (the archived value at
35 C). Outputs: results/stage4_module.csv, results/stage4_resolved_cell.csv, verification/stage4/v1.json.

Usage: python3 src/stage4/module.py
"""
import argparse, csv, json, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as CM                     # noqa: E402
import thermal_lever as TL              # noqa: E402  (Stage 1, unchanged)

CELLS = range(8)


def module_case(T_C, c, h, R):
    q = (c * TL.CAP_AH) ** 2 * R                              # Stage 1's heat_per_cell, with the ambient's resistance
    T, bal, info = TL.module_steady(q, h, T_amb=T_C + 273.15)
    i_hot, i_cold = int(T.argmax()), int(T.argmin())
    band = CM.BAND_OF_AMBIENT[T_C]
    lo, hi = CM.BAND_SETPOINTS_C[band]
    row = {"ambient_C": T_C, "band": band, "crate": c, "rate": CM.RATE_LABEL[c], "h": h, "R_ohm": R, "q_cell_W": q,
           "T_max_C": float(T.max() - 273.15), "T_min_C": float(T.min() - 273.15), "T_mean_C": float(T.mean() - 273.15),
           "spread_K": float(T.max() - T.min()), "hot_minus_mean_K": float(T.max() - T.mean()),
           "i_hot": i_hot, "hot_position": TL.POSITIONS[i_hot], "i_cold": i_cold, "cold_position": TL.POSITIONS[i_cold],
           "balance_W": bal, "newton_iterations": info["iterations"],
           "converged": bool(abs(bal) <= 1e-9 and info["iterations"] < 100),
           "cells_inside_band_setpoints": bool(lo <= float(T.min() - 273.15) and float(T.max() - 273.15) <= hi)}
    for i in CELLS:
        row[f"T_cell{i}_C"] = float(T[i] - 273.15)
    for i in CELLS:
        row[f"offset_cell{i}_K"] = float(T[i] - T.mean())
    return row


def run_module(R_used):
    rows = []
    for T_C in CM.AMBIENTS_C:
        for c in CM.RATES:
            for h in CM.H_LIST:
                rows.append(module_case(T_C, c, h, R_used[f"{T_C:g}"]))
    return rows


def v1(rows):
    ref = [r for r in csv.DictReader(open(os.path.join(CM.ROOT, "results", "module_spread.csv"))) if r["arrangement"] == "2x4 in air"]
    pairs = {"q_cell_W": "q_cell_W", "T_max_C": "T_max_C", "T_min_C": "T_min_C", "T_mean_C": "T_mean_C",
             "spread_K": "dT_module_K", "hot_minus_mean_K": "dT_hot_minus_mean_K"}
    pairs.update({f"T_cell{i}_C": f"T_cell{i}_C" for i in CELLS})
    cases, worst, bitwise, printed = [], 0.0, True, True
    for r in rows:
        if r["ambient_C"] != 35.0 or r["crate"] not in (0.5, 1.0, 2.0):
            continue
        a = [x for x in ref if float(x["crate"]) == r["crate"] and float(x["h"]) == r["h"]]
        if len(a) != 1:
            raise SystemExit("V1: no unique Stage 1 row")
        a = a[0]
        diffs = {k: abs(r[k] - float(a[v])) for k, v in pairs.items()}
        same = all(r[k] == float(a[v]) for k, v in pairs.items())
        six = all(round(r[k], 6) == round(float(a[v]), 6) for k, v in pairs.items() if k != "q_cell_W")
        ints = (r["i_hot"] == int(a["i_hot"]) and r["i_cold"] == int(a["i_cold"])
                and r["newton_iterations"] == int(a["newton_iterations"]))
        worst = max(worst, max(d for k, d in diffs.items() if k != "q_cell_W"))
        bitwise &= same
        printed &= six and ints
        cases.append({"crate": r["crate"], "h": r["h"], "spread_K": r["spread_K"], "spread_stage1_K": float(a["dT_module_K"]),
                      "max_abs_diff_K": max(d for k, d in diffs.items() if k != "q_cell_W"), "q_diff_W": diffs["q_cell_W"],
                      "bitwise_identical": same, "six_decimals_identical": six, "hot_cold_iterations_identical": ints})
    ok = len(cases) == 9 and worst <= 1e-9 and printed
    return {"check": "V1 (D62): the 35 C cases at C/2, 1C and 2C against results/module_spread.csv", "cases": cases,
            "max_abs_diff_K": worst, "all_bitwise_identical": bitwise, "pass": ok}


def run_resolved(R10):
    import resolved_cell_steady as RC                      # Stage 1, unchanged (needs scikit-fem and pyamg)
    cell = RC.SteadyCell(1.0)
    ref = [r for r in csv.DictReader(open(os.path.join(CM.ROOT, "results", "resolved_cell.csv")))]
    rows = []
    for c in (1.0, 2.0):
        q10 = (c * TL.CAP_AH) ** 2 * R10
        r = cell.solve(q10, 10.0, T_amb=10.0 + 273.15)
        rows.append({"ambient_C": 10.0, "crate": c, "rate": CM.RATE_LABEL[c], "h": 10.0, "R_ohm": R10, **_flat(r)})
    for c in (1.0, 2.0):
        r = cell.solve(TL.heat_per_cell(c), 10.0)                     # Stage 1's case at 35 C, rerun
        a = [x for x in ref if float(x["crate"]) == c and float(x["h"]) == 10.0 and x["base_h"] in ("", "nan", "NaN")][0]
        d = {"ambient_C": 35.0, "crate": c, "rate": CM.RATE_LABEL[c], "h": 10.0, "R_ohm": TL.R_DC, **_flat(r)}
        d["dT_cell_K_stage1"] = float(a["dT_cell_K"])
        d["dT_cell_diff_K"] = d["dT_cell_K"] - float(a["dT_cell_K"])
        d["T_max_diff_K"] = d["T_max_C"] - float(a["T_max_C"])
        rows.append(d)
    return rows, {"mesh_nodes_half_model": cell.n, "refine": cell.refine}


def _flat(r):
    keep = ("q_W", "T_max_C", "T_min_C", "T_vol_mean_C", "T_surface_mean_C", "dT_cell_K", "dT_hot_minus_mean_K",
            "energy_residual_rel", "newton_iterations", "cg_all_converged")
    out = {k: r[k] for k in keep}
    out.update({"hot_r_mm": r["hot"]["r_mm"], "hot_z_mm": r["hot"]["z_mm"], "cold_r_mm": r["cold"]["r_mm"], "cold_z_mm": r["cold"]["z_mm"]})
    return out


def write_csv(path, rows):
    keys = list(rows[0].keys())
    for r in rows[1:]:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main(res_path, out_dir, ver_dir):
    res = json.load(open(res_path))
    R_used = res["R_used_ohm"]
    if R_used["35"] != CM.read_r35():
        raise SystemExit("the 35 C resistance is not the archived one")
    rows = run_module(R_used)
    if not all(r["converged"] for r in rows):
        print("not converged:", [(r["ambient_C"], r["crate"], r["h"]) for r in rows if not r["converged"]])
    write_csv(os.path.join(out_dir, "stage4_module.csv"), rows)
    os.makedirs(ver_dir, exist_ok=True)
    chk = v1(rows)
    json.dump(chk, open(os.path.join(ver_dir, "v1.json"), "w"), indent=1)
    print("V1 pass:", chk["pass"], "max diff", chk["max_abs_diff_K"], "bitwise", chk["all_bitwise_identical"])
    rc, meta = run_resolved(R_used["10"])
    for r in rc:
        m = [x for x in rows if x["ambient_C"] == r["ambient_C"] and x["crate"] == r["crate"] and x["h"] == r["h"]][0]
        r["module_spread_K"] = m["spread_K"]
        r["within_over_module"] = r["dT_cell_K"] / m["spread_K"]
    write_csv(os.path.join(out_dir, "stage4_resolved_cell.csv"), rc)
    json.dump({"meta": meta, "rows": rc}, open(os.path.join(out_dir, "stage4_resolved_cell.json"), "w"), indent=1)
    for r in rows:
        if r["h"] == 10.0:
            print(f"{r['ambient_C']:>4g} C {r['rate']:>3}: spread {r['spread_K']:.6f} K, T {r['T_min_C']:.2f} to {r['T_max_C']:.2f} C")
    for r in rc:
        print(f"resolved {r['ambient_C']:g} C {r['rate']}: within-cell {r['dT_cell_K']:.6f} K, module {r['module_spread_K']:.6f} K")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--resistance", default=os.path.join(CM.ROOT, "results", "stage4_resistance.json"))
    ap.add_argument("--out", default=os.path.join(CM.ROOT, "results"))
    ap.add_argument("--verification", default=os.path.join(CM.ROOT, "verification", "stage4"))
    a = ap.parse_args()
    main(a.resistance, a.out, a.verification)
