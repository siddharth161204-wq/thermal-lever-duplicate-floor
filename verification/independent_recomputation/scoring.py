"""Implied divergence, Step 3 comparison and prediction scores, from this check's own outputs only.

Registered forms (D17, D18): implied divergence = |s| x Delta T (module hottest minus coldest; cell hottest minus
coldest); Step 3 = lowest rate on C/16..2C at which the spread exceeds |Delta T_floor| at that rate (h fixed);
scoring rule for P1 to P5 as written in predictions.md.
Secondary (D18): the law evaluated at each module cell's own temperature, exact pair difference between hottest and
coldest cells and exact worst-of-8 shortfall, where the cells lie inside the law's range (Law 1: 15 to 55 C,
Law 2: 0 to 55 C); values outside the range are also computed and flagged.
"""
import math

import numpy as np

import law1
import law2
from common import FLOORS, FLOOR_KEYS, H_LIST, RATES, RATE_NAMES, load_json, save, OUT
import os

LAWS = ["law_1_kucinskis2022", "law_2_schimpe2018"]
RANGE = {"law_1_kucinskis2022": (15.0, 55.0), "law_2_schimpe2018": (0.0, 55.0)}


def main(cell_level=None):
    L1 = load_json(os.path.join(OUT, "law1.json"))
    L2 = load_json(os.path.join(OUT, "law2.json"))
    MOD = load_json(os.path.join(OUT, "module.json"))
    RC = load_json(os.path.join(OUT, "resolved_cell.json"))
    s = {"law_1_kucinskis2022": {n: L1["C2"]["s_centred_secant"] for n in RATE_NAMES},
         "law_2_schimpe2018": {n: L2["rates"][n]["s"] for n in RATE_NAMES}}
    dTf = {law: {n: {k: FLOORS[k] / s[law][n] for k in FLOOR_KEYS} for n in RATE_NAMES} for law in LAWS}
    mod = {(r["crate"], r["h"]): r for r in MOD["module"]}
    cell = {(r["crate"], r["h"]): r for r in RC["all_cases"] if r["base_h"] is None}
    p2 = law2.params()

    def L_of(law, T_C, c):
        if law == "law_1_kucinskis2022":
            return math.exp(law1.lnL_rate(0.5, T_C))
        return law2.lifetime(273.15 + T_C, c, p2)[0]

    rows = []
    for law in LAWS:
        for c, n in zip(RATES, RATE_NAMES):
            for h in H_LIST:
                m, ce = mod[(c, h)], cell[(c, h)]
                ss = abs(s[law][n])
                Ts = m["T_cells_C"]
                lo, hi = RANGE[law]
                inside = (min(Ts) >= lo - 1e-9) and (max(Ts) <= hi + 1e-9)
                Ls = [L_of(law, T, c) for T in Ts]
                La, Lb = L_of(law, m["T_min_C"], c), L_of(law, m["T_max_C"], c)
                exact_pair = abs(La - Lb) / (0.5 * (La + Lb))
                mean_L = float(np.mean(Ls))
                exact_w8 = (mean_L - min(Ls)) / mean_L
                row = {"law": law, "rate": n, "crate": c, "h": h, "s_per_K": s[law][n],
                       "dT_module_K": m["dT_module_K"], "dT_module_hot_minus_mean_K": m["dT_hot_minus_mean_K"],
                       "dT_cell_K": ce["dT_cell_K"], "dT_cell_hot_minus_mean_K": ce["dT_hot_minus_mean_K"],
                       "implied_module": ss * m["dT_module_K"], "implied_module_hot_minus_mean": ss * m["dT_hot_minus_mean_K"],
                       "implied_cell": ss * ce["dT_cell_K"],
                       "exact_pair_hottest_coldest": exact_pair, "exact_worst_of_8_shortfall": exact_w8,
                       "cells_inside_law_range": inside}
                for k in FLOOR_KEYS:
                    row[f"module_exceeds_{k}"] = m["dT_module_K"] > abs(dTf[law][n][k])
                    row[f"cell_exceeds_{k}"] = ce["dT_cell_K"] > abs(dTf[law][n][k])
                rows.append(row)
    # Step 3
    step3 = {}
    for law in LAWS:
        step3[law] = {}
        for h in H_LIST:
            step3[law][str(h)] = {}
            for k in FLOOR_KEYS:
                res = {}
                for what, key in (("module", "dT_module_K"), ("cell", "dT_cell_K")):
                    first = "none up to 2C"
                    for c, n in zip(RATES, RATE_NAMES):
                        val = mod[(c, h)][key] if what == "module" else cell[(c, h)][key]
                        if val > abs(dTf[law][n][k]):
                            first = n
                            break
                    res[what] = first
                step3[law][str(h)][k] = res
    # dataset conditions (C/16..C/2, every h), brief's heat
    ds = {}
    for law in LAWS:
        sel = [r for r in rows if r["law"] == law and r["crate"] <= 0.5]
        ds[law] = {"max_implied_module": max(r["implied_module"] for r in sel),
                   "max_implied_cell": max(r["implied_cell"] for r in sel),
                   "max_dT_module_K": max(r["dT_module_K"] for r in sel)}
    # D19, cycle-averaged heat
    d19 = []
    for r in MOD["D19"]:
        n = r["rate"]
        d19.append({"rate": n, "h": r["h"], "q_cycle_avg_W": r["q_cell_W"], "ratio_to_brief_heat": r["ratio_to_brief_heat"],
                    "dT_module_K": r["dT_module_K"],
                    "implied_law_1_kucinskis2022": abs(s["law_1_kucinskis2022"][n]) * r["dT_module_K"],
                    "implied_law_2_schimpe2018": abs(s["law_2_schimpe2018"][n]) * r["dT_module_K"]})
    # predictions
    ds_rates = RATE_NAMES[:4]
    p1_vals = {law: {n: dTf[law][n]["p95_8.7pct"] for n in ds_rates} for law in LAWS}
    p1_per_law = {law: all(1.0 <= v <= 3.0 for v in p1_vals[law].values()) for law in LAWS}
    P1 = {"values_K": p1_vals, "held_per_law": p1_per_law, "held": all(p1_per_law.values())}
    p2v = mod[(0.5, 10.0)]["dT_module_K"]
    P2 = {"value_K": p2v, "by_h": {str(h): mod[(0.5, h)]["dT_module_K"] for h in H_LIST}, "held": p2v < 1.0}
    p3 = {str(h): cell[(0.5, h)]["dT_cell_K"] for h in H_LIST}
    P3 = {"values_K": p3, "held": all(v < 0.5 for v in p3.values())}
    order = RATE_NAMES
    P4 = {}
    for law in LAWS:
        r10 = step3[law]["10.0"]["p95_8.7pct"]["module"]
        r5 = step3[law]["5.0"]["p95_8.7pct"]["module"]
        first = (r10 == "none up to 2C") or (order.index(r10) >= order.index("1C"))
        second = (r10 != "none up to 2C") and (r5 != "none up to 2C") and (order.index(r5) == order.index(r10) - 1)
        P4[law] = {"h10": r10, "h5": r5, "first_clause": first, "second_clause": second}
    P4["held"] = all(P4[l]["first_clause"] and P4[l]["second_clause"] for l in LAWS)
    ratios = {n: abs(dTf["law_1_kucinskis2022"][n]["p95_8.7pct"]) / abs(dTf["law_2_schimpe2018"][n]["p95_8.7pct"]) for n in ds_rates}
    P5 = {"ratio": ratios, "held": all(max(v, 1 / v) < 2 for v in ratios.values())}
    out = {"s_per_K": s, "dT_floor_K": dTf, "rows": rows, "step3_lowest_rate": step3, "dataset_conditions_brief_heat": ds,
           "D19": d19, "scorecard": {"P1": P1, "P2": P2, "P3": P3, "P4": P4, "P5": P5},
           "cell_values_from": "resolved_cell.json all_cases (finite-volume grid " +
                               f"{RC['all_cases'][0]['nr']} x {RC['all_cases'][0]['nz']})"}
    save("scoring.json", out)
    for k, v in out["scorecard"].items():
        print(k, v)
    print("step3", step3)
    print("dataset", ds)


if __name__ == "__main__":
    main()
