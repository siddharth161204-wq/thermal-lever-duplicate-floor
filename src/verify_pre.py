#!/usr/bin/env python3
"""V1 and V2, run before the predictions are written (DECISIONS.md D15). Neither computes a lever at 35 C or a module
spread at a registered case.

V1  Each ageing law against what its paper reports.
    Law 1 (Kucinskis et al. 2022): the law is the paper's printed table (Fig. 7), so V1 checks the transcription:
      the cycle counts the text quotes (p. 7) must equal the transcribed cells, the text's two qualitative statements
      about where the row maxima sit must hold, and cycles must not decrease from 95 to 80 % SoH in any cell.
    Law 2 (Schimpe et al. 2018): the model, with the Table IV parameters, at 1C CC cycling (charge and discharge
      at 3 A, full range, no high-SOC term, as the paper validates its CC tests) after 2800 full equivalent cycles,
      against the text (p. A183): "less than 9%" at 25 C and "12%" at 45 C. Acceptance: within the paper's own
      model error for CC tests, 2 percentage points (p. A190), i.e. below 11 % at 25 C and 10 to 14 % at 45 C.
V2  The module network.
    (a) The energy balance closes to round-off (relative residual below 1e-9) on non-registered inputs: random heat
        per node and h = 17.3 W/(m2 K), in both arrangements; and on every registered rate and h with the strip
        conductance set to 1e3 W/K.
    (b) With the strip conductance at 1e3 W/K, all eight nodes agree within 0.01 K at the largest heat (2C) and every h.
Usage: python3 src/verify_pre.py --out verification
"""
import argparse, json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import thermal_lever as TL


def v1_law1():
    k = TL.KUC
    quotes = [  # (soh, rate, T, value quoted in the text, p. 7)
        (0.80, 0.2, 25.0, 662), (0.80, 0.6, 25.0, 115), (0.80, 0.6, 45.0, 387)]
    checks = [{"soh": s, "charge_rate": c, "T_C": t, "text": v, "table": k[(s, c)].get(t), "equal": k[(s, c)].get(t) == v}
              for s, c, t, v in quotes]
    half = k[(0.80, 0.2)][15.0] / k[(0.80, 0.2)][25.0]          # "only approximately half of the cycles ... at 15 C"
    row95 = k[(0.95, 0.4)]
    row90 = k[(0.90, 0.4)]
    statements = {
        "15C_is_about_half_of_25C_at_0.2C_to_80pct": {"ratio": half, "holds": 0.4 <= half <= 0.6},
        "0.4C_to_95pct_most_cycles_at_25C": {"argmax_T": max(row95, key=row95.get), "holds": max(row95, key=row95.get) == 25.0},
        "0.4C_to_90pct_more_cycles_at_35C": {"argmax_T": max(row90, key=row90.get), "holds": max(row90, key=row90.get) == 35.0},
    }
    mono = []
    for c in (0.2, 0.4, 0.6):
        for t in (-10.0, 0.0, 15.0, 25.0, 35.0, 45.0, 55.0):
            seq = [k[(s, c)].get(t) for s in (0.95, 0.90, 0.85, 0.80)]
            if None in seq:
                continue
            mono.append({"rate": c, "T_C": t, "cycles_95_90_85_80": seq, "non_decreasing": all(b >= a for a, b in zip(seq, seq[1:]))})
    ok = all(x["equal"] for x in checks) and all(v["holds"] for v in statements.values()) and all(m["non_decreasing"] for m in mono)
    ties = [m for m in mono if any(b == a for a, b in zip(m["cycles_95_90_85_80"], m["cycles_95_90_85_80"][1:]))]
    return {"quoted_values": checks, "statements": statements, "monotone_cells": len(mono),
            "cells_with_equal_successive_values": ties, "pass": ok}


def v1_law2():
    p = TL.schimpe_params("table")
    out = {}
    for T_C, stated, lo, hi in ((25.0, "less than 9%", 0.0, 0.11), (45.0, "12%", 0.10, 0.14)):
        a, b, parts, _ = TL.schimpe_cycle_coeffs(T_C + 273.15, 1.0, 1.0, p, terms=("cal", "high", "low"))
        N = 2800.0
        loss = a * math.sqrt(N)
        out[f"{T_C:.0f}C"] = {"stated": stated, "model_loss_pct": 100 * loss,
                             "parts_pct": {k: 100 * v * math.sqrt(N) for k, v in parts.items()},
                             "accept_range_pct": [100 * lo, 100 * hi], "pass": lo <= loss <= hi}
    pt = TL.schimpe_params("text")
    for T_C in (25.0, 45.0):
        a, _, _, _ = TL.schimpe_cycle_coeffs(T_C + 273.15, 1.0, 1.0, pt, terms=("cal", "high", "low"))
        out[f"{T_C:.0f}C"]["model_loss_pct_text_values"] = 100 * a * math.sqrt(2800.0)
    out["pass"] = all(v["pass"] for k, v in out.items() if k.endswith("C"))
    out["note"] = ("1 full equivalent cycle = 6 Ah throughput for the 3 Ah cell (p. A185); cycle time 2 h at 1C; calendar factor "
                   "averaged over a uniform SOC from 0 to 1; ambient temperatures, without the cells' self-heating")
    return out


def v2():
    rng = np.random.default_rng(20261010)
    res = {"a_random_inputs": [], "a_registered_with_strip_1e3": [], "b_strip_1e3_uniformity": []}
    for base in (False, True):
        pairs = TL.NEIGHBOURS_ROW if base else TL.NEIGHBOURS
        q = rng.uniform(0.01, 2.0, 8)
        T, bal, info = TL.module_steady(q, 17.3, pairs=pairs, base=base)
        res["a_random_inputs"].append({"arrangement": "1x8 cold plate" if base else "2x4 air", "sources_W": info["sources_W"],
                                       "balance_W": bal, "rel": bal / info["sources_W"], "iterations": info["iterations"]})
    for c in TL.CRATES:
        for h in TL.H_LIST:
            T, bal, info = TL.module_steady(TL.heat_per_cell(c), h, g_strip=1e3)
            res["a_registered_with_strip_1e3"].append({"crate": c, "h": h, "rel": bal / info["sources_W"], "iterations": info["iterations"]})
    for h in TL.H_LIST:
        T, bal, info = TL.module_steady(TL.heat_per_cell(2.0), h, g_strip=1e3)
        res["b_strip_1e3_uniformity"].append({"h": h, "q_cell_W": TL.heat_per_cell(2.0), "max_minus_min_K": float(T.max() - T.min())})
    res["a_pass"] = all(abs(x["rel"]) < 1e-9 for x in res["a_random_inputs"] + res["a_registered_with_strip_1e3"])
    res["b_pass"] = all(x["max_minus_min_K"] < 0.01 for x in res["b_strip_1e3_uniformity"])
    res["pass"] = res["a_pass"] and res["b_pass"]
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default="verification"); a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    out = {"V1_law1_kucinskis2022": v1_law1(), "V1_law2_schimpe2018": v1_law2(), "V2_module": v2()}
    with open(os.path.join(a.out, "v1_v2.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
