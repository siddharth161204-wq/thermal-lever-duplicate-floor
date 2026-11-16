#!/usr/bin/env python3
"""V3 of Stage 6 (DECISIONS.md D96 (1)): two sets by hand, in plain Python with nothing of src/stage6.

The first set of two or more cells of 2017-05-12 and the first of 2019-01-24, in the order of results/stage6_sets.csv.
For each: every cell's T (the mean of Tavg over cycles 2 to 100, from the committed extract, with the site's swap of
channels 15 and 16 where it applies) and ln L (from inputs/matr/stage6_lifetimes.csv), the set's means, deviations,
sums of squares and products, its within-set slope, its pairs with dT and d, and whether its warmest cell is the first
to die; written out, and compared with results/stage6_cells.csv to 1e-12 relative.

Usage: python3 -I verification/stage6/v3_hand.py
"""
import csv, json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
LABEL = {"2017-05-12": "b1", "2017-06-30": "b2", "2018-04-12": "b3", "2019-01-24": "b9"}


def first_set(batch):
    for r in csv.DictReader(open(os.path.join(ROOT, "results", "stage6_sets.csv"))):
        if r["batch"] == batch and int(r["cells"]) >= 2:
            return r
    return None


def main():
    cells = {r["cell"]: r for r in csv.DictReader(open(os.path.join(ROOT, "results", "stage6_cells.csv")))}
    life = {r["cell"]: float(r["cycle_life"]) for r in csv.DictReader(open(os.path.join(ROOT, "inputs", "matr",
                                                                                      "stage6_lifetimes.csv")))}
    out, ok = {"check": "V3 (D96 (1))", "sets": []}, True
    for batch in ("2017-05-12", "2019-01-24"):
        s = first_set(batch)
        lab = LABEL[batch]
        desc = {int(r["index"]): r for r in csv.DictReader(open(os.path.join(ROOT, "inputs", "matr", "a0",
                                                                              f"{lab}_a0_cells.csv")))}
        cyc = {}
        for r in csv.DictReader(open(os.path.join(ROOT, "inputs", "matr", "a0", f"{lab}_a0_cycles.csv"))):
            cyc.setdefault(int(r["index"]), {})[int(r["cycle"])] = r
        members = [c for c in cells.values() if c["set"] == s["set"]]
        rows = []
        for c in members:
            i = int(c["index"])
            ch = desc[i]["channel_id"].strip()
            src = i
            if batch in ("2017-05-12", "2017-06-30") and ch in ("15", "16"):     # the site's swap (D87 (3))
                other = "16" if ch == "15" else "15"
                src = next(k for k, d in desc.items() if d["channel_id"].strip() == other)
            tav = [float(cyc[src][n]["Tavg"]) for n in range(2, 101)]
            T = sum(tav) / 99.0
            rows.append({"cell": c["cell"], "channel": ch, "T": T, "lnL": math.log(life[c["cell"]]),
                         "T_pipeline": float(c["T_reg"])})
        n = len(rows)
        Tm = sum(r["T"] for r in rows) / n
        ym = sum(r["lnL"] for r in rows) / n
        Sxx = sum((r["T"] - Tm) ** 2 for r in rows)
        Sxy = sum((r["T"] - Tm) * (r["lnL"] - ym) for r in rows)
        pairs = []
        for a in range(n):
            for b in range(a + 1, n):
                w, c_ = (rows[a], rows[b]) if rows[a]["T"] > rows[b]["T"] else (rows[b], rows[a])
                pairs.append({"warm": w["cell"], "cool": c_["cell"], "dT": w["T"] - c_["T"], "d": w["lnL"] - c_["lnL"]})
        warm = max(rows, key=lambda r: r["T"])
        first = min(rows, key=lambda r: r["lnL"])
        agree = all(abs(r["T"] - r["T_pipeline"]) <= 1e-12 * abs(r["T_pipeline"]) for r in rows)
        ok &= agree
        out["sets"].append({"set": s["set"], "cells": rows, "T_mean": Tm, "lnL_mean": ym, "S_TT": Sxx, "S_Ty": Sxy,
                            "within_slope": Sxy / Sxx, "set_s_w_per_K": -Sxy / Sxx, "pairs": pairs,
                            "pairs_check": [abs(sum(p["dT"] * p["d"] for p in pairs) - n * Sxy),
                                            abs(sum(p["dT"] ** 2 for p in pairs) - n * Sxx)],
                            "warmest": warm["cell"], "first_to_die": first["cell"],
                            "warmest_first": warm["cell"] == first["cell"], "T_agrees_with_pipeline": agree,
                            "S_TT_pipeline_check": abs(float(s["T_reg_sd_within_K"]) ** 2 * (n - 1) - Sxx)})
        print(s["set"], n, "cells; T agrees:", agree, "; set s_w", -Sxy / Sxx, "; warmest first:", warm["cell"] == first["cell"])
    out["pass"] = bool(ok)
    json.dump(out, open(os.path.join(HERE, "v3.json"), "w"), indent=1)
    print("V3 pass:", ok)


if __name__ == "__main__":
    main()
