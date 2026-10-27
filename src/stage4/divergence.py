#!/usr/bin/env python3
"""Stage 4: the implied divergence and the map (DECISIONS.md D60), read from results/stage4_module.csv and the
archived levers and floors.

For a lever s and a case: each cell's implied lifetime offset from the module mean is -s (T_i - T_mean) in ln L;
the module's implied divergence is |s| x spread, with the interval spread x [least, greatest] |s| over the lever's
95 % interval. Per case, the median and the 95th percentile (numpy's linear interpolation) over the band's levers; a
case offers a thermal share under a floor when the median exceeds the floor's mean, a detectable one when the median
exceeds the floor's 95th percentile. The map: per ambient, cooling and floor, the smallest rate of each.

Usage: python3 src/stage4/divergence.py
"""
import argparse, csv, json, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as CM                     # noqa: E402


def abs_interval(lo, hi):
    """Least and greatest |s| over [lo, hi]."""
    if lo <= 0.0 <= hi:
        return 0.0, max(abs(lo), abs(hi))
    return min(abs(lo), abs(hi)), max(abs(lo), abs(hi))


def read_module(path):
    rows = []
    for r in csv.DictReader(open(path)):
        d = {k: (float(v) if k not in ("band", "rate", "hot_position", "cold_position", "converged",
                                       "cells_inside_band_setpoints") else v) for k, v in r.items()}
        rows.append(d)
    return rows


def case_key(r):
    return (r["ambient_C"], r["crate"], r["h"])


def divergence(module_rows, d):
    bands, floors = d["bands"], d["floors"]
    lev_rows, cases = [], []
    for m in module_rows:
        band = bands[m["band"]]
        spread = m["spread_K"]
        divs = []
        for lev in band:
            s = lev["s"]
            lo_abs, hi_abs = abs_interval(lev["lo"], lev["hi"])
            row = {"ambient_C": m["ambient_C"], "crate": m["crate"], "rate": m["rate"], "h": m["h"], "band": m["band"],
                   "row_kind": "lever", "lever": lev["label"], "source": lev["source"], "secant": lev["secant"],
                   "s_per_K": s, "s_lo": lev["lo"], "s_hi": lev["hi"], "spread_K": spread,
                   "divergence": abs(s) * spread, "divergence_lo": lo_abs * spread, "divergence_hi": hi_abs * spread}
            for i in range(8):
                row[f"offset_lnL_cell{i}"] = -s * m[f"offset_cell{i}_K"]
            lev_rows.append(row)
            divs.append(abs(s) * spread)
        divs = np.array(divs)
        med = float(np.median(divs))
        p95 = float(np.percentile(divs, 95))
        for kind, val in (("band median", med), ("band 95th percentile", p95)):
            lev_rows.append({"ambient_C": m["ambient_C"], "crate": m["crate"], "rate": m["rate"], "h": m["h"],
                             "band": m["band"], "row_kind": kind, "lever": kind, "spread_K": spread, "divergence": val})
        c = {"ambient_C": m["ambient_C"], "crate": m["crate"], "rate": m["rate"], "h": m["h"], "band": m["band"],
             "spread_K": spread, "T_min_C": m["T_min_C"], "T_max_C": m["T_max_C"], "T_mean_C": m["T_mean_C"],
             "cells_inside_band_setpoints": m["cells_inside_band_setpoints"] == "True",
             "levers": len(band), "median": med, "p95": p95, "max": float(divs.max()), "min": float(divs.min()),
             "median_lever": CM.median_lever(band)["label"]}
        for fname, f in floors.items():
            c[f"offers_{fname}"] = bool(med > f["mean"])
            c[f"detectable_{fname}"] = bool(med > f["p95"])
            c[f"levers_above_mean_{fname}"] = int((divs > f["mean"]).sum())
            c[f"levers_above_p95_{fname}"] = int((divs > f["p95"]).sum())
        cases.append(c)
    return lev_rows, cases


def the_map(cases, floors):
    out = {}
    for fname in floors:
        out[fname] = {}
        for T in CM.AMBIENTS_C:
            out[fname][f"{T:g} C"] = {}
            for h in CM.H_LIST:
                cc = sorted([c for c in cases if c["ambient_C"] == T and c["h"] == h], key=lambda c: c["crate"])
                share = [c["rate"] for c in cc if c[f"offers_{fname}"]]
                det = [c["rate"] for c in cc if c[f"detectable_{fname}"]]
                out[fname][f"{T:g} C"][f"h = {h:g}"] = {
                    "smallest_rate_offering_a_share": share[0] if share else "none up to 2C",
                    "smallest_rate_offering_a_detectable_share": det[0] if det else "none up to 2C",
                    "rates_offering": share, "rates_detectable": det}
    return out


def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main(module_path, out_dir):
    d = CM.load_all()
    mod = read_module(module_path)
    lev_rows, cases = divergence(mod, d)
    write_csv(os.path.join(out_dir, "stage4_divergence.csv"), lev_rows)
    doc = {"definition": "D60: per case, median and 95th percentile of |s| x spread over the band's levers; offers when the "
                         "median exceeds the floor's mean, detectable when it exceeds the floor's 95th percentile",
           "floors": d["floors"],
           "median_levers": {"warm": d["named"]["median_warm"]["label"], "warm_s": d["named"]["median_warm"]["s"],
                             "cold": d["named"]["median_cold"]["label"], "cold_s": d["named"]["median_cold"]["s"]},
           "cases": cases, "map": the_map(cases, d["floors"])}
    json.dump(doc, open(os.path.join(out_dir, "stage4_map.json"), "w"), indent=1)
    for c in cases:
        if c["h"] == 10.0:
            print(f"{c['ambient_C']:>4g} C {c['rate']:>3} h 10: median {100*c['median']:.3f} %, p95 {100*c['p95']:.3f} %")
    print(json.dumps(doc["map"], indent=1)[:3000])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--module", default=os.path.join(CM.ROOT, "results", "stage4_module.csv"))
    ap.add_argument("--out", default=os.path.join(CM.ROOT, "results"))
    a = ap.parse_args()
    main(a.module, a.out)
