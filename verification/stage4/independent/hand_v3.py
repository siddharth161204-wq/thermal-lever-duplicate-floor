#!/usr/bin/env python3
"""A second, plain-Python implementation (no numpy) of the two cases D62's V3 fixes: 35 C at 1C and h = 10, and
10 C at 2C and h = 5. Every lever's |s| x spread, then the band's median and 95th percentile with numpy's linear
interpolation written out: sorted values v_0..v_{k-1}, position p = 0.95 (k - 1), v_floor(p) + (p - floor(p)) x
(v_floor(p)+1 - v_floor(p)). Compared with out/divergence_cases.csv.
Output: out/hand_v3.json
"""
import csv, json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
sys.path.insert(0, HERE)
import levers as LV                                 # noqa: E402


def pct_linear(vals, p):
    v = sorted(vals)
    pos = p / 100.0 * (len(v) - 1)
    i = int(math.floor(pos))
    f = pos - i
    if i + 1 >= len(v):
        return v[-1], {"position": pos, "i": i, "fraction": f}
    return v[i] + f * (v[i + 1] - v[i]), {"position": pos, "i": i, "fraction": f, "v_i": v[i], "v_i+1": v[i + 1]}


def main():
    warm, cold, floors = LV.load()
    with open(os.path.join(OUT, "module.csv")) as f:
        mod = list(csv.DictReader(f))
    with open(os.path.join(OUT, "divergence_cases.csv")) as f:
        cases = list(csv.DictReader(f))
    out = {}
    for Ta, c, h, band in ((35.0, 1.0, 10.0, warm), (10.0, 2.0, 5.0, cold)):
        m = [r for r in mod if float(r["T_amb_C"]) == Ta and float(r["crate"]) == c and float(r["h"]) == h][0]
        spread = float(m["spread_K"])
        per = [(lv["lever"], abs(lv["s"]) * spread) for lv in band]
        vals = [x for _, x in per]
        med, med_d = pct_linear(vals, 50)
        p95, p95_d = pct_linear(vals, 95)
        reg = [r for r in cases if float(r["T_amb_C"]) == Ta and float(r["crate"]) == c and float(r["h"]) == h][0]
        out[f"{Ta:g} C, {c:g}C, h = {h:g}"] = {"spread_K": spread, "per_lever": per, "median": med, "median_detail": med_d,
                                               "p95": p95, "p95_detail": p95_d,
                                               "numpy_path_median": float(reg["median"]), "numpy_path_p95": float(reg["p95"]),
                                               "median_diff": med - float(reg["median"]), "p95_diff": p95 - float(reg["p95"])}
    with open(os.path.join(OUT, "hand_v3.json"), "w") as f:
        json.dump(out, f, indent=1)
    for k, v in out.items():
        print(k, "median", v["median"], "p95", v["p95"], "diffs", v["median_diff"], v["p95_diff"])


if __name__ == "__main__":
    main()
