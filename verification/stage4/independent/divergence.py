#!/usr/bin/env python3
"""Independent check of Stage 4, part (d): the implied divergence per D60, the per-case median and 95th percentile
over the band, offers / detectable under both floors, and the map.

Per case (ambient, rate, h) and per lever of the band the ambient picks (35 and 25 C: warm; 10 C: cold):
  implied lifetime offset of cell i from the module mean, in ln L: -s (T_i - T_mean)
  implied divergence: |s| x spread (first-order pair difference); exact pair difference 2 tanh(|s| x spread / 2)
  interval: spread x [least, greatest] |s| over the lever's 95 % interval, zero at the low end where it holds zero
Per case over the band: median and 95th percentile (numpy linear interpolation) of |s| x spread; the median lever
(median of |s|, D57); offers = median > floor mean; detectable = median > floor 95th percentile; both floors.
Map: for each ambient, cooling and floor, the smallest rate that offers, and the smallest that is detectable.

Usage: python3 verification/stage4/independent/divergence.py
Output: out/divergence.csv (case x lever), out/divergence_cases.csv (per case), out/map.json
"""
import csv, json, math, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
sys.path.insert(0, HERE)
import levers as LV                                 # noqa: E402

RATE_NAME = {0.25: "C/4", 0.5: "C/2", 1.0: "1C", 2.0: "2C"}


def abs_interval(lo, hi):
    if lo <= 0.0 <= hi:
        return 0.0, max(abs(lo), abs(hi))
    return min(abs(lo), abs(hi)), max(abs(lo), abs(hi))


def main():
    warm, cold, floors = LV.load()
    with open(os.path.join(OUT, "module.csv")) as f:
        cases = list(csv.DictReader(f))
    rows, crow = [], []
    for c in cases:
        Ta, rate, h = float(c["T_amb_C"]), float(c["crate"]), float(c["h"])
        spread = float(c["spread_K"])
        T = np.array([float(c[f"T_cell{i}_C"]) for i in range(8)])
        Tmean = float(c["T_mean_C"])
        band = warm if Ta in (35.0, 25.0) else cold
        bname = "warm" if band is warm else "cold"
        divs = []
        for lv in band:
            a = abs(lv["s"])
            lo, hi = abs_interval(lv["lo"], lv["hi"])
            div = a * spread
            divs.append(div)
            r = {"T_amb_C": Ta, "crate": rate, "h": h, "band": bname, "lever": lv["lever"], "s_per_K": lv["s"],
                 "s_ci_lo": lv["lo"], "s_ci_hi": lv["hi"], "spread_K": spread,
                 "implied_divergence": div, "implied_divergence_lo": lo * spread, "implied_divergence_hi": hi * spread,
                 "exact_pair_difference": 2.0 * math.tanh(div / 2.0),
                 "T_min_C": float(c["T_min_C"]), "T_max_C": float(c["T_max_C"]),
                 "lever_secant_T_lo_C": lv["T_lo_C"], "lever_secant_T_hi_C": lv["T_hi_C"]}
            for i in range(8):
                r[f"ln_offset_cell{i}"] = -lv["s"] * (T[i] - Tmean)
            for fk, fv in floors.items():
                r[f"above_{fk}_mean"] = bool(div > fv["mean"])
                r[f"above_{fk}_p95"] = bool(div > fv["p95"])
            rows.append(r)
        divs = np.array(divs)
        med = float(np.median(divs))
        p95 = float(np.percentile(divs, 95))          # numpy's default, linear interpolation
        abs_s = np.array([abs(lv["s"]) for lv in band])
        order = np.argsort(abs_s, kind="stable")
        k = len(band)
        med_lever = band[order[k // 2]]["lever"] if k % 2 == 1 else None
        cr = {"T_amb_C": Ta, "crate": rate, "rate": RATE_NAME[rate], "h": h, "band": bname, "n_levers": k,
              "spread_K": spread, "T_min_C": float(c["T_min_C"]), "T_max_C": float(c["T_max_C"]),
              "T_mean_C": Tmean, "median": med, "p95": p95, "max": float(divs.max()), "min": float(divs.min()),
              "median_lever": med_lever, "median_abs_s": float(np.median(abs_s)),
              "median_equals_median_abs_s_times_spread": bool(med == float(np.median(abs_s)) * spread),
              "median_minus_median_abs_s_times_spread": med - float(np.median(abs_s)) * spread}
        for fk, fv in floors.items():
            cr[f"offers_{fk}"] = bool(med > fv["mean"])
            cr[f"detectable_{fk}"] = bool(med > fv["p95"])
            cr[f"levers_above_{fk}_mean"] = int((divs > fv["mean"]).sum())
            cr[f"levers_above_{fk}_p95"] = int((divs > fv["p95"]).sum())
        crow.append(cr)
    with open(os.path.join(OUT, "divergence.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: repr(v) if isinstance(v, float) else v for k, v in r.items()})
    with open(os.path.join(OUT, "divergence_cases.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(crow[0].keys()))
        w.writeheader()
        for r in crow:
            w.writerow({k: repr(v) if isinstance(v, float) else v for k, v in r.items()})
    # the map (D60, brief B3)
    mp = {}
    for Ta in (35.0, 25.0, 10.0):
        for h in (5.0, 10.0, 25.0):
            for fk in floors:
                sub = sorted([r for r in crow if r["T_amb_C"] == Ta and r["h"] == h], key=lambda r: r["crate"])
                off = [r["rate"] for r in sub if r[f"offers_{fk}"]]
                det = [r["rate"] for r in sub if r[f"detectable_{fk}"]]
                mono_off = [r[f"offers_{fk}"] for r in sub] == sorted([r[f"offers_{fk}"] for r in sub])
                mono_det = [r[f"detectable_{fk}"] for r in sub] == sorted([r[f"detectable_{fk}"] for r in sub])
                mp[f"{Ta:g} C, h = {h:g}, {fk}"] = {"ambient_C": Ta, "h": h, "floor": fk,
                                                     "smallest_rate_offers": off[0] if off else "none up to 2C",
                                                     "smallest_rate_detectable": det[0] if det else "none up to 2C",
                                                     "monotone_in_rate": bool(mono_off and mono_det)}
    with open(os.path.join(OUT, "map.json"), "w") as f:
        json.dump({"floors": floors, "map": mp}, f, indent=1)
    for r in crow:
        print(f"{r['T_amb_C']:4.0f} C {r['rate']:>3s} h={r['h']:<4g} spread={r['spread_K']:9.6f} K  median={100*r['median']:8.4f} %  "
              f"p95={100*r['p95']:8.4f} %  offers G/K={int(r['offers_Geslin'])}/{int(r['offers_KIT'])}  "
              f"detect G/K={int(r['detectable_Geslin'])}/{int(r['detectable_KIT'])}  med lever={r['median_lever']}")
    for k, v in mp.items():
        print(k, v["smallest_rate_offers"], v["smallest_rate_detectable"], v["monotone_in_rate"])


if __name__ == "__main__":
    main()
