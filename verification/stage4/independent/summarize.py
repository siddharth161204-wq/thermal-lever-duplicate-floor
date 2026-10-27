#!/usr/bin/env python3
"""Collects the checker's frozen numbers into out/summary.json and writes out/manifest.json (SHA-256 of every output).
Run last, before the freeze commit.
"""
import csv, glob, hashlib, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
sys.path.insert(0, HERE)
import levers as LV                                 # noqa: E402


def rd(name):
    with open(os.path.join(OUT, name)) as f:
        return list(csv.DictReader(f))


def main():
    warm, cold, floors = LV.load()
    res = json.load(open(os.path.join(OUT, "resistance.json")))
    chk = json.load(open(os.path.join(OUT, "resistance_check.json")))
    cases = rd("divergence_cases.csv")
    span = {"warm": (min(lv["T_lo_C"] for lv in warm), max(lv["T_hi_C"] for lv in warm)),
            "cold": (min(lv["T_lo_C"] for lv in cold), max(lv["T_hi_C"] for lv in cold))}
    table = []
    for r in cases:
        lo, hi = span[r["band"]]
        table.append({"case": f"{float(r['T_amb_C']):g} C {r['rate']} h={float(r['h']):g}", "spread_K": float(r["spread_K"]),
                      "T_min_C": float(r["T_min_C"]), "T_max_C": float(r["T_max_C"]),
                      "median_pct": 100 * float(r["median"]), "p95_pct": 100 * float(r["p95"]),
                      "offers_Geslin": r["offers_Geslin"] == "True", "detectable_Geslin": r["detectable_Geslin"] == "True",
                      "offers_KIT": r["offers_KIT"] == "True", "detectable_KIT": r["detectable_KIT"] == "True",
                      "cells_inside_band_secant_span": float(r["T_min_C"]) >= lo and float(r["T_max_C"]) <= hi,
                      "band_secant_span_C": [lo, hi]})
    exp = rd("experiment.csv")
    sim = {(int(x["case_index"])): x for x in rd("simulation_lognormal_all.csv")}
    simn = {(int(x["case_index"])): x for x in rd("simulation_normal_all.csv")}
    e_rows = []
    for i, r in enumerate(exp):
        s = sim.get(i) if r["shape"] == "lognormal" else simn.get(i)
        e_rows.append({"group": r["group"], "lever": r["lever"], "dT_K": float(r["dT_K"]), "floor": r["floor"], "shape": r["shape"],
                       "n": int(r["n"]), "n_sim_10k": int(s["n_sim"]) if s else None,
                       "V4_pass_5pct": (s["V4_pass_5pct"] == "True") if s else None, "central": r["central"] == "True"})
    summ = {"freeze_note": "written before opening any Stage 4 pipeline file",
            "resistance_ohm": {k: v["R_registered_ohm"] for k, v in res["ambients"].items()},
            "resistance_ratios": {"25/35": res["used"]["ratio_25_over_35_archived"], "10/35": res["used"]["ratio_10_over_35_archived"]},
            "resistance_failed_calls": {k: v["failed_calls"] for k, v in res["ambients"].items()},
            "resistance_separate_code_spike_screened_ohm": {str(r["T_C"]): r["R_ohm_spike_screened"] for r in chk["rows"]},
            "resistance_separate_code_unscreened_ohm": {str(r["T_C"]): r["R_ohm"] for r in chk["rows"]},
            "cases": table, "map": json.load(open(os.path.join(OUT, "map.json")))["map"],
            "experiment": e_rows, "scorecard": json.load(open(os.path.join(OUT, "scorecard.json"))),
            "resolved_cell": rd("resolved_cell.csv"),
            "sensitivity_own_temperature_verdicts": json.load(open(os.path.join(OUT, "sensitivity_own_temperature.json")))["verdicts_under_this_heat"]}
    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(summ, f, indent=1)
    man = {}
    for p in sorted(glob.glob(os.path.join(OUT, "*")) + glob.glob(os.path.join(HERE, "*.py"))):
        if os.path.basename(p) == "manifest.json":
            continue
        man[os.path.relpath(p, HERE)] = hashlib.sha256(open(p, "rb").read()).hexdigest()
    with open(os.path.join(OUT, "manifest.json"), "w") as f:
        json.dump(man, f, indent=1)
    print(len(man), "files hashed")
    print(json.dumps({k: summ[k] for k in ("resistance_ohm", "resistance_ratios", "resistance_failed_calls",
                                           "resistance_separate_code_spike_screened_ohm")}, indent=1))
    print({k: v["held"] for k, v in summ["scorecard"].items()})


if __name__ == "__main__":
    main()
