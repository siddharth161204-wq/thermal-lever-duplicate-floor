#!/usr/bin/env python3
"""Stage 4: every number the note and the README quote, read from the results files (results/stage4_paper_numbers.json).

Usage: python3 src/stage4/paper_numbers.py
"""
import csv, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as CM                     # noqa: E402

RES = os.path.join(CM.ROOT, "results")
COMMITS = {"predictions": {"commit": "c1e58a19d1653bb477de485c4fc9f68e956ee070", "pushed_utc": "2026-10-27T09:50:05Z"},
           "code": {"commit": "e5ecf58a971ce153a0ea7a94ffedfb1282db8d77"},
           "power_fix_D64": {"commit": "aa74702183cb1ed324cbf1a9933a1a539ab50401"},
           "results": {"commit": "6f658f1fe534c8d9a158308dc1ff90eb12fd0318", "pushed_utc": "2026-10-27T10:07:16Z"}}


def rng(vals, k=3):
    return [round(min(vals), k), round(max(vals), k)]


def main():
    res = json.load(open(os.path.join(RES, "stage4_resistance.json")))
    mod = list(csv.DictReader(open(os.path.join(RES, "stage4_module.csv"))))
    rc = list(csv.DictReader(open(os.path.join(RES, "stage4_resolved_cell.csv"))))
    mp = json.load(open(os.path.join(RES, "stage4_map.json")))
    ex = list(csv.DictReader(open(os.path.join(RES, "stage4_experiment.csv"))))
    sc = json.load(open(os.path.join(RES, "stage4_scorecard.json")))
    ph = json.load(open(os.path.join(RES, "stage4_posthoc.json")))
    v = {f"v{i}": json.load(open(os.path.join(CM.ROOT, "verification", "stage4", f"v{i}.json"))) for i in (1, 2, 3, 4)}

    out = {"commits": COMMITS, "scorecard": sc["_summary"]}
    out["resistance_mOhm"] = {k: 1e3 * x for k, x in res["R_used_ohm"].items()}
    out["resistance_ratio"] = res["ratio_to_35C"]
    sp = {}
    for T in ("35", "25", "10"):
        for rate in ("C/4", "C/2", "1C", "2C"):
            rows = [m for m in mod if f"{float(m['ambient_C']):g}" == T and m["rate"] == rate]
            sp[f"{T} C {rate}"] = {"h10_K": float([m for m in rows if float(m["h"]) == 10.0][0]["spread_K"]),
                                   "range_over_h_K": rng([float(m["spread_K"]) for m in rows]),
                                   "cells_C": rng([float(m["T_min_C"]) for m in rows] + [float(m["T_max_C"]) for m in rows], 1)}
    out["spread"] = sp
    out["resolved_10C"] = {r["rate"]: {"within_cell_K": float(r["dT_cell_K"]), "module_K": float(r["module_spread_K"]),
                                       "ratio": float(r["within_over_module"])} for r in rc if float(r["ambient_C"]) == 10.0}
    med = {}
    for T in (35.0, 25.0, 10.0):
        for rate in ("C/4", "C/2", "1C", "2C"):
            cc = [c for c in mp["cases"] if c["ambient_C"] == T and c["rate"] == rate]
            med[f"{T:g} C {rate}"] = {"median_pct_range_over_h": rng([100 * c["median"] for c in cc], 2),
                                      "p95_pct_range_over_h": rng([100 * c["p95"] for c in cc], 2),
                                      "max_lever_pct_range_over_h": rng([100 * c["max"] for c in cc], 2),
                                      "levers_above_Geslin_p95": [c["levers_above_p95_Geslin"] for c in sorted(cc, key=lambda c: c["h"])]}
    out["median_divergence"] = med
    out["median_levers_pct_per_K"] = {"warm": [mp["median_levers"]["warm"], 100 * mp["median_levers"]["warm_s"]],
                                      "cold": [mp["median_levers"]["cold"], 100 * mp["median_levers"]["cold_s"]]}
    out["map"] = {f: {T: {h: [x["smallest_rate_offering_a_share"], x["smallest_rate_offering_a_detectable_share"]] for h, x in hh.items()}
                      for T, hh in tt.items()} for f, tt in mp["map"].items()}

    def ex_get(part, lever, dT, floor, chamber=""):
        r = [e for e in ex if e["part"] == part and e["lever"] == lever and float(e["dT_K"]) == dT and e["floor"] == floor
             and e["chamber"].startswith(chamber)]
        if len(r) != 1:
            raise SystemExit(f"no unique experiment row: {part} {chamber} {lever} {dT} {floor}")
        return {"n": int(r[0]["n_pairs"]), "n_sim_lognormal": int(r[0]["n_sim_lognormal"]), "n_normal_shape": int(r[0]["n_pairs_normal_shape"])}
    w, c = mp["median_levers"]["warm"], mp["median_levers"]["cold"]
    out["experiment"] = {
        "warm_central": {f"{dT:g} K {f}": ex_get("C2", w, dT, f, "warm") for dT in (2.0, 5.0) for f in ("Geslin", "KIT")},
        "cold_central": {f"{dT:g} K {f}": ex_get("C2", c, dT, f, "cold") for dT in (2.0, 5.0) for f in ("Geslin", "KIT")},
        "c1": {lev: {f"{dT:g} K {f}": ex_get("C1", lev, dT, f) for dT in (1.0, 2.0, 5.0, 10.0) for f in ("Geslin", "KIT")}
               for lev in ("KIT 10-100 1.00/1.00C", "Sandia NMC 25 to 35 C", "KIT 10 to 25 C median", "KIT 10-90 1.67/1.00C")},
        "warm_band_range_2K_Geslin": rng([int(e["n_pairs"]) for e in ex if e["part"] == "C2" and e["chamber"].startswith("warm")
                                          and float(e["dT_K"]) == 2.0 and e["floor"] == "Geslin"], 0),
        "cold_band_range_5K_Geslin": rng([int(e["n_pairs"]) for e in ex if e["part"] == "C2" and e["chamber"].startswith("cold")
                                          and float(e["dT_K"]) == 5.0 and e["floor"] == "Geslin"], 0),
        "cold_band_range_2K_Geslin": rng([int(e["n_pairs"]) for e in ex if e["part"] == "C2" and e["chamber"].startswith("cold")
                                          and float(e["dT_K"]) == 2.0 and e["floor"] == "Geslin"], 0)}
    out["checks"] = {"V1": {"pass": v["v1"]["pass"], "bitwise": v["v1"]["all_bitwise_identical"]},
                     "V2": {"pass": v["v2"]["pass"], "max_abs_rel": max(abs(a["rel_diff"]) for a in v["v2"]["ambients"].values())},
                     "V3": {"pass": v["v3"]["pass"]},
                     "V4": ph["v4"]}
    out["posthoc"] = {"cold_bracket_p3_ratios": ph["cold_bracket"]["p3_ratio_with"], "margins": ph["margins"],
                      "signed_median": {k: ph["signed_median"][k] for k in ("signed_median_s", "lever", "warm_chamber_dT2_Geslin_n",
                                                                             "warm_chamber_dT5_Geslin_n")}}
    cb = [x for x in ph["cold_bracket"]["cases"] if x["ambient_C"] == 10.0]
    out["posthoc"]["cold_bracket_1C"] = {f"R({k} C)": {"spread_K": rng([x["spread_K"] for x in cb if x["resistance_of_C"] == float(k) and x["rate"] == "1C"]),
                                                       "median_pct": rng([100 * x["median_divergence"] for x in cb if x["resistance_of_C"] == float(k) and x["rate"] == "1C"], 2),
                                                       "cells_C": rng([x["T_min_C"] for x in cb if x["resistance_of_C"] == float(k) and x["rate"] == "1C"]
                                                                      + [x["T_max_C"] for x in cb if x["resistance_of_C"] == float(k) and x["rate"] == "1C"], 1)}
                                         for k in ("25", "35")}
    own = json.load(open(os.path.join(RES, "stage4_posthoc_own_temperature.json")))
    oc = own["cases"]
    pick = lambda T, rate, key: rng([x[key] for x in oc if x["ambient_C"] == T and x["rate"] == rate], 4)
    out["posthoc"]["own_temperature"] = {
        "R_mOhm_at": {k: 1e3 * v for k, v in own["R_of_T_ohm"].items() if k in ("10", "35", "50", "90")},
        "spread_K": {"35 C 1C": pick(35.0, "1C", "own_spread_K"), "35 C 2C": pick(35.0, "2C", "own_spread_K"),
                     "10 C 1C": pick(10.0, "1C", "own_spread_K"), "10 C 2C": pick(10.0, "2C", "own_spread_K")},
        "median": {"35 C 2C": pick(35.0, "2C", "own_median"), "10 C 1C": pick(10.0, "1C", "own_median")},
        "max_lever_35C_2C_h5": [x["own_max_lever"] for x in oc if x["ambient_C"] == 35.0 and x["rate"] == "2C" and x["h"] == 5.0][0],
        "verdicts": own["verdicts"]}
    ind = json.load(open(os.path.join(CM.ROOT, "verification", "stage4", "independent", "record.json")))
    out["independent"] = {"freeze_commit": "fd5d1e73f6fa87e5b60cc9c2689b26e072a6b1d0", "frozen_utc": "2026-10-27T10:33:27Z",
                          "pushed_utc": "2026-10-27T10:33:36Z", "comparison_commit": "2c17f4d",
                          "largest_deterministic_difference": "2e-12 (power at n 1.7e-12; divergences 8.8e-13 relative; module 1.7e-13 K)",
                          "verdicts_identical": True, "record": "verification/stage4/independent/record.json", "freeze_in_record": ind["freeze"]["commit"]}
    json.dump(out, open(os.path.join(RES, "stage4_paper_numbers.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
