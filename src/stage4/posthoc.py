#!/usr/bin/env python3
"""Stage 4, after the registered results (DECISIONS.md, Stage 4 after the results): conditional readings that change
no registered number or verdict. Output: results/stage4_posthoc.json.

1. The cold heat as a bracket: the module at 10 C ambient with the resistance of 25 C and of 35 C (both from the
   procedure of D58) in place of the 10 C value, and at 25 C ambient with the 35 C value; the cells' temperatures under
   the registered heat; what the cold verdicts would read with each.
2. The reading the author did not choose: the signed median of the warm band (-0.42 % per K) in place of the median
   of |s|, for the warm divergences and the warm chamber's pairs.
3. The resistance procedure's rates other than C/2: failures and the C/16 charge values.
4. Margins of the registered verdicts.

Usage: python3 src/stage4/posthoc.py
"""
import csv, json, math, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as CM                     # noqa: E402
import module as MD                     # noqa: E402
import experiment as EX                 # noqa: E402

RES = os.path.join(CM.ROOT, "results")


def cold_bracket(R, d):
    floors = d["floors"]
    s_cold = abs(d["named"]["median_cold"]["s"])
    s_warm = abs(d["named"]["median_warm"]["s"])
    out = {"note": "conditional: the registered 10 C heat uses R(10 C); these use R(25 C) or R(35 C) at 10 C ambient, "
                   "and R(35 C) at 25 C ambient", "cases": []}
    for T_C, key, s in ((10.0, "25", s_cold), (10.0, "35", s_cold), (25.0, "35", s_warm)):
        for c in CM.RATES:
            for h in CM.H_LIST:
                r = MD.module_case(T_C, c, h, R[key])
                med = s * r["spread_K"]
                out["cases"].append({"ambient_C": T_C, "resistance_of_C": float(key), "R_ohm": R[key], "rate": CM.RATE_LABEL[c],
                                     "h": h, "spread_K": r["spread_K"], "T_min_C": r["T_min_C"], "T_max_C": r["T_max_C"],
                                     "median_divergence": med,
                                     "offers_Geslin": med > floors["Geslin"]["mean"], "detectable_Geslin": med > floors["Geslin"]["p95"],
                                     "offers_KIT": med > floors["KIT"]["mean"], "detectable_KIT": med > floors["KIT"]["p95"]})
    reg = {(float(m["ambient_C"]), m["rate"], float(m["h"])): m for m in csv.DictReader(open(os.path.join(RES, "stage4_module.csv")))}
    out["p3_ratio_with"] = {}
    for key in ("25", "35"):
        out["p3_ratio_with"][f"R({key} C)"] = {
            f"h {h:g}": [x for x in out["cases"] if x["ambient_C"] == 10.0 and x["resistance_of_C"] == float(key)
                         and x["rate"] == "1C" and x["h"] == h][0]["spread_K"] / float(reg[(35.0, "1C", h)]["spread_K"]) for h in CM.H_LIST}
    out["registered_cell_temperatures_10C"] = {f"{r} h {h:g}": [float(reg[(10.0, r, h)]["T_min_C"]), float(reg[(10.0, r, h)]["T_max_C"])]
                                               for r in ("1C", "2C") for h in CM.H_LIST}
    return out


def signed_median(d, mp):
    warm = d["bands"]["warm"]
    s_signed = float(np.median([l["s"] for l in warm]))
    lev = [l for l in warm if l["s"] == s_signed]
    rows = []
    for c in mp["cases"]:
        if c["band"] != "warm":
            continue
        med = abs(s_signed) * c["spread_K"]
        rows.append({"ambient_C": c["ambient_C"], "rate": c["rate"], "h": c["h"], "median_registered": c["median"],
                     "median_signed_reading": med, "offers_Geslin": med > d["floors"]["Geslin"]["mean"],
                     "offers_KIT": med > d["floors"]["KIT"]["mean"]})
    sd = EX.lognormal_shape(d["floors"]["Geslin"]["mean"], d["floors"]["Geslin"]["p95"])["sd_e"]
    n2 = EX.n_required(abs(s_signed) * 2.0, sd)
    return {"note": "conditional: the reading not chosen on 27 October 2026", "signed_median_s": s_signed,
            "lever": lev[0]["label"] if lev else None, "cases": rows,
            "warm_chamber_dT2_Geslin_n": n2, "warm_chamber_dT5_Geslin_n": EX.n_required(abs(s_signed) * 5.0, sd)}


def procedure_notes(res):
    out = {}
    for T, run in res["runs"].items():
        failed = [f"{x['c_rate']} {x['direction']}" for x in run["cycle_heat"] if "failed" in x]
        failed += [f"pulse {x['c_rate']} {x['pulse_s']} s" for x in run["pulses"] if "failed" in x]
        c16 = [x for x in run["cycle_heat"] if abs(x["c_rate"] - 0.0625) < 1e-9 and "failed" not in x]
        c2 = [x for x in run["cycle_heat"] if abs(x["c_rate"] - 0.5) < 1e-9]
        out[T] = {"failed": failed,
                  "C16": {x["direction"]: {"R_eff_ohm": x["R_eff_irreversible_ohm"], "mean_total_W": x["mean Total heating [W]"],
                                           "mean_reversible_W": x["mean Reversible heating [W]"],
                                           "mean_ohmic_plus_irreversible_W": x["mean Ohmic heating [W]"] + x["mean Irreversible electrochemical heating [W]"]}
                          for x in c16},
                  "C2": {x["direction"]: {"duration_s": x["duration_s"], "charge_Ah": x["charge_Ah"], "R_eff_ohm": x["R_eff_irreversible_ohm"]}
                         for x in c2}}
    return out


def v4_summary():
    rows = list(csv.DictReader(open(os.path.join(RES, "stage4_experiment.csv"))))
    n = np.array([int(r["n_pairs"]) for r in rows])
    ns = np.array([int(r["n_sim_lognormal"]) for r in rows])
    nn = np.array([int(r["n_sim_normal_same_sd"]) for r in rows])
    rel = (ns - n) / n
    reln = (nn - n) / n
    big, small = n > 100, n <= 50
    pw = np.array([float(r["sim_power_at_n_calc_lognormal"]) for r in rows])
    nrm = np.array([int(r["n_pairs_normal_shape"]) for r in rows])
    return {"cases": int(n.size), "v4_within_5pct": int((np.abs(rel) <= 0.05).sum()),
            "lognormal_rel_min": float(rel.min()), "lognormal_rel_max": float(rel.max()),
            "lognormal_abs_rel_max_where_n_above_100": float(np.abs(rel[big]).max()),
            "lognormal_largest_shortfall_pairs": int((n - ns).max()),
            "sim_power_at_calculated_n_min": float(pw.min()),
            "normal_check_abs_pairs_max_where_n_at_most_50": int(np.abs(nn - n)[small].max()),
            "normal_check_abs_rel_max_where_n_above_100": float(np.abs(reln[big]).max()),
            "normal_check_within_5pct": int((np.abs(reln) <= 0.05).sum()),
            "normal_shape_over_registered_n": {"min": float((nrm / n).min()), "max": float((nrm / n).max()),
                                               "where_n_above_100_mean": float((nrm / n)[big].mean())}}


def main():
    d = CM.load_all()
    res = json.load(open(os.path.join(RES, "stage4_resistance.json")))
    mp = json.load(open(os.path.join(RES, "stage4_map.json")))
    sc = json.load(open(os.path.join(RES, "stage4_scorecard.json")))
    R = res["R_used_ohm"]
    out = {"cold_bracket": cold_bracket(R, d), "signed_median": signed_median(d, mp), "procedure": procedure_notes(res),
           "v4": v4_summary()}
    warm35 = [c for c in mp["cases"] if c["ambient_C"] == 35.0]
    out["margins"] = {"P1_largest_median_to_1C": max(c["median"] for c in warm35 if c["crate"] <= 1.0),
                      "P5_largest_35C_median": max(c["median"] for c in warm35),
                      "P5_KIT_mean": d["floors"]["KIT"]["mean"],
                      "cases_with_cells_outside_band_setpoints": sum(not c["cells_inside_band_setpoints"] for c in mp["cases"]),
                      "of_cases": len(mp["cases"])}
    out["scorecard"] = sc.get("_summary")
    json.dump(out, open(os.path.join(RES, "stage4_posthoc.json"), "w"), indent=1)
    print(json.dumps({"p3_ratio_with": out["cold_bracket"]["p3_ratio_with"], "signed": {k: out["signed_median"][k] for k in
                      ("signed_median_s", "lever", "warm_chamber_dT2_Geslin_n", "warm_chamber_dT5_Geslin_n")},
                      "margins": out["margins"]}, indent=1))
    for x in out["cold_bracket"]["cases"]:
        if x["h"] == 10.0:
            print(f"{x['ambient_C']:g} C with R({x['resistance_of_C']:g} C) {x['rate']}: spread {x['spread_K']:.3f} K, "
                  f"T {x['T_min_C']:.1f} to {x['T_max_C']:.1f} C, median {100 * x['median_divergence']:.2f} %")


if __name__ == "__main__":
    main()
