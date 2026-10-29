#!/usr/bin/env python3
"""Independent Stage 5 check: the scorecard P1 to P8 by D78, from out/B.json and out/C.json, with each verdict's
margin and, where an interval decides it, the share of the 200 replicate bootstraps (out/B_crosscheck.json) and the
1,000,000-draw interval that give the same verdict.

Outputs: out/scorecard.json, out/scorecard.csv.
Usage: PYTHONDONTWRITEBYTECODE=1 python3 verification/stage5/independent/score.py
"""
import json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ind_common as C   # noqa: E402


def main():
    B = json.load(open(os.path.join(C.OUT, "B.json")))
    Cj = json.load(open(os.path.join(C.OUT, "C.json")))
    X = json.load(open(os.path.join(C.OUT, "B_crosscheck.json")))["mc_spread"]
    A = B["analyses"]
    med = B["band_medians"]
    sm_w, sm_c = abs(med["warm"]["s"]) * 100.0, abs(med["cold"]["s"]) * 100.0
    sc = {}

    lo, hi = A["KIT"]["beta_ci95"]
    sc["P1"] = {"rule": "the 95 % interval of beta for KIT (B1) contains zero, its ends included",
                "held": lo <= 0.0 <= hi, "value": [lo, hi], "margin": min(-lo, hi),
                "mc_share_same": X["KIT"]["beta"]["contains_share"]["0.0"],
                "big_draw_interval": X["KIT"]["beta"]["big_ci95"]}
    sc["P2"] = {"rule": "the KIT beta interval contains one", "held": lo <= 1.0 <= hi, "value": [lo, hi],
                "margin": min(1.0 - lo, hi - 1.0), "mc_share_same": X["KIT"]["beta"]["contains_share"]["1.0"],
                "big_draw_interval": X["KIT"]["beta"]["big_ci95"]}
    w = A["KIT warm"]["s_w_pct_per_K"]
    sc["P3"] = {"rule": "s_w of the KIT warm band (B2) is above zero", "held": w > 0.0, "value_pct_per_K": w,
                "interval_pct_per_K": A["KIT warm"]["s_w_ci95_pct_per_K"], "note": "point estimate, no Monte Carlo"}
    c = A["KIT cold"]["s_w_pct_per_K"]
    sc["P4"] = {"rule": "s_w of the KIT cold band is below zero", "held": c < 0.0, "value_pct_per_K": c,
                "interval_pct_per_K": A["KIT cold"]["s_w_ci95_pct_per_K"], "note": "point estimate, no Monte Carlo"}
    sg = A["KIT"]["sign"]
    sc["P5"] = {"rule": "the exact two-sided binomial p of the KIT sign agreement (B3) exceeds 0.05",
                "held": sg["p_two_sided"] > 0.05, "agree": sg["agree"], "pairs": sg["pairs_tested"],
                "left_out": sg["left_out"], "share": sg["share"], "p_two_sided": sg["p_two_sided"]}
    row = [r for r in Cj["C1"] if r["case"] == "KIT" and r["floor"] == "KIT"][0]
    calc, sim = row["detectable_s_calc_pct_per_K"], row["detectable_s_sim_lognormal_pct_per_K"]
    held_calc = calc > sm_w and calc > sm_c
    held_sim = sim > sm_w and sim > sm_c
    scored = row["v4_pass"] or held_sim == held_calc
    sc["P6"] = {"rule": "the KIT detectable |s| of C1 (calculation, KIT floor, registered shape, 106 pairs) exceeds both "
                        "band medians; not scored if V4 fails there and the simulated value gives the other verdict",
                "held": held_calc if scored else None, "scored": scored, "detectable_calc_pct_per_K": calc,
                "detectable_sim_lognormal_pct_per_K": sim, "v4_rel_diff": row["v4_rel_diff"], "v4_pass": row["v4_pass"],
                "band_medians_abs_pct_per_K": [sm_w, sm_c], "margin_over_cold_pct_per_K": calc - sm_c,
                "under_Geslin_floor_beside": [r["detectable_s_calc_pct_per_K"] for r in Cj["C1"]
                                              if r["case"] == "KIT" and r["floor"] == "Geslin"][0]}
    ws = A["Sandia"]["beta_ci95_width"]
    wk = A["KIT"]["beta_ci95_width"]
    sc["P7"] = {"rule": "the width of the Sandia beta interval is at least twice that of the KIT interval",
                "held": ws >= 2.0 * wk, "width_sandia": ws, "width_kit": wk, "ratio": ws / wk,
                "mc_mean_widths": [X["Sandia"]["beta"]["width_mean"], X["KIT"]["beta"]["width_mean"]],
                "big_draw_ratio": (X["Sandia"]["beta"]["big_ci95"][1] - X["Sandia"]["beta"]["big_ci95"][0]) /
                                  (X["KIT"]["beta"]["big_ci95"][1] - X["KIT"]["beta"]["big_ci95"][0])}
    plo, phi = A["pooled"]["beta_ci95"]
    sc["P8"] = {"rule": "the upper end of the pooled beta interval is at most 3", "held": phi <= 3.0, "value": [plo, phi],
                "margin": 3.0 - phi, "mc_share_same": 1.0 - X["pooled"]["beta"]["contains_share"]["3.0"],
                "big_draw_interval": X["pooled"]["beta"]["big_ci95"]}
    summary = {k: v["held"] for k, v in sc.items()}
    doc = {"scoring": "D78, read by the independent checker", "verdicts": summary, "details": sc}
    C.write_json("scorecard.json", doc)
    C.write_csv("scorecard.csv", [{"prediction": k, "held": v["held"], "rule": v["rule"],
                                   "numbers": json.dumps({kk: vv for kk, vv in v.items() if kk not in ("rule", "held")})}
                                  for k, v in sc.items()])
    print(json.dumps(summary))
    for k, v in sc.items():
        print(k, json.dumps({kk: vv for kk, vv in v.items() if kk != "rule"}))


if __name__ == "__main__":
    main()
