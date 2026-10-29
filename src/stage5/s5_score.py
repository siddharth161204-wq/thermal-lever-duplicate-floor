#!/usr/bin/env python3
"""Stage 5: scoring of predictions_stage5.md by the rules of DECISIONS.md D78.

Reads results/stage5_beta.json and results/stage5_power.json; writes results/stage5_scorecard.json.

Usage: python3 src/stage5/s5_score.py
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s5_inputs as A                   # noqa: E402


def score(beta, power):
    K, S, Pd = beta["dataset"]["KIT"], beta["dataset"]["Sandia"], beta["dataset"]["pooled"]
    kw, kc = beta["band"]["KIT"]["warm"], beta["band"]["KIT"]["cold"]
    med = beta["medians"]
    lo, hi = K["beta_boot"]["ci95"]
    sc = {}
    sc["P1"] = {"held": bool(lo <= 0.0 <= hi), "beta": K["beta"], "ci95": [lo, hi]}
    sc["P2"] = {"held": bool(lo <= 1.0 <= hi), "beta": K["beta"], "ci95": [lo, hi]}
    sc["P3"] = {"held": bool(kw["s_w_per_K"] > 0.0), "s_w_pct_per_K": kw["s_w_pct_per_K"],
                "ci95_pct": kw["s_w_boot_pct"]["ci95"], "band_median_pct_per_K": 100.0 * med["warm"]["s_per_K"]}
    sc["P4"] = {"held": bool(kc["s_w_per_K"] < 0.0), "s_w_pct_per_K": kc["s_w_pct_per_K"],
                "ci95_pct": kc["s_w_boot_pct"]["ci95"], "band_median_pct_per_K": 100.0 * med["cold"]["s_per_K"]}
    sg = K["sign"]
    sc["P5"] = {"held": bool(sg["p_two_sided"] > 0.05), "agree": sg["agree"], "pairs": sg["pairs"], "share": sg["share"],
                "p_two_sided": sg["p_two_sided"], "left_out_zero": sg["left_out_zero"]}
    c1 = [x for x in power["C1"] if x["registered"] and x["case"] == "KIT" and x["floor"] == "KIT"]
    if len(c1) != 1:
        raise SystemExit("C1 KIT case under the KIT floor not found")
    c1 = c1[0]
    thr = [med["warm"]["abs_s_per_K"], med["cold"]["abs_s_per_K"]]
    calc_held = all(c1["detectable_s_per_K"] > t for t in thr)
    sim_held = all(c1["sim_lognormal_detectable_s_per_K"] > t for t in thr)
    scored = c1["v4_pass"] or (calc_held == sim_held)
    sc["P6"] = {"held": bool(calc_held) if scored else None, "scored": bool(scored),
                "detectable_s_pct_per_K": 100.0 * c1["detectable_s_per_K"],
                "simulated_pct_per_K": 100.0 * c1["sim_lognormal_detectable_s_per_K"], "v4_pass": c1["v4_pass"],
                "band_median_abs_s_pct_per_K": [100.0 * t for t in thr], "pairs": c1["pairs"]}
    wK, wS = K["beta_boot"]["width"], S["beta_boot"]["width"]
    sc["P7"] = {"held": bool(wS >= 2.0 * wK), "width_Sandia": wS, "width_KIT": wK, "ratio": wS / wK,
                "ci95_Sandia": S["beta_boot"]["ci95"], "ci95_KIT": [lo, hi]}
    phi = Pd["beta_boot"]["ci95"][1]
    sc["P8"] = {"held": bool(phi <= 3.0), "pooled_beta": Pd["beta"], "ci95": Pd["beta_boot"]["ci95"]}
    return sc


if __name__ == "__main__":
    R = os.path.join(A.ROOT, "results")
    beta = json.load(open(os.path.join(R, "stage5_beta.json")))
    power = json.load(open(os.path.join(R, "stage5_power.json")))
    sc = score(beta, power)
    doc = {"rules": "DECISIONS.md D78", "scorecard": sc,
           "held": [k for k, v in sc.items() if v["held"] is True],
           "failed": [k for k, v in sc.items() if v["held"] is False],
           "not_scored": [k for k, v in sc.items() if v["held"] is None]}
    json.dump(doc, open(os.path.join(R, "stage5_scorecard.json"), "w"), indent=1)
    print(json.dumps({k: doc[k] for k in ("held", "failed", "not_scored")}))
