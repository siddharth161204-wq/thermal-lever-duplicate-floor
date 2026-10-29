#!/usr/bin/env python3
"""Stage 5: the numbers the note, the README and the state document quote, gathered from the results files.

Reads results/stage5_beta.json, stage5_power.json, stage5_scorecard.json, stage5_posthoc.json, verification/stage5/
(V1 to V4 and, once it exists, independent/record.json); writes results/stage5_paper_numbers.json.

Usage: python3 src/stage5/s5_paper.py
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s5_inputs as A                   # noqa: E402

COMMITS = {
    "predictions": {"commit": "da613a1488838fd8a73032882d4f9fc5d6a0750d", "pushed_utc": "2026-10-29T04:09:13Z",
                    "predictions_stage5_md_sha256": "db5ba3ff88b721d12fe5f556e2673acf132b394bfc412b6fd28b268ee4f87799"},
    "code": {"commit": "ee4a403", "pushed_utc": "2026-10-29T04:16:27Z"},
    "results": {"commit": "101ca54a199def64ab7227a41277e2408b7f5e91", "pushed_utc": "2026-10-29T04:22:03Z"},
}


def pct(x):
    return 100.0 * x


def main():
    R = os.path.join(A.ROOT, "results")
    V = os.path.join(A.ROOT, "verification", "stage5")
    b = json.load(open(os.path.join(R, "stage5_beta.json")))
    p = json.load(open(os.path.join(R, "stage5_power.json")))
    s = json.load(open(os.path.join(R, "stage5_scorecard.json")))
    ph = json.load(open(os.path.join(R, "stage5_posthoc.json")))
    v = {k: json.load(open(os.path.join(V, f"{k}.json"))) for k in ("v1", "v2", "v3", "v4")}
    out = {"commits": COMMITS, "counts": b["counts"]}
    out["beta"] = {k: {"beta": b["dataset"][k]["beta"], "ci95": b["dataset"][k]["beta_boot"]["ci95"],
                       "width": b["dataset"][k]["beta_boot"]["width"], "pairs": b["dataset"][k]["pairs"],
                       "sets": b["dataset"][k]["sets"]} for k in ("KIT", "Sandia", "pooled")}
    out["within_set_lever_pct_per_K"] = {
        k: {band: {"s_w": b["band"][k][band]["s_w_pct_per_K"], "ci95": b["band"][k][band]["s_w_boot_pct"]["ci95"],
                   "pairs": b["band"][k][band]["pairs"], "band_median": pct(b["band"][k][band]["band_median"]["s_per_K"])}
            for band in ("warm", "cold")} for k in ("KIT", "Sandia", "pooled")}
    out["within_set_lever_pct_per_K"]["KIT 0 C (unscored)"] = {"s_w": b["kit_0C_unscored"]["s_w_pct_per_K"],
                                                               "ci95": b["kit_0C_unscored"]["s_w_boot_pct"]["ci95"],
                                                               "pairs": b["kit_0C_unscored"]["pairs"]}
    out["sign"] = {k: {"all": b["dataset"][k]["sign"], **{band: b["band"][k][band]["sign"] for band in ("warm", "cold")}}
                   for k in ("KIT", "Sandia", "pooled")}
    out["residual"] = {k: {"sd": b["dataset"][k]["resid_sd"], "ci95": b["dataset"][k]["resid_sd_boot"]["ci95"],
                           "floors": b["dataset"][k]["floor_inside_interval"]} for k in ("KIT", "Sandia", "pooled")}
    br = b["brackets"]
    out["brackets"] = {
        "KIT untied": {"pairs": br["KIT untied"]["pairs"], "beta": br["KIT untied"]["beta"], "ci95": br["KIT untied"]["beta_boot"]["ci95"]},
        "Sandia untied": br["Sandia untied"], "pooled untied": br["pooled untied"],
        "Sandia measure (b)": {"beta": br["Sandia measure (b)"]["beta"], "ci95": br["Sandia measure (b)"]["beta_boot"]["ci95"],
                               "sign": br["Sandia measure (b)"]["sign"]},
        "Sandia without the NCA 15 C 1C pair": {"beta": br["Sandia without the NCA 15 C 1C pair"]["beta"],
                                                "ci95": br["Sandia without the NCA 15 C 1C pair"]["beta_boot"]["ci95"],
                                                "sign": br["Sandia without the NCA 15 C 1C pair"]["sign"]}}
    out["C1"] = [{k: x[k] for k in ("case", "floor", "registered", "pairs", "dT_rms_K", "detectable_s_per_K",
                                    "sim_lognormal_detectable_s_per_K", "detectable_s_shared_cells_per_K", "shared_cell_F",
                                    "detectable_s_normal_shape_per_K", "v4_rel_diff", "v4_pass")} for x in p["C1"]]
    out["C2"] = [{k: x[k] for k in ("dataset", "band", "floor", "registered", "pairs", "dT_rms_K", "abs_s_per_K", "power_calc",
                                    "power_sim_lognormal", "pairs_needed_calc", "pairs_needed_sim_lognormal_resampled",
                                    "pairs_needed_shared_cells", "pairs_needed_normal_shape")} for x in p["C2"]]
    kit_c1 = [x for x in p["C1"] if x["registered"] and x["case"] == "KIT" and x["floor"] == "KIT"][0]
    out["headline"] = {
        "KIT detectable over warm median": kit_c1["detectable_s_per_K"] / b["medians"]["warm"]["abs_s_per_K"],
        "KIT detectable over cold median": kit_c1["detectable_s_per_K"] / b["medians"]["cold"]["abs_s_per_K"],
        "Sandia interval width over KIT": out["beta"]["Sandia"]["width"] / out["beta"]["KIT"]["width"],
        "Sandia dT rms over KIT": [x for x in p["C1"] if x["registered"] and x["case"] == "Sandia" and x["floor"] == "KIT"][0]["dT_rms_K"] / kit_c1["dT_rms_K"],
    }
    out["scorecard"] = {k: v_["held"] for k, v_ in s["scorecard"].items()}
    out["posthoc"] = ph
    out["checks"] = {"V1": v["v1"]["pass"], "V2": v["v2"]["pass"],
                     "V2_max_point_diff": max(r["abs_diff"] for r in v["v2"]["rows"]),
                     "V2_max_end_diff_in_mc_sd": max(max(r["ends_in_mc_sd"]) for r in v["v2"]["rows"]),
                     "V3": v["v3"]["pass"], "V4_registered": v["v4"]["pass_registered"], "V4_all": v["v4"]["pass_all"],
                     "V4_failures": v["v4"]["failures"],
                     "V4_registered_rel_range": [min(r["v4_rel_diff"] for r in v["v4"]["rows"] if r["registered"]),
                                                 max(r["v4_rel_diff"] for r in v["v4"]["rows"] if r["registered"])],
                     "normal_check_max_abs_rel": max(abs(r["normal_check_rel_diff"]) for r in v["v4"]["rows"])}
    rec = os.path.join(V, "independent", "record.json")
    if os.path.exists(rec):
        out["independent"] = json.load(open(rec))
    json.dump(out, open(os.path.join(R, "stage5_paper_numbers.json"), "w"), indent=1)
    print(json.dumps({k: out[k] for k in ("beta", "headline", "scorecard", "checks")}, indent=1))


if __name__ == "__main__":
    main()
