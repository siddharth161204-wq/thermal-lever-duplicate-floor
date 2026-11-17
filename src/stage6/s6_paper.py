#!/usr/bin/env python3
"""Stage 6, every number the note, the README and the state document quote, and Part D's table (DECISIONS.md D95 (4)).

Reads the Stage 6 results and the archived Stage 3 and Stage 5 values and writes results/stage6_paper_numbers.json:
Part A0 (files, counts, contact flags, reliability, resistance, channel, power), Part B, Part C, the brackets and
unscored lines, the scorecard, the checks, and Part D: the within-set and between-policy levers beside Sandia's levers
for the same cell and beside what Stage 5 and this stage could detect, with the protocols' differences.

Usage: python3 -I src/stage6/s6_paper.py
"""
import csv, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s6_common as C                   # noqa: E402

PROTOCOLS = {
    "Sandia (Stage 3)": "charged at 0.5C, discharged at 0.5C to 3C (the levers at 1C and 2C), cells in 15, 25 and 35 C air, "
                        "a K-type thermocouple on the skin, periodic 0.5C capacity checks",
    "these cells (Stage 6)": "charged at 3.6C to 8C to 80 % state of charge in one, two or four steps, then 1C CC-CV to 3.6 V, "
                             "discharged at 4C CC-CV to 2.0 V, in a 30 C chamber, a type T thermocouple on the can, "
                             "lifetime counted cycle by cycle to 0.88 Ah",
}


def j(rel):
    return json.load(open(C.path(rel)))


def pct(x):
    return {"value": 100 * x["value"], "ci95": [100 * v for v in x["ci95"]]}


def main():
    arch = C.archive()
    a0, pw = j("results/stage6_a0.json"), j("results/stage6_power.json")["cases"]
    W, B = j("results/stage6_within.json"), j("results/stage6_between.json")
    sc = j("results/stage6_scorecard.json")
    cells = list(csv.DictReader(open(C.path("results/stage6_cells.csv"))))
    ent = [r for r in cells if r["enters"] == "True"]
    files = j("inputs/matr/files.json")["files"]
    out = {"files": {k: {"bytes": v["bytes"], "sha256": v["sha256"]} for k, v in files.items()},
           "counts": a0["counts"],
           "contact": {"flagged": sum(r["contact_flag"] == "True" for r in ent),
                       "flagged_in_sets_of_two_or_more": sum(r["contact_flag"] == "True" and int(r["set_size"]) >= 2 for r in ent),
                       "by_batch": {b: sum(1 for r in ent if r["batch"] == b and r["contact_flag"] == "True")
                                    for b in sorted(set(r["batch"] for r in ent))}},
           "reliability": a0["reliability"], "resistance": a0["resistance"],
           "channel": {k: v for k, v in a0["channel"].items() if k != "instrument"},
           "power": [{k: v for k, v in p.items() if k in ("scope", "floor", "cells", "sets", "T_dev_rms_K", "detectable_s_per_K",
                                                           "sim_registered_shape_detectable_s_per_K", "v5_rel_diff", "v5_pass",
                                                           "pairs", "pairs_form_detectable_s_per_K", "shared_cell_F",
                                                           "pairs_form_shared_cells_detectable_s_per_K", "without_flagged")}
                     for p in pw]}
    reg = {}
    for s in ("pooled", "Severson", "Attia"):
        w, b = W["registered"][s], B["registered"][s]
        r = {"cells_within": w["cells_within"], "sets_within": w["sets_within"], "cells_all": b["cells_all"],
             "sets_all": b["sets_all"], "s_w_pct_per_K": pct(w["B1"]["s_w_per_K"]),
             "pairs_form_pct_per_K": pct(w["B1"]["pairs_form_s_per_K"]), "pairs": w["B1"]["pairs"], "B3": w["B3"],
             "sd_cell": w["B4"]["sd_cell"], "sqrt2_sd_cell": w["B4"]["sqrt2_sd_cell"],
             "sqrt2_sd_cell_ci95": w["B4"]["sqrt2_sd_cell_ci95"], "s_b_pct_per_K": pct(b["C1"]["s_b_per_K"]),
             "s_b_minus_s_w_pct_per_K": pct(b["C2"]["s_b_minus_s_w_per_K"]), "C3_share": b["C3"]["share"],
             "C3_share_residual_form": b["C3"]["share_residual_form"], "reliability_bracket": w["reliability_bracket"]}
        if "s_w_R_held_per_K" in w["B2"]:
            r["B2"] = {"cells": w["B2"]["cells"], "sets": w["B2"]["sets"],
                       "s_w_R_held_pct_per_K": pct(w["B2"]["s_w_R_held_per_K"]),
                       "s_w_same_cells_pct_per_K": pct(w["B2"]["s_w_same_cells_per_K"]),
                       "lnL_per_lnR": w["B2"]["lnL_per_lnR"], "lnR_per_K_within": w["B2"]["lnR_per_K_within"]}
        else:
            r["B2"] = w["B2"]
        reg[s] = r
    out["registered"] = reg
    out["brackets"] = {br: {s: {"cells_within": W["brackets"][br][s]["cells_within"],
                                "s_w_pct_per_K": pct(W["brackets"][br][s]["B1"]["s_w_per_K"]),
                                "s_b_pct_per_K": pct(B["brackets"][br][s]["C1"]["s_b_per_K"]),
                                "B3": W["brackets"][br][s]["B3"],
                                "s_w_R_held_pct_per_K": (pct(W["brackets"][br][s]["B2"]["s_w_R_held_per_K"])
                                                         if "s_w_R_held_per_K" in W["brackets"][br][s].get("B2", {}) else None)}
                            for s in ("pooled", "Severson", "Attia")} for br in W["brackets"]}
    out["unscored"] = {"T_rise": {s: pct(W["unscored"]["T_rise"][s]["B1"]["s_w_per_K"]) for s in ("pooled", "Severson", "Attia")},
                       "R_cycles_2_10": {s: ({"s_w_R_held_pct_per_K": pct(v["s_w_R_held_per_K"]), "lnL_per_lnR": v["lnL_per_lnR"]}
                                             if "s_w_R_held_per_K" in v else v)
                                         for s, v in W["unscored"]["R_cycles_2_10"].items()},
                       "C1_sets_of_two_or_more": {s: pct(B["unscored"]["C1_sets_of_two_or_more"][s]["C1"]["s_b_per_K"])
                                                  for s in ("pooled", "Severson", "Attia")}}
    lev = arch["levers"]
    s5 = arch["stage5"]
    pk = {(p["scope"], p["floor"]): p for p in pw}
    out["part_D"] = {
        "rows": [
            {"what": "Stage 6, within a set, pooled (s_w)", **pct(W["registered"]["pooled"]["B1"]["s_w_per_K"])},
            {"what": "Stage 6, within a set, without flagged cells", **pct(W["brackets"]["without_flagged"]["pooled"]["B1"]["s_w_per_K"])},
            {"what": "Stage 6, within a set, the first paper", **pct(W["registered"]["Severson"]["B1"]["s_w_per_K"])},
            {"what": "Stage 6, within a set, the second paper", **pct(W["registered"]["Attia"]["B1"]["s_w_per_K"])},
            {"what": "Stage 6, between policies, pooled (s_b)", **pct(B["registered"]["pooled"]["C1"]["s_b_per_K"])},
            {"what": "Sandia, 15 to 25 C at 1C (Stage 3, registered)",
             "value": 100 * lev[("registered 1C", "15 to 25 C")][0], "ci95": [100 * lev[("registered 1C", "15 to 25 C")][1],
                                                                             100 * lev[("registered 1C", "15 to 25 C")][2]]},
            {"what": "Sandia, 15 to 25 C at 2C (Stage 3)",
             "value": 100 * lev[("secondary 2C", "15 to 25 C")][0], "ci95": [100 * lev[("secondary 2C", "15 to 25 C")][1],
                                                                            100 * lev[("secondary 2C", "15 to 25 C")][2]]},
            {"what": "Sandia, 25 to 35 C at 2C (Stage 3)",
             "value": 100 * lev[("secondary 2C", "25 to 35 C")][0], "ci95": [100 * lev[("secondary 2C", "25 to 35 C")][1],
                                                                            100 * lev[("secondary 2C", "25 to 35 C")][2]]}],
        "detectable_pct_per_K": {
            "Stage 6, pooled": {"KIT floor": 100 * pk[("pooled", "KIT")]["detectable_s_per_K"],
                                "Geslin floor": 100 * pk[("pooled", "Geslin")]["detectable_s_per_K"]},
            "Stage 6, pooled, without flagged cells": {"KIT floor": 100 * pk[("pooled", "KIT")]["without_flagged"]["detectable_s_per_K"],
                                                       "Geslin floor": 100 * pk[("pooled", "Geslin")]["without_flagged"]["detectable_s_per_K"]},
            "Stage 5, KIT replicates": {"KIT floor": 100 * s5[("KIT", "KIT")]["detectable_s_per_K"],
                                        "Geslin floor": 100 * s5[("KIT", "Geslin")]["detectable_s_per_K"]},
            "Stage 5, Sandia replicates": {"KIT floor": 100 * s5[("Sandia", "KIT")]["detectable_s_per_K"],
                                           "Geslin floor": 100 * s5[("Sandia", "Geslin")]["detectable_s_per_K"]},
            "Stage 5, pooled": {"KIT floor": 100 * s5[("pooled", "KIT")]["detectable_s_per_K"],
                                "Geslin floor": 100 * s5[("pooled", "Geslin")]["detectable_s_per_K"]}},
        "protocols": PROTOCOLS}
    out["scorecard"] = {k: sc[k].get("holds") for k in ("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8")}
    checks = {}
    for v in ("v1", "v2", "v3", "v4", "v5"):
        p = C.path(f"verification/stage6/{v}.json")
        checks[v.upper()] = j(f"verification/stage6/{v}.json").get("pass") if os.path.exists(p) else None
    out["checks"] = checks
    json.dump(out, open(C.path("results/stage6_paper_numbers.json"), "w"), indent=1)
    print(json.dumps(out["part_D"], indent=1))
    print(out["scorecard"], out["checks"])


if __name__ == "__main__":
    main()
