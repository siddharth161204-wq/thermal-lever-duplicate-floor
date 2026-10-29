#!/usr/bin/env python3
"""Stage 5 independent check, written after the freeze (ea3a5b4): the checker's own recomputation of the post-hoc
numbers that D80, the README and the note quote, from the archived inputs and the checker's frozen pair table (nothing
of src/stage5 is imported). The pipeline's values are read from results/stage5_posthoc.json at 157827f only to set
them beside.

(1) The ties: per dataset, the median and largest lifetime difference over the pair's resolution, and the median
    resolution over the lifetime (KIT: each cell's res90_efc over its L90_efc, over the distinct cells of the scored
    pairs and, beside it, over the two cells of every scored pair; Sandia: the pair's resolution over its mean lifetime).
(2) The sign test with one pair per set: each set's warmest and coolest cell (KIT by T_C, Sandia by measure (a)).
(3) Facts the text quotes: Stage 4's pairs needed at 2 K in the warm chamber under the central lever; the SHA-256 of
    predictions_stage5.md at da613a1; the dT spreads and ratios.

Output: out_compare/postfreeze_d80_checks.json.
Usage: PYTHONDONTWRITEBYTECODE=1 python3 verification/stage5/independent/postfreeze_d80_checks.py
"""
import csv, json, os, statistics, subprocess, sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ind_common as C   # noqa: E402

POST = "157827f1a77ef43e1f4545fc1c6fb24c2fb42270"
OUTC = os.path.join(C.HERE, "out_compare")


def git_json(rev, rel):
    return json.loads(subprocess.run(["git", "-C", C.ROOT, "show", f"{rev}:{rel}"], capture_output=True, check=True).stdout)


def binom(k, n):
    return float(stats.binomtest(k, n, 0.5, alternative="two-sided").pvalue)


def main():
    rows, _ = C.build_pairs()
    cells = {r["cell"]: r for r in csv.DictReader(open(C.path("results/stage2_cells.csv")))}
    out = {}
    # (1) ties
    K = C.select(rows, "KIT")
    S = C.select(rows, "Sandia")
    kd = [abs(x["L_warm"] - x["L_cool"]) / x["resolution_efc"] for x in K]
    kcells = sorted({x["warm"] for x in K} | {x["cool"] for x in K})
    k_res_over_L_cells = [float(cells[c]["res90_efc"]) / float(cells[c]["L90_efc"]) for c in kcells]
    k_res_over_L_entries = [float(cells[c]["res90_efc"]) / float(cells[c]["L90_efc"]) for x in K for c in (x["warm"], x["cool"])]
    k_pairres_over_meanL = [x["resolution_efc"] / (0.5 * (x["L_warm"] + x["L_cool"])) for x in K]
    sd_ = [abs(x["L_warm"] - x["L_cool"]) / x["resolution_efc"] for x in S]
    s_res_over_L = [x["resolution_efc"] / (0.5 * (x["L_warm"] + x["L_cool"])) for x in S]
    out["ties"] = {
        "KIT": {"pairs": len(K), "ties": sum(x["tie"] for x in K), "median_dL_over_res": statistics.median(kd), "max_dL_over_res": max(kd),
                "median_res_over_L_distinct_cells": statistics.median(k_res_over_L_cells), "distinct_cells": len(kcells),
                "median_res_over_L_cell_entries": statistics.median(k_res_over_L_entries),
                "median_pair_res_over_mean_L": statistics.median(k_pairres_over_meanL)},
        "Sandia": {"pairs": len(S), "ties": sum(x["tie"] for x in S), "median_dL_over_res": statistics.median(sd_), "max_dL_over_res": max(sd_),
                   "median_res_over_mean_L": statistics.median(s_res_over_L)},
    }
    # (2) one pair per set: each set's warmest and coolest cell
    def one_pair(ps, temp_of):
        by = {}
        for x in ps:
            by.setdefault((x["dataset"], x["set"]), []).append(x)
        chosen = []
        for k, xs in by.items():
            cs = {c for x in xs for c in (x["warm"], x["cool"])}
            hot = max(cs, key=temp_of)
            cold = min(cs, key=temp_of)
            m = [x for x in xs if {x["warm"], x["cool"]} == {hot, cold}]
            if len(m) != 1:
                C.stop("one pair per set: no pair of the warmest and coolest cell in " + str(k))
            chosen.append(m[0])
        return chosen

    rise = {}
    for x in S:
        rise[x["cell_a"]] = x["rise_a_K"]
        rise[x["cell_b"]] = x["rise_b_K"]
    temp = lambda c: float(cells[c]["T_C"]) if c in cells else rise[c]
    out["one_pair_per_set"] = {}
    for ds, ps in (("KIT", K), ("Sandia", S), ("pooled", K + S)):
        ch = one_pair(ps, temp)
        rec = {}
        for b in ("all", "warm", "cold"):
            xs = [x for x in ch if b == "all" or x["band"] == b]
            agree = sum(int(np.sign(x["d"]) == np.sign(x["p"])) for x in xs)
            rec[b] = {"pairs": len(xs), "agree": agree, "p_two_sided": binom(agree, len(xs))}
        out["one_pair_per_set"][ds] = rec
    # (3) quoted facts
    s4 = [r for r in csv.DictReader(open(C.path("results/stage4_experiment.csv")))
          if r["part"] == "C2" and r["chamber"].startswith("warm") and r["central"] == "True" and float(r["dT_K"]) == 2.0]
    out["stage4_warm_chamber_central_2K"] = {r["floor"]: {"n_pairs": int(r["n_pairs"]), "n_sim_lognormal": int(r["n_sim_lognormal"])}
                                             for r in s4}
    out["predictions_stage5_sha256_at_da613a1"] = C.sha256_at(C.BASE_COMMIT, "predictions_stage5.md")
    rms = lambda ps: float(np.sqrt(np.mean([x["dT_K"] ** 2 for x in ps])))
    out["dT_rms_K"] = {"KIT": rms(K), "KIT warm": rms(C.select(rows, "KIT", "warm")), "KIT cold": rms(C.select(rows, "KIT", "cold")),
                       "Sandia": rms(S), "Sandia over KIT": rms(S) / rms(K)}
    out["dT_mean_K"] = {"KIT warm": float(np.mean([x["dT_K"] for x in C.select(rows, "KIT", "warm")]))}
    out["sum_p2"] = {"KIT": float(sum(x["p"] ** 2 for x in K)), "Sandia": float(sum(x["p"] ** 2 for x in S))}
    # the pipeline's values beside
    out["pipeline_posthoc_157827f"] = git_json(POST, "results/stage5_posthoc.json")
    os.makedirs(OUTC, exist_ok=True)
    with open(os.path.join(OUTC, "postfreeze_d80_checks.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "pipeline_posthoc_157827f"}, indent=1))


if __name__ == "__main__":
    main()
