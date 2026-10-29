#!/usr/bin/env python3
"""Independent Stage 5 check, Part A: the pairs (D71, D72, D73), the input hashes and the counts.

Outputs: out/pairs.csv (every KIT and Sandia pair), out/pairs.json (counts, band choices, hashes, shapes, medians).
Usage: PYTHONDONTWRITEBYTECODE=1 python3 verification/stage5/independent/pairs.py
"""
import collections, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ind_common as C   # noqa: E402


def counts(rows):
    out = {}
    for ds in ("KIT", "Sandia"):
        allp = [x for x in rows if x["dataset"] == ds]
        sc = [x for x in allp if x["scored"]]
        by_band = {}
        for b in ("warm", "cold"):
            xs = [x for x in sc if x["band"] == b]
            by_band[b] = {"pairs": len(xs), "sets": len({x["set"] for x in xs}),
                          "ties": sum(x["tie"] for x in xs), "untied": sum(not x["tie"] for x in xs)}
        reasons = collections.Counter(x["reason_not_scored"] for x in allp if not x["scored"])
        out[ds] = {"pairs_in_file": len(allp), "sets_in_file": len({x["set"] for x in allp}),
                   "scored_pairs": len(sc), "scored_sets": len({x["set"] for x in sc}),
                   "scored_ties": sum(x["tie"] for x in sc), "scored_untied": sum(not x["tie"] for x in sc),
                   "scored_censored": sum(x["censored"] for x in sc),
                   "by_band": by_band, "not_scored_by_reason": dict(reasons),
                   "by_temperature_scored": dict(collections.Counter(str(x["temperature_C"]) for x in sc)),
                   "by_temperature_in_file": dict(collections.Counter(str(x["temperature_C"]) for x in allp))}
        if ds == "Sandia":
            with_a = [x for x in allp if x["dT_K"] is not None]
            out[ds]["pairs_with_measure_a_both"] = len(with_a)
            out[ds]["ties_among_measure_a_pairs"] = sum(x["tie"] for x in with_a)
            out[ds]["untied_measure_a_pairs"] = [x["set"] for x in with_a if not x["tie"]]
            out[ds]["scored_with_measure_b"] = sum(x["dT_b_K"] is not None for x in sc)
            out[ds]["scored_orientation_b_differs"] = sum((x["warm_b"] != x["warm"]) for x in sc)
            out[ds]["by_chemistry_band"] = dict(collections.Counter(f"{x['protocol_or_chemistry']} {x['band']}" for x in sc))
        if ds == "KIT":
            out[ds]["untied_pairs"] = [f"{x['set']} {x['warm']}-{x['cool']} ({x['band']})" for x in allp if not x["tie"]]
            out[ds]["ties_all_142"] = sum(x["tie"] for x in allp)
    sc = [x for x in rows if x["scored"]]
    out["pooled"] = {"scored_pairs": len(sc), "scored_sets": len({(x["dataset"], x["set"]) for x in sc})}
    return out


def main():
    rows, bands = C.build_pairs()
    v1o = C.v1_oriented_check(rows)
    sh = C.shapes()
    med = C.band_medians()
    hashes = []
    for rel in C.FILES_D71:
        here = C.sha256_file(C.path(rel))
        at_tag = C.sha256_at(C.ARCHIVE_TAG, rel)
        first = None
        for tag in ("v1.0.0", "v1.1.0", "v1.2.0", "v1.3.0"):
            hsh = C.sha256_at(tag, rel)
            if hsh is not None:
                first = {"tag": tag, "sha256": hsh}
                break
        hashes.append({"file": rel, "sha256_worktree": here, "sha256_at_" + C.ARCHIVE_TAG: at_tag,
                       "equal_to_archive": here == at_tag, "first_archived": first})
    if not all(h["equal_to_archive"] for h in hashes):
        C.stop("an input differs from v1.3.0")
    cnt = counts(rows)
    keys = ["dataset", "set", "temperature_C", "protocol_or_chemistry", "window", "rates", "band", "set_band",
            "reason_not_scored", "warm", "cool", "dT_K", "L_warm", "L_cool", "d", "lever", "s_per_K", "s_ci95_lo",
            "s_ci95_hi", "p", "tie", "resolution_efc", "res_warm_efc", "res_cool_efc", "censored", "scored",
            "warm_b", "cool_b", "dT_b_K", "d_b", "p_b", "warm_cell_T_C", "cool_cell_T_C",
            "cells_in_group", "cell_a", "cell_b", "L_a_efc", "L_b_efc", "rise_a_K", "rise_b_K", "rest_a_C", "rest_b_C"]
    C.write_csv("pairs.csv", rows, keys)
    doc = {"counts": cnt, "band_of_25C_sets": bands, "v1_oriented_measure_a": v1o, "input_hashes": hashes,
           "shapes_rebuilt": {k: {kk: v[kk] for kk in ("sigma", "mu", "roots", "sd_e", "sd_e_normal_shape", "median",
                                                         "mean_check", "p95_check")} | {"floor": v["floor"]}
                              for k, v in sh.items()},
           "band_medians": med}
    C.write_json("pairs.json", doc)
    import json
    print(json.dumps(cnt, indent=1))
    print(json.dumps(v1o))
    print(json.dumps({k: v["sd_e"] for k, v in sh.items()}))


if __name__ == "__main__":
    main()
