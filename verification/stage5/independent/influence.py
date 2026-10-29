#!/usr/bin/env python3
"""Independent Stage 5 check: leave-one-set-out influence on beta (the checker's own sensitivity; not registered).
For KIT, Sandia and pooled, beta with each set left out in turn, each set's share of sum p^2, and for the most
influential set the beta interval without it (10,000 draws, own seeds [BASE_SEED, 990 + i]).

Output: out/influence.json.
Usage: PYTHONDONTWRITEBYTECODE=1 python3 verification/stage5/independent/influence.py
"""
import json, os, sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ind_common as C   # noqa: E402
import partb as B        # noqa: E402


def beta_of(ps):
    p = np.array([x["p"] for x in ps])
    d = np.array([x["d"] for x in ps])
    return float(np.sum(p * d) / np.sum(p * p))


def main():
    rows, _ = C.build_pairs()
    out = {}
    for i, (name, ds) in enumerate((("KIT", "KIT"), ("Sandia", "Sandia"), ("pooled", None))):
        ps = C.select(rows, ds)
        tot = sum(x["p"] ** 2 for x in ps)
        sets = []
        for k in dict.fromkeys((x["dataset"], x["set"]) for x in ps):
            rest = [x for x in ps if (x["dataset"], x["set"]) != k]
            share = sum(x["p"] ** 2 for x in ps if (x["dataset"], x["set"]) == k) / tot
            sets.append({"set": f"{k[0]} {k[1]}", "p2_share": share, "beta_without": beta_of(rest)})
        sets.sort(key=lambda r: -r["p2_share"])
        top = sets[0]["set"].split(" ", 1)
        rest = [x for x in ps if not (x["dataset"] == top[0] and x["set"] == top[1])]
        strata = [[x for x in rest if x["dataset"] == "KIT"], [x for x in rest if x["dataset"] == "Sandia"]]
        res = B.analysis(f"{name} without {sets[0]['set']}", [s for s in strata if s], (990 + i,), registered=False)
        swing = sorted(sets, key=lambda r: r["beta_without"])
        out[name] = {"beta": beta_of(ps), "top_sets_by_p2_share": sets[:5],
                     "beta_without_range": [swing[0]["beta_without"], swing[-1]["beta_without"]],
                     "lowest_without": swing[0], "highest_without": swing[-1],
                     "without_top_set": {"set": sets[0]["set"], "pairs": res["n_pairs"], "beta": res["beta"],
                                         "beta_ci95": res["bootstrap"]["beta_ci95"], "seed": res["bootstrap"]["seed"]}}
        print(name, json.dumps(out[name]["without_top_set"]), "range without one set:", out[name]["beta_without_range"])
        for s in sets[:4]:
            print("   ", s)
    C.write_json("influence.json", out)


if __name__ == "__main__":
    main()
