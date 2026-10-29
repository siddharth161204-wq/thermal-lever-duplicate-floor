#!/usr/bin/env python3
"""Stage 5, after the independent check (DECISIONS.md D81): its flags recomputed with this pipeline's code.

None of these changes a registered number or verdict. (1) How much of each fit one set carries: each set's share of
sum(p^2), and beta without the set of largest share. (2) beta per band, with the bootstrap of sets (unscored). (3) The
brief's "resampling protocols": beta's interval when whole protocols (KIT: a protocol's sets at 10, 25 and 40 C; Sandia:
a chemistry, window and rates at every ambient) are drawn instead of sets. (4) Every registered interval again with
1,000,000 draws, to show how far 10,000 draws move an end. (5) The Sandia pairs that measure (b) orients the other way,
and the median dT of the scored pairs.

Writes results/stage5_posthoc_flags.json.

Usage: python3 src/stage5/s5_posthoc_flags.py
"""
import json, os, statistics, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s5_inputs as A                   # noqa: E402
import s5_within as W                   # noqa: E402

BIG = 1_000_000
CHUNK = 100_000


def shares(rows):
    t = W.set_table(rows)
    tot = float(t["pp"].sum())
    order = np.argsort(-t["pp"])
    top = [{"set": t["keys"][j][1], "dataset": t["keys"][j][0], "share_of_sum_p2": float(t["pp"][j] / tot)} for j in order[:3]]
    big = t["keys"][order[0]]
    rest = [r for r in rows if (r["dataset"], r["set"]) != big]
    return {"largest": top, "beta_all": W.beta_point(rows), "beta_without_largest": W.beta_point(rest),
            "pairs_without_largest": len(rest)}


def cluster_key(r, level):
    if level == "set":
        return (r["dataset"], r["set"])
    if r["dataset"] == "KIT":
        return ("KIT", r["protocol"])
    chem, amb, rest = r["set"].split("_", 2)          # e.g. NMC_25C_0-100_0.5/1C -> chemistry, ambient, window and rates
    return ("Sandia", chem + "_" + rest)


def boot_ends(rows, seed, level="set", n_boot=10000, what="beta"):
    """Percentile ends of beta (or s_w) with draws of clusters at `level`, stratified by dataset, in chunks."""
    keys, sums = [], {}
    for r in rows:
        k = cluster_key(r, level)
        if k not in sums:
            keys.append(k)
            sums[k] = [0.0, 0.0, 0.0, 0.0]
        s = sums[k]
        s[0] += r["p"] * r["d"]
        s[1] += r["p"] ** 2
        s[2] += r["dT_K"] * r["d"]
        s[3] += r["dT_K"] ** 2
    S = np.array([sums[k] for k in keys])
    strata = [k[0] for k in keys]
    rng = np.random.default_rng(seed)
    vals = []
    done = 0
    while done < n_boot:
        b = min(CHUNK, n_boot - done)
        C = W.boot_counts(strata, b, rng)
        if what == "beta":
            vals.append((C @ S[:, 0]) / (C @ S[:, 1]))
        else:
            vals.append(-(C @ S[:, 2]) / (C @ S[:, 3]))
        done += b
    v = np.concatenate(vals)
    lo, hi = np.percentile(v, [2.5, 97.5])
    return {"ci95": [float(lo), float(hi)], "clusters": len(keys), "draws": int(n_boot), "seed": list(seed) if isinstance(seed, (list, tuple)) else seed}


def main():
    import csv
    P = A.load_all()["pairs"]
    S = W.subsets(P)
    kit, snl, pool, band = S["kit"], S["snl"], S["pool"], S["band"]
    out = {"note": "after the independent check; no registered number or verdict changes (D81)"}
    out["influence"] = {name: shares(rows) for name, rows in (("KIT", kit), ("Sandia", snl), ("pooled", pool))}
    out["beta_by_band"] = {}
    k = 960
    for name, rows in (("KIT", kit), ("Sandia", snl), ("pooled", pool)):
        for b in ("warm", "cold"):
            sub = band(rows, b)
            out["beta_by_band"][f"{name} {b}"] = {"beta": W.beta_point(sub), **boot_ends(sub, [A.SEED, k])}
            k += 1
    out["resampling_protocols"] = {name: {"beta": W.beta_point(rows), **boot_ends(rows, [A.SEED, 970 + i], level="protocol")}
                                   for i, (name, rows) in enumerate((("KIT", kit), ("Sandia", snl), ("pooled", pool)))}
    out["million_draws"] = {}
    k = 980
    for name, rows in (("KIT", kit), ("Sandia", snl), ("pooled", pool)):
        out["million_draws"][f"beta {name}"] = boot_ends(rows, [A.SEED, k], n_boot=BIG)
        k += 1
        for b in ("warm", "cold"):
            out["million_draws"][f"s_w {name} {b} (per K)"] = boot_ends(band(rows, b), [A.SEED, k], n_boot=BIG, what="sw")
            k += 1
    flips = [r for r in snl if r.get("warm_b") and r["warm_b"] != r["warm"]]
    out["sandia_measure_b_orients_other_way"] = {"pairs": len(flips), "of": len(snl)}
    out["median_dT_K"] = {"KIT scored": statistics.median(r["dT_K"] for r in kit),
                          "KIT all 142": statistics.median(r["dT_K"] for r in P if r["dataset"] == "KIT"),
                          "Sandia scored": statistics.median(r["dT_K"] for r in snl)}
    json.dump(out, open(os.path.join(A.ROOT, "results", "stage5_posthoc_flags.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
