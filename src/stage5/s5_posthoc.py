#!/usr/bin/env python3
"""Stage 5, after the results (DECISIONS.md D80): numbers that change no registered value or verdict.

(1) The sign test of B3 with one pair per set (each set's warmest and coolest cell, the pair of largest dT), since the
three pairs of a set share cells and the registered test counts them as independent. (2) How deep the ties are: each
scored pair's lifetime difference over its resolution, and the resolution over the lifetime. (3) The untied KIT pairs.

Writes results/stage5_posthoc.json.

Usage: python3 src/stage5/s5_posthoc.py
"""
import csv, json, os, statistics, sys

from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s5_inputs as A                   # noqa: E402


def f(x):
    return float(x)


def one_per_set(rows):
    best = {}
    for r in rows:
        k = (r["dataset"], r["set"])
        if k not in best or f(r["dT_K"]) > f(best[k]["dT_K"]):
            best[k] = r
    return list(best.values())


def sign(rows):
    agree = sum((f(r["d"]) > 0) == (f(r["p"]) > 0) for r in rows)
    return {"pairs": len(rows), "agree": agree, "share": agree / len(rows),
            "p_two_sided": float(stats.binomtest(agree, len(rows), 0.5).pvalue)}


def main():
    P = list(csv.DictReader(open(os.path.join(A.ROOT, "results", "stage5_pairs.csv"))))
    sc = [r for r in P if r["scored"] == "True"]
    groups = {"KIT": [r for r in sc if r["dataset"] == "KIT"], "Sandia": [r for r in sc if r["dataset"] == "Sandia"]}
    groups["pooled"] = groups["KIT"] + groups["Sandia"]
    out = {"note": "after the results; no registered number or verdict changes (D80)", "one_pair_per_set": {}, "ties": {}}
    for name, rows in groups.items():
        out["one_pair_per_set"][name] = {"all": sign(one_per_set(rows))}
        for b in ("warm", "cold"):
            out["one_pair_per_set"][name][b] = sign(one_per_set([r for r in rows if r["band"] == b]))
    cells = {r["cell"]: r for r in csv.DictReader(open(os.path.join(A.ROOT, "results", "stage2_cells.csv")))}
    for name in ("KIT", "Sandia"):
        rows = groups[name]
        depth = [abs(f(r["L_warm_efc"]) - f(r["L_cool_efc"])) / f(r["resolution_efc"]) for r in rows]
        if name == "KIT":
            share = [f(cells[c]["res90_efc"]) / f(cells[c]["L90_efc"]) for r in rows for c in (r["warm"], r["cool"])]
        else:
            share = [f(r["resolution_efc"]) / (0.5 * (f(r["L_warm_efc"]) + f(r["L_cool_efc"]))) for r in rows]
        out["ties"][name] = {"pairs": len(rows), "ties": sum(r["tie"] == "True" for r in rows),
                             "median_lifetime_difference_over_resolution": statistics.median(depth),
                             "max_lifetime_difference_over_resolution": max(depth),
                             "median_resolution_over_lifetime": statistics.median(share),
                             "resolution_basis": "each cell's res90_efc over its L90 (KIT); the pair's resolution over the pair's mean lifetime (Sandia)"}
    out["kit_untied_pairs"] = [{k: r[k] for k in ("set", "temp_C", "protocol", "band", "warm", "cool", "dT_K", "d", "p",
                                                  "L_warm_efc", "L_cool_efc", "resolution_efc")}
                               for r in groups["KIT"] if r["tie"] != "True"]
    json.dump(out, open(os.path.join(A.ROOT, "results", "stage5_posthoc.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
