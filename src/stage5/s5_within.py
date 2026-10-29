#!/usr/bin/env python3
"""Stage 5, Part B: the within-set lever (DECISIONS.md D73 and D74).

B1 beta = sum(p d) / sum(p^2) over the scored pairs of a dataset; B2 the within-set lever s_w = -sum(dT d) / sum(dT^2)
per band; B3 the sign of d against the sign of p with the exact two-sided binomial test; B4 the residual standard
deviation of d - beta p beside Stage 4's sd(e). Intervals: 10,000 bootstrap draws of sets (KIT parameter sets, Sandia
groups) with replacement, the pooled draws stratified by dataset, percentile intervals; numpy's default_rng seeded
with [20261029, k]. The brackets of D74 run the same code on subsets.

Writes results/stage5_beta.json.

Usage: python3 src/stage5/s5_within.py
"""
import json, math, os, sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s5_inputs as A                   # noqa: E402

N_BOOT = 10000
Q = (0.025, 0.975)


# --------------------------------------------------------------------------------------------------------------- helpers
def set_table(rows, dkey="d", pkey="p", tkey="dT_K"):
    """Per set, in order of first appearance: stratum (dataset), n, and the sums the fits need."""
    order, sums = [], {}
    for r in rows:
        k = (r["dataset"], r["set"])
        if k not in sums:
            order.append(k)
            sums[k] = {"n": 0, "pd": 0.0, "pp": 0.0, "dd": 0.0, "Td": 0.0, "TT": 0.0}
        s = sums[k]
        d, T = r[dkey], r[tkey]
        p = r[pkey] if pkey else None
        s["n"] += 1
        s["dd"] += d * d
        s["Td"] += T * d
        s["TT"] += T * T
        if p is not None:
            s["pd"] += p * d
            s["pp"] += p * p
    arr = {f: np.array([sums[k][f] for k in order], dtype=float) for f in ("n", "pd", "pp", "dd", "Td", "TT")}
    arr["strata"] = [k[0] for k in order]
    arr["keys"] = order
    return arr


def boot_counts(strata, n_boot, rng):
    """Counts of each set in each draw: every stratum's sets drawn with replacement, as many as it has."""
    J = len(strata)
    C = np.zeros((n_boot, J))
    for st in ("KIT", "Sandia"):                         # fixed order of the strata
        idx = np.array([j for j in range(J) if strata[j] == st], dtype=int)
        m = len(idx)
        if m == 0:
            continue
        draw = rng.integers(0, m, size=(n_boot, m))
        cs = np.zeros((n_boot, m))
        np.add.at(cs, (np.arange(n_boot)[:, None], draw), 1.0)
        C[:, idx] = cs
    return C


def interval(x):
    """Percentile interval (numpy's linear interpolation) and the Monte Carlo spread of each end (D77's V2 rule)."""
    x = np.asarray(x, dtype=float)
    lo, hi = np.percentile(x, [100 * Q[0], 100 * Q[1]])
    xs = np.sort(x)
    B = len(xs)
    mc = []
    for q in Q:
        r = q * (B - 1)
        dr = math.sqrt(B * q * (1.0 - q))
        a, b = int(math.floor(r - dr)), int(math.ceil(r + dr))
        mc.append(0.5 * (xs[min(b, B - 1)] - xs[max(a, 0)]))
    return {"ci95": [float(lo), float(hi)], "width": float(hi - lo), "boot_mean": float(x.mean()),
            "boot_sd": float(x.std(ddof=1)), "mc_sd_of_ends": [float(v) for v in mc], "draws_used": int(B)}


def beta_point(rows, dkey="d", pkey="p"):
    num = sum(r[pkey] * r[dkey] for r in rows)
    den = sum(r[pkey] ** 2 for r in rows)
    return num / den


def sw_point(rows, dkey="d", tkey="dT_K"):
    return -sum(r[tkey] * r[dkey] for r in rows) / sum(r[tkey] ** 2 for r in rows)


def resid_sd_point(rows, beta, dkey="d", pkey="p"):
    n = len(rows)
    return math.sqrt(sum((r[dkey] - beta * r[pkey]) ** 2 for r in rows) / (n - 1))


def sign_test(rows, dkey="d", pkey="p"):
    use = [r for r in rows if r[dkey] != 0 and r[pkey] != 0]
    agree = sum((r[dkey] > 0) == (r[pkey] > 0) for r in use)
    n = len(use)
    p = float(stats.binomtest(agree, n, 0.5, alternative="two-sided").pvalue) if n else None
    return {"pairs": n, "left_out_zero": len(rows) - n, "agree": int(agree), "share": agree / n if n else None,
            "p_two_sided": p}


# -------------------------------------------------------------------------------------------------------------- analyses
def analyse(rows, k, what, dkey="d", pkey="p", tkey="dT_K", noise=None):
    """Point values and the bootstrap over sets for one subset; `what` lists the quantities D74 asks for."""
    t = set_table(rows, dkey, pkey if ("beta" in what or "resid" in what) else None, tkey)
    rng = np.random.default_rng([A.SEED, k])
    C = boot_counts(t["strata"], N_BOOT, rng)
    out = {"pairs": len(rows), "sets": len(t["keys"]),
           "sets_by_dataset": {st: t["strata"].count(st) for st in ("KIT", "Sandia") if st in t["strata"]},
           "ties": int(sum(bool(r["tie"]) for r in rows)), "seed": [A.SEED, k], "draws": N_BOOT}
    if "beta" in what or "resid" in what:
        num, den = C @ t["pd"], C @ t["pp"]
        keep = den > 0
        bb = num[keep] / den[keep]
        beta = beta_point(rows, dkey, pkey)
        out["draws_left_out_zero_sum_p2"] = int((~keep).sum())
        if "beta" in what:
            out["beta"] = beta
            out["beta_boot"] = interval(bb)
        if "resid" in what:
            n_star = (C @ t["n"])[keep]
            ssr = (C @ t["dd"])[keep] - 2.0 * bb * num[keep] + bb * bb * den[keep]
            sd_b = np.sqrt(np.maximum(ssr, 0.0) / (n_star - 1.0))
            sd = resid_sd_point(rows, beta, dkey, pkey)
            iv = interval(sd_b)
            out["resid_sd"] = sd
            out["resid_sd_boot"] = iv
            out["floor_inside_interval"] = {name: {"sd_e": f["sd_e"], "inside": bool(iv["ci95"][0] <= f["sd_e"] <= iv["ci95"][1]),
                                                   "ratio_resid_to_sd_e": sd / f["sd_e"]} for name, f in noise.items()}
    if "sw" in what:
        sw = sw_point(rows, dkey, tkey)
        swb = -(C @ t["Td"]) / (C @ t["TT"])
        out["s_w_per_K"] = sw
        out["s_w_pct_per_K"] = 100.0 * sw
        out["s_w_boot"] = interval(swb)
        out["s_w_boot_pct"] = {"ci95": [100.0 * v for v in out["s_w_boot"]["ci95"]]}
    if "sign" in what:
        out["sign"] = sign_test(rows, dkey, pkey)
    return out


def subsets(P):
    sc = [r for r in P if r["scored"]]
    kit = [r for r in sc if r["dataset"] == "KIT"]
    snl = [r for r in sc if r["dataset"] == "Sandia"]
    band = lambda rows, b: [r for r in rows if r["band"] == b]
    kit0 = [r for r in P if r["dataset"] == "KIT" and r["temp_C"] == 0.0]
    untied = lambda rows: [r for r in rows if not r["tie"]]
    snl_b = [r for r in snl if r.get("p_b") is not None]
    snl_x = [r for r in snl if r["set"] != A.NCA_15C_1C_GROUP]
    return {"kit": kit, "snl": snl, "pool": kit + snl, "kit0": kit0, "band": band, "untied": untied,
            "snl_b": snl_b, "snl_x": snl_x}


def run(P, noise, medians):
    S = subsets(P)
    band, untied = S["band"], S["untied"]
    res = {"method": "D73, D74: beta = sum(p d)/sum(p^2), s_w = -sum(dT d)/sum(dT^2), bootstrap of sets with replacement "
                     "(pooled stratified by dataset), 10,000 draws, percentile 95 % intervals; sign test exact two-sided; "
                     "residual sd sqrt(sum r^2/(n - 1))",
           "seed_base": A.SEED, "draws": N_BOOT, "medians": medians,
           "noise": {k: {"sd_e": v["sd_e"], "mean": v["mean"], "p95": v["p95"]} for k, v in noise.items()}}
    # B1, B3 and B4 per dataset
    res["dataset"] = {
        "KIT": analyse(S["kit"], 1, ("beta", "resid", "sign"), noise=noise),
        "Sandia": analyse(S["snl"], 2, ("beta", "resid", "sign"), noise=noise),
        "pooled": analyse(S["pool"], 3, ("beta", "resid", "sign"), noise=noise),
    }
    # B2 and B3 per band
    res["band"] = {}
    for name, rows, base in (("KIT", S["kit"], 10), ("Sandia", S["snl"], 20), ("pooled", S["pool"], 30)):
        res["band"][name] = {}
        for i, b in enumerate(("warm", "cold")):
            r = analyse(band(rows, b), base + 1 + i, ("sw", "sign"))
            r["band_median"] = medians[b]
            res["band"][name][b] = r
    # unscored: the 0 C KIT sets (no between-set lever)
    res["kit_0C_unscored"] = analyse(S["kit0"], 13, ("sw",))
    # brackets, unscored
    br = {}
    ku = untied(S["kit"])
    br["KIT untied"] = analyse(ku, 41, ("beta", "resid", "sign"), noise=noise)
    br["KIT untied, by band"] = {b: analyse(band(ku, b), 44 + i, ("sw", "sign")) for i, b in enumerate(("warm", "cold"))}
    su = untied(S["snl"])
    br["Sandia untied"] = {"pairs": len(su), "note": "no scored Sandia pair is untied (D72)"}
    if su:
        raise SystemExit("D72 records no untied scored Sandia pair")
    br["pooled untied"] = {"same_pairs_as": "KIT untied", "pairs": len(ku),
                           "note": "no scored Sandia pair is untied, so the pooled untied pairs are the KIT untied pairs and their numbers are KIT's"}
    br["Sandia measure (b)"] = analyse(S["snl_b"], 51, ("beta", "sign"), dkey="d_b", pkey="p_b", tkey="dT_b_K")
    br["Sandia measure (b)"]["pairs_without_measure_b"] = len(S["snl"]) - len(S["snl_b"])
    br["Sandia without the NCA 15 C 1C pair"] = analyse(S["snl_x"], 52, ("beta", "sign"))
    res["brackets"] = br
    res["counts"] = {
        "KIT pairs in file": sum(r["dataset"] == "KIT" for r in P),
        "KIT scored": len(S["kit"]), "KIT scored sets": len({r["set"] for r in S["kit"]}),
        "KIT 0 C pairs, not scored": len(S["kit0"]),
        "KIT warm": len(band(S["kit"], "warm")), "KIT cold": len(band(S["kit"], "cold")),
        "KIT ties among scored": sum(bool(r["tie"]) for r in S["kit"]),
        "KIT ties among all 142": sum(bool(r["tie"]) for r in P if r["dataset"] == "KIT"),
        "Sandia pairs in file": sum(r["dataset"] == "Sandia" for r in P),
        "Sandia with measure (a)": sum(r["dataset"] == "Sandia" and r["dT_K"] is not None for r in P),
        "Sandia scored": len(S["snl"]), "Sandia scored groups": len({r["set"] for r in S["snl"]}),
        "Sandia warm": len(band(S["snl"], "warm")), "Sandia cold": len(band(S["snl"], "cold")),
        "Sandia ties among scored": sum(bool(r["tie"]) for r in S["snl"]),
        "Sandia ties among the 60 with measure (a)": sum(bool(r["tie"]) for r in P if r["dataset"] == "Sandia" and r["dT_K"] is not None),
        "Sandia with measure (a) but no lever": sum(r["dataset"] == "Sandia" and r["dT_K"] is not None and not r["scored"] for r in P),
        "pooled scored": len(S["pool"]),
    }
    return res


if __name__ == "__main__":
    d = A.load_all()
    res = run(d["pairs"], d["noise"], d["medians"])
    json.dump(res, open(os.path.join(A.ROOT, "results", "stage5_beta.json"), "w"), indent=1)
    for k in ("KIT", "Sandia", "pooled"):
        x = res["dataset"][k]
        print(k, "beta", round(x["beta"], 4), [round(v, 4) for v in x["beta_boot"]["ci95"]], "sign", x["sign"]["agree"], "/",
              x["sign"]["pairs"], "p", round(x["sign"]["p_two_sided"], 4), "resid sd", round(x["resid_sd"], 5),
              [round(v, 5) for v in x["resid_sd_boot"]["ci95"]])
    for k in ("KIT", "Sandia", "pooled"):
        for b in ("warm", "cold"):
            x = res["band"][k][b]
            print(k, b, "s_w %", round(x["s_w_pct_per_K"], 3), [round(v, 3) for v in x["s_w_boot_pct"]["ci95"]], "sign",
                  x["sign"]["agree"], "/", x["sign"]["pairs"], "p", round(x["sign"]["p_two_sided"], 4))
    print(json.dumps(res["counts"], indent=1))
