#!/usr/bin/env python3
"""Stage 5 independent check, written after the freeze (ea3a5b4): the numbers behind the checker's flags 5, 7, 9 and 17,
which D81 is to record. Frozen values are read from out/; beside each bootstrap interval of a single 10,000-draw run,
an interval from 1,000,000 draws (own SFC64 seeds [BASE_SEED, 1400 + i]) gives the value free of Monte Carlo noise to
the digits quoted. The post-freeze additions are labelled; none is a registered quantity.

Flag 5: the KIT warm band against the between-set lever.
Flag 7: the single KIT set and the single Sandia group that carry the fits.
Flag 9: the 1C lever on the Sandia 2C and 3C groups.
Flag 17: resampling protocols instead of sets.

Output: out_compare/postfreeze_flag_numbers.json.
Usage: PYTHONDONTWRITEBYTECODE=1 python3 verification/stage5/independent/postfreeze_flag_numbers.py
"""
import csv, json, os, sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ind_common as C   # noqa: E402
import partb as B        # noqa: E402

BIG = 1_000_000
OUTC = os.path.join(C.HERE, "out_compare")
_k = [1400]


def big(strata_pairs, key=None):
    xys = [B.as_xy(ps, "a", key) for ps in strata_pairs if ps]
    _k[0] += 1
    b = B.bootstrap([B.stratum(xy) for xy in xys], (_k[0],), n_draws=BIG)
    out = {"seed": b["seed"], "draws": BIG, "s_w_ci95_pct_per_K": [100 * v for v in b["s_w_ci95"]]}
    if "beta_ci95" in b:
        out["beta_ci95"] = b["beta_ci95"]
    return out


def beta_of(ps):
    p = np.array([x["p"] for x in ps])
    d = np.array([x["d"] for x in ps])
    return float(np.sum(p * d) / np.sum(p * p))


def main():
    rows, _ = C.build_pairs()
    Bf = json.load(open(os.path.join(C.OUT, "B.json")))["analyses"]
    X = json.load(open(os.path.join(C.OUT, "B_crosscheck.json")))["mc_spread"]
    INF = json.load(open(os.path.join(C.OUT, "influence.json")))
    D80 = json.load(open(os.path.join(OUTC, "postfreeze_d80_checks.json")))
    K, S = C.select(rows, "KIT"), C.select(rows, "Sandia")
    Kw, Kc = C.select(rows, "KIT", "warm"), C.select(rows, "KIT", "cold")
    proto = lambda x: (x["dataset"], x["protocol_or_chemistry"] if x["dataset"] == "KIT" else
                       f"{x['protocol_or_chemistry']} {x['window']} {x['rates']}")
    out = {"note": "post-freeze; frozen values from out/, 1,000,000-draw intervals added beside them"}

    # ---- flag 5
    loo = []
    for s in dict.fromkeys(x["set"] for x in Kw):
        loo.append((s, beta_of([x for x in Kw if x["set"] != s])))
    loo.sort(key=lambda t: t[1])
    out["flag5"] = {
        "KIT warm s_w_pct_per_K": Bf["KIT warm"]["s_w_pct_per_K"],
        "KIT warm s_w ci95 (frozen 10,000)": Bf["KIT warm"]["s_w_ci95_pct_per_K"],
        "KIT warm s_w ci95 (1,000,000)": [100 * v for v in X["KIT warm"]["s_w"]["big_ci95"]],
        "KIT warm sign (registered per-band B3)": Bf["KIT warm"]["sign"],
        "KIT warm sign, one pair per set (post hoc)": D80["one_pair_per_set"]["KIT"]["warm"],
        "KIT warm beta (checker extra)": Bf["KIT warm"]["beta"],
        "KIT warm beta ci95 (frozen 10,000)": Bf["KIT warm"]["beta_ci95"],
        "KIT warm beta ci95 (1,000,000)": X["KIT warm"]["beta"]["big_ci95"],
        "KIT warm beta ci95, 12 protocol clusters (frozen 10,000)": Bf["sens KIT warm, protocol clusters"]["beta_ci95"],
        "KIT warm beta ci95, 12 protocol clusters (1,000,000)": big([Kw], proto)["beta_ci95"],
        "KIT warm beta, leaving out one set: range": [loo[0][1], loo[-1][1]], "lowest_without": loo[0], "highest_without": loo[-1],
        "KIT cold beta (checker extra)": Bf["KIT cold"]["beta"],
        "KIT cold beta ci95 (frozen 10,000)": Bf["KIT cold"]["beta_ci95"],
        "KIT cold beta ci95 (1,000,000)": X["KIT cold"]["beta"]["big_ci95"],
        "KIT beta": Bf["KIT"]["beta"],
    }
    # ---- flag 7
    def without(ps, ds_set):
        return [x for x in ps if (x["dataset"], x["set"]) != ds_set]

    k7 = ("KIT", "P036")
    s7 = ("Sandia", "NMC_25C_0-100_0.5/3C")
    g = [x for x in S if x["set"] == s7[1]]
    out["flag7"] = {
        "KIT P036": {"set": "P036, 10 C, 10-100 1.67/1.00C, lever -11.59 % per K",
                     "p2_share_of_KIT": INF["KIT"]["top_sets_by_p2_share"][0]["p2_share"],
                     "KIT beta without it": INF["KIT"]["without_top_set"]["beta"],
                     "ci95 (frozen 10,000)": INF["KIT"]["without_top_set"]["beta_ci95"],
                     "ci95 (1,000,000)": big([without(K, k7)])["beta_ci95"],
                     "KIT beta leaving out any one set: range": INF["KIT"]["beta_without_range"]},
        "Sandia NMC 25 C 0-100 0.5/3C": {
            "pairs": len(g), "cells": {x: None for x in ()},
            "p2_share_of_Sandia": INF["Sandia"]["top_sets_by_p2_share"][0]["p2_share"],
            "p2_share_of_pooled": INF["pooled"]["top_sets_by_p2_share"][0]["p2_share"],
            "Sandia share of pooled sum p^2": D80["sum_p2"]["Sandia"] / (D80["sum_p2"]["Sandia"] + D80["sum_p2"]["KIT"]),
            "Sandia beta without it": INF["Sandia"]["without_top_set"]["beta"],
            "Sandia ci95 (frozen 10,000)": INF["Sandia"]["without_top_set"]["beta_ci95"],
            "Sandia ci95 (1,000,000)": big([without(S, s7)])["beta_ci95"],
            "pooled beta without it": INF["pooled"]["without_top_set"]["beta"],
            "pooled ci95 (frozen 10,000)": INF["pooled"]["without_top_set"]["beta_ci95"],
            "pooled ci95 (1,000,000)": big([K, without(S, s7)])["beta_ci95"],
            "Sandia beta leaving out any one group: range": INF["Sandia"]["beta_without_range"],
            "pooled beta leaving out any one set: range": INF["pooled"]["beta_without_range"]},
    }
    cells = {}
    for x in g:
        cells[x["cell_a"][-1]] = (x["rise_a_K"], x["L_a_efc"])
        cells[x["cell_b"][-1]] = (x["rise_b_K"], x["L_b_efc"])
    out["flag7"]["Sandia NMC 25 C 0-100 0.5/3C"]["cells"] = {c: {"swing_K": v[0], "L_efc": v[1]} for c, v in sorted(cells.items())}
    wS = out["flag7"]["Sandia NMC 25 C 0-100 0.5/3C"]["Sandia ci95 (1,000,000)"]
    wK = X["KIT"]["beta"]["big_ci95"]
    out["flag7"]["Sandia NMC 25 C 0-100 0.5/3C"]["P7 width ratio without it (1,000,000 draws each)"] = (wS[1] - wS[0]) / (wK[1] - wK[0])
    # ---- flag 9
    hi = [x for x in S if x["rates"].endswith(("/2.0C", "/3.0C"))]
    lo = [x for x in S if not x["rates"].endswith(("/2.0C", "/3.0C"))]
    lev = {(r["series"], r["chemistry"], r["secant"]): r for r in csv.DictReader(open(C.path("results/stage3_lever.csv")))}
    out["flag9"] = {
        "Sandia scored pairs at 2C or 3C discharge": len(hi), "groups": len({x["set"] for x in hi}),
        "pairs at 1C": len(lo), "groups at 1C": len({x["set"] for x in lo}),
        "share of Sandia sum p^2 at 2C or 3C": float(sum(x["p"] ** 2 for x in hi) / sum(x["p"] ** 2 for x in S)),
        "levers_pct_per_K": {f"{c} {sec}": {"registered 1C": (100 * float(lev[("registered 1C", c, sec)]["s_per_K"])
                                                              if lev[("registered 1C", c, sec)]["determinable"] == "True" else None),
                                            "secondary 2C": (100 * float(lev[("secondary 2C", c, sec)]["s_per_K"])
                                                             if lev[("secondary 2C", c, sec)]["determinable"] == "True" else None)}
                             for c in ("LFP", "NCA", "NMC") for sec in ("25 to 35 C", "15 to 25 C")},
        "sensitivity (post-freeze, not registered): Sandia beta over the 1C groups only":
            {"beta": beta_of(lo), "ci95 (1,000,000)": big([lo])["beta_ci95"]},
        "sensitivity (post-freeze, not registered): Sandia beta over the 2C and 3C groups only":
            {"beta": beta_of(hi), "ci95 (1,000,000)": big([hi])["beta_ci95"]},
        "note": "no 3C lever exists; D72 gives every Sandia pair its chemistry's registered 1C lever",
    }
    # ---- flag 17
    out["flag17"] = {}
    for name, strata, fz in (("KIT beta", [K], "sens KIT, protocol clusters"), ("Sandia beta", [S], "sens Sandia, protocol clusters"),
                             ("pooled beta", [K, S], "sens pooled, protocol clusters"),
                             ("KIT warm s_w", [Kw], "sens KIT warm, protocol clusters"),
                             ("KIT cold s_w", [Kc], "sens KIT cold, protocol clusters")):
        b = big(strata, proto)
        rec = {"clusters": Bf[fz]["n_sets"], "frozen 10,000": (Bf[fz]["beta_ci95"] if "beta" in name else Bf[fz]["s_w_ci95_pct_per_K"]),
               "1,000,000": (b["beta_ci95"] if "beta" in name else b["s_w_ci95_pct_per_K"]),
               "registered sets, 1,000,000 (checker)": (X[name.split()[0]]["beta"]["big_ci95"] if "beta" in name else
                                                       [100 * v for v in X[" ".join(name.split()[:2])]["s_w"]["big_ci95"]])}
        out["flag17"][name] = rec
    kb = out["flag17"]["KIT beta"]["1,000,000"]
    sb = out["flag17"]["Sandia beta"]["1,000,000"]
    pb = out["flag17"]["pooled beta"]["1,000,000"]
    out["flag17"]["verdicts under protocol clusters (1,000,000 draws)"] = {
        "P1 (KIT contains 0)": kb[0] <= 0 <= kb[1], "P2 (KIT contains 1)": kb[0] <= 1 <= kb[1],
        "P7 ratio": (sb[1] - sb[0]) / (kb[1] - kb[0]), "P7": (sb[1] - sb[0]) >= 2 * (kb[1] - kb[0]),
        "P8 (pooled upper <= 3)": pb[1] <= 3.0, "P3, P4, P5, P6": "unchanged: point values, the sign test and C1 use no bootstrap"}
    os.makedirs(OUTC, exist_ok=True)
    with open(os.path.join(OUTC, "postfreeze_flag_numbers.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
