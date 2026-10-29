#!/usr/bin/env python3
"""Independent Stage 5 check: second computation of Part B's point values and the Monte Carlo spread of the
bootstrap intervals (the checker's own robustness study; not a registered quantity).

(1) beta and s_w by numpy.linalg.lstsq and by math.fsum, against partb.py's closed form.
(2) For the registered analyses and the near-threshold unscored ones, 200 independent replicate bootstraps of
    10,000 draws (own SFC64 seeds [BASE_SEED, 900 + i, r]) and one of 1,000,000 draws: the spread of every interval
    end and the share of replicates in which a verdict read from it would change.

Output: out/B_crosscheck.json.
Usage: PYTHONDONTWRITEBYTECODE=1 python3 verification/stage5/independent/crosscheck_b.py
"""
import json, math, os, sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ind_common as C   # noqa: E402
import partb as B        # noqa: E402

REPS, BIG = 200, 1_000_000


def lstsq_slope(x, y):
    sol, *_ = np.linalg.lstsq(np.asarray(x)[:, None], np.asarray(y), rcond=None)
    return float(sol[0])


def main():
    rows, _ = C.build_pairs()
    sh = C.shapes()
    Bj = json.load(open(os.path.join(C.OUT, "B.json")))["analyses"]
    K, S = C.select(rows, "KIT"), C.select(rows, "Sandia")
    sets = {
        "KIT": [K], "Sandia": [S], "pooled": [K, S],
        "KIT warm": [C.select(rows, "KIT", "warm")], "KIT cold": [C.select(rows, "KIT", "cold")],
        "Sandia warm": [C.select(rows, "Sandia", "warm")], "Sandia cold": [C.select(rows, "Sandia", "cold")],
        "pooled warm": [C.select(rows, "KIT", "warm"), C.select(rows, "Sandia", "warm")],
        "pooled cold": [C.select(rows, "KIT", "cold"), C.select(rows, "Sandia", "cold")],
        "KIT 0 C": [C.select(rows, "KIT", scored=False, temperature=0.0)],
        "Sandia without NCA 15C 1C": [C.select(rows, "Sandia", exclude_set={"NCA_15C_0-100_0.5/1C"})],
    }
    out = {"point_second_method": {}, "mc_spread": {}}
    worst = 0.0
    for name, strata in sets.items():
        ps = [x for st in strata for x in st]
        dT = [x["dT_K"] for x in ps]
        d = [x["d"] for x in ps]
        sw_l = -lstsq_slope(dT, d)
        sw_f = -math.fsum(a * b for a, b in zip(dT, d)) / math.fsum(a * a for a in dT)
        rec = {"s_w_lstsq": sw_l, "s_w_fsum": sw_f, "s_w_partb": Bj[name]["s_w_per_K"]}
        worst = max(worst, abs(sw_l - Bj[name]["s_w_per_K"]) / max(1, abs(sw_l)), abs(sw_f - Bj[name]["s_w_per_K"]))
        if ps[0]["p"] is not None:
            p = [x["p"] for x in ps]
            bl = lstsq_slope(p, d)
            bf = math.fsum(a * b for a, b in zip(p, d)) / math.fsum(a * a for a in p)
            r = [b_ - bf * a for a, b_ in zip(p, d)]
            sdf = math.sqrt(math.fsum(v * v for v in r) / (len(r) - 1))
            rec.update({"beta_lstsq": bl, "beta_fsum": bf, "beta_partb": Bj[name]["beta"], "resid_sd_fsum": sdf,
                        "resid_sd_partb": Bj[name]["resid_sd"]})
            worst = max(worst, abs(bl - Bj[name]["beta"]) / max(1, abs(bl)), abs(bf - Bj[name]["beta"]) / max(1, abs(bf)),
                        abs(sdf - Bj[name]["resid_sd"]))
        out["point_second_method"][name] = rec
    out["point_second_method_max_diff"] = worst
    print("second method, max difference:", worst, flush=True)

    # Monte Carlo spread of the interval ends
    thresholds = {"beta": {"KIT": [0.0, 1.0], "pooled": [3.0], "Sandia without NCA 15C 1C": [0.0]},
                  "s_w": {"pooled warm": [0.0], "Sandia warm": [0.0], "KIT warm": [0.0], "KIT cold": [0.0]},
                  "sd": {"KIT": [sh["KIT"]["sd_e"], sh["Geslin"]["sd_e"]], "Sandia": [sh["KIT"]["sd_e"], sh["Geslin"]["sd_e"]],
                         "pooled": [sh["KIT"]["sd_e"], sh["Geslin"]["sd_e"]]}}
    for i, (name, strata) in enumerate(sets.items()):
        xys = [B.as_xy(st) for st in strata]
        St = [B.stratum(xy) for xy in xys if len(xy)]
        reps = []
        for r in range(REPS):
            reps.append(B.bootstrap(St, (900 + i, r)))
        big = B.bootstrap(St, (950 + i,), n_draws=BIG)
        rec = {"replicates": REPS, "draws_each": B.N_DRAWS, "big_draws": BIG}
        for q, key in (("s_w", "s_w_ci95"), ("beta", "beta_ci95"), ("sd", "resid_sd_ci95")):
            if key not in reps[0] or reps[0][key] is None:
                continue
            arr = np.array([x[key] for x in reps])
            rec[q] = {"lo_mean": float(arr[:, 0].mean()), "lo_sd": float(arr[:, 0].std(ddof=1)),
                      "hi_mean": float(arr[:, 1].mean()), "hi_sd": float(arr[:, 1].std(ddof=1)),
                      "big_ci95": big[key], "registered_run_ci95": Bj[name]["bootstrap"][key]}
            if name in thresholds.get(q, {}):
                rec[q]["contains_share"] = {repr(t): float(np.mean((arr[:, 0] <= t) & (t <= arr[:, 1])))
                                            for t in thresholds[q][name]}
                rec[q]["big_contains"] = {repr(t): bool(big[key][0] <= t <= big[key][1]) for t in thresholds[q][name]}
        if "beta" in rec:
            w = np.array([x["beta_ci95"][1] - x["beta_ci95"][0] for x in reps])
            rec["beta"]["width_mean"], rec["beta"]["width_sd"] = float(w.mean()), float(w.std(ddof=1))
        out["mc_spread"][name] = rec
        print(name, json.dumps({k: v for k, v in rec.items() if k in ("beta", "s_w", "sd")})[:900], flush=True)
    C.write_json("B_crosscheck.json", out)


if __name__ == "__main__":
    main()
