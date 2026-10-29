#!/usr/bin/env python3
"""V2 of Stage 5 (DECISIONS.md D77): beta and the within-set levers by separate code.

Imports nothing of src/stage5 or src/stage4. Reads the pairs from the archived files itself (D71, D72), forms p and d
(D73), fits beta and s_w by least squares through the origin with numpy.linalg.lstsq, and draws its own bootstrap of
sets with Python's random.Random (Mersenne Twister, seed 20261030), the pooled draws stratified by dataset. Passes
where beta and every s_w agree with results/stage5_beta.json to 1e-10 (|a - b| <= 1e-10 max(1, |a|)) and every interval
end lies within 4 sqrt(2) times the Monte Carlo standard deviation of one 10,000-draw percentile that the pipeline
records for it.

Usage: python3 verification/stage5/v2_separate.py
"""
import csv, json, math, os, random

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
B = 10000


def kit():
    lev = json.load(open(os.path.join(ROOT, "results", "stage2_lever.json")))
    warm = {r["protocol"]: r for r in lev["lever_25_40"]}
    cold = {r["protocol"]: r for r in lev["lever_10_25"]}
    out = []
    for r in csv.DictReader(open(os.path.join(ROOT, "results", "stage2_warmer_pairs.csv"))):
        T = int(r["age_temp"])
        pr = r["protocol"]
        if T == 40:
            band = "warm"
        elif T == 10:
            band = "cold"
        elif T == 25:
            wl, cl = warm[pr]["ci95"], cold[pr]["ci95"]
            band = "warm" if (wl[1] - wl[0]) < (cl[1] - cl[0]) else "cold"
        else:
            continue                                       # 0 C: no lever reaches it (D72)
        s = (warm if band == "warm" else cold)[pr]["s_per_K"]
        dT = float(r["dT_K"])
        out.append({"ds": "KIT", "set": r["parameter_id"], "band": band, "dT": dT, "p": -s * dT,
                    "d": math.log(float(r["L_warm"]) / float(r["L_cool"]))})
    return out


def sandia():
    lev = {}
    for r in csv.DictReader(open(os.path.join(ROOT, "results", "stage3_lever.csv"))):
        if r["series"] == "registered 1C" and r["determinable"] == "True":
            lev[(r["chemistry"], r["secant"])] = (float(r["s_per_K"]), float(r["ci95_hi"]) - float(r["ci95_lo"]))
    out = []
    for r in csv.DictReader(open(os.path.join(ROOT, "results", "stage3_pairs.csv"))):
        if not r["rise_a_K"] or not r["rise_b_K"]:
            continue
        ch, amb = r["chemistry"], float(r["ambient_C"])
        w, c = (ch, "25 to 35 C"), (ch, "15 to 25 C")
        if amb == 35.0:
            key = w if w in lev else None
        elif amb == 15.0:
            key = c if c in lev else None
        else:
            if w in lev and c in lev:
                key = w if lev[w][1] < lev[c][1] else c
            else:
                key = w if w in lev else (c if c in lev else None)
        if key is None:
            continue
        a, b = float(r["rise_a_K"]), float(r["rise_b_K"])
        La, Lb = float(r["L_a_efc"]), float(r["L_b_efc"])
        dT, d = (a - b, math.log(La / Lb)) if a > b else (b - a, math.log(Lb / La))
        out.append({"ds": "Sandia", "set": r["group"], "band": "warm" if key[1] == "25 to 35 C" else "cold", "dT": dT,
                    "p": -lev[key][0] * dT, "d": d})
    return out


def lstsq_slope(x, y):
    sol, *_ = np.linalg.lstsq(np.asarray(x, dtype=float)[:, None], np.asarray(y, dtype=float), rcond=None)
    return float(sol[0])


def boot(rows, xkey, sign, rnd):
    """Draws of sets with replacement within each dataset; returns the slope of each draw."""
    sets = {}
    for r in rows:
        s = sets.setdefault((r["ds"], r["set"]), [0.0, 0.0])
        s[0] += r[xkey] * r["d"]
        s[1] += r[xkey] * r[xkey]
    strata = {}
    for k, v in sets.items():
        strata.setdefault(k[0], []).append(v)
    vals = []
    for _ in range(B):
        num = den = 0.0
        for st in ("KIT", "Sandia"):
            g = strata.get(st, [])
            m = len(g)
            for _ in range(m):
                v = g[rnd.randrange(m)]
                num += v[0]
                den += v[1]
        vals.append(sign * num / den)
    return vals


def main():
    pipe = json.load(open(os.path.join(ROOT, "results", "stage5_beta.json")))
    K, S = kit(), sandia()
    data = {"KIT": K, "Sandia": S, "pooled": K + S}
    rnd = random.Random(20261030)
    out = {"check": "V2 (D77)", "rows": []}
    ok = True

    def compare(label, point, pipe_point, vals, pipe_boot, scale=1.0):
        nonlocal ok
        lo, hi = np.percentile(vals, [2.5, 97.5])
        mc = pipe_boot["mc_sd_of_ends"]
        tol = [4.0 * math.sqrt(2.0) * m for m in mc]
        d_pt = abs(point - pipe_point)
        ends = [abs(lo - pipe_boot["ci95"][0]), abs(hi - pipe_boot["ci95"][1])]
        good = d_pt <= 1e-10 * max(1.0, abs(pipe_point)) and ends[0] <= tol[0] and ends[1] <= tol[1]
        ok &= good
        out["rows"].append({"quantity": label, "separate": point, "pipeline": pipe_point, "abs_diff": d_pt,
                            "separate_ci95": [float(lo), float(hi)], "pipeline_ci95": pipe_boot["ci95"],
                            "end_diffs": [float(e) for e in ends], "tolerances": tol,
                            "ends_in_mc_sd": [float(e / m) if m > 0 else None for e, m in zip(ends, mc)], "pass": bool(good)})
        print(f"{label:22s} separate {point:+.12f} pipeline {pipe_point:+.12f} diff {d_pt:.1e}  ends "
              f"[{lo:+.4f}, {hi:+.4f}] against [{pipe_boot['ci95'][0]:+.4f}, {pipe_boot['ci95'][1]:+.4f}] "
              f"({ends[0] / mc[0]:.2f} and {ends[1] / mc[1]:.2f} MC sd)  {good}")

    for name in ("KIT", "Sandia", "pooled"):
        rows = data[name]
        beta = lstsq_slope([r["p"] for r in rows], [r["d"] for r in rows])
        compare(f"beta {name}", beta, pipe["dataset"][name]["beta"], boot(rows, "p", 1.0, rnd), pipe["dataset"][name]["beta_boot"])
    for name in ("KIT", "Sandia", "pooled"):
        for b in ("warm", "cold"):
            rows = [r for r in data[name] if r["band"] == b]
            sw = -lstsq_slope([r["dT"] for r in rows], [r["d"] for r in rows])
            compare(f"s_w {name} {b}", sw, pipe["band"][name][b]["s_w_per_K"], boot(rows, "dT", -1.0, rnd),
                    pipe["band"][name][b]["s_w_boot"])
    out["counts"] = {"KIT": len(K), "Sandia": len(S), "pooled": len(K) + len(S)}
    out["pass"] = bool(ok)
    json.dump(out, open(os.path.join(HERE, "v2.json"), "w"), indent=1)
    print(out["counts"], "V2 pass:", ok)


if __name__ == "__main__":
    main()
