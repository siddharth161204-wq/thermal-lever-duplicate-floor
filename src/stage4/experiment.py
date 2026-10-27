#!/usr/bin/env python3
"""Stage 4: the experiment that would settle it (DECISIONS.md D61), with V4 (D62).

n split pairs, each two duplicate cells with one held dT warmer; each pair gives d = ln L_warm - ln L_cool = -s dT + e,
e the duplicate difference, symmetric, |e| distributed as the floor's pair difference. Registered shape: |e|
log-normal with the floor's mean m and 95th percentile q (the smaller root). Test: the paired t-test of the d against
zero, two-sided at 5 %. Calculation: the noncentral t power with noncentrality sqrt(n) |s| dT / sd(e); n is the
smallest integer from 2 with power >= 0.80. Normal shape (one line): sd(e) = m sqrt(pi / 2).
V4: 10,000 synthetic experiments at each n tried, drawing e from the registered shape (and, beside it, normal e of the
same sd), the same test; the simulated n by bisection between 2 and four times the calculated n.

Outputs: results/stage4_experiment.csv, results/stage4_experiment.json (shapes, method, V4 summary),
results/stage4_experiment_curves.csv (pairs needed against dT for the figure), verification/stage4/v4.json.

Usage: python3 src/stage4/experiment.py [--no-sim]
"""
import argparse, csv, json, math, os, sys, time

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as CM                     # noqa: E402

ALPHA, POWER, N_EXP, SEED = 0.05, 0.80, 10000, 20261027
Z95 = float(stats.norm.ppf(0.95))
C1_DT = (1.0, 2.0, 5.0, 10.0)
C2_DT = (2.0, 5.0)
CHAMBERS = {"warm": 35.0, "cold": 10.0}


def lognormal_shape(m, q):
    """|e| log-normal with mean m and 95th percentile q: sigma^2/2 - z sigma + ln(q/m) = 0, the smaller root."""
    L = math.log(q / m)
    disc = Z95 * Z95 - 2.0 * L
    if disc < 0:
        raise SystemExit("no log-normal has this mean and 95th percentile")
    roots = (Z95 - math.sqrt(disc), Z95 + math.sqrt(disc))
    sig = roots[0]
    mu = math.log(m) - 0.5 * sig * sig
    mu_big = math.log(m) - 0.5 * roots[1] ** 2
    return {"sigma": sig, "mu": mu, "roots": list(roots), "median": math.exp(mu), "median_with_larger_root": math.exp(mu_big),
            "mean_check": math.exp(mu + 0.5 * sig * sig), "p95_check": math.exp(mu + Z95 * sig),
            "sd_e": math.exp(mu + sig * sig)}                            # sqrt(E[e^2]) = sqrt(exp(2 mu + 2 sigma^2))


def power_t(n, delta, sd):
    df = n - 1
    tc = stats.t.ppf(1.0 - ALPHA / 2.0, df)
    nc = math.sqrt(n) * delta / sd
    return float(stats.nct.sf(tc, df, nc) + stats.nct.cdf(-tc, df, nc))


def n_required(delta, sd):
    """Smallest integer n >= 2 with power >= 0.80 (power rises with n)."""
    if delta <= 0:
        return None
    if power_t(2, delta, sd) >= POWER:
        return 2
    lo, hi = 2, 4
    while power_t(hi, delta, sd) < POWER:
        lo, hi = hi, hi * 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if power_t(mid, delta, sd) >= POWER:
            hi = mid
        else:
            lo = mid
    return hi


def sim_power(n, delta, shape, seed, n_exp=N_EXP, chunk=5_000_000):
    rng = np.random.default_rng(seed)
    tc = stats.t.ppf(1.0 - ALPHA / 2.0, n - 1)
    rej, done = 0, 0
    per = max(1, int(chunk // n))
    while done < n_exp:
        k = min(per, n_exp - done)
        if shape["kind"] == "lognormal":
            mag = np.exp(shape["mu"] + shape["sigma"] * rng.standard_normal((k, n)))
            sign = rng.integers(0, 2, size=(k, n), dtype=np.int8).astype(np.float64) * 2.0 - 1.0
            d = delta + sign * mag
        else:
            d = delta + shape["sd_e"] * rng.standard_normal((k, n))
        mean = d.mean(axis=1)
        sdv = d.std(axis=1, ddof=1)
        tstat = mean / (sdv / math.sqrt(n))
        rej += int((np.abs(tstat) > tc).sum())
        done += k
    return rej / n_exp


def n_simulated(delta, shape, n_calc, case_id, code):
    """Bisection between 2 and 4 n_calc (doubling beyond if needed) on the simulated power, fresh experiments per n."""
    tried = {}

    def ok(n):
        if n not in tried:
            tried[n] = sim_power(n, delta, shape, [SEED, case_id, n, code])
        return tried[n] >= POWER

    if ok(2):
        return 2, tried
    lo, hi = 2, max(4 * n_calc, 3)
    while not ok(hi):
        lo, hi = hi, hi * 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ok(mid):
            hi = mid
        else:
            lo = mid
    return hi, tried


def cases(d):
    named, bands = d["named"], d["bands"]
    out = []
    for lev in named["c1"]:
        for dT in C1_DT:
            out.append({"part": "C1", "chamber": "", "lever_name": lev["name"], "lever": lev["label"], "source": lev["source"],
                        "s_per_K": lev["s"], "conditional": lev["conditional"], "central": False, "dT_K": dT})
    central = {"warm": named["median_warm"]["label"], "cold": named["median_cold"]["label"]}
    for ch in ("warm", "cold"):
        for lev in bands[ch]:
            for dT in C2_DT:
                out.append({"part": "C2", "chamber": f"{ch} ({CHAMBERS[ch]:g} C)", "lever_name": "central lever" if lev["label"] == central[ch] else "",
                            "lever": lev["label"], "source": lev["source"], "s_per_K": lev["s"], "conditional": False,
                            "central": lev["label"] == central[ch], "dT_K": dT})
    return out


def main(out_dir, ver_dir, simulate):
    d = CM.load_all()
    shapes = {}
    for fname, f in d["floors"].items():
        ln = lognormal_shape(f["mean"], f["p95"])
        shapes[fname] = {"lognormal": {"kind": "lognormal", **ln},
                         "normal": {"kind": "normal", "sd_e": f["mean"] * math.sqrt(math.pi / 2.0)},
                         "floor": f}
    rows = []
    t_start = time.time()
    for cid, c in enumerate(cases(d)):
        delta = abs(c["s_per_K"]) * c["dT_K"]
        for fname, sh in shapes.items():
            r = dict(c)
            r.update({"floor": fname, "delta_lnL": delta, "sd_e_registered": sh["lognormal"]["sd_e"]})
            n = n_required(delta, sh["lognormal"]["sd_e"])
            r["n_pairs"] = n
            r["power_at_n"] = power_t(n, delta, sh["lognormal"]["sd_e"])
            r["sd_e_normal_shape"] = sh["normal"]["sd_e"]
            r["n_pairs_normal_shape"] = n_required(delta, sh["normal"]["sd_e"])
            if simulate:
                key = cid * 2 + (0 if fname == "Geslin" else 1)
                ns, tried = n_simulated(delta, sh["lognormal"], n, key, 1)
                nn, tried_n = n_simulated(delta, {"kind": "normal", "sd_e": sh["lognormal"]["sd_e"]}, n, key, 2)
                r["n_sim_lognormal"] = ns
                r["v4_rel_diff"] = (ns - n) / n
                r["v4_pass"] = abs(ns - n) <= 0.05 * n
                r["n_sim_normal_same_sd"] = nn
                r["normal_check_rel_diff"] = (nn - n) / n
                r["normal_check_pass"] = abs(nn - n) <= 0.05 * n
                r["sim_power_at_n_sim_lognormal"] = tried[ns]
                r["sim_power_at_n_calc_lognormal"] = tried[n] if n in tried else sim_power(n, delta, sh["lognormal"], [SEED, key, n, 1])
            rows.append(r)
        if simulate and cid % 10 == 0:
            print(f"case {cid}: {time.time() - t_start:.0f} s", flush=True)
    # curves for the figure: the two chambers' central levers, dT 0.5 to 10 K, both floors, registered shape
    curves = []
    for ch in ("warm", "cold"):
        lev = d["named"]["median_" + ch]
        for fname, sh in shapes.items():
            for dT in np.round(np.arange(0.5, 10.0001, 0.1), 1):
                curves.append({"chamber": ch, "lever": lev["label"], "s_per_K": lev["s"], "floor": fname, "dT_K": float(dT),
                               "n_pairs": n_required(abs(lev["s"]) * float(dT), sh["lognormal"]["sd_e"]),
                               "n_pairs_normal_shape": n_required(abs(lev["s"]) * float(dT), sh["normal"]["sd_e"])})
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(out_dir, "stage4_experiment.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(out_dir, "stage4_experiment_curves.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(curves[0].keys()))
        w.writeheader()
        w.writerows(curves)
    doc = {"method": "D61: paired t-test of ln L_warm - ln L_cool against zero, two-sided 5 %, power 0.80 from the noncentral t "
                     "with sd(e) of the registered log-normal shape; normal shape sd(e) = m sqrt(pi/2)",
           "z_0.95": Z95,
           "shapes": {k: {"lognormal": v["lognormal"], "normal": v["normal"], "floor": v["floor"]} for k, v in shapes.items()},
           "named": {"median_warm": d["named"]["median_warm"]["label"], "median_cold": d["named"]["median_cold"]["label"],
                     "c1": [{"name": x["name"], "lever": x["label"], "s": x["s"], "conditional": x["conditional"]} for x in d["named"]["c1"]]},
           "cases": len(rows), "simulated": simulate}
    if simulate:
        doc["v4"] = {"pass_all": all(r["v4_pass"] for r in rows), "failures": sum(not r["v4_pass"] for r in rows),
                     "max_abs_rel_diff": max(abs(r["v4_rel_diff"]) for r in rows),
                     "normal_check_pass_all": all(r["normal_check_pass"] for r in rows),
                     "normal_check_max_abs_rel_diff": max(abs(r["normal_check_rel_diff"]) for r in rows)}
        os.makedirs(ver_dir, exist_ok=True)
        json.dump({"check": "V4 (D62)", **doc["v4"],
                   "rows": [{k: r[k] for k in ("part", "chamber", "lever", "s_per_K", "dT_K", "floor", "n_pairs", "n_sim_lognormal",
                                               "v4_rel_diff", "v4_pass", "n_sim_normal_same_sd", "normal_check_rel_diff",
                                               "normal_check_pass")} for r in rows]},
                  open(os.path.join(ver_dir, "v4.json"), "w"), indent=1)
    json.dump(doc, open(os.path.join(out_dir, "stage4_experiment.json"), "w"), indent=1)
    print(json.dumps({k: doc[k] for k in doc if k in ("shapes", "v4")}, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(CM.ROOT, "results"))
    ap.add_argument("--verification", default=os.path.join(CM.ROOT, "verification", "stage4"))
    ap.add_argument("--no-sim", action="store_true")
    a = ap.parse_args()
    main(a.out, a.verification, not a.no_sim)
