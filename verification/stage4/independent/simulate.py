#!/usr/bin/env python3
"""Independent check of Stage 4, part (e), the simulation in the spirit of V4 (D62).

For a case of out/experiment.csv (lever s, dT, floor, shape): each synthetic experiment draws n differences
d_j = -s dT + e_j, with e_j = (random sign) x exp(mu + sigma Z_j) under the registered log-normal shape, or
e_j = sd_e Z_j under the normal shape of the same sd(e), and applies the two-sided paired t-test at 5 %
(|t| > t_0.975,n-1, t = mean(d) / (sd(d) / sqrt(n))). The simulated power at n is the share of 10,000 experiments
that reject. The simulated n is the smallest n whose simulated power reaches 0.80, found by bisection between 2 and
four times the calculated n (D62). Seeds: numpy SeedSequence([20261027, case index, n, shape code]); the pipeline's
own seed scheme is not known to this checker, so the two simulated n can differ by Monte Carlo noise.

Usage:
  python3 simulate.py --shape lognormal [--subset required|all] [--nexp 10000]
  python3 simulate.py --precise   (P7's two cases, 200,000 experiments per n over a range of n)
Output: out/simulation_<shape>_<subset>.csv, out/simulation_precise.csv
"""
import argparse, csv, json, math, os, sys, time

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
sys.path.insert(0, HERE)
import levers as LV                                 # noqa: E402
import power as PW                                  # noqa: E402

SHAPE_CODE = {"lognormal": 1, "normal": 2}


def shapes():
    _, _, floors = LV.load()
    out = {}
    for fk, fv in floors.items():
        ln = PW.lognormal_shape(fv["mean"], fv["p95"])["small"]
        out[fk] = {"mu": ln["mu"], "sigma": ln["sigma"], "sd_ln": ln["sd_e"], "sd_norm": PW.normal_sd(fv["mean"])}
    return out


def sim_power(n, shift, shape, par, seed, nexp=10000, cells_per_chunk=4_000_000):
    """shift = -s dT; returns (power, wrong-sign share of all experiments)."""
    rng = np.random.default_rng(np.random.SeedSequence(seed))
    tc = stats.t.ppf(0.975, n - 1)
    rows = max(1, cells_per_chunk // n)
    rej = wrong = 0
    done = 0
    while done < nexp:
        k = min(rows, nexp - done)
        z = rng.standard_normal((k, n))
        if shape == "lognormal":
            e = np.exp(par["mu"] + par["sigma"] * z)
            sgn = rng.integers(0, 2, size=(k, n), dtype=np.int8)
            e = np.where(sgn == 1, e, -e)
        else:
            e = par["sd_norm"] * z
        d = shift + e
        m = d.mean(axis=1)
        v = d.var(axis=1, ddof=1)
        t = m / np.sqrt(v / n)
        r = np.abs(t) > tc
        rej += int(r.sum())
        wrong += int((r & (np.sign(m) != np.sign(shift))).sum())
        done += k
    return rej / nexp, wrong / nexp


def bisect_n(n_calc, shift, shape, par, case_idx, nexp):
    code = SHAPE_CODE[shape]
    cache = {}

    def P(n):
        if n not in cache:
            cache[n] = sim_power(n, shift, shape, par, [20261027, case_idx, n, code], nexp)
        return cache[n][0]

    if P(2) >= 0.8:
        return 2, cache
    lo, hi = 2, max(3, 4 * n_calc)
    while P(hi) < 0.8:                       # not expected; widen if it happens and say so
        lo, hi = hi, 2 * hi
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if P(mid) >= 0.8:
            hi = mid
        else:
            lo = mid
    return hi, cache


def run(shape, subset, nexp):
    par_all = shapes()
    with open(os.path.join(OUT, "experiment.csv")) as f:
        rows = list(csv.DictReader(f))
    outp = os.path.join(OUT, f"simulation_{shape}_{subset}.csv")
    res = []
    t_start = time.time()
    for idx, r in enumerate(rows):
        if r["shape"] != shape:
            continue
        dT = float(r["dT_K"])
        central = r["central"] == "True"
        if subset == "required" and not (central or (r["group"] == "C1" and dT in (2.0, 5.0))):
            continue
        n_calc = int(r["n"])
        s = float(r["s_per_K"])
        par = par_all[r["floor"]]
        t0 = time.time()
        n_sim, cache = bisect_n(n_calc, -s * dT, shape, par, idx, nexp)
        p_at = cache[n_sim]
        p_prev = cache.get(n_sim - 1)
        if p_prev is None and n_sim > 2:
            p_prev = sim_power(n_sim - 1, -s * dT, shape, par, [20261027, idx, n_sim - 1, SHAPE_CODE[shape]], nexp)
        rec = {"case_index": idx, "group": r["group"], "lever": r["lever"], "s_per_K": s, "dT_K": dT, "floor": r["floor"],
               "shape": shape, "central": central, "n_calc": n_calc, "n_sim": n_sim,
               "rel_diff": (n_sim - n_calc) / n_calc, "V4_pass_5pct": abs(n_sim - n_calc) <= 0.05 * n_calc,
               "sim_power_at_n_sim": p_at[0], "sim_wrong_sign_at_n_sim": p_at[1],
               "sim_power_at_n_sim_minus_1": None if p_prev is None else p_prev[0],
               "sim_power_at_n_calc": cache[n_calc][0] if n_calc in cache else None,
               "calc_power_at_n_calc": float(r["power_at_n"]), "n_tried": len(cache), "seconds": time.time() - t0}
        res.append(rec)
        print(f"[{time.time()-t_start:7.1f} s] {r['group']:22s} {r['lever']:30.30s} dT={dT:4g} {r['floor']:6s} n_calc={n_calc:5d} "
              f"n_sim={n_sim:5d} ({100*rec['rel_diff']:+.1f} %) pass={rec['V4_pass_5pct']}", flush=True)
        with open(outp, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(res[0].keys()))
            w.writeheader()
            for x in res:
                w.writerow({k: repr(v) if isinstance(v, float) else v for k, v in x.items()})


def precise():
    """P7's two cases under the registered shape, Geslin floor: 200,000 experiments per n."""
    par = shapes()["Geslin"]
    warm, cold, _ = LV.load()
    cw = sorted(warm, key=lambda lv: abs(lv["s"]))[len(warm) // 2]
    cc = sorted(cold, key=lambda lv: abs(lv["s"]))[len(cold) // 2]
    out = []
    for name, lv, dT, ns in (("warm chamber central, dT = 2 K", cw, 2.0, range(18, 41)),
                             ("cold chamber central, dT = 5 K", cc, 5.0, range(2, 8))):
        for n in ns:
            for shape in ("lognormal", "normal"):
                p, wrong = sim_power(n, -lv["s"] * dT, shape, par, [20261027, 999, n, SHAPE_CODE[shape], 7], 200000)
                eff = abs(lv["s"]) * dT / (par["sd_ln"] if shape == "lognormal" else par["sd_norm"])
                out.append({"case": name, "lever": lv["lever"], "s_per_K": lv["s"], "dT_K": dT, "shape": shape, "n": n,
                            "sim_power_200k": p, "mc_se": math.sqrt(p * (1 - p) / 200000), "wrong_sign": wrong,
                            "calc_power": PW.power_quad(n, eff)[0]})
                print(out[-1], flush=True)
    with open(os.path.join(OUT, "simulation_precise.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader()
        for x in out:
            w.writerow(x)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--shape", default="lognormal", choices=["lognormal", "normal"])
    ap.add_argument("--subset", default="required", choices=["required", "all"])
    ap.add_argument("--nexp", type=int, default=10000)
    ap.add_argument("--precise", action="store_true")
    a = ap.parse_args()
    if a.precise:
        precise()
    else:
        run(a.shape, a.subset, a.nexp)
