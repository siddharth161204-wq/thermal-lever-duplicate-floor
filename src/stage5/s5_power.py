#!/usr/bin/env python3
"""Stage 5, Part C: could the data have shown it (DECISIONS.md D75), with V4 (D77).

The test is B2's: the least-squares slope b of d on dT through the origin, t = b sqrt(sum dT^2) / s_r with
s_r^2 = sum (d - b dT)^2 / (n - 1), two-sided at 5 % against t with n - 1 degrees of freedom.
C1, calculation: Stage 4's power_t(n, |s| dT_rms, sd(e)), the noncentral t with noncentrality |s| sqrt(sum dT^2) / sd(e);
the detectable |s| is where that power is 0.80 (Brent, 1e-12 relative). Simulation: 10,000 synthetic datasets with the
observed dT values and e from Stage 4's registered shape (log-normal magnitude, random sign), d = |s| dT + e; the
residuals do not depend on s, so one set of draws serves every s, and the smallest s whose simulated power reaches
0.80 is found by bisection (1e-6 relative). Beside: normal e of the same sd(e), the normal shape of D61, the shared-cell
line (sd(e) sqrt(F)). C2: per band, the power at the band's median between-set |s| and the pairs needed at the band's
dT spread, by the calculation and by simulation (each synthetic dataset drawing its n dT values from the band's).

Writes results/stage5_power.json and verification/stage5/v4.json.

Usage: python3 src/stage5/s5_power.py
"""
import json, math, os, sys, time

import numpy as np
from scipy import optimize, stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s5_inputs as A                   # noqa: E402
import s5_within as W                   # noqa: E402

EX4 = A.EX4
ALPHA, POWER, N_EXP = 0.05, 0.80, 10000
CHUNK = 5_000_000


def dT_of(rows, tkey="dT_K"):
    return np.array([r[tkey] for r in rows], dtype=float)


def shared_cell_F(rows, tkey="dT_K"):
    """F = sum over cells of w^2 / (2 sum over pairs of dT^2); w = the cell's dT as the warm cell minus as the cool one."""
    w = {}
    for r in rows:
        kw, kc = (r["dataset"], r["warm"]), (r["dataset"], r["cool"])
        w[kw] = w.get(kw, 0.0) + r[tkey]
        w[kc] = w.get(kc, 0.0) - r[tkey]
    return sum(v * v for v in w.values()) / (2.0 * sum(r[tkey] ** 2 for r in rows))


def detectable_calc(n, dT_rms, sd):
    """The |s| at which Stage 4's power_t(n, |s| dT_rms, sd) equals 0.80."""
    f = lambda s: EX4.power_t(n, s * dT_rms, sd) - POWER
    hi = 1e-4
    while f(hi) < 0:
        hi *= 2.0
    return optimize.brentq(f, 0.0, hi, xtol=1e-18, rtol=1e-12, maxiter=500)


def draws(shape, size, rng):
    if shape["kind"] == "lognormal":
        mag = np.exp(shape["mu"] + shape["sigma"] * rng.standard_normal(size))
        sign = rng.integers(0, 2, size=size, dtype=np.int8).astype(np.float64) * 2.0 - 1.0
        return sign * mag
    return shape["sd_e"] * rng.standard_normal(size)


class FixedDesign:
    """10,000 synthetic datasets at the observed dT values; the statistic's pieces that do not depend on s."""

    def __init__(self, dT, shape, seed):
        rng = np.random.default_rng(seed)
        n = len(dT)
        e = draws(shape, (N_EXP, n), rng)
        S = float(dT @ dT)
        self.b_e = e @ dT / S
        r = e - self.b_e[:, None] * dT[None, :]
        s_r = np.sqrt((r * r).sum(axis=1) / (n - 1))
        tc = stats.t.ppf(1.0 - ALPHA / 2.0, n - 1)
        self.c = tc * s_r / math.sqrt(S)

    def power(self, s):
        b = s + self.b_e
        return float(np.mean((b > self.c) | (b < -self.c)))

    def detectable(self, s_start):
        lo, hi = 0.0, max(s_start, 1e-6)
        while self.power(hi) < POWER:
            lo, hi = hi, hi * 2.0
        while hi - lo > 1e-6 * hi:
            mid = 0.5 * (lo + hi)
            if self.power(mid) >= POWER:
                hi = mid
            else:
                lo = mid
        return hi, self.power(hi)


def sim_power_resampled(n, s, pool, shape, seed):
    """Simulated power at n pairs, each dataset drawing its n dT values with replacement from `pool`."""
    rng = np.random.default_rng(seed)
    tc = stats.t.ppf(1.0 - ALPHA / 2.0, n - 1)
    rej, done = 0, 0
    per = max(1, CHUNK // n)
    m = len(pool)
    while done < N_EXP:
        k = min(per, N_EXP - done)
        x = pool[rng.integers(0, m, size=(k, n))]
        d = s * x + draws(shape, (k, n), rng)
        Sxx = (x * x).sum(axis=1)
        b = (x * d).sum(axis=1) / Sxx
        r = d - b[:, None] * x
        s_r = np.sqrt((r * r).sum(axis=1) / (n - 1))
        t = b * np.sqrt(Sxx) / s_r
        rej += int((np.abs(t) > tc).sum())
        done += k
    return rej / N_EXP


def n_simulated_resampled(s, pool, shape, n_calc, case):
    """Stage 4's bisection: from 2, between 2 and four times the calculated n, doubling beyond it if needed."""
    tried = {}

    def ok(n):
        if n not in tried:
            tried[n] = sim_power_resampled(n, s, pool, shape, [A.SEED, 200 + case, n])
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


def shapes_of(noise):
    out = {}
    for name, f in noise.items():
        out[name] = {"lognormal": {"kind": "lognormal", "mu": f["mu"], "sigma": f["sigma"], "sd_e": f["sd_e"]},
                     "normal_same_sd": {"kind": "normal", "sd_e": f["sd_e"]},
                     "normal_shape_sd_e": f["sd_e_normal_shape"]}
    return out


def c1_case(name, rows, floor, sh, case, registered, tkey="dT_K"):
    dT = dT_of(rows, tkey)
    n = len(dT)
    rms = math.sqrt(float(dT @ dT) / n)
    sd = sh["lognormal"]["sd_e"]
    F = shared_cell_F(rows, tkey)
    calc = detectable_calc(n, rms, sd)
    t0 = time.time()
    s_ln, pw_ln = FixedDesign(dT, sh["lognormal"], [A.SEED, 100 + case, 1]).detectable(calc)
    s_nm, pw_nm = FixedDesign(dT, sh["normal_same_sd"], [A.SEED, 100 + case, 2]).detectable(calc)
    out = {"case": name, "registered": registered, "floor": floor, "pairs": n, "sets": len({(r["dataset"], r["set"]) for r in rows}),
           "dT_rms_K": rms, "sum_dT2_K2": float(dT @ dT), "dT_median_K": float(np.median(dT)), "sd_e": sd,
           "detectable_s_per_K": calc, "power_check_at_detectable": EX4.power_t(n, calc * rms, sd),
           "detectable_s_normal_shape_per_K": detectable_calc(n, rms, sh["normal_shape_sd_e"]),
           "shared_cell_F": F, "detectable_s_shared_cells_per_K": detectable_calc(n, rms, sd * math.sqrt(F)),
           "sim_lognormal_detectable_s_per_K": s_ln, "sim_lognormal_power_there": pw_ln,
           "sim_normal_same_sd_detectable_s_per_K": s_nm, "sim_normal_same_sd_power_there": pw_nm,
           "v4_rel_diff": (s_ln - calc) / calc, "v4_pass": abs(s_ln - calc) <= 0.05 * calc,
           "normal_check_rel_diff": (s_nm - calc) / calc, "normal_check_pass": abs(s_nm - calc) <= 0.05 * calc,
           "seeds": {"lognormal": [A.SEED, 100 + case, 1], "normal_same_sd": [A.SEED, 100 + case, 2]},
           "seconds": round(time.time() - t0, 1)}
    return out


def c2_case(dataset, band, rows, floor, sh, med, case, registered):
    dT = dT_of(rows)
    n = len(dT)
    rms = math.sqrt(float(dT @ dT) / n)
    sd = sh["lognormal"]["sd_e"]
    s = med["abs_s_per_K"]
    F = shared_cell_F(rows)
    t0 = time.time()
    n_calc = EX4.n_required(s * rms, sd)
    n_sim, tried = n_simulated_resampled(s, dT, sh["lognormal"], n_calc, case)
    fd = FixedDesign(dT, sh["lognormal"], [A.SEED, 300 + case, 1])
    return {"dataset": dataset, "band": band, "registered": registered, "floor": floor, "pairs": n, "dT_rms_K": rms,
            "band_median": med["label"], "abs_s_per_K": s, "sd_e": sd,
            "power_calc": EX4.power_t(n, s * rms, sd), "power_sim_lognormal": fd.power(s),
            "pairs_needed_calc": n_calc, "power_calc_at_needed": EX4.power_t(n_calc, s * rms, sd),
            "pairs_needed_normal_shape": EX4.n_required(s * rms, sh["normal_shape_sd_e"]),
            "shared_cell_F": F, "pairs_needed_shared_cells": EX4.n_required(s * rms, sd * math.sqrt(F)),
            "pairs_needed_sim_lognormal_resampled": n_sim, "sim_power_at_needed": tried[n_sim],
            "sim_rel_diff": (n_sim - n_calc) / n_calc,
            "seeds": {"power_sim": [A.SEED, 300 + case, 1], "n_sim": f"[{A.SEED}, {200 + case}, n]"},
            "seconds": round(time.time() - t0, 1)}


def curves(sets, noise):
    grid = sorted({int(round(x)) for x in np.logspace(math.log10(3), math.log10(3000), 90)})
    out = []
    for name, rows in sets:
        dT = dT_of(rows)
        rms = math.sqrt(float(dT @ dT) / len(dT))
        for floor, f in noise.items():
            for n in grid:
                out.append({"dataset": name, "floor": floor, "pairs": n, "dT_rms_K": rms,
                            "detectable_s_per_K": detectable_calc(n, rms, f["sd_e"])})
    return out


def run(P, noise, medians):
    S = W.subsets(P)
    band = S["band"]
    sh = shapes_of(noise)
    c1 = []
    case = 0
    for name, rows in (("KIT", S["kit"]), ("Sandia", S["snl"]), ("pooled", S["pool"])):
        for floor in ("KIT", "Geslin"):
            c1.append(c1_case(name, rows, floor, sh[floor], case, True))
            case += 1
    for name, rows in (("KIT", S["kit"]), ("Sandia", S["snl"]), ("pooled", S["pool"])):
        for b in ("warm", "cold"):
            for floor in ("KIT", "Geslin"):
                c1.append(c1_case(f"{name} {b}", band(rows, b), floor, sh[floor], case, False))
                case += 1
    allkit = [r for r in P if r["dataset"] == "KIT"]
    for floor in ("KIT", "Geslin"):
        c1.append(c1_case("KIT all 142 pairs", allkit, floor, sh[floor], case, False))
        case += 1
    print("C1 done", flush=True)
    c2 = []
    case = 0
    for name, rows in (("KIT", S["kit"]), ("Sandia", S["snl"]), ("pooled", S["pool"])):
        for b in ("warm", "cold"):
            for floor in ("KIT", "Geslin"):
                reg = name == "KIT" and floor == "KIT"
                c2.append(c2_case(name, b, band(rows, b), floor, sh[floor], medians[b], case, reg))
                print("C2", name, b, floor, c2[-1]["pairs_needed_calc"], c2[-1]["pairs_needed_sim_lognormal_resampled"],
                      c2[-1]["seconds"], "s", flush=True)
                case += 1
    out = {"method": "D75", "alpha": ALPHA, "power": POWER, "synthetic_datasets": N_EXP, "seed_base": A.SEED,
           "shapes": {k: {"mu": v["mu"], "sigma": v["sigma"], "sd_e": v["sd_e"], "sd_e_normal_shape": v["sd_e_normal_shape"]}
                      for k, v in noise.items()},
           "band_medians": medians, "C1": c1, "C2": c2,
           "curves": curves((("KIT", S["kit"]), ("Sandia", S["snl"])), noise)}
    v4 = {"check": "V4 (D77): simulated detectable |s| under the registered shape within 5 % of the calculated one",
          "cases": len(c1), "pass_all": all(x["v4_pass"] for x in c1),
          "pass_registered": all(x["v4_pass"] for x in c1 if x["registered"]),
          "failures": [x["case"] + " / " + x["floor"] for x in c1 if not x["v4_pass"]],
          "normal_check_pass_all": all(x["normal_check_pass"] for x in c1),
          "rows": [{k: x[k] for k in ("case", "floor", "registered", "pairs", "detectable_s_per_K", "sim_lognormal_detectable_s_per_K",
                                      "v4_rel_diff", "v4_pass", "sim_normal_same_sd_detectable_s_per_K", "normal_check_rel_diff",
                                      "normal_check_pass")} for x in c1]}
    out["V4"] = {k: v4[k] for k in ("cases", "pass_all", "pass_registered", "failures", "normal_check_pass_all")}
    return out, v4


if __name__ == "__main__":
    d = A.load_all()
    out, v4 = run(d["pairs"], d["noise"], d["medians"])
    json.dump(out, open(os.path.join(A.ROOT, "results", "stage5_power.json"), "w"), indent=1)
    vdir = os.path.join(A.ROOT, "verification", "stage5")
    os.makedirs(vdir, exist_ok=True)
    json.dump(v4, open(os.path.join(vdir, "v4.json"), "w"), indent=1)
    for x in out["C1"]:
        print(f"{x['case']:24s} {x['floor']:6s} n {x['pairs']:4d} rms {x['dT_rms_K']:.3f} calc {100 * x['detectable_s_per_K']:.3f} "
              f"sim {100 * x['sim_lognormal_detectable_s_per_K']:.3f} normal {100 * x['sim_normal_same_sd_detectable_s_per_K']:.3f} "
              f"shared {100 * x['detectable_s_shared_cells_per_K']:.3f} (F {x['shared_cell_F']:.3f}) V4 {x['v4_pass']}")
    print(json.dumps(out["V4"], indent=1))
