#!/usr/bin/env python3
"""Checks of the Stage 4 code on made-up inputs (no archived lever, floor or Stage 4 result is used).

1. The paired t-test sample size of src/stage4/experiment.py against published values for a one-sample t-test,
   two-sided 5 %, power 0.80 (Cohen's d = 0.2, 0.5, 0.8 and 1.0 give 199, 34, 15 and 10 in G*Power's matched-pairs
   test), and against a simulation with normal noise.
2. The log-normal fit returns the mean and 95th percentile it was given, for made-up pairs (m, q).
3. The implied divergence's median, 95th percentile, interval and verdicts on a made-up module and made-up levers,
   against values written out by hand.
4. Stage 1's network on a made-up heat: the four corner cells share one temperature and the four inner cells another,
   and the energy balance closes.

Usage: python3 tests/stage4_synthetic.py
"""
import math, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src", "stage4"))
sys.path.insert(0, os.path.join(ROOT, "src"))
import experiment as EX                 # noqa: E402
import divergence as DV                 # noqa: E402
import thermal_lever as TL              # noqa: E402

fails = []


def check(name, cond, detail=""):
    print(("ok   " if cond else "FAIL ") + name + (f": {detail}" if detail else ""))
    if not cond:
        fails.append(name)


# 1. sample sizes
for d, want in ((0.2, 199), (0.5, 34), (0.8, 15), (1.0, 10)):
    n = EX.n_required(d, 1.0)
    check(f"n for d = {d}", n == want, f"{n} against {want}")
pw = [EX.power_t(n, d, 1.0) for d in (4.0, 6.0, 8.0, 10.0) for n in (2, 3, 4)]
check("power defined and below 1 at large noncentrality", all(math.isfinite(p) and 0.0 < p <= 1.0 for p in pw))
ns_d = [EX.n_required(d, 1.0) for d in np.linspace(0.05, 12.0, 240)]
check("n never rises with the effect", all(a >= b for a, b in zip(ns_d, ns_d[1:])), f"{ns_d[-5:]}")
check("n at d = 8.1 (nc 14 at n = 3) is 3", EX.n_required(8.1, 1.0) == 3, f"{EX.n_required(8.1, 1.0)}")
sim = EX.sim_power(34, 0.5, {"kind": "normal", "sd_e": 1.0}, [1, 2, 3])
check("simulated power at d = 0.5, n = 34 (normal noise)", abs(sim - EX.power_t(34, 0.5, 1.0)) < 0.015,
      f"{sim:.4f} against {EX.power_t(34, 0.5, 1.0):.4f}")
ns, _ = EX.n_simulated(0.5, {"kind": "normal", "sd_e": 1.0}, 34, 99, 2)
check("simulated n at d = 0.5 (normal noise)", abs(ns - 34) <= 2, f"{ns}")

# 2. log-normal fit
try:                                   # a log-normal's 95th percentile is at most exp(z^2 / 2) = 3.87 times its mean
    EX.lognormal_shape(0.01, 0.05)
    check("no log-normal for q / m = 5", False)
except SystemExit:
    check("no log-normal for q / m = 5", True)
for m, q in ((0.03, 0.09), (0.07, 0.2), (0.01, 0.035)):
    s = EX.lognormal_shape(m, q)
    check(f"log-normal fit m = {m}, q = {q}", abs(s["mean_check"] - m) < 1e-12 and abs(s["p95_check"] - q) < 1e-12,
          f"sigma {s['sigma']:.4f}")
    z = np.random.default_rng(5).standard_normal(400000)
    x = np.exp(s["mu"] + s["sigma"] * z)
    check(f"  sample mean, 95th percentile and rms", abs(x.mean() / m - 1) < 0.01 and abs(np.percentile(x, 95) / q - 1) < 0.01
          and abs(math.sqrt((x ** 2).mean()) / s["sd_e"] - 1) < 0.02)

# 3. divergence on made-up numbers: spread 2 K, levers -0.03, -0.01, +0.02 (warm), floors mean 0.03, p95 0.05
mod = [{"ambient_C": 35.0, "crate": 1.0, "rate": "1C", "h": 10.0, "band": "warm", "spread_K": 2.0, "T_min_C": 40.0,
        "T_max_C": 42.0, "T_mean_C": 41.0, "cells_inside_band_setpoints": "False",
        **{f"offset_cell{i}_K": (1.0 if i in (1, 2, 5, 6) else -1.0) for i in range(8)}}]
levs = [{"label": "a", "source": "x", "secant": "", "s": -0.03, "lo": -0.04, "hi": -0.02},
        {"label": "b", "source": "x", "secant": "", "s": -0.01, "lo": -0.02, "hi": 0.005},
        {"label": "c", "source": "x", "secant": "", "s": 0.02, "lo": 0.01, "hi": 0.03}]
dd = {"bands": {"warm": levs}, "floors": {"F": {"mean": 0.03, "p95": 0.05}}}
rows, cases = DV.divergence(mod, dd)
c = cases[0]
check("median of |s| x spread", abs(c["median"] - 0.04) < 1e-15, f"{c['median']}")      # sorted 0.02, 0.04, 0.06
check("95th percentile", abs(c["p95"] - (0.04 + 0.9 * 0.02)) < 1e-15, f"{c['p95']}")    # position 1.9
check("offers and detectable", c["offers_F"] is True and c["detectable_F"] is False)
b = [r for r in rows if r.get("lever") == "b"][0]
check("interval holding zero", b["divergence_lo"] == 0.0 and abs(b["divergence_hi"] - 0.04) < 1e-15)
a = [r for r in rows if r.get("lever") == "a"][0]
check("cell offsets in ln L", abs(a["offset_lnL_cell1"] - 0.03) < 1e-15 and abs(a["offset_lnL_cell0"] + 0.03) < 1e-15,
      "a hot cell under a negative lever lives longer")

# 4. the network's symmetry on a made-up heat
T, bal, info = TL.module_steady(0.3, 10.0, T_amb=290.0)
check("corner cells equal", np.ptp(T[[0, 3, 4, 7]]) < 1e-12, f"{np.ptp(T[[0, 3, 4, 7]]):.1e}")
check("inner cells equal", np.ptp(T[[1, 2, 5, 6]]) < 1e-12, f"{np.ptp(T[[1, 2, 5, 6]]):.1e}")
check("energy balance", abs(bal) < 1e-9, f"{bal:.1e} W")

print("all passed" if not fails else f"{len(fails)} failed: {fails}")
sys.exit(1 if fails else 0)
