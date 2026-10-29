#!/usr/bin/env python3
"""Stage 5: the code on made-up pairs with a planted within-set lever (DECISIONS.md D76), before the archived pairs.

(1) Without noise, the pairwise slope recovers a planted lever exactly and beta is exactly one; with lifetimes that do
not depend on temperature beta is zero. (2) For complete three-cell sets the pairwise slope through the origin equals
the within-set (fixed-effects) slope, and F is 1.5; for two-cell sets F is 1 (D73, D75). (3) With noise, the bootstrap
over sets covers the planted lever in most of 20 datasets, and beta's interval covers one. (4) The calculated
detectable |s| agrees with a simulation of independent normal pair differences within 3 %. (5) The sign test.

Usage: python3 tests/stage5_synthetic.py
"""
import math, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "src", "stage5"))
import s5_power as PW                   # noqa: E402
import s5_within as W                   # noqa: E402


def make_pairs(rng, n_sets, cells, s_true, s_between, t_sd, cell_sd, dataset="KIT"):
    rows, fe_num, fe_den = [], 0.0, 0.0
    for j in range(n_sets):
        T = 25.0 + t_sd * rng.standard_normal(cells)
        lnL = 6.0 - s_true * (T - T.mean()) + cell_sd * rng.standard_normal(cells)
        fe_num += float(((T - T.mean()) * (lnL - lnL.mean())).sum())
        fe_den += float(((T - T.mean()) ** 2).sum())
        for a in range(cells):
            for b in range(a + 1, cells):
                w, c = (a, b) if T[a] > T[b] else (b, a)
                dT = float(T[w] - T[c])
                d = float(lnL[w] - lnL[c])
                rows.append({"dataset": dataset, "set": f"S{j:03d}", "warm": f"S{j:03d}-{w}", "cool": f"S{j:03d}-{c}",
                             "dT_K": dT, "d": d, "p": -s_between * dT, "tie": False, "band": "warm"})
    return rows, -fe_num / fe_den


def main():
    ok = True
    rng = np.random.default_rng(5)
    # (1) exact recovery without noise
    rows, fe = make_pairs(rng, 30, 3, 0.04, 0.04, 0.6, 0.0)
    sw, beta = W.sw_point(rows), W.beta_point(rows)
    t1 = abs(sw - 0.04) < 1e-12 and abs(beta - 1.0) < 1e-12
    rows0, _ = make_pairs(rng, 30, 3, 0.0, 0.04, 0.6, 0.0)
    t1b = abs(W.beta_point(rows0)) < 1e-12 and abs(W.sw_point(rows0)) < 1e-12
    print("(1) exact recovery:", t1, t1b, sw, beta)
    ok &= t1 and t1b
    # (2) the fixed-effects identity and F
    rows, fe = make_pairs(rng, 25, 3, 0.03, 0.03, 0.5, 0.05)
    t2 = abs(W.sw_point(rows) - fe) < 1e-12 * max(1.0, abs(fe))
    F3 = PW.shared_cell_F(rows)
    rows2, _ = make_pairs(rng, 25, 2, 0.03, 0.03, 0.5, 0.05)
    F2 = PW.shared_cell_F(rows2)
    rows4, _ = make_pairs(rng, 25, 4, 0.03, 0.03, 0.5, 0.05)
    F4 = PW.shared_cell_F(rows4)
    t2b = abs(F3 - 1.5) < 1e-12 and abs(F2 - 1.0) < 1e-12 and abs(F4 - 2.0) < 1e-12
    print("(2) fixed-effects identity:", t2, "F (3, 2, 4 cells):", F3, F2, F4, t2b)
    ok &= t2 and t2b
    # (3) coverage of the bootstrap over sets
    cover_sw, cover_beta = 0, 0
    for k in range(20):
        r, _ = make_pairs(np.random.default_rng(100 + k), 40, 3, 0.03, 0.03, 0.5, 0.05)
        a = W.analyse(r, 900 + k, ("beta", "sw", "sign"))
        lo, hi = a["s_w_boot"]["ci95"]
        cover_sw += lo <= 0.03 <= hi
        lo, hi = a["beta_boot"]["ci95"]
        cover_beta += lo <= 1.0 <= hi
    t3 = cover_sw >= 16 and cover_beta >= 16
    print("(3) bootstrap coverage of 20:", cover_sw, cover_beta, t3)
    ok &= t3
    # (4) calculation against a simulation with independent normal pair differences
    dT = np.abs(np.random.default_rng(7).normal(0.0, 0.6, 120))
    rms = math.sqrt(float(dT @ dT) / len(dT))
    calc = PW.detectable_calc(len(dT), rms, 0.1)
    sim, _ = PW.FixedDesign(dT, {"kind": "normal", "sd_e": 0.1}, [1, 2, 3]).detectable(calc)
    t4 = abs(sim - calc) <= 0.03 * calc
    print("(4) detectable |s| calculated and simulated:", calc, sim, t4)
    ok &= t4
    # (5) the sign test
    agree = [{"d": 1.0, "p": 1.0}] * 30
    half = [{"d": 1.0, "p": 1.0}] * 15 + [{"d": -1.0, "p": 1.0}] * 15
    s1, s2 = W.sign_test(agree), W.sign_test(half)
    t5 = s1["p_two_sided"] < 1e-8 and abs(s2["p_two_sided"] - 1.0) < 1e-12 and s1["share"] == 1.0 and s2["share"] == 0.5
    print("(5) sign test:", s1["p_two_sided"], s2["p_two_sided"], t5)
    ok &= t5
    print("stage5_synthetic pass:", bool(ok))
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
