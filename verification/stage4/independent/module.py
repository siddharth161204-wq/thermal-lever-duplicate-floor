#!/usr/bin/env python3
"""Independent check of Stage 4, part (b): the 36 module cases (D58), with the checker's own solver.

The network is written here from the physics of DECISIONS.md D5 to D12 and the constants that src/thermal_lever.py
states (D9 to D11), not by calling its solver:
  two rows of four LG M50 21700 cells (D = 21.1 mm, L = 70.1 mm) at 23.1 mm pitch, a 2 mm gap;
  per neighbour pair (horizontal and vertical, ten pairs, no diagonals): one nickel strip 90 W/mK x 0.2 mm x 8 mm over
  the pitch, air-gap conduction 0.03 W/mK x (A_lat / 4) / gap, and grey exchange eps^2 sigma F A_lat (Ti^4 - Tj^4)
  with F the view factor of parallel equal cylinders per whole lateral area and eps = 0.8;
  to the ambient: convection h over A_outer minus A_lat / 4 per neighbour, radiation eps sigma over A_outer minus
  n eps F A_lat (n neighbours), the ambient being the temperature of the air and of the surroundings (D58);
  heat per cell (c x 5 A)^2 x R(ambient).

Primary solution: the symmetric reduction (four corner cells with one corner and one inner neighbour, four inner cells
with one corner and two inner neighbours) solved by nested bracketing root finding (scipy brentq), which shares no
algorithm with module_steady's Newton iteration. Reference: a Newton solve of the full eight-node system in mpmath at
50 significant digits, which does not assume the symmetry. Cross-check: thermal_lever.module_steady called as the
pipeline is said to call it (D58).

Usage: python3 verification/stage4/independent/module.py
Output: verification/stage4/independent/out/module.csv and module.json
"""
import csv, json, math, os, sys

import numpy as np
from scipy.optimize import brentq
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = os.path.join(HERE, "out")

# ---------------- constants, written from D6, D7, D9 to D11 ----------------
SIGMA = 5.670374419e-8
D, L = 0.0211, 0.0701
CAP = 5.0
PITCH = 0.0231
GAP = PITCH - D
K_AIR = 0.03
EPS = 0.8
A_LAT = math.pi * D * L
A_END = math.pi * (D / 2) ** 2
A_OUT = A_LAT + 2 * A_END
A_FACE = A_LAT / 4
G_STRIP = 90.0 * 0.2e-3 * 8e-3 / PITCH          # one strip per neighbour pair (D9)
G_GAP = K_AIR * A_FACE / GAP
X = 1.0 + GAP / D
F = (math.sqrt(X * X - 1.0) + math.asin(1.0 / X) - X) / math.pi
G_PAIR = G_STRIP + G_GAP
R_PAIR = EPS ** 2 * SIGMA * F * A_LAT

PAIRS = [(0, 1), (1, 2), (2, 3), (4, 5), (5, 6), (6, 7), (0, 4), (1, 5), (2, 6), (3, 7)]
NNB = [sum(1 for p in PAIRS if i in p) for i in range(8)]
CORNERS, INNERS = [0, 3, 4, 7], [1, 2, 5, 6]
assert all(NNB[i] == 2 for i in CORNERS) and all(NNB[i] == 3 for i in INNERS)
A_CONV = [A_OUT - n * A_FACE for n in NNB]
A_RAD = [A_OUT - n * EPS * F * A_LAT for n in NNB]

AMBIENTS_C = [35.0, 25.0, 10.0]
RATES = [0.25, 0.5, 1.0, 2.0]
H_LIST = [5.0, 10.0, 25.0]
R35_ARCHIVED = 0.01871279570954531


def loss_amb(T, Ta, h, i):
    return h * A_CONV[i] * (T - Ta) + EPS * SIGMA * A_RAD[i] * (T ** 4 - Ta ** 4)


def pair_flow(Ta_, Tb_):
    return G_PAIR * (Ta_ - Tb_) + R_PAIR * (Ta_ ** 4 - Tb_ ** 4)


def solve_symmetric(q, h, Ta):
    """Corner temperature Tc and inner temperature Ti by nested brentq.
    Corner balance: q = loss(Tc) + flow(Tc -> Ti); inner balance: q = loss(Ti) + flow(Ti -> Tc)."""
    ic, ii = CORNERS[0], INNERS[0]

    def Ti_of_Tc(Tc):
        # flow(Tc, Ti) = q - loss_c(Tc); flow is decreasing in Ti, so bracket Ti above Tc
        target = q - loss_amb(Tc, Ta, h, ic)
        g = lambda Ti: pair_flow(Tc, Ti) - target
        lo, hi = Tc - 1e-9, Tc + 1.0
        while g(hi) > 0:
            hi = Tc + 2 * (hi - Tc)
        while g(lo) < 0:
            lo = Tc - 2 * (Tc - lo) - 1e-6
        return brentq(g, lo, hi, xtol=1e-15, rtol=4 * np.finfo(float).eps, maxiter=500)

    def resid(Tc):
        Ti = Ti_of_Tc(Tc)
        return loss_amb(Ti, Ta, h, ii) + pair_flow(Ti, Tc) - q

    # the total balance 2q = loss_c(Tc) + loss_i(Ti) puts Tc between Ta and the single-cell temperature of an exposed cell
    lo = Ta
    hi = Ta + 1.0
    while loss_amb(hi, Ta, h, ic) < 2 * q:
        hi = Ta + 2 * (hi - Ta)
    Tc = brentq(resid, lo, hi, xtol=1e-15, rtol=4 * np.finfo(float).eps, maxiter=500)
    Ti = Ti_of_Tc(Tc)
    T = np.empty(8)
    T[CORNERS] = Tc
    T[INNERS] = Ti
    return T


def solve_mp(q, h, Ta, T0, dps=50):
    """Full eight-node Newton in mpmath (no symmetry assumed)."""
    mp.mp.dps = dps
    qm, hm, Tam = mp.mpf(q), mp.mpf(h), mp.mpf(Ta)
    sig, eps = mp.mpf(SIGMA), mp.mpf(EPS)
    Dm, Lm, pitch = mp.mpf("0.0211"), mp.mpf("0.0701"), mp.mpf("0.0231")
    gap = pitch - Dm
    alat = mp.pi * Dm * Lm
    aend = mp.pi * (Dm / 2) ** 2
    aout = alat + 2 * aend
    aface = alat / 4
    Xm = 1 + gap / Dm
    Fm = (mp.sqrt(Xm * Xm - 1) + mp.asin(1 / Xm) - Xm) / mp.pi
    gpair = mp.mpf(90) * mp.mpf("0.2e-3") * mp.mpf("8e-3") / pitch + mp.mpf("0.03") * aface / gap
    rpair = eps ** 2 * sig * Fm * alat
    aconv = [aout - n * aface for n in NNB]
    arad = [aout - n * eps * Fm * alat for n in NNB]
    T = [mp.mpf(float(t)) for t in T0]
    for it in range(60):
        Fv = [qm - hm * aconv[i] * (T[i] - Tam) - eps * sig * arad[i] * (T[i] ** 4 - Tam ** 4) for i in range(8)]
        J = mp.zeros(8, 8)
        for i in range(8):
            J[i, i] = -hm * aconv[i] - 4 * eps * sig * arad[i] * T[i] ** 3
        for a, b in PAIRS:
            fl = gpair * (T[a] - T[b]) + rpair * (T[a] ** 4 - T[b] ** 4)
            Fv[a] -= fl
            Fv[b] += fl
            da = gpair + 4 * rpair * T[a] ** 3
            db = gpair + 4 * rpair * T[b] ** 3
            J[a, a] -= da; J[a, b] += db; J[b, a] += da; J[b, b] -= db
        dT = mp.lu_solve(J, mp.matrix([-f for f in Fv]))
        T = [T[i] + dT[i] for i in range(8)]
        if max(abs(dT[i]) for i in range(8)) < mp.mpf(10) ** (-(dps - 8)):
            break
    return T, it + 1


def balance(T, q, h, Ta):
    return float(8 * q - sum(loss_amb(T[i], Ta, h, i) for i in range(8)))


def main():
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "resistance.json")) as f:
        RES = json.load(f)
    Rmap = {35.0: R35_ARCHIVED, 25.0: RES["ambients"]["25"]["R_registered_ohm"], 10.0: RES["ambients"]["10"]["R_registered_ohm"]}
    sys.path.insert(0, os.path.join(ROOT, "src"))
    import thermal_lever as TL                       # cross-check only
    const_check = {"G_strip": (G_STRIP, TL.G_STRIP), "G_gap": (G_GAP, TL.G_GAP), "F_view": (F, TL.F_VIEW),
                   "A_lat": (A_LAT, TL.A_LAT), "A_outer": (A_OUT, TL.A_OUTER), "A_face": (A_FACE, TL.A_FACE),
                   "R_DC_35": (R35_ARCHIVED, TL.R_DC)}
    rows = []
    for Tc_amb in AMBIENTS_C:
        Ta = Tc_amb + 273.15
        R = Rmap[Tc_amb]
        for c in RATES:
            q = (c * CAP) ** 2 * R
            for h in H_LIST:
                T = solve_symmetric(q, h, Ta)
                Tm, itm = solve_mp(q, h, Ta, T)
                Tm_f = np.array([float(x) for x in Tm])
                spread_mp = float(max(Tm) - min(Tm))
                mean_mp = float(sum(Tm) / 8)
                # cross-check: Stage 1's solver
                Tp, bal_p, info_p = TL.module_steady(q, h, T_amb=Ta)
                spread = float(T.max() - T.min())
                mean = float(T.mean())
                r = {"T_amb_C": Tc_amb, "crate": c, "h": h, "R_ohm": R, "q_cell_W": q,
                     "T_max_C": float(T.max() - 273.15), "T_min_C": float(T.min() - 273.15), "T_mean_C": mean - 273.15,
                     "spread_K": spread, "hot_minus_mean_K": float(T.max() - mean),
                     "i_hot": int(np.argmax(T)), "i_cold": int(np.argmin(T)),
                     "balance_W": balance(T, q, h, Ta)}
                for i in range(8):
                    r[f"T_cell{i}_C"] = float(T[i] - 273.15)
                for i in range(8):
                    r[f"offset{i}_K"] = float(T[i] - mean)
                r.update({"mp_spread_K": spread_mp, "mp_newton_iterations": itm,
                          "own_minus_mp_max_abs_K": float(np.abs(T - Tm_f).max()),
                          "own_spread_minus_mp_K": spread - spread_mp, "own_mean_minus_mp_K": mean - mean_mp,
                          "ms_spread_K": float(Tp.max() - Tp.min()), "ms_balance_W": bal_p,
                          "ms_newton_iterations": info_p["iterations"],
                          "ms_converged_D58": bool(abs(bal_p) <= 1e-9 and info_p["iterations"] < 100),
                          "ms_minus_mp_max_abs_K": float(np.abs(Tp - Tm_f).max()),
                          "ms_minus_own_max_abs_K": float(np.abs(Tp - T).max()),
                          "ms_spread_minus_mp_K": float(Tp.max() - Tp.min()) - spread_mp})
                for i in range(8):
                    r[f"ms_T_cell{i}_C"] = float(Tp[i] - 273.15)
                rows.append(r)
                print(f"{Tc_amb:4.0f} C c={c:<4} h={h:<4} q={q:.6f} W  spread={spread:.6f} K mean={mean-273.15:.4f} C "
                      f"Tmax={T.max()-273.15:.3f} C  |own-mp|={r['own_minus_mp_max_abs_K']:.1e} |ms-mp|={r['ms_minus_mp_max_abs_K']:.1e} "
                      f"ms_it={info_p['iterations']} ms_bal={bal_p:.1e}", flush=True)
    keys = list(rows[0].keys())
    with open(os.path.join(OUT, "module.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: repr(v) if isinstance(v, float) else v for k, v in r.items()})

    # V1-style comparison of the 35 C cases with Stage 1's archived results/module_spread.csv
    v1 = []
    with open(os.path.join(ROOT, "results", "module_spread.csv")) as f:
        arch = [x for x in csv.DictReader(f) if x["arrangement"] == "2x4 in air"]
    for x in arch:
        c, h = float(x["crate"]), float(x["h"])
        mine = [r for r in rows if r["T_amb_C"] == 35.0 and r["crate"] == c and r["h"] == h]
        if not mine:
            continue
        r = mine[0]
        cells_a = np.array([float(x[f"T_cell{i}_C"]) for i in range(8)])
        cells_own = np.array([r[f"T_cell{i}_C"] for i in range(8)])
        cells_ms = np.array([r[f"ms_T_cell{i}_C"] for i in range(8)])
        v1.append({"crate": c, "h": h,
                   "q_equal_bitwise": float(x["q_cell_W"]) == r["q_cell_W"],
                   "own_max_abs_K": float(np.abs(cells_own - cells_a).max()),
                   "own_spread_diff_K": r["spread_K"] - float(x["dT_module_K"]),
                   "own_mean_diff_K": r["T_mean_C"] - float(x["T_mean_C"]),
                   "own_six_decimals_equal": bool(all(f"{a:.6f}" == f"{b:.6f}" for a, b in zip(cells_own, cells_a))
                                                  and f"{r['spread_K']:.6f}" == f"{float(x['dT_module_K']):.6f}"),
                   "ms_bitwise_cells": bool(np.array_equal(cells_ms, cells_a)),
                   "ms_bitwise_spread": r["ms_spread_K"] == float(x["dT_module_K"]),
                   "ms_iterations_equal": r["ms_newton_iterations"] == int(x["newton_iterations"]),
                   "hot_cold_index_equal": (r["i_hot"], r["i_cold"]) == (int(x["i_hot"]), int(x["i_cold"]))})
    summary = {"constants_own_vs_thermal_lever": {k: {"own": a, "thermal_lever": b, "diff": a - b} for k, (a, b) in const_check.items()},
               "R_used_ohm": {f"{k:g}": v for k, v in Rmap.items()},
               "V1_35C_against_module_spread_csv": v1,
               "max_abs_own_minus_mp_K": max(r["own_minus_mp_max_abs_K"] for r in rows),
               "max_abs_ms_minus_mp_K": max(r["ms_minus_mp_max_abs_K"] for r in rows),
               "max_abs_ms_minus_own_K": max(r["ms_minus_own_max_abs_K"] for r in rows),
               "all_ms_converged_D58": all(r["ms_converged_D58"] for r in rows),
               "max_abs_ms_balance_W": max(abs(r["ms_balance_W"]) for r in rows),
               "max_ms_iterations": max(r["ms_newton_iterations"] for r in rows),
               "max_abs_own_balance_W": max(abs(r["balance_W"]) for r in rows)}
    with open(os.path.join(OUT, "module.json"), "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps({k: v for k, v in summary.items() if k != "V1_35C_against_module_spread_csv"}, indent=1))
    for v in v1:
        print(v)


if __name__ == "__main__":
    main()
