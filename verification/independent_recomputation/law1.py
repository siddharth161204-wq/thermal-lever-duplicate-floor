"""Law 1 (Kucinskis et al. 2022, D3): lever from the printed cycles to 90 % SoH, Fig. 7b.

Definitions used (DECISIONS.md D3, D16, D17):
  ln L at a C/2 charge = mean of ln L at 0.4C and 0.6C at each printed temperature,
  ln L interpolated linearly in T between the printed temperatures,
  s(35 C) = -(ln L(45) - ln L(25)) / 20 K   (equal to the 0.5 K centred difference of the piecewise-linear curve),
  Delta T_floor = f / s,  exact pair form 2 atanh(f/2) / s.
"""
import csv
import math
import os

import numpy as np

from common import FLOORS, FLOOR_KEYS, REPO, save

CSV = os.path.join(REPO, "inputs", "kucinskis2022_fig7_cycles.csv")


def load_table(soh=0.9):
    tab = {}
    with open(CSV) as f:
        for row in csv.DictReader(f):
            if abs(float(row["soh"]) - soh) > 1e-9:
                continue
            if row["cycles_to_soh"].strip() == "":
                continue
            tab[(float(row["charge_c_rate"]), float(row["ambient_T_C"]))] = float(row["cycles_to_soh"])
    return tab


TAB = load_table(0.9)
T_GRID = [15.0, 25.0, 35.0, 45.0, 55.0]


def lnL_rate(rate, T):
    """ln L at a printed charge rate (0.2, 0.4, 0.6) or at C/2 (mean of 0.4 and 0.6 in ln L), linear in T."""
    if abs(rate - 0.5) < 1e-12:
        return 0.5 * (lnL_rate(0.4, T) + lnL_rate(0.6, T))
    Ts = [t for t in T_GRID if (rate, t) in TAB]
    ys = [math.log(TAB[(rate, t)]) for t in Ts]
    if T < Ts[0] - 1e-12 or T > Ts[-1] + 1e-12:
        return float("nan")
    return float(np.interp(T, Ts, ys))


def s_centred_secant(rate):
    return -(lnL_rate(rate, 45.0) - lnL_rate(rate, 25.0)) / 20.0


def s_centred_difference(rate, T=35.0, d=0.25):
    return -(lnL_rate(rate, T + d) - lnL_rate(rate, T - d)) / (2 * d)


def floors_for(s):
    return ({k: FLOORS[k] / s for k in FLOOR_KEYS},
            {k: 2 * math.atanh(FLOORS[k] / 2) / s for k in FLOOR_KEYS})


def main():
    out = {}
    s = s_centred_secant(0.5)
    s_cd_quarter = s_centred_difference(0.5, d=0.25)
    s_cd_half = s_centred_difference(0.5, d=0.5)
    lin, exact = floors_for(s)
    out["C2"] = {
        "L_25": math.exp(lnL_rate(0.5, 25)), "L_35": math.exp(lnL_rate(0.5, 35)), "L_45": math.exp(lnL_rate(0.5, 45)),
        "s_centred_secant": s, "s_centred_diff_pm0.25K": s_cd_quarter, "s_centred_diff_pm0.5K": s_cd_half,
        "dT_floor_K": lin, "dT_floor_exact_pair_K": exact,
    }
    # printed charge rates
    out["by_charge_rate"] = {}
    for r in (0.2, 0.4, 0.6, 0.5):
        sr = s_centred_secant(r)
        lin_r, exact_r = floors_for(sr)
        out["by_charge_rate"][str(r)] = {
            "L_25_35_45": [math.exp(lnL_rate(r, T)) for T in (25, 35, 45)],
            "s_centred": sr,
            "s_25_to_35": -(lnL_rate(r, 35) - lnL_rate(r, 25)) / 10,
            "s_35_to_45": -(lnL_rate(r, 45) - lnL_rate(r, 35)) / 10,
            "dT_floor_K": lin_r, "dT_floor_exact_pair_K": exact_r,
        }
    # lifetime curve at C/2, 15 to 55 C
    out["curve_C2"] = [{"T_C": float(T), "L": math.exp(lnL_rate(0.5, float(T)))} for T in range(15, 56)]

    # D16: assumed noise per printed value, sigma(ln L) = 0.0265 / sqrt(2)
    sig = 0.0265 / math.sqrt(2)
    # analytic propagation: s = -[(a45 + b45)/2 - (a25 + b25)/2] / 20, four independent printed values
    se_analytic = math.sqrt(4 * (0.5 / 20) ** 2) * sig            # = sig / 20
    se_single_rate = math.sqrt(2) * sig / 20                       # one printed rate (the D16 hand form)
    res = {}
    for n_draws, seed in ((20000, 20261010), (2_000_000, 12345)):
        rng = np.random.default_rng(seed)
        a25 = math.log(TAB[(0.4, 25.0)]) + sig * rng.standard_normal(n_draws)
        b25 = math.log(TAB[(0.6, 25.0)]) + sig * rng.standard_normal(n_draws)
        a45 = math.log(TAB[(0.4, 45.0)]) + sig * rng.standard_normal(n_draws)
        b45 = math.log(TAB[(0.6, 45.0)]) + sig * rng.standard_normal(n_draws)
        s_mc = -((a45 + b45) / 2 - (a25 + b25) / 2) / 20
        dT = np.abs(FLOORS["p95_8.7pct"] / s_mc)
        res[str(n_draws)] = {
            "seed": seed,
            "s_2.5_50_97.5": np.percentile(s_mc, [2.5, 50, 97.5]).tolist(),
            "s_std": float(s_mc.std(ddof=1)),
            "fraction_s_positive": float((s_mc > 0).mean()),
            "abs_dT_floor_p95_2.5_50_97.5_K": np.percentile(dT, [2.5, 50, 97.5]).tolist(),
        }
    out["assumed_noise_mc"] = {
        "sigma_lnL_per_printed_value": sig,
        "se_s_analytic_four_values": se_analytic,
        "se_s_single_printed_rate_two_values": se_single_rate,
        "D16_hand_figure": 0.0013,
        "draws": res,
        "normal_approx_s_2.5_97.5": [s - 1.959964 * se_analytic, s + 1.959964 * se_analytic],
    }
    save("law1.json", out)
    print("Law 1, C/2 charge: L(25,35,45) =", [round(x, 4) for x in (out["C2"]["L_25"], out["C2"]["L_35"], out["C2"]["L_45"])])
    print("  s =", s, " (centred diff +-0.25 K:", s_cd_quarter, ", +-0.5 K:", s_cd_half, ")")
    print("  dT_floor =", lin)
    print("  exact pair =", exact)
    for r in ("0.2", "0.4", "0.6"):
        b = out["by_charge_rate"][r]
        print("  rate", r, "s =", b["s_centred"], "dT_floor(8.7) =", b["dT_floor_K"]["p95_8.7pct"])
    print("  MC:", res, "\n  SE analytic (C/2, four values) =", se_analytic, " single printed rate =", se_single_rate)


if __name__ == "__main__":
    main()
