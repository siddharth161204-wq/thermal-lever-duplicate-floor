#!/usr/bin/env python3
"""Independent check of Stage 4, part (e): the experiment of D61, calculated (no simulation here).

n split pairs; d = ln L_warm - ln L_cool = -s dT + e; |e| log-normal with the floor's mean and 95th percentile (the
registered shape, smaller root), sd(e) = sqrt(E[e^2]) = exp(mu + sigma^2); the normal shape sd(e) = m sqrt(pi/2) for
the sensitivity. Paired t-test against zero, two-sided 5 %, 80 % power, noncentral t (power.py, own quadrature);
n is the smallest integer from 2 with power >= 0.80.

C1: s = the central KIT value at 1C (10-100 1.00/1.00C, D57), the Sandia 1C value (NMC 25 to 35 C), the central cold
KIT value (the mean of the sixth and seventh |s| of the twelve KIT 10 to 25 C levers, negative), and, conditional, the
nine-protocol reading of the central KIT value at 1C (10-90 1.67/1.00C); dT = 1, 2, 5, 10 K.
C2: the warm chamber (35 C, every lever of the warm band) and the cold chamber (10 C, every lever of the cold band),
dT = 2 and 5 K, the central (median |s|) lever named. Both floors throughout. Figure grid: dT 0.5 to 10 K in 0.5 K
steps under the two chambers' central levers, both floors.

Usage: python3 verification/stage4/independent/experiment.py
Output: out/experiment.csv, out/experiment_grid.csv, out/experiment.json
"""
import csv, json, math, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
sys.path.insert(0, HERE)
import levers as LV                                 # noqa: E402
import power as PW                                  # noqa: E402


def central(band):
    a = sorted(band, key=lambda lv: abs(lv["s"]))
    assert len(a) % 2 == 1
    return a[len(a) // 2]


def main():
    warm, cold, floors = LV.load()
    shapes = {}
    for fk, fv in floors.items():
        ln = PW.lognormal_shape(fv["mean"], fv["p95"])
        shapes[fk] = {"lognormal": ln["small"], "lognormal_larger_root_not_used": ln["large"],
                      "measured_median": fv["median"], "normal_sd_e": PW.normal_sd(fv["mean"]),
                      "m": fv["mean"], "q": fv["p95"]}
    kit_warm = [lv for lv in warm if lv["source"] == "KIT"]
    kit_cold = [lv for lv in cold if lv["source"] == "KIT"]
    one_one = [lv for lv in kit_warm if lv["protocol"].endswith(" 1.00/1.00C")]
    assert len(one_one) == 3
    c1_kit = sorted(one_one, key=lambda lv: abs(lv["s"]))[1]
    nine = [lv for lv in kit_warm if lv["protocol"].endswith("/1.00C")]
    assert len(nine) == 9
    c1_kit_nine = sorted(nine, key=lambda lv: abs(lv["s"]))[4]
    c1_sandia = [lv for lv in warm if lv["source"] == "Sandia"][0]
    kc = sorted(kit_cold, key=lambda lv: abs(lv["s"]))
    c1_cold_abs = 0.5 * (abs(kc[5]["s"]) + abs(kc[6]["s"]))
    assert kc[5]["s"] < 0 and kc[6]["s"] < 0
    c1_cold = {"lever": f"KIT 10 to 25 C median of |s| ({kc[5]['protocol']} and {kc[6]['protocol']})", "s": -c1_cold_abs}
    cw, cc = central(warm), central(cold)
    rows = []

    def add(group, label, lv, dT, central_flag=False, conditional=False):
        for fk in floors:
            for shape in ("lognormal", "normal"):
                sd = shapes[fk]["lognormal"]["sd_e"] if shape == "lognormal" else shapes[fk]["normal_sd_e"]
                eff = abs(lv["s"]) * dT / sd
                n = PW.n_required(eff)
                p_n, up, low = PW.power_quad(n, eff)
                p_nm1 = PW.power_quad(n - 1, eff)[0] if n > 2 else None
                n_upper_only = PW.n_required(eff, lower_tail=False)
                try:
                    p_nct = PW.power_nct(n, eff)
                except Exception:
                    p_nct = float("nan")
                rows.append({"group": group, "label": label, "lever": lv["lever"], "s_per_K": lv["s"], "dT_K": dT,
                             "floor": fk, "shape": shape, "registered_shape": shape == "lognormal", "sd_e": sd,
                             "effect_size": eff, "noncentrality_at_n": math.sqrt(n) * eff, "n": n,
                             "power_at_n": p_n, "power_at_n_minus_1": p_nm1, "lower_tail_at_n": low,
                             "n_upper_tail_only": n_upper_only, "power_nct_at_n": p_nct,
                             "nct_minus_quad": (p_nct - p_n) if p_nct == p_nct else float("nan"),
                             "nct_is_nan": not (p_nct == p_nct),
                             "central": central_flag, "conditional": conditional})

    for dT in (1.0, 2.0, 5.0, 10.0):
        add("C1", "central KIT value at 1C (1.00/1.00C protocols)", c1_kit, dT)
        add("C1", "Sandia 1C value", c1_sandia, dT)
        add("C1", "central cold KIT value", c1_cold, dT)
        add("C1", "central KIT value at 1C, nine-protocol reading (conditional)", c1_kit_nine, dT, conditional=True)
    for dT in (2.0, 5.0):
        for lv in warm:
            add("C2 warm chamber 35 C", "warm band lever", lv, dT, central_flag=lv is cw)
        for lv in cold:
            add("C2 cold chamber 10 C", "cold band lever", lv, dT, central_flag=lv is cc)
    with open(os.path.join(OUT, "experiment.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: repr(v) if isinstance(v, float) else v for k, v in r.items()})
    grid = []
    for dT in np.arange(0.5, 10.0001, 0.5):
        for nm, lv in (("warm chamber central", cw), ("cold chamber central", cc)):
            for fk in floors:
                for shape in ("lognormal", "normal"):
                    sd = shapes[fk]["lognormal"]["sd_e"] if shape == "lognormal" else shapes[fk]["normal_sd_e"]
                    eff = abs(lv["s"]) * float(dT) / sd
                    grid.append({"chamber": nm, "lever": lv["lever"], "s_per_K": lv["s"], "dT_K": float(dT), "floor": fk,
                                 "shape": shape, "n": PW.n_required(eff)})
    with open(os.path.join(OUT, "experiment_grid.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(grid[0].keys()))
        w.writeheader()
        for r in grid:
            w.writerow(r)
    summ = {"shapes": shapes, "c1_levers": {"central KIT 1C": c1_kit["lever"], "central KIT 1C s": c1_kit["s"],
                                            "nine-protocol reading": c1_kit_nine["lever"], "nine s": c1_kit_nine["s"],
                                            "Sandia 1C": c1_sandia["lever"], "Sandia s": c1_sandia["s"],
                                            "central cold KIT s": c1_cold["s"]},
            "central_warm_lever": cw["lever"], "central_warm_s": cw["s"],
            "central_cold_lever": cc["lever"], "central_cold_s": cc["s"],
            "max_abs_nct_minus_quad": float(np.nanmax([abs(r["nct_minus_quad"]) for r in rows])),
            "nct_nan_count": int(sum(r["nct_is_nan"] for r in rows)),
            "cases_where_upper_tail_only_changes_n": [(r["lever"], r["dT_K"], r["floor"], r["shape"], r["n"], r["n_upper_tail_only"])
                                                      for r in rows if r["n"] != r["n_upper_tail_only"]]}
    with open(os.path.join(OUT, "experiment.json"), "w") as f:
        json.dump(summ, f, indent=1)
    print(json.dumps(summ, indent=1))
    for r in rows:
        if r["group"] == "C1" or r["central"]:
            print(f"{r['group']:22s} {r['lever']:44s} s={100*r['s_per_K']:+.4f} dT={r['dT_K']:4g} {r['floor']:6s} {r['shape']:9s} "
                  f"sd={r['sd_e']:.6f} n={r['n']:6d} pow={r['power_at_n']:.4f} prev={r['power_at_n_minus_1']}")


if __name__ == "__main__":
    main()
