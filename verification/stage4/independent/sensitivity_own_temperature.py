#!/usr/bin/env python3
"""The checker's own sensitivity, NOT a registered quantity: the module with the heat at each cell's own temperature.

D58 takes the resistance at the ambient, not at each cell's temperature, and calls the cold heat a bracket. This
script shows which way and how far: R(T) is the same C/2 quantity (DFN, ORegan2022, isothermal, mean of the C/2
discharge and charge time-averaged irreversible heat over I^2) computed by resistance_check.py --sweep from 10 to 90 C
with the sub-second heat spike of the cold charges screened out; ln R is interpolated linearly in T and held at the
90 C value above 90 C. Two variants: (i) each cell's heat at its own temperature; (ii) one heat for all cells at the
module mean temperature. Each is a fixed point, iterated with the module network (thermal_lever.module_steady, which
this checker's own solver reproduces to 1e-13 K) until no temperature moves by more than 1e-10 K.

Output: out/sensitivity_own_temperature.csv, out/sensitivity_own_temperature.json
"""
import csv, json, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = os.path.join(HERE, "out")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)
import thermal_lever as TL                        # noqa: E402
import levers as LV                               # noqa: E402


def main():
    with open(os.path.join(OUT, "resistance_sweep.json")) as f:
        sw = json.load(f)["rows"]
    Ts = np.array([r["T_C"] for r in sw])
    lnR = np.log(np.array([r["R_ohm_spike_screened"] for r in sw]))
    R_of = lambda T_C: float(np.exp(np.interp(T_C, Ts, lnR)))       # held at the end values outside 10 to 90 C
    warm, cold, floors = LV.load()
    med_abs = {"warm": float(np.median([abs(lv["s"]) for lv in warm])), "cold": float(np.median([abs(lv["s"]) for lv in cold]))}
    with open(os.path.join(OUT, "module.csv")) as f:
        reg = {(float(r["T_amb_C"]), float(r["crate"]), float(r["h"])): r for r in csv.DictReader(f)}
    rows = []
    for (Ta_C, c, h), r in sorted(reg.items()):
        Ta = Ta_C + 273.15
        out = {"T_amb_C": Ta_C, "crate": c, "h": h, "registered_spread_K": float(r["spread_K"]),
               "registered_T_mean_C": float(r["T_mean_C"])}
        for variant in ("own", "mean"):
            T = np.array([float(r[f"T_cell{i}_C"]) for i in range(8)]) + 273.15
            for it in range(500):
                if variant == "own":
                    q = np.array([(c * 5.0) ** 2 * R_of(t - 273.15) for t in T])
                else:
                    q = np.full(8, (c * 5.0) ** 2 * R_of(T.mean() - 273.15))
                Tn, bal, info = TL.module_steady(q, h, T_amb=Ta)
                d = np.abs(Tn - T).max()
                T = 0.5 * (T + Tn) if it > 50 else Tn                # damped late, if needed
                if d < 1e-10:
                    break
            band = "warm" if Ta_C in (35.0, 25.0) else "cold"
            spread = float(T.max() - T.min())
            out.update({f"{variant}_spread_K": spread, f"{variant}_T_mean_C": float(T.mean() - 273.15),
                        f"{variant}_T_max_C": float(T.max() - 273.15), f"{variant}_iterations": it + 1,
                        f"{variant}_q_mean_W": float(q.mean()), f"{variant}_median_divergence": med_abs[band] * spread,
                        f"{variant}_beyond_90C": bool(T.max() - 273.15 > 90.0)})
        rows.append(out)
    with open(os.path.join(OUT, "sensitivity_own_temperature.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for x in rows:
            w.writerow({k: repr(v) if isinstance(v, float) else v for k, v in x.items()})
    G, K = floors["Geslin"], floors["KIT"]
    get = lambda Ta, c, h: [x for x in rows if x["T_amb_C"] == Ta and x["crate"] == c and x["h"] == h][0]
    verdicts = {}
    for v in ("own", "mean"):
        p3 = {h: get(10.0, 1.0, h)[f"{v}_spread_K"] / get(35.0, 1.0, h)[f"{v}_spread_K"] for h in (5.0, 10.0, 25.0)}
        p4a = all(get(10.0, 1.0, h)[f"{v}_median_divergence"] > G["mean"] for h in (5.0, 10.0, 25.0))
        p4b = get(10.0, 1.0, 5.0)[f"{v}_median_divergence"] > G["p95"]
        p5a = all(x[f"{v}_median_divergence"] <= K["mean"] for x in rows if x["T_amb_C"] == 35.0)
        p5b = all(get(10.0, 2.0, h)[f"{v}_median_divergence"] > K["mean"] for h in (5.0, 10.0, 25.0))
        p5c = all(x[f"{v}_median_divergence"] <= K["mean"] for x in rows if x["T_amb_C"] == 10.0 and x["crate"] < 2.0)
        p1 = all(x[f"{v}_median_divergence"] <= G["mean"] for x in rows if x["T_amb_C"] == 35.0 and x["crate"] < 2.0)
        p2a = get(35.0, 2.0, 5.0)[f"{v}_median_divergence"] > G["mean"]
        p2b = max(abs(lv["s"]) for lv in warm) * get(35.0, 2.0, 5.0)[f"{v}_spread_K"] <= G["p95"]
        verdicts[v] = {"P1": p1, "P2": p2a and p2b, "P3_ratios": p3, "P3": all(x >= 1.30 for x in p3.values()),
                       "P4": p4a and p4b, "P4_first": p4a, "P4_second": p4b, "P5": p5a and p5b and p5c,
                       "P5_clauses": [p5a, p5b, p5c]}
    with open(os.path.join(OUT, "sensitivity_own_temperature.json"), "w") as f:
        json.dump({"note": "checker's sensitivity, not registered", "R_of_T": dict(zip([float(x) for x in Ts], [float(np.exp(x)) for x in lnR])),
                   "verdicts_under_this_heat": verdicts}, f, indent=1)
    for x in rows:
        print(f"{x['T_amb_C']:4.0f} C c={x['crate']:<4} h={x['h']:<4} reg {x['registered_spread_K']:8.4f} K ({x['registered_T_mean_C']:6.1f} C) | "
              f"own {x['own_spread_K']:8.4f} K ({x['own_T_mean_C']:6.1f} C, med {100*x['own_median_divergence']:6.2f} %) | "
              f"mean {x['mean_spread_K']:8.4f} K ({x['mean_T_mean_C']:6.1f} C)")
    print(json.dumps(verdicts, indent=1))


if __name__ == "__main__":
    main()
