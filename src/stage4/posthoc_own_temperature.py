#!/usr/bin/env python3
"""Stage 4, after the results (DECISIONS.md D69): the module with each cell's heat at its own temperature, a
conditional reading that changes no registered number or verdict.

The registered heat (D58) takes the resistance at the ambient. Here R(T) is the same C/2 quantity, computed with the
separate code of V2 (verification/stage4/v2_resistance.py: tight solver tolerances, 1 s grid, no solver spikes) at 10
to 90 C; ln R is interpolated linearly in T and held at the end values outside that range. Two fixed points with
Stage 1's network (thermal_lever.module_steady): (own) each cell's heat at its own temperature; (mean) one heat for
all cells at the module's mean temperature. This repeats, with this study's own code, a sensitivity the independent
checker computed (verification/stage4/independent/sensitivity_own_temperature.py).

Output: results/stage4_posthoc_own_temperature.json
Usage: python3 src/stage4/posthoc_own_temperature.py
"""
import csv, json, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as CM                     # noqa: E402
import thermal_lever as TL              # noqa: E402
sys.path.insert(0, os.path.join(CM.ROOT, "verification", "stage4"))
import v2_resistance as V2              # noqa: E402

TEMPS = [10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 70.0, 80.0, 90.0]


def r_curve():
    rows = []
    for T in TEMPS:
        dis, ch = V2.one_direction(T, +1), V2.one_direction(T, -1)
        rows.append({"T_C": T, "R_ohm": 0.5 * (dis["R_ohm"] + ch["R_ohm"]), "charge_Ah_discharge": dis["current_A"] * dis["duration_s"] / 3600,
                     "charge_Ah_charge": -ch["current_A"] * ch["duration_s"] / 3600, "terminations": [dis["termination"], ch["termination"]]})
        print(f"{T:g} C: {rows[-1]['R_ohm']:.6g} ohm", flush=True)
    return rows


def fixed_point(T_amb_C, c, h, R_of, variant, T0):
    T = np.array(T0, dtype=float) + 273.15
    for it in range(500):
        if variant == "own":
            q = np.array([(c * TL.CAP_AH) ** 2 * R_of(t - 273.15) for t in T])
        else:
            q = np.full(8, (c * TL.CAP_AH) ** 2 * R_of(T.mean() - 273.15))
        Tn, bal, info = TL.module_steady(q, h, T_amb=T_amb_C + 273.15)
        d = float(np.abs(Tn - T).max())
        T = 0.5 * (T + Tn) if it > 50 else Tn
        if d < 1e-10:
            break
    return T, it + 1, d


def main():
    d = CM.load_all()
    curve = r_curve()
    Ts = np.array([r["T_C"] for r in curve])
    lnR = np.log([r["R_ohm"] for r in curve])
    R_of = lambda T_C: float(np.exp(np.interp(T_C, Ts, lnR)))
    s_med = {"warm": abs(d["named"]["median_warm"]["s"]), "cold": abs(d["named"]["median_cold"]["s"])}
    s_max = {b: max(abs(l["s"]) for l in d["bands"][b]) for b in ("warm", "cold")}
    reg = list(csv.DictReader(open(os.path.join(CM.ROOT, "results", "stage4_module.csv"))))
    cases = []
    for r in reg:
        T_amb, c, h = float(r["ambient_C"]), float(r["crate"]), float(r["h"])
        band = CM.BAND_OF_AMBIENT[T_amb]
        T0 = [float(r[f"T_cell{i}_C"]) for i in range(8)]
        x = {"ambient_C": T_amb, "rate": r["rate"], "h": h, "registered_spread_K": float(r["spread_K"]),
             "registered_median": s_med[band] * float(r["spread_K"])}
        for variant in ("own", "mean"):
            T, it, dlast = fixed_point(T_amb, c, h, R_of, variant, T0)
            sp = float(T.max() - T.min())
            x.update({f"{variant}_spread_K": sp, f"{variant}_T_min_C": float(T.min() - 273.15), f"{variant}_T_max_C": float(T.max() - 273.15),
                      f"{variant}_median": s_med[band] * sp, f"{variant}_max_lever": s_max[band] * sp, f"{variant}_iterations": it})
        cases.append(x)
    F = d["floors"]
    get = lambda T, rate, h: [x for x in cases if x["ambient_C"] == T and x["rate"] == rate and x["h"] == h][0]
    verdicts = {}
    for v in ("own", "mean"):
        p1 = not any(x[f"{v}_median"] > F["Geslin"]["mean"] for x in cases if x["ambient_C"] == 35.0 and x["rate"] != "2C")
        c2 = get(35.0, "2C", 5.0)
        p2 = c2[f"{v}_median"] > F["Geslin"]["mean"] and c2[f"{v}_max_lever"] <= F["Geslin"]["p95"]
        p3r = {f"h {h:g}": get(10.0, "1C", h)[f"{v}_spread_K"] / get(35.0, "1C", h)[f"{v}_spread_K"] for h in CM.H_LIST}
        p4 = all(get(10.0, "1C", h)[f"{v}_median"] > F["Geslin"]["mean"] for h in CM.H_LIST) and get(10.0, "1C", 5.0)[f"{v}_median"] > F["Geslin"]["p95"]
        p5 = (not any(x[f"{v}_median"] > F["KIT"]["mean"] for x in cases if x["ambient_C"] == 35.0)
              and all(get(10.0, "2C", h)[f"{v}_median"] > F["KIT"]["mean"] for h in CM.H_LIST)
              and not any(x[f"{v}_median"] > F["KIT"]["mean"] for x in cases if x["ambient_C"] == 10.0 and x["rate"] != "2C"))
        verdicts[v] = {"P1": p1, "P2": p2, "P3": all(r >= 1.30 for r in p3r.values()), "P3_ratios": p3r, "P4": p4, "P5": p5}
    out = {"note": "conditional, after the results: the heat at each cell's own temperature (own) or at the module mean (mean)",
           "R_of_T_ohm": {f"{r['T_C']:g}": r["R_ohm"] for r in curve}, "R_curve": curve, "cases": cases, "verdicts": verdicts}
    json.dump(out, open(os.path.join(CM.ROOT, "results", "stage4_posthoc_own_temperature.json"), "w"), indent=1)
    for x in cases:
        if x["rate"] in ("1C", "2C"):
            print(f"{x['ambient_C']:g} C {x['rate']} h {x['h']:g}: registered {x['registered_spread_K']:.3f} K, own {x['own_spread_K']:.3f} K "
                  f"({x['own_T_min_C']:.1f} to {x['own_T_max_C']:.1f} C, median {100 * x['own_median']:.2f} %, max lever {100 * x['own_max_lever']:.2f} %), "
                  f"mean {x['mean_spread_K']:.3f} K ({100 * x['mean_median']:.2f} %)")
    print(json.dumps(verdicts, indent=1))


if __name__ == "__main__":
    main()
