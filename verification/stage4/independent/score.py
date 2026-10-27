#!/usr/bin/env python3
"""Independent check of Stage 4, part (f): P1 to P7 scored by D63 from the checker's own outputs.

Usage: python3 verification/stage4/independent/score.py
Output: out/scorecard.json
"""
import csv, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
sys.path.insert(0, HERE)
import levers as LV                                 # noqa: E402


def rd(name):
    with open(os.path.join(OUT, name)) as f:
        return list(csv.DictReader(f))


def main():
    warm, cold, floors = LV.load()
    G, K = floors["Geslin"], floors["KIT"]
    cases = rd("divergence_cases.csv")
    divs = rd("divergence.csv")
    res = rd("resolved_cell.csv")
    exp = rd("experiment.csv")
    sim_path = os.path.join(OUT, "simulation_lognormal_all.csv")
    if not os.path.exists(sim_path):
        sim_path = os.path.join(OUT, "simulation_lognormal_required.csv")
    sim = rd(os.path.basename(sim_path))

    def case(Ta, c, h):
        x = [r for r in cases if float(r["T_amb_C"]) == Ta and float(r["crate"]) == c and float(r["h"]) == h]
        assert len(x) == 1
        return x[0]

    sc = {}
    # P1
    p1 = {f"{c:g}C h={h:g}": float(case(35.0, c, h)["median"]) for c in (0.25, 0.5, 1.0) for h in (5.0, 10.0, 25.0)}
    sc["P1"] = {"held": all(v <= G["mean"] for v in p1.values()), "median_warm": p1, "max": max(p1.values()),
                "floor": G["mean"]}
    # P2
    c2 = case(35.0, 2.0, 5.0)
    lev = [r for r in divs if float(r["T_amb_C"]) == 35.0 and float(r["crate"]) == 2.0 and float(r["h"]) == 5.0]
    above = [(r["lever"], float(r["implied_divergence"])) for r in lev if float(r["implied_divergence"]) > G["p95"]]
    beside = {}
    for h in (10.0, 25.0):
        lv_h = [r for r in divs if float(r["T_amb_C"]) == 35.0 and float(r["crate"]) == 2.0 and float(r["h"]) == h]
        beside[f"h={h:g}"] = {"median": float(case(35.0, 2.0, h)["median"]),
                              "levers_above_p95": sum(float(r["implied_divergence"]) > G["p95"] for r in lv_h)}
    sc["P2"] = {"held": float(c2["median"]) > G["mean"] and not above, "median_h5": float(c2["median"]),
                "first_clause": float(c2["median"]) > G["mean"], "second_clause": not above,
                "levers_above_Geslin_p95": above, "max_lever_divergence": max(float(r["implied_divergence"]) for r in lev),
                "other_coolings": beside}
    # P3
    p3 = {}
    for h in (5.0, 10.0, 25.0):
        a, b = float(case(10.0, 1.0, h)["spread_K"]), float(case(35.0, 1.0, h)["spread_K"])
        p3[f"h={h:g}"] = {"spread_10C": a, "spread_35C": b, "ratio": a / b}
    sc["P3"] = {"held": all(v["ratio"] >= 1.30 for v in p3.values()), "cases": p3}
    # P4
    p4 = {f"h={h:g}": float(case(10.0, 1.0, h)["median"]) for h in (5.0, 10.0, 25.0)}
    sc["P4"] = {"held": all(v > G["mean"] for v in p4.values()) and p4["h=5"] > G["p95"],
                "first_clause": all(v > G["mean"] for v in p4.values()), "second_clause": p4["h=5"] > G["p95"],
                "median_cold": p4}
    # P5
    warm35 = {f"{float(r['crate']):g}C h={float(r['h']):g}": float(r["median"]) for r in cases if float(r["T_amb_C"]) == 35.0}
    cold10 = {f"{float(r['crate']):g}C h={float(r['h']):g}": float(r["median"]) for r in cases if float(r["T_amb_C"]) == 10.0}
    cl1 = all(v <= K["mean"] for v in warm35.values())
    cl2a = all(cold10[f"2C h={h:g}"] > K["mean"] for h in (5.0, 10.0, 25.0))
    below = {k: v for k, v in cold10.items() if not k.startswith("2C")}
    cl2b = all(v <= K["mean"] for v in below.values())
    sc["P5"] = {"held": cl1 and cl2a and cl2b, "clause_35C_none_above_KIT_mean": cl1, "max_35C": max(warm35.values()),
                "clause_10C_every_2C_above": cl2a, "clause_10C_none_below_2C_above": cl2b,
                "cases_below_2C_above_KIT_mean": {k: v for k, v in below.items() if v > K["mean"]},
                "floor": K["mean"]}
    # P6
    r6 = [r for r in res if float(r["T_amb_C"]) == 10.0 and float(r["crate"]) == 1.0][0]
    sc["P6"] = {"held": float(r6["dT_cell_K"]) < float(r6["module_spread_K"]) / 3.0, "dT_cell_K": float(r6["dT_cell_K"]),
                "module_spread_K": float(r6["module_spread_K"]), "ratio": float(r6["ratio_cell_over_module"])}
    # P7
    cw = sorted(warm, key=lambda lv: abs(lv["s"]))[len(warm) // 2]["lever"]
    cc = sorted(cold, key=lambda lv: abs(lv["s"]))[len(cold) // 2]["lever"]

    def pick(group, lever, dT):
        x = [r for r in exp if r["group"] == group and r["lever"] == lever and float(r["dT_K"]) == dT
             and r["floor"] == "Geslin" and r["shape"] == "lognormal"]
        assert len(x) == 1
        y = [r for r in sim if r["group"] == group and r["lever"] == lever and float(r["dT_K"]) == dT
             and r["floor"] == "Geslin"]
        return int(x[0]["n"]), (int(y[0]["n_sim"]) if y else None), (y[0]["V4_pass_5pct"] == "True" if y else None)

    nc, nc_sim, v4c = pick("C2 cold chamber 10 C", cc, 5.0)
    nw, nw_sim, v4w = pick("C2 warm chamber 35 C", cw, 2.0)
    held = nc < 10 and nw > 30
    other_c = (not v4c) and ((nc_sim < 10) != (nc < 10))
    other_w = (not v4w) and ((nw_sim > 30) != (nw > 30))
    sc["P7"] = {"held": None if (other_c or other_w) else held, "scored": not (other_c or other_w),
                "cold_central_lever": cc, "n_cold_dT5": nc, "n_cold_sim": nc_sim, "V4_cold_pass": v4c,
                "warm_central_lever": cw, "n_warm_dT2": nw, "n_warm_sim": nw_sim, "V4_warm_pass": v4w,
                "cold_clause": nc < 10, "warm_clause": nw > 30, "simulation_file": os.path.basename(sim_path)}
    with open(os.path.join(OUT, "scorecard.json"), "w") as f:
        json.dump(sc, f, indent=1)
    for k, v in sc.items():
        print(k, "held" if v["held"] else ("NOT SCORED" if v["held"] is None else "did not hold"))
    print(json.dumps(sc, indent=1))


if __name__ == "__main__":
    main()
