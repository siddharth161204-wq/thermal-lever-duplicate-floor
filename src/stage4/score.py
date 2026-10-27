#!/usr/bin/env python3
"""Stage 4: scoring of predictions_stage4.md by the rules of DECISIONS.md D63.

Reads results/stage4_module.csv, results/stage4_map.json, results/stage4_resolved_cell.csv and
results/stage4_experiment.csv; writes results/stage4_scorecard.json.

Usage: python3 src/stage4/score.py
"""
import argparse, csv, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as CM                     # noqa: E402


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return x


def load(results):
    mod = [{k: _f(v) for k, v in r.items()} for r in csv.DictReader(open(os.path.join(results, "stage4_module.csv")))]
    mp = json.load(open(os.path.join(results, "stage4_map.json")))
    rc = [{k: _f(v) for k, v in r.items()} for r in csv.DictReader(open(os.path.join(results, "stage4_resolved_cell.csv")))]
    ex = [{k: _f(v) for k, v in r.items()} for r in csv.DictReader(open(os.path.join(results, "stage4_experiment.csv")))]
    return mod, mp, rc, ex


def score(mod, mp, rc, ex, floors):
    G, K = floors["Geslin"], floors["KIT"]
    case = lambda T, c, h: [x for x in mp["cases"] if x["ambient_C"] == T and x["crate"] == c and x["h"] == h][0]
    spread = lambda T, c, h: [x for x in mod if x["ambient_C"] == T and x["crate"] == c and x["h"] == h][0]["spread_K"]
    H = CM.H_LIST
    sc = {}

    p1 = [case(35.0, c, h) for c in (0.25, 0.5, 1.0) for h in H]
    sc["P1"] = {"held": not any(x["median"] > G["mean"] for x in p1), "floor": G["mean"],
                "medians": {f"{x['rate']} h {x['h']:g}": x["median"] for x in p1},
                "largest": max(x["median"] for x in p1)}

    x = case(35.0, 2.0, 5.0)
    a, b = x["median"] > G["mean"], x["max"] <= G["p95"]
    sc["P2"] = {"held": a and b, "median_above_mean_floor": a, "no_lever_above_p95": b, "median": x["median"],
                "max_lever": x["max"], "levers_above_p95": x["levers_above_p95_Geslin"],
                "beside_other_coolings": {f"h {h:g}": {"median": case(35.0, 2.0, h)["median"], "max_lever": case(35.0, 2.0, h)["max"],
                                                       "levers_above_p95": case(35.0, 2.0, h)["levers_above_p95_Geslin"]} for h in (10.0, 25.0)}}

    ratios = {f"h {h:g}": spread(10.0, 1.0, h) / spread(35.0, 1.0, h) for h in H}
    sc["P3"] = {"held": all(r >= 1.30 for r in ratios.values()), "ratios": ratios,
                "spreads_10C": {f"h {h:g}": spread(10.0, 1.0, h) for h in H}, "spreads_35C": {f"h {h:g}": spread(35.0, 1.0, h) for h in H}}

    p4 = {f"h {h:g}": case(10.0, 1.0, h)["median"] for h in H}
    a = all(v > G["mean"] for v in p4.values())
    b = case(10.0, 1.0, 5.0)["median"] > G["p95"]
    sc["P4"] = {"held": a and b, "share_at_every_cooling": a, "detectable_at_h5": b, "medians": p4}

    warm35 = [case(35.0, c, h) for c in CM.RATES for h in H]
    a = not any(x["median"] > K["mean"] for x in warm35)
    cold = {f"{CM.RATE_LABEL[c]} h {h:g}": case(10.0, c, h)["median"] > K["mean"] for c in CM.RATES for h in H}
    b = all(cold[f"2C h {h:g}"] for h in H) and not any(v for k, v in cold.items() if not k.startswith("2C"))
    sc["P5"] = {"held": a and b, "none_at_35C": a, "only_2C_at_10C": b, "largest_35C_median": max(x["median"] for x in warm35),
                "offers_at_10C": cold, "medians_10C": {f"{CM.RATE_LABEL[c]} h {h:g}": case(10.0, c, h)["median"] for c in CM.RATES for h in H}}

    r = [x for x in rc if x["ambient_C"] == 10.0 and x["crate"] == 1.0 and x["h"] == 10.0][0]
    sc["P6"] = {"held": r["dT_cell_K"] < spread(10.0, 1.0, 10.0) / 3.0, "within_cell_K": r["dT_cell_K"],
                "module_spread_K": spread(10.0, 1.0, 10.0), "ratio": r["dT_cell_K"] / spread(10.0, 1.0, 10.0)}

    def ex_row(ch, dT):
        rr = [e for e in ex if e["part"] == "C2" and str(e["chamber"]).startswith(ch) and e["central"] in (True, "True")
              and e["dT_K"] == dT and e["floor"] == "Geslin"]
        if len(rr) != 1:
            raise SystemExit("P7: no unique experiment row")
        return rr[0]
    cr, wr = ex_row("cold", 5.0), ex_row("warm", 2.0)
    a, b = cr["n_pairs"] < 10, wr["n_pairs"] > 30
    v4_ok = all(str(e.get("v4_pass")) == "True" for e in (cr, wr))
    held = a and b
    note = None
    if not v4_ok:
        a2 = (cr["n_sim_lognormal"] < 10) if str(cr.get("v4_pass")) != "True" else a
        b2 = (wr["n_sim_lognormal"] > 30) if str(wr.get("v4_pass")) != "True" else b
        if (a2 and b2) != held:
            held, note = None, "not scored: V4 fails in a P7 case and the simulated n gives the other verdict (D63)"
    sc["P7"] = {"held": held, "note": note, "cold_n": cr["n_pairs"], "cold_lever": cr["lever"], "warm_n": wr["n_pairs"],
                "warm_lever": wr["lever"], "cold_n_sim": cr.get("n_sim_lognormal"), "warm_n_sim": wr.get("n_sim_lognormal"),
                "v4_pass_both": v4_ok, "cold_n_normal_shape": cr["n_pairs_normal_shape"], "warm_n_normal_shape": wr["n_pairs_normal_shape"]}
    return sc


def main(results):
    floors = CM.read_floors()
    sc = score(*load(results), floors)
    sc["_summary"] = {k: v["held"] for k, v in sc.items() if k.startswith("P")}
    json.dump(sc, open(os.path.join(results, "stage4_scorecard.json"), "w"), indent=1)
    print(json.dumps(sc, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=os.path.join(CM.ROOT, "results"))
    main(ap.parse_args().results)
