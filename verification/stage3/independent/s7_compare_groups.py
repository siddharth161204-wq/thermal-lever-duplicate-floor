#!/usr/bin/env python3
"""Stage 3 independent check, step 7 (after the freeze): every group-level number and verdict of mine against the
pipeline's results/stage3_lever.json, stage3_floor.csv, stage3_pairs.csv, stage3_pairs_oriented_measure_[ab].csv,
stage3_thermal_share.json and stage3_scorecard.json.

Two runs of my code are compared:
  frozen    my primary readings, the outputs frozen in frozen/ (cells_full.json, groups.json);
  pipeline  my code run after the freeze with the two readings the per-cell comparison showed the pipeline uses
            (check voltage at or below cut-off + 0.02 V; a unit's discharge directly follows its charge), and with
            the pipeline's thermal-share layout (slopes only where every pair has a prediction; one pair per group =
            warmest and coolest), to test reproduction to rounding, bootstrap intervals included.
The frozen scripts s0 to s5 and lib.py are imported, not changed.

Usage: python3 -I s7_compare_groups.py <out_json>
"""
import copy, csv, itertools, json, math, os, sys
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib, s1_cells, s3_groups

REPO = "/home/claude/thermal-lever-duplicate-floor"
RES = os.path.join(REPO, "results")
CHEMS = ["LFP", "NCA", "NMC"]
NB = 10000
SEED = 20261024


# ----------------------------------------------------------------------------------------------- pipeline files
def load_pipeline():
    P = {"lever": json.load(open(os.path.join(RES, "stage3_lever.json"))),
         "share": json.load(open(os.path.join(RES, "stage3_thermal_share.json"))),
         "score": json.load(open(os.path.join(RES, "stage3_scorecard.json"))),
         "floor": list(csv.DictReader(open(os.path.join(RES, "stage3_floor.csv")))),
         "pairs": list(csv.DictReader(open(os.path.join(RES, "stage3_pairs.csv")))),
         "or_a": list(csv.DictReader(open(os.path.join(RES, "stage3_pairs_oriented_measure_a.csv")))),
         "or_b": list(csv.DictReader(open(os.path.join(RES, "stage3_pairs_oriented_measure_b.csv"))))}
    return P


def f(x):
    if x in ("", None):
        return float("nan")
    return float(x)


# ----------------------------------------------------------------------------------------------- comparison helpers
class Cmp:
    def __init__(self):
        self.rows = []

    def num(self, name, mine, pipe):
        mine = float("nan") if mine is None else float(mine)
        pipe = float("nan") if pipe is None else float(pipe)
        if math.isnan(mine) and math.isnan(pipe):
            self.rows.append((name, mine, pipe, 0.0, 0.0, "both missing"))
            return
        if math.isnan(mine) or math.isnan(pipe):
            self.rows.append((name, mine, pipe, float("inf"), float("inf"), "MISSING ON ONE SIDE"))
            return
        d = abs(mine - pipe)
        r = d / abs(pipe) if pipe != 0 else d
        self.rows.append((name, mine, pipe, d, r, ""))

    def eq(self, name, mine, pipe):
        self.rows.append((name, mine, pipe, 0.0 if mine == pipe else float("inf"), 0.0 if mine == pipe else float("inf"),
                          "" if mine == pipe else "DIFFERS"))

    def worst(self, prefix=""):
        rr = [r for r in self.rows if r[0].startswith(prefix)]
        fin = [r for r in rr if r[5] != "both missing"]
        ma = max(fin, key=lambda r: r[3]) if fin else None
        mr = max(fin, key=lambda r: r[4]) if fin else None
        return {"n": len(rr), "n_issues": sum(1 for r in rr if r[5] not in ("", "both missing")),
                "max_abs": None if ma is None else {"field": ma[0], "value": ma[3]},
                "max_rel": None if mr is None else {"field": mr[0], "value": mr[4]}}


SEC = {"25-35": "s_25_35", "15-25": "s_15_25"}
SER = {"1C": "registered_1C", "2C": "secondary_2C"}


def compare_levers(c, levers, P):
    for key, m in levers.items():
        chem, ser, sec = key.split("_")
        p = P["lever"][SER[ser]][chem][SEC[sec]]
        base = f"lever/{ser}/{chem}/{sec}"
        if p is None:
            c.eq(base + "/determinable", bool(m.get("group_temperatures_defined", True) and not math.isnan(m.get("s_per_K", float("nan")))), False)
            # the pipeline still reports the group temperatures and mean ln L of the series
            pl = P["lever"][SER[ser]][chem]
            lo, hi = sec.split("-")
            c.num(base + "/T_lo_C (series table)", m["T_lo_C"], pl["T_group_C"][lo + "C"])
            c.num(base + "/T_hi_C over cells with T_cyc (series table)", m["alt_T_over_cells_with_one"]["T_hi_C"], pl["T_group_C"][hi + "C"])
            c.num(base + "/mean_lnL_lo (series table)", m["mean_lnL_lo"], pl["mean_lnL"][lo + "C"])
            c.num(base + "/mean_lnL_hi (series table)", m["mean_lnL_hi"], pl["mean_lnL"][hi + "C"])
            continue
        c.eq(base + "/determinable", True, True)
        for mk, pk in (("s_per_K", "s"), ("T_lo_C", "T_lo_C"), ("T_hi_C", "T_hi_C"), ("dT_K", "dT_K"), ("mean_lnL_lo", "mean_lnL_lo"),
                       ("mean_lnL_hi", "mean_lnL_hi"), ("s_setpoint_per_K", "s_setpoint"), ("n_lo", "n_lo"), ("n_hi", "n_hi")):
            c.num(f"{base}/{mk}", m[mk], p[pk])
        c.num(base + "/ci95_lo", m["ci95"][0], p["ci95"][0])
        c.num(base + "/ci95_hi", m["ci95"][1], p["ci95"][1])
        c.num(base + "/law1", m["law1_s_per_K"], p["laws"]["law_1_kucinskis2022"])
        c.num(base + "/law2", m["law2_s_per_K"], p["laws"]["law_2_schimpe2018"])
        c.num(base + "/law1_setpoint", m["law1_s_setpoint_per_K"], p["laws_setpoint"]["law_1_kucinskis2022"])
        c.num(base + "/law2_setpoint", m["law2_s_setpoint_per_K"], p["laws_setpoint"]["law_2_schimpe2018"])


def compare_floor(c, floor, P):
    pf = {r["chemistry"]: r for r in P["floor"] if r["source"].startswith("Sandia")}
    for ch in CHEMS:
        m, p = floor[ch], pf[ch]
        for mk, pk in (("pairs_total", "pairs_total"), ("pairs_observed", "pairs_with_difference"), ("pairs_censored", "pairs_censored"),
                       ("groups_with_pairs", "groups"), ("mean", "mean"), ("median", "median"), ("p95", "p95"),
                       ("pairs_below_resolution", "below_resolution")):
            c.num(f"floor/{ch}/{mk}", m[mk], f(p[pk]))


def compare_pairs(c, pairs, P):
    pp = {(r["cell_a"], r["cell_b"]): r for r in P["pairs"]}
    c.eq("pairs/count", len(pairs), len(pp))
    for q in pairs:
        p = pp.get((q["a"], q["b"]))
        if p is None:
            c.eq(f"pairs/{q['a']}|{q['b']}/present", True, False)
            continue
        c.num(f"pairs/{q['a']}|{q['b']}/rel_diff", q.get("rel_diff"), f(p["rel_diff"]))
        c.num(f"pairs/{q['a']}|{q['b']}/resolution", q.get("resolution_efc"), f(p["resolution_efc"]))
        c.eq(f"pairs/{q['a']}|{q['b']}/below_resolution", bool(q.get("below_resolution")), p["below_resolution"] == "True")


def compare_oriented(c, ori, P, meas):
    pp = {frozenset((r["warm"], r["cool"])): r for r in P["or_" + meas]}
    mine = [q for q in ori if q["status"] != "no_measure"]
    c.eq(f"oriented_{meas}/count", len(mine), len(pp))
    for q in mine:
        p = pp.get(frozenset((q["warmer"], q["cooler"])))
        tag = f"oriented_{meas}/{q['warmer']}|{q['cooler']}"
        if p is None:
            c.eq(tag + "/present", True, False)
            continue
        c.eq(tag + "/warmer", q["warmer"], p["warm"])
        c.num(tag + "/dT", q["dT"], f(p["dT_K"]))
        c.num(tag + "/pred", q["predicted"], f(p["pred_ln_ratio"]))
        c.num(tag + "/obs", q.get("observed"), f(p["obs_ln_ratio"]))
        outcome = "tie" if q.get("tie") else ("warmer shorter" if q.get("warmer_shorter") else "warmer longer")
        c.eq(tag + "/outcome", outcome, p["outcome"])


def subset_stats(sub):
    ent = [q for q in sub if q.get("enters_sign_test")]
    k = sum(q["warmer_shorter"] for q in ent)
    dts = np.array([q["dT"] for q in sub]) if sub else np.array([])
    out = {"pairs": len(sub), "warmer_shorter": k, "warmer_longer": len(ent) - k, "ties": sum(1 for q in sub if q.get("tie")),
           "share": k / len(ent) if ent else float("nan"),
           "p1": float(stats.binomtest(k, len(ent), 0.5, alternative="greater").pvalue) if ent else float("nan"),
           "p2": float(stats.binomtest(k, len(ent), 0.5).pvalue) if ent else float("nan"),
           "groups": len({q["group"] for q in sub}),
           "dT_median": float(np.median(dts)) if len(dts) else float("nan"), "dT_mean": float(np.mean(dts)) if len(dts) else float("nan"),
           "dT_max": float(np.max(dts)) if len(dts) else float("nan"),
           "dT_below_1K": float(np.mean(dts < 1.0)) if len(dts) else float("nan")}
    items = [q for q in sub if q.get("status") == "both_observed"]
    if items and all(not math.isnan(q["predicted"]) for q in items):
        P_ = np.array([q["predicted"] for q in items])
        O_ = np.array([q["observed"] for q in items])
        out["slope"] = float(np.sum(P_ * O_) / np.sum(P_ * P_))
        out["rms_pred"] = float(np.sqrt(np.mean(P_ * P_)))
        out["rms_obs"] = float(np.sqrt(np.mean(O_ * O_)))
        out["ratio"] = out["rms_pred"] / out["rms_obs"]
        out["_items"] = items
    return out


def boot(items, rng, order):
    groups = [g for g in order if any(q["group"] == g for q in items)]
    PP = np.array([sum(q["predicted"] ** 2 for q in items if q["group"] == g) for g in groups])
    PO = np.array([sum(q["predicted"] * q["observed"] for q in items if q["group"] == g) for g in groups])
    sb = np.empty(NB)
    for it in range(NB):
        ix = rng.integers(0, len(groups), len(groups))
        sb[it] = PO[ix].sum() / PP[ix].sum()
    return [float(np.percentile(sb, 2.5)), float(np.percentile(sb, 97.5))]


def pipeline_layout_share(oriented, group_order, rng, draw_all=False):
    res = {}
    for meas in ("a", "b"):
        ori = [q for q in oriented[meas] if q["status"] != "no_measure"]
        subs = [("all", ori), ("1C or below", [q for q in ori if q["c_discharge"] <= 1]),
                ("2C and 3C", [q for q in ori if q["c_discharge"] >= 2])]
        subs += [(f"chemistry {ch}", [q for q in ori if q["chemistry"] == ch]) for ch in CHEMS]
        # one pair per group: the warmest and coolest cells of the group by the measure
        one = []
        for g in group_order:
            gq = [q for q in ori if q["group"] == g]
            if not gq:
                continue
            best = max(gq, key=lambda q: q["dT"])
            one.append(best)
        subs.append(("one pair per group (warmest and coolest)", one))
        r = {}
        for name, sub in subs:
            st = subset_stats(sub)
            if "slope" in st:
                st["ci"] = boot(st.pop("_items"), rng, group_order)
            elif draw_all and sub:
                items = [q for q in sub if q.get("status") == "both_observed"]
                gs = [g for g in group_order if any(q["group"] == g for q in items)]
                for it in range(NB):
                    rng.integers(0, len(gs), len(gs))
            st.pop("_items", None)
            r[name] = st
        res[meas] = r
    return res


def compare_share(c, share, P, tag):
    for meas in ("a", "b"):
        ps = P["share"]["measure_" + meas]
        for name, m in share[meas].items():
            p = ps[name]
            base = f"{tag}share_{meas}/{name}"
            for mk, pk in (("pairs", "pairs"), ("warmer_shorter", "warmer_shorter"), ("warmer_longer", "warmer_longer"),
                           ("ties", "ties"), ("share", "share"), ("p1", "sign_test_p_one_sided"), ("p2", "sign_test_p_two_sided"),
                           ("groups", "groups"), ("dT_median", "dT_K_median"), ("dT_mean", "dT_K_mean"), ("dT_max", "dT_K_max"),
                           ("dT_below_1K", "dT_below_1K_fraction")):
                c.num(f"{base}/{mk}", m.get(mk), p.get(pk))
            if "slope" in p or "slope" in m:
                c.num(base + "/slope", m.get("slope"), p.get("slope"))
                c.num(base + "/rms_pred", m.get("rms_pred"), p.get("rms_pred"))
                c.num(base + "/rms_obs", m.get("rms_obs"), p.get("rms_obs"))
                c.num(base + "/ratio", m.get("ratio"), p.get("sd_ratio_pred_over_obs"))
                if "ci" in m:
                    c.num(base + "/ci_lo", m["ci"][0], p["slope_ci95"][0])
                    c.num(base + "/ci_hi", m["ci"][1], p["slope_ci95"][1])


def compare_module(c, g, share_mine, P, tag):
    pm = P["share"]["module"]
    m = g["module"]
    c.num(tag + "module/s_nca", m["s_NCA"], pm["s_nca_per_K"])
    c.num(tag + "module/offer_K", m["offer_K"], pm["module_largest_dT_K"])
    c.num(tag + "module/implied", m["divergence"], pm["implied_divergence"])
    c.num(tag + "module/implied_ci_lo", m["divergence_ci95_from_s_ci"][0], pm["implied_divergence_ci95"][0])
    c.num(tag + "module/implied_ci_hi", m["divergence_ci95_from_s_ci"][1], pm["implied_divergence_ci95"][1])
    for mk, pk in (("NCA_pairs", "pairs_nca"), ("all_pairs", "pairs_all")):
        c.num(f"{tag}module/{mk}/n", m[mk]["n"], pm[pk]["n"])
        c.num(f"{tag}module/{mk}/mean", m[mk]["mean"], pm[pk]["implied_mean"])
        c.num(f"{tag}module/{mk}/max", m[mk]["max"], pm[pk]["implied_max"])


def compare_score(c, g, P, tag):
    S = P["score"]
    v = g["verdicts"]
    for k in [f"P{i}" for i in range(1, 11)]:
        p = S[k]
        pv = "not scored" if not p["scored"] else ("holds" if p["held"] else "did not hold")
        c.eq(f"{tag}verdict/{k}", v[k], pv)
    L = g["levers"]
    c.num(tag + "score/P1/s", L["NMC_1C_25-35"]["s_per_K"], S["P1"]["s_per_K"])
    c.num(tag + "score/P1/ci_lo", L["NMC_1C_25-35"]["ci95"][0], S["P1"]["ci95"][0])
    c.num(tag + "score/P1/ci_hi", L["NMC_1C_25-35"]["ci95"][1], S["P1"]["ci95"][1])
    c.num(tag + "score/P2/s_15_25", L["NMC_1C_15-25"]["s_per_K"], S["P2"]["s_15_25_per_K"])
    c.num(tag + "score/P3/s", L["NCA_1C_25-35"]["s_per_K"], S["P3"]["s_per_K"])
    for ch in CHEMS:
        c.num(f"{tag}score/P5/{ch}/mean", g["floor"][ch]["mean"], S["P5"]["chemistries"][ch]["mean"])
        c.num(f"{tag}score/P5/{ch}/pairs", g["floor"][ch]["pairs_observed"], S["P5"]["chemistries"][ch]["pairs"])
        c.num(f"{tag}score/P10/{ch}", L[f"{ch}_1C_15-25"]["s_per_K"], S["P10"]["s_15_25_per_K"][ch])
    c.num(tag + "score/P6/fraction", g["P6"]["share_below_1K"], S["P6"]["fraction_below_1K"])
    c.num(tag + "score/P6/median", g["P6"]["median_dT_K"], S["P6"]["median_dT_K"])
    c.num(tag + "score/P6/pairs", g["P6"]["n_pairs"], S["P6"]["pairs"])
    sa = g["thermal_share"]["a"]["all"]
    c.num(tag + "score/P7/share", sa["share"], S["P7"]["share"])
    c.num(tag + "score/P7/ties", sa["n_ties"], S["P7"]["ties"])
    c.num(tag + "score/P7/warmer_shorter", sa["k_warmer_shorter"], S["P7"]["warmer_shorter"])
    c.num(tag + "score/P7/p1", sa["p_one_sided"], S["P7"]["p_one_sided"])
    c.num(tag + "score/P7/p2", sa["p_two_sided"], S["P7"]["p_two_sided"])
    c.num(tag + "score/P8/n_2C3C", g["thermal_share"]["a"]["rate_2C_and_3C"]["n_entering"], S["P8"]["pairs_entering_sign_test"]["2C and 3C"])
    c.num(tag + "score/P8/n_1C", g["thermal_share"]["a"]["rate_1C_or_below"]["n_entering"], S["P8"]["pairs_entering_sign_test"]["1C or below"])
    c.num(tag + "score/P9/implied", g["module"]["divergence"], S["P9"]["implied_divergence"])
    c.num(tag + "score/P9/ci_lo", g["module"]["divergence_ci95_from_s_ci"][0], S["P9"]["implied_divergence_ci95"][0])
    c.num(tag + "score/P9/ci_hi", g["module"]["divergence_ci95_from_s_ci"][1], S["P9"]["implied_divergence_ci95"][1])


def advance_lever_draws(levers):
    """A generator positioned after the lever bootstraps, in the order of D44 and the brief, drawing nothing for a
    secant that is not determinable."""
    rng = np.random.default_rng(SEED)
    for ser in ("1C", "2C"):
        for chem in CHEMS:
            for sec in ("25-35", "15-25"):
                m = levers[f"{chem}_{ser}_{sec}"]
                if math.isnan(m.get("s_per_K", float("nan"))):
                    continue
                for it in range(NB):
                    rng.integers(0, m["n_lo"], m["n_lo"])
                    rng.integers(0, m["n_hi"], m["n_hi"])
    return rng


def run(tag, cells, P, out):
    g = s3_groups.analyse(copy.deepcopy(cells))
    c = Cmp()
    compare_levers(c, g["levers"], P)
    compare_floor(c, g["floor"], P)
    compare_pairs(c, g["pairs_all"], P)
    compare_oriented(c, g["oriented_pairs"]["a"], P, "a")
    compare_oriented(c, g["oriented_pairs"]["b"], P, "b")
    order = []
    for q in P["or_a"]:
        if q["group"] not in order:
            order.append(q["group"])
    rng = advance_lever_draws(g["levers"])
    sh = pipeline_layout_share(g["oriented_pairs"], order, rng)
    rng2 = advance_lever_draws(g["levers"])
    sh_all = pipeline_layout_share(g["oriented_pairs"], order, rng2, draw_all=True)
    compare_share(c, sh, P, "")
    c2 = Cmp()
    compare_share(c2, sh_all, P, "drawall_")
    compare_module(c, g, sh, P, "")
    compare_score(c, g, P, "")
    rec = {"summary": {k: c.worst(k) for k in ("lever/", "floor/", "pairs/", "oriented_a/", "oriented_b/", "share_a/", "share_b/",
                                                  "module/", "verdict/", "score/")},
           "overall": c.worst(""),
           "issues": [r for r in c.rows if r[5] not in ("", "both missing")],
           "share_ci_if_every_subset_draws": {"overall": c2.worst(""), "ci_rows": [r for r in c2.rows if "/ci_" in r[0]]},
           "rows": c.rows}
    out[tag] = rec
    return g, sh


def main(out_path):
    P = load_pipeline()
    out = {}
    frozen = json.load(open(os.path.join(HERE, "frozen", "cells_full.json")))
    run("frozen_primary", frozen, P, out)
    idx = lib.index()
    cells = [s1_cells.process(cid, idx, {"check_v": "one_sided", "units": "U3"}) for cid in sorted(idx)]
    g, sh = run("pipeline_readings", cells, P, out)
    out["pipeline_readings_share"] = sh
    for tag in ("frozen_primary", "pipeline_readings"):
        r = out[tag]
        print("=====", tag, "overall", r["overall"])
        for k, v in r["summary"].items():
            print(f"   {k:14s} n={v['n']:4d} issues={v['n_issues']:3d} max_abs={v['max_abs']} max_rel={v['max_rel']}")
        print("   issues (first 40):")
        for row in r["issues"][:40]:
            print("     ", row)
        print("   share CIs if every subset draws:", r["share_ci_if_every_subset_draws"]["overall"])
    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=1, default=float)


if __name__ == "__main__":
    main(sys.argv[1])
