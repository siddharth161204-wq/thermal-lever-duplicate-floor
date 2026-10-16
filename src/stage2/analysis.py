#!/usr/bin/env python3
"""Stage 2, computations 1 to 4 and the scoring of predictions_stage2.md.

Reads the per-cell table built by process_cells.py (conditions, lifetime, resolution, censoring and the two
temperature measures) and the Stage 1 results; writes results/stage2_lever.csv and .json, results/stage2_pairs.csv,
results/stage2_thermal_share.json, results/stage2_module.json and results/stage2_scorecard.json.

1. The measured lever per chemistry: s = -(mean ln L at 35 C - mean ln L at 25 C) / 10 K on the registered 1C
   temperature series (0 to 100 %, 0.5C charge), with a bootstrap over cells within each temperature; the 15 to 25 C
   secant and the curvature; the 2C series and a pooled common-slope fit as secondary values; both laws over the same
   secants at the same conditions and threshold.
2. The replicate floor per chemistry: |L_a - L_b| / mean(L_a, L_b) over every pair of every replicate group.
3. The thermal share: per pair, the self-heating difference, the predicted ln ratio -s dT and the observed ln ratio;
   the sign share with the sign test, the slope through the origin with a bootstrap over replicate groups, and the
   ratio of root-mean-square predicted to observed ln ratios; by discharge-rate group; measure (b) as secondary.
4. The measured NCA lever on the module's largest offer at the Geslin rates (Stage 1) and on the pairs' measured
   self-heating differences, against the Geslin floor.

Usage: python3 src/stage2/analysis.py --cells results/stage2_cells.csv --stage1 results --out results
"""
import argparse, itertools, json, math, os, sys

import numpy as np
import pandas as pd
from scipy import stats as SS

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import laws as LW  # noqa: E402

SEED = 20261016
N_BOOT = 10000
TEMPS = (15.0, 25.0, 35.0)
CHEMS = ("LFP", "NCA", "NMC")
REGISTERED_RATE = 1.0
SECONDARY_RATE = 2.0


def temperature_series(cells, chem, c_dis):
    return cells[(cells.chemistry == chem) & (cells.soc_min == 0) & (cells.soc_max == 100) & (cells.c_charge == 0.5)
                 & (cells.c_discharge == c_dis) & (cells.ambient_C.isin(TEMPS))]


def _secants(lnL):
    m = {T: float(np.mean(v)) for T, v in lnL.items()}
    s2535 = -(m[35.0] - m[25.0]) / 10.0
    s1525 = -(m[25.0] - m[15.0]) / 10.0 if 15.0 in m else np.nan
    return s2535, s1525, m


def lever(cells, chem, c_dis, frac, rng):
    ts = temperature_series(cells, chem, c_dis)
    out = {"chemistry": chem, "c_discharge": c_dis, "soh_threshold": frac,
           "n_cells": {f"{int(T)}C": int((ts.ambient_C == T).sum()) for T in TEMPS},
           "censored_cells": sorted(ts.loc[ts.censored, "cell_id"].tolist())}
    if ts.censored.any() or any((ts.ambient_C == T).sum() == 0 for T in (25.0, 35.0)):
        out["determinable"] = False
        return out
    lnL = {T: np.log(ts.loc[ts.ambient_C == T, "L_efc"].to_numpy()) for T in TEMPS if (ts.ambient_C == T).any()}
    s2535, s1525, m = _secants(lnL)
    boots = np.empty((N_BOOT, 2))
    for b in range(N_BOOT):
        rs = {T: rng.choice(v, size=v.size, replace=True) for T, v in lnL.items()}
        boots[b] = _secants(rs)[:2]
    lo, hi = np.percentile(boots[:, 0], [2.5, 97.5])
    out.update({"determinable": True, "s_25_35_per_K": s2535, "s_25_35_ci95": [float(lo), float(hi)],
                "s_15_25_per_K": s1525,
                "s_15_25_ci95": [float(x) for x in np.percentile(boots[:, 1], [2.5, 97.5])] if np.isfinite(s1525) else None,
                "curvature_per_K2": (s2535 - s1525) / 10.0 if np.isfinite(s1525) else np.nan,
                "mean_lnL": {f"{int(T)}C": v for T, v in m.items()},
                "L_efc": {f"{int(T)}C": ts.loc[ts.ambient_C == T, "L_efc"].round(3).tolist() for T in TEMPS},
                "cells": {f"{int(T)}C": ts.loc[ts.ambient_C == T, "cell_id"].tolist() for T in TEMPS}})
    sec = LW.law_secants(25.0, 35.0, c_dis, frac)
    sec15 = LW.law_secants(15.0, 25.0, c_dis, frac)
    out["laws_25_35_per_K"] = sec
    out["laws_15_25_per_K"] = sec15
    out["law_curves"] = LW.law_curves([15.0, 20.0, 25.0, 30.0, 35.0], c_dis, frac)
    return out


def pooled_lever(cells, chem, frac):
    """Secondary: one temperature slope with an intercept per discharge rate, 25 and 35 C, 1C and 2C series."""
    rows = []
    for c in (REGISTERED_RATE, SECONDARY_RATE):
        ts = temperature_series(cells, chem, c)
        ts = ts[ts.ambient_C.isin([25.0, 35.0])]
        if ts.censored.any() or ts.empty:
            return {"determinable": False}
        rows += [(c, T, math.log(L)) for T, L in zip(ts.ambient_C, ts.L_efc)]
    r = np.array(rows)
    X = np.column_stack([r[:, 1], (r[:, 0] == REGISTERED_RATE).astype(float), (r[:, 0] == SECONDARY_RATE).astype(float)])
    beta, *_ = np.linalg.lstsq(X, r[:, 2], rcond=None)
    return {"determinable": True, "s_per_K": float(-beta[0]), "n_cells": int(len(r))}


def pair_rows(cells, s_by_chem):
    rows = []
    for g, sub in cells.groupby("group", sort=True):
        sub = sub.sort_values("replicate")
        n = len(sub)
        for (_, a), (_, b) in itertools.combinations(sub.iterrows(), 2):
            row = {"group": g, "chemistry": a.chemistry, "ambient_C": a.ambient_C, "window": f"{int(a.soc_min)}-{int(a.soc_max)}",
                   "c_charge": a.c_charge, "c_discharge": a.c_discharge,
                   "rate_group": "1C or below" if a.c_discharge <= 1.0 else "2C and 3C",
                   "cells_in_group": n, "cell_a": a.cell_id, "cell_b": b.cell_id,
                   "L_a_efc": a.L_efc, "L_b_efc": b.L_efc, "censored_a": bool(a.censored), "censored_b": bool(b.censored),
                   "efc_last_check_a": a.efc_last_check, "efc_last_check_b": b.efc_last_check,
                   "resolution_efc": float(np.nanmax([a.resolution_efc, b.resolution_efc])) if not (a.censored or b.censored) else np.nan,
                   "rise_a_K": a.get("rise_K", np.nan), "rise_b_K": b.get("rise_K", np.nan),
                   "rest_a_C": a.get("rest_T_C", np.nan), "rest_b_C": b.get("rest_T_C", np.nan),
                   "s_chem_per_K": s_by_chem.get(a.chemistry, np.nan)}
            if not (a.censored or b.censored):
                row["rel_diff"] = abs(a.L_efc - b.L_efc) / (0.5 * (a.L_efc + b.L_efc))
                row["below_resolution"] = abs(a.L_efc - b.L_efc) < row["resolution_efc"]
            else:
                row["rel_diff"] = np.nan
                row["below_resolution"] = np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def floor_stats(pairs, floor):
    out = {}
    for chem in CHEMS:
        p = pairs[pairs.chemistry == chem]
        ok = p[p.rel_diff.notna()]
        out[chem] = {"pairs_total": int(len(p)), "pairs_with_difference": int(len(ok)),
                     "pairs_censored": int(len(p) - len(ok)),
                     "groups": int(p.group.nunique()), "groups_with_3plus_cells": int((p.groupby("group").cells_in_group.first() >= 3).sum()),
                     "mean": float(ok.rel_diff.mean()) if len(ok) else np.nan,
                     "median": float(ok.rel_diff.median()) if len(ok) else np.nan,
                     "p95": float(np.percentile(ok.rel_diff, 95)) if len(ok) else np.nan,
                     "below_resolution": int(ok.below_resolution.sum()) if len(ok) else 0}
    out["geslin_floor"] = {"mean": floor["mean"], "p95": floor["p95"], "n_pairs": floor["n_pairs"]}
    return out


def oriented(pairs, measure):
    """Orient each pair so that cell w is the warmer by the measure ('rise' for (a), 'rest' for (b))."""
    ka, kb = ("rise_a_K", "rise_b_K") if measure == "rise" else ("rest_a_C", "rest_b_C")
    rows = []
    for _, p in pairs.iterrows():
        ta, tb = p[ka], p[kb]
        if not (np.isfinite(ta) and np.isfinite(tb)) or ta == tb:
            continue
        w, c = ("a", "b") if ta > tb else ("b", "a")
        dT = abs(ta - tb)
        Lw, Lc = p[f"L_{w}_efc"], p[f"L_{c}_efc"]
        cw, cc = p[f"censored_{w}"], p[f"censored_{c}"]
        r = {"group": p.group, "chemistry": p.chemistry, "rate_group": p.rate_group, "c_discharge": p.c_discharge,
             "warm": p[f"cell_{w}"], "cool": p[f"cell_{c}"], "dT_K": dT, "pred_ln_ratio": -p.s_chem_per_K * dT}
        if not (cw or cc):
            r["obs_ln_ratio"] = math.log(Lw / Lc)
            if abs(Lw - Lc) < p.resolution_efc:
                r["outcome"] = "tie"
            else:
                r["outcome"] = "warmer shorter" if Lw < Lc else "warmer longer"
        elif cw and cc:
            r["obs_ln_ratio"] = np.nan; r["outcome"] = "both censored"
        elif cw:                                   # warm censored: it lived longer if its last check passed the other's L
            r["obs_ln_ratio"] = np.nan
            r["outcome"] = "warmer longer" if p[f"efc_last_check_{w}"] > Lc else "censored, order unknown"
        else:
            r["obs_ln_ratio"] = np.nan
            r["outcome"] = "warmer shorter" if p[f"efc_last_check_{c}"] > Lw else "censored, order unknown"
        rows.append(r)
    return pd.DataFrame(rows)


def share_stats(o, rng):
    n_ws = int((o.outcome == "warmer shorter").sum())
    n_wl = int((o.outcome == "warmer longer").sum())
    n_tie = int((o.outcome == "tie").sum())
    n = n_ws + n_wl
    res = {"pairs": int(len(o)), "warmer_shorter": n_ws, "warmer_longer": n_wl, "ties": n_tie,
           "excluded_censored": int(o.outcome.isin(["both censored", "censored, order unknown"]).sum()),
           "share": n_ws / n if n else np.nan,
           "sign_test_p_one_sided": float(SS.binomtest(n_ws, n, 0.5, alternative="greater").pvalue) if n else np.nan,
           "sign_test_p_two_sided": float(SS.binomtest(n_ws, n, 0.5).pvalue) if n else np.nan,
           "groups": int(o.group.nunique())}
    u = o[o.obs_ln_ratio.notna()]
    x, y = u.pred_ln_ratio.to_numpy(), u.obs_ln_ratio.to_numpy()
    res["pairs_for_slope"] = int(len(u))
    if len(u) and np.sum(x * x) > 0:
        res["slope"] = float(np.sum(x * y) / np.sum(x * x))
        gids = u.group.to_numpy()
        ug = np.unique(gids)
        idx = {g: np.flatnonzero(gids == g) for g in ug}
        bs = []
        for _ in range(N_BOOT):
            pick = np.concatenate([idx[g] for g in rng.choice(ug, size=ug.size, replace=True)])
            sx = np.sum(x[pick] ** 2)
            if sx > 0:
                bs.append(np.sum(x[pick] * y[pick]) / sx)
        res["slope_ci95"] = [float(v) for v in np.percentile(bs, [2.5, 97.5])]
        res["rms_pred"] = float(np.sqrt(np.mean(x * x)))
        res["rms_obs"] = float(np.sqrt(np.mean(y * y)))
        res["sd_ratio_pred_over_obs"] = res["rms_pred"] / res["rms_obs"] if res["rms_obs"] > 0 else np.nan
    res["dT_K_median"] = float(o.dT_K.median()) if len(o) else np.nan
    res["dT_K_mean"] = float(o.dT_K.mean()) if len(o) else np.nan
    res["dT_K_max"] = float(o.dT_K.max()) if len(o) else np.nan
    res["dT_below_1K_fraction"] = float((o.dT_K < 1.0).mean()) if len(o) else np.nan
    return res


def independent_pairs(o_all, pairs_all, measure):
    """Secondary: one pair per replicate group, its warmest and its coolest cell by the measure."""
    ka, kb = ("rise_a_K", "rise_b_K") if measure == "rise" else ("rest_a_C", "rest_b_C")
    keep = []
    for g, sub in o_all.groupby("group"):
        temps = {}
        for _, p in pairs_all[pairs_all.group == g].iterrows():
            temps[p.cell_a] = p[ka]
            temps[p.cell_b] = p[kb]
        temps = {k: v for k, v in temps.items() if np.isfinite(v)}
        if len(temps) < 2:
            continue
        hot, cold = max(temps, key=temps.get), min(temps, key=temps.get)
        keep += list(sub[(sub.warm == hot) & (sub.cool == cold)].index)
    return o_all.loc[keep]


def thermal_share(pairs, measure, rng):
    o = oriented(pairs, measure)
    out = {"measure": "self-heating rise (a)" if measure == "rise" else "rest-period skin temperature (b)",
           "all": share_stats(o, rng)}
    for rg in ("1C or below", "2C and 3C"):
        out[rg] = share_stats(o[o.rate_group == rg], rng)
    for chem in CHEMS:
        out[f"chemistry {chem}"] = share_stats(o[o.chemistry == chem], rng)
    ind = independent_pairs(o, pairs, measure)
    out["one pair per group (warmest and coolest)"] = share_stats(ind, rng)
    return out, o


def module_application(lev_nca, stage1_dir, pairs, floor):
    ms = pd.read_csv(os.path.join(stage1_dir, "module_spread.csv"))
    geslin = ms[ms.crate <= 0.5 + 1e-9]
    dT_mod = float(geslin.dT_module_K.max())
    row = geslin.loc[geslin.dT_module_K.idxmax()]
    s = lev_nca["s_25_35_per_K"]
    lo, hi = lev_nca["s_25_35_ci95"]
    nca = pairs[(pairs.chemistry == "NCA")]
    dTn = (nca.rise_a_K - nca.rise_b_K).abs().dropna()
    dTa = (pairs.rise_a_K - pairs.rise_b_K).abs().dropna()
    return {"s_nca_per_K": s, "s_nca_ci95": [lo, hi],
            "module_largest_dT_K": dT_mod, "module_case": {"crate": float(row.crate), "h": float(row.h)},
            "implied_divergence": abs(s) * dT_mod,
            "implied_divergence_ci95": sorted([abs(lo) * dT_mod, abs(hi) * dT_mod]) if lo * hi > 0 else [0.0, max(abs(lo), abs(hi)) * dT_mod],
            "pairs_nca": {"n": int(len(dTn)), "implied_mean": float(abs(s) * dTn.mean()) if len(dTn) else np.nan,
                          "implied_max": float(abs(s) * dTn.max()) if len(dTn) else np.nan},
            "pairs_all": {"n": int(len(dTa)), "implied_mean": float(abs(s) * dTa.mean()) if len(dTa) else np.nan,
                          "implied_max": float(abs(s) * dTa.max()) if len(dTa) else np.nan},
            "geslin_floor": {"mean": floor["mean"], "p95": floor["p95"], "half_mean_registered": 0.015}}


def score(lev, flo, share_a, mod):
    s = {c: lev[c]["s_25_35_per_K"] if lev[c].get("determinable") else np.nan for c in CHEMS}
    law2_lfp = lev["LFP"]["laws_25_35_per_K"]["law_2_schimpe2018"] if lev["LFP"].get("determinable") else np.nan
    sc = {}
    sc["P1"] = {"held": bool(0.01 <= s["NCA"] <= 0.06), "computed": s["NCA"], "band": [0.01, 0.06], "central": 0.03}
    sc["P2"] = {"held": bool(s["LFP"] > s["NCA"] and s["LFP"] > s["NMC"]), "computed": s}
    ratio = s["LFP"] / law2_lfp if law2_lfp else np.nan
    sc["P3"] = {"held": bool(np.sign(s["LFP"]) == np.sign(law2_lfp) and 0.5 <= ratio <= 2.0),
                "computed": {"s_LFP": s["LFP"], "law_2": law2_lfp, "ratio": ratio}}
    p4 = {}
    for c in CHEMS:
        if flo[c]["pairs_with_difference"] >= 5:
            p4[c] = {"mean": flo[c]["mean"], "within": bool(0.015 <= flo[c]["mean"] <= 0.06), "pairs": flo[c]["pairs_with_difference"]}
    sc["P4"] = {"held": bool(p4) and all(v["within"] for v in p4.values()), "computed": p4, "band": [0.015, 0.06]}
    low = share_a["1C or below"]
    sc["P5"] = {"held": bool(low["dT_below_1K_fraction"] >= 0.80), "computed": {"fraction_below_1K": low["dT_below_1K_fraction"],
                "median_dT_K": low["dT_K_median"], "pairs": low["pairs"]}, "central_K": 0.3}
    al = share_a["all"]
    sc["P6"] = {"held": bool(0.5 < al["share"] < 0.75), "computed": {"share": al["share"], "warmer_shorter": al["warmer_shorter"],
                "warmer_longer": al["warmer_longer"], "ties": al["ties"]}, "central": 0.60}
    hi, lo = share_a["2C and 3C"], share_a["1C or below"]
    n_hi, n_lo = hi["warmer_shorter"] + hi["warmer_longer"], lo["warmer_shorter"] + lo["warmer_longer"]
    if n_hi >= 5 and n_lo >= 5:
        sc["P7"] = {"scored": True, "held": bool(hi["share"] > lo["share"]), "computed": {"2C and 3C": hi["share"], "1C or below": lo["share"],
                    "pairs": {"2C and 3C": n_hi, "1C or below": n_lo}}}
    else:
        sc["P7"] = {"scored": False, "held": None, "computed": {"pairs": {"2C and 3C": n_hi, "1C or below": n_lo}}}
    sc["P8"] = {"held": bool(mod["implied_divergence"] < 0.015), "computed": mod["implied_divergence"], "limit": 0.015}
    return sc


def main(cells_csv, stage1_dir, out_dir, thresholds):
    rng = np.random.default_rng(SEED)
    cells = pd.read_csv(cells_csv)
    floor = json.load(open(os.path.join(os.path.dirname(os.path.dirname(HERE)), "inputs", "floor.json")))
    lev = {c: lever(cells, c, REGISTERED_RATE, thresholds[c], rng) for c in CHEMS}
    lev2 = {c: lever(cells, c, SECONDARY_RATE, thresholds[c], rng) for c in CHEMS}
    pooled = {c: pooled_lever(cells, c, thresholds[c]) for c in CHEMS}
    s_by_chem = {c: lev[c]["s_25_35_per_K"] for c in CHEMS if lev[c].get("determinable")}
    pairs = pair_rows(cells, s_by_chem)
    flo = floor_stats(pairs, floor)
    share_a, o_a = thermal_share(pairs, "rise", rng)
    share_b, o_b = thermal_share(pairs, "rest", rng)
    mod = module_application(lev["NCA"], stage1_dir, pairs, floor) if lev["NCA"].get("determinable") else {}
    sc = score(lev, flo, share_a, mod) if mod else {}
    os.makedirs(out_dir, exist_ok=True)
    lrows = []
    for tag, L in (("registered 1C", lev), ("secondary 2C", lev2)):
        for c in CHEMS:
            e = L[c]
            if not e.get("determinable"):
                lrows.append({"series": tag, "chemistry": c, "determinable": False}); continue
            lrows.append({"series": tag, "chemistry": c, "determinable": True, "soh_threshold": e["soh_threshold"],
                          "s_25_35_per_K": e["s_25_35_per_K"], "ci95_lo": e["s_25_35_ci95"][0], "ci95_hi": e["s_25_35_ci95"][1],
                          "s_15_25_per_K": e["s_15_25_per_K"], "curvature_per_K2": e["curvature_per_K2"],
                          "law1_25_35_per_K": e["laws_25_35_per_K"]["law_1_kucinskis2022"],
                          "law2_25_35_per_K": e["laws_25_35_per_K"]["law_2_schimpe2018"],
                          "law1_15_25_per_K": e["laws_15_25_per_K"]["law_1_kucinskis2022"],
                          "law2_15_25_per_K": e["laws_15_25_per_K"]["law_2_schimpe2018"],
                          "n_15C": e["n_cells"]["15C"], "n_25C": e["n_cells"]["25C"], "n_35C": e["n_cells"]["35C"]})
    pd.DataFrame(lrows).to_csv(os.path.join(out_dir, "stage2_lever.csv"), index=False)
    json.dump({"registered_rate_C": REGISTERED_RATE, "thresholds": thresholds, "seed": SEED, "n_boot": N_BOOT,
               "registered": lev, "secondary_2C": lev2, "pooled_1C_2C": pooled},
              open(os.path.join(out_dir, "stage2_lever.json"), "w"), indent=1, default=float)
    pairs.to_csv(os.path.join(out_dir, "stage2_pairs.csv"), index=False)
    o_a.to_csv(os.path.join(out_dir, "stage2_pairs_oriented_measure_a.csv"), index=False)
    o_b.to_csv(os.path.join(out_dir, "stage2_pairs_oriented_measure_b.csv"), index=False)
    json.dump({"floor": flo, "measure_a": share_a, "measure_b": share_b}, open(os.path.join(out_dir, "stage2_thermal_share.json"), "w"),
              indent=1, default=float)
    json.dump(mod, open(os.path.join(out_dir, "stage2_module.json"), "w"), indent=1, default=float)
    json.dump(sc, open(os.path.join(out_dir, "stage2_scorecard.json"), "w"), indent=1, default=float)
    return lev, flo, share_a, share_b, mod, sc


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", required=True)
    ap.add_argument("--stage1", default="results")
    ap.add_argument("--out", default="results")
    ap.add_argument("--thresholds", required=True, help='JSON, e.g. {"LFP": 0.95, "NCA": 0.90, "NMC": 0.90}')
    a = ap.parse_args()
    main(a.cells, a.stage1, a.out, json.loads(a.thresholds))
