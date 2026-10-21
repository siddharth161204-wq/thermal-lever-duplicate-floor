#!/usr/bin/env python3
"""Stage 2 on the KIT LG INR18650HG2 dataset (Luh and Blank, doi:10.35097/1969): lifetimes, cell temperatures, the
lever per protocol, the replicate floor, the warmer-replicate test, the module translation and the scoring, as
DECISIONS.md D27 to D34 define them.

Usage: python3 src/stage2/kit.py --eoc DIR_WITH_cell_eocv2_CSVs --out results
"""
import argparse, glob, json, math, os, sys
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats as SS

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
Q_NOM, R_GAS = 3.0, 8.314462618
PELTIER_UTC = pd.Timestamp("2023-04-23T00:00:00Z").timestamp()
NOV_EXC = (pd.Timestamp("2022-11-27T00:00:00Z").timestamp(), pd.Timestamp("2022-11-30T00:00:00Z").timestamp())
OFFER_REG, OFFER_DYN = 0.391279, 1.17          # K: results/module_spread.csv (C/2, h = 10); D24 post hoc
LAW1, LAW2 = -0.0099, 0.052                    # per K, Stage 1 (README Table 1)
SEED_W8 = 20261021
COLS = ["timestamp_s", "timestamp_origin", "cyc_condition", "cyc_charged", "v_min_target_V", "i_dischg_min_A",
        "cap_aged_est_Ah", "delta_q_dischg_Ah", "t_start_degC", "t_end_degC", "cyc_duration_s", "num_cycles_checkup",
        "total_q_dischg_sum_Ah"]


def load(path):
    e = pd.read_csv(path, sep=";", usecols=COLS)
    return e[e.timestamp_origin != 5].sort_values("timestamp_s").reset_index(drop=True)


def checkups(e):
    c = e[(e.cyc_condition == 2) & (e.cyc_charged == 0) & ((e.v_min_target_V - 2.5).abs() < 1e-3)
          & ((e.i_dischg_min_A + 1.0).abs() < 1e-3)].copy()
    c["q"] = c.delta_q_dischg_Ah.abs()
    c = c.loc[c.groupby("num_cycles_checkup").q.idxmax()].sort_values("timestamp_s")
    c["efc"] = (c.total_q_dischg_sum_Ah - c.q) / Q_NOM
    return c[["num_cycles_checkup", "timestamp_s", "efc", "q"]].reset_index(drop=True)


def crossing(x, y, thr, t=None):
    for k in range(1, len(y)):
        if y[k] < thr:
            if y[k - 1] < thr:
                return np.nan, np.nan, np.nan
            w = (y[k - 1] - thr) / (y[k - 1] - y[k])
            tt = t[k - 1] + w * (t[k] - t[k - 1]) if t is not None else np.nan
            return x[k - 1] + w * (x[k] - x[k - 1]), x[k] - x[k - 1], tt
    return np.nan, np.nan, np.nan


def wmean_T(r):
    w = r.cyc_duration_s.to_numpy()
    return float(np.sum(0.5 * (r.t_start_degC + r.t_end_degC) * w) / w.sum()) if w.sum() > 0 else np.nan


def cell_record(meta, e):
    cu = checkups(e)
    q0 = float(cu.q.iloc[0])
    out = dict(meta, n_checkups=len(cu), q_first_Ah=q0, efc_last_checkup=float(cu.efc.iloc[-1]),
               t_first_checkup=float(cu.timestamp_s.iloc[0]), t_last_checkup=float(cu.timestamp_s.iloc[-1]))
    for f in (0.90, 0.80):
        L, res, tL = crossing(cu.efc.to_numpy(), cu.q.to_numpy(), f * q0, cu.timestamp_s.to_numpy())
        tag = int(f * 100)
        out[f"L{tag}_efc"], out[f"res{tag}_efc"], out[f"t_cross{tag}"] = L, res, tL
        out[f"censored{tag}"] = bool(not np.isfinite(L))
    reg = e[(e.cyc_condition == 1) & (e.cyc_duration_s > 0)].copy()
    reg["efc"] = reg.total_q_dischg_sum_Ah / Q_NOM
    reg = reg[reg.timestamp_s >= out["t_first_checkup"]]
    upto = reg if out["censored90"] else reg[reg.efc <= out["L90_efc"]]
    out["T_C"] = wmean_T(upto)
    out["T_pre_peltier_C"] = wmean_T(upto[upto.timestamp_s < PELTIER_UTC])
    out["T_post_peltier_C"] = wmean_T(upto[upto.timestamp_s >= PELTIER_UTC])
    out["runs_in_T"] = int(len(upto))
    rm = reg.cap_aged_est_Ah.rolling(11, center=True, min_periods=1).median().to_numpy()
    out["L90_est_efc"] = crossing(reg.efc.to_numpy(), rm, 0.90 * q0)[0]
    end = out["t_cross90"] if not out["censored90"] else out["t_last_checkup"]
    out["nov2022_inside_window"] = bool(out["t_first_checkup"] <= NOV_EXC[1] and end >= NOV_EXC[0])
    return out, cu


def set_stats(cells, Lcol, ccol):
    g = []
    for pid, s in cells.groupby("parameter_id"):
        u = s[~s[ccol]]
        g.append({"parameter_id": pid, "protocol": s.protocol.iloc[0], "age_temp": int(s.age_temp.iloc[0]),
                  "n_uncensored": len(u), "enters": len(u) >= 2,
                  "mean_lnL": float(np.log(u[Lcol]).mean()) if len(u) else np.nan,
                  "var_lnL": float(np.log(u[Lcol]).var(ddof=1)) if len(u) >= 2 else np.nan,
                  "T_group_C": float(u.T_C.mean()) if len(u) else np.nan,
                  "T_group_pre_peltier_C": float(u.T_pre_peltier_C.mean()) if len(u) else np.nan,
                  "lower_bound_efc": float(s.loc[s[ccol], "efc_last_checkup"].min()) if s[ccol].any() else np.nan})
    g = pd.DataFrame(g)
    e = g[g.enters]
    sigma2 = float(((e.n_uncensored - 1) * e.var_lnL).sum() / (e.n_uncensored - 1).sum())
    return g, sigma2


def levers(g, sigma2, lo, hi):
    rows = []
    for prot, s in g.groupby("protocol"):
        a, b = s[s.age_temp == lo].iloc[0], s[s.age_temp == hi].iloc[0]
        r = {"protocol": prot, "T_lo_set": lo, "T_hi_set": hi, "determinable": bool(a.enters and b.enters),
             "n_lo": int(a.n_uncensored), "n_hi": int(b.n_uncensored)}
        if r["determinable"]:
            dT = b.T_group_C - a.T_group_C
            dl = b.mean_lnL - a.mean_lnL
            se = math.sqrt(sigma2 * (1 / a.n_uncensored + 1 / b.n_uncensored)) / dT
            Tm = 0.5 * (a.T_group_C + b.T_group_C) + 273.15
            r.update({"T_lo_C": a.T_group_C, "T_hi_C": b.T_group_C, "dT_K": dT, "s_per_K": -dl / dT, "se": se,
                      "ci95": [-dl / dT - 1.96 * se, -dl / dT + 1.96 * se], "s_setpoint_per_K": -dl / (hi - lo),
                      "Ea_kJ_per_mol": -dl / dT * R_GAS * Tm * Tm / 1000.0})
            if lo == 25:
                r["s_25_without_after_peltier"] = -dl / (b.T_group_C - a.T_group_pre_peltier_C)
            r["module_reg"] = abs(r["s_per_K"]) * OFFER_REG
            r["module_dyn"] = abs(r["s_per_K"]) * OFFER_DYN
        else:
            r["lower_bounds_efc"] = {str(lo): a.lower_bound_efc, str(hi): b.lower_bound_efc}
        rows.append(r)
    return rows


def floor(cells, sigma2, Lcol, ccol):
    pre, allp = [], []
    for pid, s in cells.groupby("parameter_id"):
        s = s.set_index("parameter_nr")
        if 1 in s.index and 2 in s.index and not s.loc[1, ccol] and not s.loc[2, ccol]:
            L1, L2 = s.loc[1, Lcol], s.loc[2, Lcol]
            pre.append({"parameter_id": pid, "protocol": s.protocol.iloc[0], "age_temp": int(s.age_temp.iloc[0]),
                        "L1": L1, "L2": L2, "rel_diff": abs(L1 - L2) / (0.5 * (L1 + L2))})
        u = s[~s[ccol]]
        for i, j in combinations(u.index, 2):
            allp.append(abs(u.loc[i, Lcol] - u.loc[j, Lcol]) / (0.5 * (u.loc[i, Lcol] + u.loc[j, Lcol])))
    pre = pd.DataFrame(pre)
    rng = np.random.default_rng(SEED_W8)
    L = np.exp(rng.normal(0.0, math.sqrt(sigma2), size=(200_000, 8)))
    w8 = float(np.median((L.mean(1) - L.min(1)) / L.mean(1)))
    return pre, {"pairs": int(len(pre)), "mean": float(pre.rel_diff.mean()), "median": float(pre.rel_diff.median()),
                 "p95": float(np.percentile(pre.rel_diff, 95)), "worst_of_8": w8, "sigma_lnL": math.sqrt(sigma2),
                 "all_pairs": {"pairs": len(allp), "mean": float(np.mean(allp)), "p95": float(np.percentile(allp, 95))}}


def warmer(cells, Lcol, ccol):
    rows = []
    for pid, s in cells.groupby("parameter_id"):
        u = s[~s[ccol]]
        for (_, a), (_, b) in combinations(u.iterrows(), 2):
            if a.T_C == b.T_C:
                continue
            w, c = (a, b) if a.T_C > b.T_C else (b, a)
            rows.append({"parameter_id": pid, "age_temp": int(a.age_temp), "protocol": a.protocol, "warm": w.cell, "cool": c.cell,
                         "dT_K": w.T_C - c.T_C, "L_warm": w[Lcol], "L_cool": c[Lcol], "warmer_shorter": bool(w[Lcol] < c[Lcol])})
    p = pd.DataFrame(rows)
    def summ(x):
        k, n = int(x.warmer_shorter.sum()), len(x)
        return {"pairs": n, "warmer_shorter": k, "share": k / n if n else np.nan,
                "p_one_sided": float(SS.binomtest(k, n, 0.5, alternative="greater").pvalue) if n else np.nan,
                "dT_median_K": float(x.dT_K.median()) if n else np.nan}
    res = {"all": summ(p)}
    for t, x in p.groupby("age_temp"):
        res[f"{t} C"] = summ(x)
    return p, res


def score(lev, lev_lo, flo, war):
    s = {r["protocol"]: r for r in lev}
    sl = {r["protocol"]: r for r in lev_lo}
    low = [p for p in s if "0.33/" in p]
    sc = {}
    if all(s[p]["determinable"] for p in low):
        m = float(np.mean([s[p]["s_per_K"] for p in low]))
        sc["P1"] = {"held": bool(all(s[p]["s_per_K"] > 0 for p in low) and 0.01 <= m <= 0.06), "mean": m,
                    "values": {p: s[p]["s_per_K"] for p in low}}
        sc["P7"] = {"held": bool(LAW1 < m < LAW2), "mean": m}
    else:
        sc["P1"] = sc["P7"] = {"held": None, "scorable": False, "not_determinable": [p for p in low if not s[p]["determinable"]]}
    p2 = {}
    for w in ("0-100", "10-100", "10-90"):
        hi, lo = s[f"{w} 1.67/1.00C"], s[f"{w} 0.33/1.00C"]
        p2[w] = {"held": bool(hi["s_per_K"] < lo["s_per_K"]), "s_1.67/1.00": hi["s_per_K"], "s_0.33/1.00": lo["s_per_K"],
                 "s_0.33/0.33 (unscored)": s[f"{w} 0.33/0.33C"].get("s_per_K")} if hi["determinable"] and lo["determinable"] else {"held": None}
    sc["P2"] = {"held": None if any(v["held"] is None for v in p2.values()) else all(v["held"] for v in p2.values()), "windows": p2}
    det = [r["s_per_K"] for r in lev if r["determinable"]]
    if not det:
        return dict(sc, P3={"held": None, "scorable": False}, P4=None, P5=None, P6={"held": None, "scorable": False})
    sc["P3"] = {"held": bool(min(det) < 0 < max(det)), "n_determinable": len(det), "min": min(det), "max": max(det)}
    sc["P4"] = {"held": bool(0.015 <= flo["mean"] <= 0.06 and 0.044 <= flo["p95"] <= 0.174), "mean": flo["mean"], "p95": flo["p95"], "pairs": flo["pairs"]}
    sc["P5"] = {"held": bool(0.45 <= war["all"]["share"] <= 0.65), **war["all"]}
    sc["P6"] = {"held": bool(max(det) <= 0.077), "max": max(det)}
    fast = [p for p in sl if ("1.00/" in p.split()[1] or "1.67/" in p.split()[1])]
    neg = [p for p in fast if sl[p]["determinable"] and sl[p]["s_per_K"] < 0]
    sc["P8"] = {"held": bool(len(neg) >= 3), "negative": neg, "of": fast,
                "values": {p: sl[p].get("s_per_K") for p in fast}}
    return sc


def main(eoc_dir, out):
    cfg = pd.read_csv(os.path.join(ROOT, "inputs", "kit", "cyclic_cells.csv"))
    recs, cus = [], []
    for _, m in cfg.iterrows():
        f = glob.glob(os.path.join(eoc_dir, f"cell_eocv2_P{m.parameter_id:03d}_{m.parameter_nr}_*.csv"))
        assert len(f) == 1, (m.cell, f)
        r, cu = cell_record(m[["cell", "parameter_id", "parameter_nr", "slave_channel", "age_temp", "window", "age_chg_rate",
                               "age_dischg_rate", "protocol"]].to_dict(), load(f[0]))
        recs.append(r)
        cus.append(cu.assign(cell=m.cell))
    cells = pd.DataFrame(recs)
    res = {"cells": len(cells), "censored90": int(cells.censored90.sum()), "censored80": int(cells.censored80.sum())}
    g90, s2 = set_stats(cells, "L90_efc", "censored90")
    g80, s2_80 = set_stats(cells, "L80_efc", "censored80")
    lev = levers(g90, s2, 25, 40)
    lev_lo = levers(g90, s2, 10, 25)
    lev80 = levers(g80, s2_80, 25, 40)
    pre, flo = floor(cells, s2, "L90_efc", "censored90")
    wp, war = warmer(cells, "L90_efc", "censored90")
    pre80, flo80 = floor(cells, s2_80, "L80_efc", "censored80")          # secondary threshold (D28)
    res["resolution_over_L90"] = {"median": float((cells.res90_efc / cells.L90_efc).median()),
                                  "p95": float((cells.res90_efc / cells.L90_efc).quantile(0.95))}
    res["resolution_over_L80"] = {"median": float((cells.res80_efc / cells.L80_efc).median()),
                                  "p95": float((cells.res80_efc / cells.L80_efc).quantile(0.95))}
    sc = score(lev, lev_lo, flo, war)
    res.update({"sigma2_lnL": s2, "sets_entering": int(g90.enters.sum()), "lever_25_40": lev, "lever_10_25": lev_lo,
                "lever_25_40_at_80pct": lev80, "floor": flo, "floor_at_80pct": flo80, "warmer": war, "scorecard": sc,
                "module_offers_K": {"registered": OFFER_REG, "dynamic_post_hoc": OFFER_DYN}})
    os.makedirs(out, exist_ok=True)
    cells.to_csv(os.path.join(out, "stage2_cells.csv"), index=False)
    pd.concat(cus).to_csv(os.path.join(out, "stage2_checkups.csv"), index=False)
    g90.to_csv(os.path.join(out, "stage2_sets.csv"), index=False)
    pd.DataFrame([{k: v for k, v in r.items() if not isinstance(v, (list, dict))} for r in lev + lev_lo]).to_csv(
        os.path.join(out, "stage2_lever.csv"), index=False)
    pre.to_csv(os.path.join(out, "stage2_pairs.csv"), index=False)
    wp.to_csv(os.path.join(out, "stage2_warmer_pairs.csv"), index=False)
    json.dump(res, open(os.path.join(out, "stage2_lever.json"), "w"), indent=1, default=float)
    json.dump(sc, open(os.path.join(out, "stage2_scorecard.json"), "w"), indent=1, default=float)
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--eoc", required=True)
    ap.add_argument("--out", default="results")
    a = ap.parse_args()
    r = main(a.eoc, a.out)
    print(json.dumps({k: r[k] for k in ("cells", "censored90", "censored80", "sets_entering", "sigma2_lnL")}))
