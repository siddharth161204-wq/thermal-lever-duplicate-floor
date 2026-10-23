#!/usr/bin/env python3
"""Stage 3, post hoc (DECISIONS.md D51 to D54), written after the registered results were committed at 9c6bb84.
Nothing here changes a registered number or verdict; every value is labelled post hoc in results/stage3_posthoc.json.

1. P4 under two readings of the LFP 35 C group temperature, which D43 and D44 leave open when a cell has no T_cyc.
2. The floor without the groups discharged at 0.5C, and the NMC 20 to 80 % 0.5C cells with their low checks taken out.
3. The thermal-share pairs by the raw sign of the lifetime difference, ignoring the resolution; the slope over the NCA
   and NMC pairs together.
4. The lifetime resolution relative to the lifetime, and the pairs below resolution.
5. D41's last-voltage test read as two-sided (within 0.02 V of the cut-off) in the four NMC cells where it matters.
6. V1: what the check lifetime and the sensitivity lifetime measure in the cells that disagree.

Usage: python3 -I src/stage3/posthoc.py --work WORK_DIR --results results
"""
import argparse, json, os, sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "src", "stage2"))
import lifetimes as LT     # noqa: E402
import laws as LW          # noqa: E402
import stage3 as S3        # noqa: E402


def group(c, chem, T, cd):
    return c[(c.chemistry == chem) & (c.ambient_C == T) & (c.soc_min == 0) & (c.soc_max == 100) & (c.c_discharge == cd)].sort_values("cell_id")


def p4_readings(c, sc_reg):
    a, b = group(c, "LFP", 25.0, 1.0), group(c, "LFP", 35.0, 1.0)
    la, lb = np.log(a.L_efc.to_numpy()), np.log(b.L_efc.to_numpy())
    out = {"mean_lnL_25C": float(la.mean()), "mean_lnL_35C": float(lb.mean()),
           "cells_35C_without_T_cyc": b.loc[b.T_cyc_C.isna(), "cell_id"].tolist()}
    s_nca, s_nmc = sc_reg["P4_inputs"]["NCA"], sc_reg["P4_inputs"]["NMC"]
    for name, T1, T2 in (("group temperature over the cells that have one", a.T_cyc_C.mean(), b.T_cyc_C.dropna().mean()),
                         ("setpoint, 10 K", 25.0, 35.0)):
        s = -(lb.mean() - la.mean()) / (T2 - T1)
        l2 = LW.law_secants(T1, T2, 1.0, 0.95)["law_2_schimpe2018"]
        parts = {"exceeds_NCA": bool(s > s_nca), "exceeds_NMC": bool(s > s_nmc), "sign_of_law_2": bool(np.sign(s) == np.sign(l2)),
                 "ratio_within_factor_3": bool(1 / 3 <= s / l2 <= 3)}
        out[name] = {"T_lo_C": float(T1), "T_hi_C": float(T2), "s_per_K": float(s), "law_2_per_K": float(l2),
                     "ratio": float(s / l2), "parts": parts, "P4_would_hold": all(parts.values())}
    b3 = b[b.T_cyc_C.notna()]
    out["without the cell lacking T_cyc"] = {"s_per_K": float(-(np.log(b3.L_efc).mean() - la.mean()) / (b3.T_cyc_C.mean() - a.T_cyc_C.mean()))}
    return out


def floor_variants(pairs):
    q = pairs[pairs.rel_diff.notna()]

    def st(x):
        return {"pairs": int(len(x)), "groups": int(x.group.nunique()), "mean": float(x.rel_diff.mean()),
                "median": float(x.rel_diff.median()), "p95": float(np.percentile(x.rel_diff, 95)),
                "below_resolution": int(x.below_resolution.astype(bool).sum())}
    return {"registered, every pair": {ch: st(x) for ch, x in q.groupby("chemistry")},
            "groups discharged faster than 0.5C": {ch: st(x) for ch, x in q[q.c_discharge > 0.5].groupby("chemistry")},
            "groups discharged at 0.5C": {ch: st(x) for ch, x in q[q.c_discharge == 0.5].groupby("chemistry")}}


def low_checks(checks, cells):
    """Checks followed by a check more than 5 points of the first check higher: not a capacity fade."""
    ch = checks[checks.kept].copy()
    ch["frac"] = ch.q_Ah / ch.groupby("cell_id").q_Ah.transform("first")
    out, relife = [], {}
    for cid, x in ch.groupby("cell_id", sort=True):
        f, e = x.frac.to_numpy(), x.efc.to_numpy()
        dip = np.flatnonzero(np.r_[f[1:], np.nan] - f > 0.05)
        if not dip.size:
            continue
        thr = float(cells.loc[cid, "soh_threshold"])
        k = np.flatnonzero(f < thr)
        out.append({"cell_id": cid, "low_checks": [{"efc": float(e[i]), "frac": float(f[i]), "next_frac": float(f[i + 1])} for i in dip],
                    "registered_crossing_at_a_low_check": bool(k.size and k[0] in dip)})
        if k.size and k[0] in dip:
            y = x.drop(index=x.index[dip]).reset_index(drop=True)
            relife[cid] = {"registered_L_efc": float(cells.loc[cid, "L_efc"]),
                           **{kk: v for kk, v in LT.check_lifetime(y[["efc", "q_Ah"]], thr).items() if kk in ("L_efc", "resolution_efc", "censored")}}
    return out, relife


def nmc_partial_floor(pairs, relife):
    g = pairs[(pairs.group == "NMC_25C_20-80_0.5/0.5C")].copy()
    for side in ("a", "b"):
        g[f"L_{side}_ph"] = [relife.get(cid, {}).get("L_efc", L) for cid, L in zip(g[f"cell_{side}"], g[f"L_{side}_efc"])]
    g["rel_ph"] = (g.L_a_ph - g.L_b_ph).abs() / (0.5 * (g.L_a_ph + g.L_b_ph))
    other = pairs[(pairs.chemistry == "NMC") & (pairs.group != "NMC_25C_20-80_0.5/0.5C") & pairs.rel_diff.notna()].rel_diff
    allv = np.r_[other.to_numpy(), g.rel_ph.to_numpy()]
    return {"group_registered_mean": float(g.rel_diff.mean()), "group_without_low_checks_mean": float(g.rel_ph.mean()),
            "NMC_floor_without_low_checks": {"mean": float(allv.mean()), "median": float(np.median(allv)),
                                             "p95": float(np.percentile(allv, 95)), "pairs": int(allv.size)}}


def raw_sign_share(o):
    u = o[o.obs_ln_ratio.notna()]
    ws, wl = int((u.obs_ln_ratio < 0).sum()), int((u.obs_ln_ratio > 0).sum())
    from scipy import stats as SS
    v = u[u.pred_ln_ratio.notna()]
    x, y = v.pred_ln_ratio.to_numpy(), v.obs_ln_ratio.to_numpy()
    return {"pairs": int(len(u)), "warmer_shorter_by_raw_sign": ws, "warmer_longer_by_raw_sign": wl,
            "share": ws / (ws + wl) if ws + wl else np.nan,
            "p_two_sided": float(SS.binomtest(ws, ws + wl, 0.5).pvalue) if ws + wl else np.nan,
            "slope_NCA_and_NMC_pairs": float(np.sum(x * y) / np.sum(x * x)), "pairs_with_prediction": int(len(v)),
            "signs_agree_pairs_with_prediction": int((np.sign(x) == np.sign(y)).sum())}


def two_sided_v41(work, cells):
    out = {}
    for cid in ("SNL_18650_NMC_25C_0-100_0.5/0.5C_a", "SNL_18650_NMC_25C_20-80_0.5/0.5C_a", "SNL_18650_NMC_25C_20-80_0.5/0.5C_b",
                "SNL_18650_NMC_25C_20-80_0.5/0.5C_c"):
        m = cells.loc[cid]
        st = pd.read_csv(os.path.join(work, "steps3", cid.replace("/", "-") + "_steps.csv"))
        st["ah"] = st.ah_counter
        every = m.c_discharge == 0.5 and m.soc_min == 0 and m.soc_max == 100
        st2 = st.copy()
        low = (st2.kind == "discharge") & (st2.v_last < 2.0 - 0.02)
        st2.loc[low, "v_last"] = 9.9                          # outside the two-sided band: not a check
        d = LT.discharges(st2, m.rated_Ah, 2.0, 0.02, 0.05, every)
        ch = LT.checks_table(d)
        keep = ch.q_Ah >= 0.5 * ch.q_Ah.iloc[0]
        life = LT.check_lifetime(ch[keep].reset_index(drop=True), 0.90)
        out[cid] = {"discharges_ending_below_1.98_V": int(low.sum()), "registered_L_efc": float(m.L_efc), "two_sided_L_efc": life["L_efc"]}
    return out


def v1_notes(cells):
    ids = ["SNL_18650_LFP_15C_0-100_0.5/2C_a", "SNL_18650_LFP_15C_0-100_0.5/2C_b", "SNL_18650_LFP_25C_0-100_0.5/0.5C_a"]
    return {cid: {k: float(cells.loc[cid, k]) for k in ("L_efc", "L_sens_efc", "resolution_efc", "q_ref_sens_Ah", "q_first_Ah")} for cid in ids}


def main(work, res):
    cells = pd.read_csv(os.path.join(res, "stage3_cells.csv"))
    ci = cells.set_index("cell_id")
    pairs = pd.read_csv(os.path.join(res, "stage3_pairs.csv"))
    checks = pd.read_csv(os.path.join(res, "stage3_checks.csv"))
    o = pd.read_csv(os.path.join(res, "stage3_pairs_oriented_measure_a.csv"))
    lev = json.load(open(os.path.join(res, "stage3_lever.json")))
    sc_reg = {"P4_inputs": {ch: lev["registered_1C"][ch]["s_25_35"]["s"] for ch in ("NCA", "NMC")}}
    lows, relife = low_checks(checks, ci)
    out = {"note": "Post hoc, after the registered results were committed at 9c6bb84; changes no registered number or verdict.",
           "P4_readings": p4_readings(cells, sc_reg),
           "floor_variants": floor_variants(pairs),
           "low_checks": lows, "lifetimes_without_low_checks": relife,
           "NMC_20_80_0.5C_without_low_checks": nmc_partial_floor(pairs, relife),
           "thermal_share_raw_sign": raw_sign_share(o),
           "resolution_over_L": {ch: {"median": float(x.median()), "min": float(x.min()), "max": float(x.max())}
                                 for ch, x in (cells.resolution_efc / cells.L_efc).groupby(cells.chemistry)},
           "check_rule_two_sided": two_sided_v41(work, ci),
           "v1_cells": v1_notes(ci)}
    json.dump(out, open(os.path.join(res, "stage3_posthoc.json"), "w"), indent=1, default=S3.jsonable)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--results", default="results")
    a = ap.parse_args()
    r = main(a.work, a.results)
    print(json.dumps({k: r[k] for k in ("P4_readings", "NMC_20_80_0.5C_without_low_checks", "thermal_share_raw_sign", "check_rule_two_sided")},
                     indent=1, default=S3.jsonable))
