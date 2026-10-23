#!/usr/bin/env python3
"""Stage 3 on the Sandia cells of Preger et al. 2020 (DECISIONS.md D38 to D50).

1. A temperature pass over the rows of every time series in cycle order (src/stage3/rows.py), labelled exactly as the
   structure pass of src/stage2/steps.py and checked against the step tables of src/stage3/prepare.py step by step:
   per step the mean cell and environment temperatures over the samples that hold a value, and the last such sample
   (D43, D50 (3)).
2. Per cell, with the ampere-hours of D50 (2): capacity checks with the validity rule of D41, the lifetime at the
   chemistry's threshold (D42), the sensitivity lifetime, and the temperature measures of D43 over the units up to
   the lifetime.
3. The levers (D44), the floor (D45), the thermal share (D46), the module (D47) and the scorecard (D48).

The Stage 2 reader and pipeline are imported unchanged; nothing in src/stage2 is edited.

Usage: python3 -I src/stage3/stage3.py --raw RAW_DIR --work WORK_DIR --out results
(RAW_DIR holds the unpacked *_timeseries.csv files; WORK_DIR the steps3/ tables written by prepare.py)
"""
import argparse, glob, json, os, sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "src", "stage2"))
import rows as RW             # noqa: E402
import sandia_io as IO        # noqa: E402
import lifetimes as LT        # noqa: E402
import laws as LW             # noqa: E402
import analysis as AN         # noqa: E402

CHEMS = ("LFP", "NCA", "NMC")
V_CUT = {"LFP": 2.0, "NCA": 2.5, "NMC": 2.0}           # D41
THRESH = {"LFP": 0.95, "NCA": 0.90, "NMC": 0.90}       # D42
TOL_V, TOL_C, I_REST = 0.02, 0.05, 0.005               # D41, D40
SEED, N_BOOT = 20261024, 10000                         # D44
OFFER_K = 0.391279                                      # Stage 1's largest module offer at the Geslin rates (D47)
REGISTERED, SECONDARY = 1.0, 2.0                        # D39


def jsonable(o):
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


# ---------------------------------------------------------------------------------------------------------------
# 1. Temperature pass
# ---------------------------------------------------------------------------------------------------------------
def temperature_pass(arr, rated, st, name):
    """Per-step temperature statistics over valid samples, on the rows in cycle order (D50 (1), (3)).

    The samples are labelled as in steps.segment (rest when |I| < 0.005C, charge or discharge by the sign of the
    current, discharge negative) and the runs must coincide with the step table in number, kind and start time;
    otherwise the computation stops."""
    t = arr["time"]
    temps = {k: arr[k] for k in ("t_cell", "t_env")}
    lab = RW.labels(arr["current"], rated, I_REST)
    starts, ends = RW.runs(lab)
    kinds = pd.Series(lab[starts]).map({1: "charge", -1: "discharge", 0: "rest"}).to_numpy()
    if not (starts.size == len(st) and np.array_equal(kinds, st.kind.to_numpy())
            and np.allclose(t[starts], st.t_start_s.to_numpy(), rtol=0, atol=1e-6)):
        raise SystemExit(f"temperature pass does not reproduce the step table: {name}")
    rid = np.repeat(np.arange(starts.size), ends - starts)
    dt = np.diff(t)
    out = {}
    for k, v in temps.items():
        val = np.isfinite(v)
        pair = (rid[:-1] == rid[1:]) & val[:-1] & val[1:]
        num = np.bincount(rid[:-1][pair], weights=(0.5 * (v[:-1] + v[1:]) * dt)[pair], minlength=starts.size)
        den = np.bincount(rid[:-1][pair], weights=dt[pair], minlength=starts.size)
        cnt = np.bincount(rid[val], minlength=starts.size)
        smp = np.bincount(rid[val], weights=v[val], minlength=starts.size)
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = np.where(den > 0, num / np.where(den > 0, den, 1.0), np.where(cnt > 0, smp / np.maximum(cnt, 1), np.nan))
        idx = np.where(val, np.arange(v.size), -1)
        last = np.maximum.reduceat(idx, starts)
        out[k + "_mean"] = mean
        out[k + "_last"] = np.where(last >= 0, v[np.maximum(last, 0)], np.nan)
        out[k + "_valid"] = cnt
    out["n_samples"] = ends - starts
    return out


# ---------------------------------------------------------------------------------------------------------------
# 2. Per cell
# ---------------------------------------------------------------------------------------------------------------
def units(st, d, limit):
    """D43 with D50 (3) to (6): per cycling unit (a cycling discharge that is not part of a check and the charge step
    before it) the swing, the adjacent rest records and the literal measure; the cell's values over its units."""
    kinds = st.kind.to_numpy()
    tm, tl, te = st.t_cell_mean.to_numpy(), st.t_cell_last.to_numpy(), st.t_env_mean.to_numpy()
    dur = st.duration_s.to_numpy()
    pos = {s: i for i, s in enumerate(st.step.to_numpy())}
    rows, skipped = [], 0
    for _, r in d[d.check_id < 0].iterrows():
        if r.efc_before > limit:
            break
        i_d = pos[r.step]
        j = i_d - 1
        while j >= 0 and kinds[j] == "rest":
            j -= 1
        if j < 0 or kinds[j] != "charge":
            continue
        if not (np.isfinite(tm[j]) and np.isfinite(tm[i_d])):
            skipped += 1
            continue
        rest = list(range(j + 1, i_d))
        k = j - 1
        while k >= 0 and kinds[k] == "rest":
            rest.append(k); k -= 1
        k = i_d + 1
        while k < len(kinds) and kinds[k] == "rest":
            rest.append(k); k += 1
        rr = [x for x in rest if np.isfinite(tm[x])]
        w_c, w_d = dur[j], dur[i_d]
        unit_cyc = (tm[j] * w_c + tm[i_d] * w_d) / (w_c + w_d) if w_c + w_d > 0 else 0.5 * (tm[j] + tm[i_d])
        rows.append({"efc": r.efc_before, "swing": tm[i_d] - tm[j],
                     "num_T": tm[j] * w_c + tm[i_d] * w_d, "den_T": w_c + w_d,
                     "num_env": (tm[j] - te[j]) * w_c + (tm[i_d] - te[i_d]) * w_d if np.isfinite(te[j]) and np.isfinite(te[i_d]) else np.nan,
                     "rest_T": float(np.mean(tm[rr])) if rr else np.nan,
                     "literal": unit_cyc - float(np.mean(tl[rr])) if rr else np.nan,
                     "n_rest": len(rest)})
    if not rows:
        return {"n_units": 0, "units_skipped_missing_T": skipped}
    u = pd.DataFrame(rows)
    env_ok = u.num_env.notna()
    return {"n_units": len(u), "units_skipped_missing_T": skipped, "swing_K": float(u.swing.mean()),
            "swing_sd_K": float(u.swing.std(ddof=1)) if len(u) > 1 else np.nan,
            "T_cyc_C": float(u.num_T.sum() / u.den_T.sum()) if u.den_T.sum() > 0 else float(np.mean(u.num_T)),
            "rest_T_C": float(u.rest_T.mean()), "units_with_rest_record": int(u.rest_T.notna().sum()),
            "literal_rise_K": float(u.literal.mean()),
            "skin_minus_env_K": float(u.num_env[env_ok].sum() / u.den_T[env_ok].sum()) if env_ok.any() and u.den_T[env_ok].sum() > 0 else np.nan,
            "efc_last_unit": float(u.efc.iloc[-1])}


def cell_record(meta, st):
    chem = meta["chemistry"]
    st["ah"] = st["ah_counter"]                                              # D40 with D50 (2)
    every = meta["c_discharge"] == 0.5 and meta["soc_min"] == 0 and meta["soc_max"] == 100
    d = LT.discharges(st, meta["rated_Ah"], V_CUT[chem], TOL_V, TOL_C, every)
    ch = LT.checks_table(d)
    if ch.empty:
        raise SystemExit(f"no capacity check found: {meta['cell_id']}")
    keep = (ch.q_Ah >= 0.5 * ch.q_Ah.iloc[0]).to_numpy()                       # D41, D50 (7)
    life = LT.check_lifetime(ch[keep].reset_index(drop=True), THRESH[chem])   # D42
    rec = dict(meta, soh_threshold=THRESH[chem], checks_found=int(len(ch)), checks_dropped=int((~keep).sum()),
               **{k: v for k, v in life.items() if k != "bracket"})
    rec["bracket_efc"] = json.dumps(life.get("bracket"))
    if meta["soc_min"] == 0 and meta["soc_max"] == 100:
        rec.update(LT.sensitivity_lifetime(d, THRESH[chem], meta["c_discharge"], TOL_C))
    limit = life["L_efc"] if not life["censored"] else np.inf               # D43, D50 (6)
    rec.update(units(st, d, limit))
    # measure (a) and (b) as analysis.py reads them; cells discharged at 0.5C carry no swing (D43, D50 (5))
    used = meta["c_discharge"] > 0.5
    rec["rise_K"] = rec.get("swing_K", np.nan) if used else np.nan
    rec["measures_used"] = bool(used and np.isfinite(rec.get("swing_K", np.nan)))
    rec["efc_total"] = float(d["cum_ah_before"].iloc[-1] + d["ah"].iloc[-1]) / meta["rated_Ah"] if len(d) else np.nan
    return rec, ch.assign(cell_id=meta["cell_id"], kept=keep)


# ---------------------------------------------------------------------------------------------------------------
# 3. Lever, floor, thermal share, module, scorecard
# ---------------------------------------------------------------------------------------------------------------
def temperature_series(cells, chem, c_dis):
    return cells[(cells.chemistry == chem) & (cells.soc_min == 0) & (cells.soc_max == 100) & (cells.c_charge == 0.5)
                 & (cells.c_discharge == c_dis)]


def lever(cells, chem, c_dis, rng):
    """D44: secants between the logged group temperatures, with the bootstrap over cells within each ambient."""
    ts = temperature_series(cells, chem, c_dis)
    g = {T: ts[ts.ambient_C == T].sort_values("cell_id") for T in (15.0, 25.0, 35.0)}
    out = {"chemistry": chem, "c_discharge": c_dis, "soh_threshold": THRESH[chem],
           "n_cells": {f"{int(T)}C": int(len(x)) for T, x in g.items()},
           "censored_cells": sorted(ts.loc[ts.censored.astype(bool), "cell_id"].tolist()),
           "cells": {f"{int(T)}C": x[["cell_id", "L_efc", "resolution_efc", "T_cyc_C", "censored"]].to_dict("records") for T, x in g.items()},
           "mean_lnL": {}, "T_group_C": {}}
    for T, x in g.items():
        if len(x) and not x.censored.astype(bool).any():
            out["mean_lnL"][f"{int(T)}C"] = float(np.log(x.L_efc.to_numpy()).mean())
            out["T_group_C"][f"{int(T)}C"] = float(x.T_cyc_C.mean())

    def sec(a, b, lo_set, hi_set):
        if len(a) == 0 or len(b) == 0 or a.censored.astype(bool).any() or b.censored.astype(bool).any():
            return None
        la, lb = np.log(a.L_efc.to_numpy()), np.log(b.L_efc.to_numpy())
        ta, tb = a.T_cyc_C.to_numpy(), b.T_cyc_C.to_numpy()
        if not (np.all(np.isfinite(ta)) and np.all(np.isfinite(tb))):
            return None
        r = {"s": float(-(lb.mean() - la.mean()) / (tb.mean() - ta.mean())), "T_lo_C": float(ta.mean()), "T_hi_C": float(tb.mean()),
             "dT_K": float(tb.mean() - ta.mean()), "mean_lnL_lo": float(la.mean()), "mean_lnL_hi": float(lb.mean()),
             "s_setpoint": float(-(lb.mean() - la.mean()) / (hi_set - lo_set)), "n_lo": int(la.size), "n_hi": int(lb.size)}
        bs = np.empty(N_BOOT)
        for k in range(N_BOOT):
            i, j = rng.integers(0, la.size, la.size), rng.integers(0, lb.size, lb.size)
            bs[k] = -(lb[j].mean() - la[i].mean()) / (tb[j].mean() - ta[i].mean())
        r["ci95"] = [float(x) for x in np.percentile(bs, [2.5, 97.5])]
        r["laws"] = LW.law_secants(r["T_lo_C"], r["T_hi_C"], c_dis, THRESH[chem])
        r["laws_setpoint"] = LW.law_secants(lo_set, hi_set, c_dis, THRESH[chem])
        return r

    out["s_25_35"] = sec(g[25.0], g[35.0], 25.0, 35.0)
    out["s_15_25"] = sec(g[15.0], g[25.0], 15.0, 25.0)
    out["law_curves"] = LW.law_curves([float(x) for x in np.arange(15.0, 40.01, 0.5)], c_dis, THRESH[chem])
    return out


def kit_values():
    """D44: the KIT levers beside NMC, read unchanged from results/stage2_lever.csv."""
    k = pd.read_csv(os.path.join(ROOT, "results", "stage2_lever.csv"))
    out = []
    for p in ("0-100 0.33/1.00C", "0-100 1.00/1.00C"):
        for lo in (25, 10):
            r = k[(k.protocol == p) & (k.T_lo_set == lo)].iloc[0]
            out.append({"protocol": p, "secant": f"{int(r.T_lo_set)} to {int(r.T_hi_set)} C", "s_per_K": float(r.s_per_K),
                        "T_lo_C": float(r.T_lo_C), "T_hi_C": float(r.T_hi_C)})
    return out


def kit_floor():
    f = json.load(open(os.path.join(ROOT, "results", "stage2_lever.json")))["floor"]
    return {"mean": f["mean"], "median": f["median"], "p95": f["p95"], "pairs": f["pairs"]}


def not_scored(why):
    return {"scored": False, "held": None, "reason": why}


def score(lev, flo, share, mod):
    """D48."""
    reg = {c: lev[c]["s_25_35"] for c in CHEMS}
    low = {c: lev[c]["s_15_25"] for c in CHEMS}
    sc = {}
    x = reg["NMC"]
    sc["P1"] = ({"scored": True, "held": bool(-0.035 <= x["s"] <= 0.005), "s_per_K": x["s"], "ci95": x["ci95"],
                 "band": [-0.035, 0.005], "central": -0.015} if x else not_scored("NMC 25 to 35 C lever not determinable"))
    y = low["NMC"]
    sc["P2"] = ({"scored": True, "held": bool(y["s"] < 0 and abs(y["s"]) > abs(x["s"])), "s_15_25_per_K": y["s"],
                 "s_25_35_per_K": x["s"], "ci95_15_25": y["ci95"]} if x and y else not_scored("an NMC lever not determinable"))
    z = reg["NCA"]
    sc["P3"] = ({"scored": True, "held": bool(0.0 <= z["s"] <= 0.05), "s_per_K": z["s"], "ci95": z["ci95"],
                 "band": [0.0, 0.05], "central": 0.02} if z else not_scored("NCA 25 to 35 C lever not determinable"))
    if all(reg.values()):
        l2 = reg["LFP"]["laws"]["law_2_schimpe2018"]
        ratio = reg["LFP"]["s"] / l2 if l2 else np.nan
        parts = {"exceeds_NCA": bool(reg["LFP"]["s"] > reg["NCA"]["s"]), "exceeds_NMC": bool(reg["LFP"]["s"] > reg["NMC"]["s"]),
                 "sign_of_law_2": bool(np.sign(reg["LFP"]["s"]) == np.sign(l2)), "ratio_within_factor_3": bool(1 / 3 <= ratio <= 3)}
        sc["P4"] = {"scored": True, "held": all(parts.values()), "parts": parts, "s_per_K": {c: reg[c]["s"] for c in CHEMS},
                    "law_2_per_K": l2, "ratio_s_over_law_2": ratio,
                    "secant_C": [reg["LFP"]["T_lo_C"], reg["LFP"]["T_hi_C"]]}
    else:
        sc["P4"] = not_scored("a 25 to 35 C lever not determinable")
    p5 = {c: {"mean": flo[c]["mean"], "pairs": flo[c]["pairs_with_difference"], "groups": flo[c]["groups"],
              "within": bool(0.015 <= flo[c]["mean"] <= 0.10)} for c in CHEMS if flo[c]["pairs_with_difference"] >= 5}
    sc["P5"] = ({"scored": True, "held": all(v["within"] for v in p5.values()), "chemistries": p5, "band": [0.015, 0.10],
                 "central": 0.04} if p5 else not_scored("no chemistry with five pairs entering the floor"))
    lo1 = share["1C or below"]
    sc["P6"] = ({"scored": True, "held": bool(lo1["dT_below_1K_fraction"] >= 0.80), "fraction_below_1K": lo1["dT_below_1K_fraction"],
                 "median_dT_K": lo1["dT_K_median"], "central_K": 0.3, "pairs": lo1["pairs"]}
                if lo1["pairs"] else not_scored("no pair at 1C or below with measure (a) for both cells"))
    a = share["all"]
    n_a = a["warmer_shorter"] + a["warmer_longer"]
    sc["P7"] = ({"scored": True, "held": bool(0.35 <= a["share"] <= 0.65), "share": a["share"], "warmer_shorter": a["warmer_shorter"],
                 "warmer_longer": a["warmer_longer"], "ties": a["ties"], "excluded_censored": a["excluded_censored"],
                 "p_one_sided": a["sign_test_p_one_sided"], "p_two_sided": a["sign_test_p_two_sided"], "band": [0.35, 0.65]}
                if n_a else not_scored("no pair enters the sign test"))
    hi = share["2C and 3C"]
    n_hi, n_lo = hi["warmer_shorter"] + hi["warmer_longer"], lo1["warmer_shorter"] + lo1["warmer_longer"]
    sc["P8"] = ({"scored": True, "held": bool(hi["share"] > lo1["share"]), "share_2C_3C": hi["share"], "share_1C_or_below": lo1["share"]}
                if n_hi >= 5 and n_lo >= 5 else not_scored("a rate group holds fewer than five pairs entering the sign test"))
    sc["P8"]["pairs_entering_sign_test"] = {"2C and 3C": n_hi, "1C or below": n_lo}
    sc["P9"] = ({"scored": True, "held": bool(mod["implied_divergence"] < 0.015), "implied_divergence": mod["implied_divergence"],
                 "implied_divergence_ci95": mod["implied_divergence_ci95"], "limit": 0.015, "offer_K": mod["module_largest_dT_K"]}
                if mod else not_scored("NCA 25 to 35 C lever not determinable"))
    sc["P10"] = ({"scored": True, "held": all(low[c]["s"] < 0 for c in CHEMS), "s_15_25_per_K": {c: low[c]["s"] for c in CHEMS},
                  "ci95": {c: low[c]["ci95"] for c in CHEMS}}
                 if all(low.values()) else not_scored("a 15 to 25 C lever not determinable"))
    sc["summary"] = {"held": sorted((k for k, v in sc.items() if v.get("held") is True), key=lambda s: int(s[1:])),
                     "did_not_hold": sorted((k for k, v in sc.items() if v.get("held") is False), key=lambda s: int(s[1:])),
                     "not_scored": sorted((k for k, v in sc.items() if v.get("scored") is False), key=lambda s: int(s[1:]))}
    return sc


def main(raw, work, out, expect_cells=86):
    rng = np.random.default_rng(SEED)
    recs, checks, tstats = [], [], {}
    paths = sorted(glob.glob(os.path.join(raw, "**", "*_timeseries.csv"), recursive=True))
    keys = [IO.parse_cell_id(os.path.basename(p))["key"] for p in paths]
    if len(set(keys)) != len(keys):
        raise SystemExit("a cell appears twice under the raw folder")
    for p in paths:
        meta = IO.parse_cell_id(os.path.basename(p))
        st = pd.read_csv(os.path.join(work, "steps3", meta["key"] + "_steps.csv"))
        _, arr, _ = RW.load_rows(p)
        tp = temperature_pass(arr, meta["rated_Ah"], st, meta["key"])
        for k in ("t_cell_mean", "t_cell_last", "t_env_mean"):
            st[k] = tp[k]
        tstats[meta["key"]] = {"steps": int(len(st)), "steps_without_cell_temperature": int((tp["t_cell_valid"] == 0).sum()),
                               "samples_without_cell_temperature": int((tp["n_samples"] - tp["t_cell_valid"]).sum()),
                               "samples": int(tp["n_samples"].sum())}
        r, ch = cell_record(meta, st)
        r.update({k: tstats[meta["key"]][k] for k in ("samples", "samples_without_cell_temperature")})
        recs.append(r); checks.append(ch)
    cells = pd.DataFrame(recs).sort_values("cell_id").reset_index(drop=True)
    if expect_cells and len(cells) != expect_cells:
        raise SystemExit(f"expected {expect_cells} cells, found {len(cells)}")
    floor = json.load(open(os.path.join(ROOT, "inputs", "floor.json")))
    lev = {c: lever(cells, c, REGISTERED, rng) for c in CHEMS}
    lev2 = {c: lever(cells, c, SECONDARY, rng) for c in CHEMS}
    s_by = {c: lev[c]["s_25_35"]["s"] for c in CHEMS if lev[c]["s_25_35"]}
    pairs = AN.pair_rows(cells, s_by)                                   # D45: every pair
    flo = AN.floor_stats(pairs, floor)
    pts = pairs[pairs.c_discharge > 0.5].reset_index(drop=True)         # D43, D50 (5): the thermal share's pairs
    share_a, o_a = AN.thermal_share(pts, "rise", rng)                  # D46, measure (a)
    share_b, o_b = AN.thermal_share(pts, "rest", rng)                  # D46, measure (b)
    nca = lev["NCA"]["s_25_35"]
    mod = AN.module_application({"s_25_35_per_K": nca["s"], "s_25_35_ci95": nca["ci95"]}, os.path.join(ROOT, "results"),
                                pts, floor) if nca else {}
    if mod and abs(mod["module_largest_dT_K"] - OFFER_K) > 5e-7:
        raise SystemExit(f"module offer {mod['module_largest_dT_K']} differs from {OFFER_K}")
    sc = score(lev, flo, share_a, mod)
    kf = kit_floor()
    os.makedirs(out, exist_ok=True)
    cells.to_csv(os.path.join(out, "stage3_cells.csv"), index=False)
    pd.concat(checks, ignore_index=True).to_csv(os.path.join(out, "stage3_checks.csv"), index=False)
    rows = []
    for tag, L in (("registered 1C", lev), ("secondary 2C", lev2)):
        for c in CHEMS:
            for key, name in (("s_25_35", "25 to 35 C"), ("s_15_25", "15 to 25 C")):
                x = L[c][key]
                row = {"series": tag, "chemistry": c, "secant": name, "soh_threshold": THRESH[c], "determinable": bool(x)}
                if x:
                    row.update({"s_per_K": x["s"], "ci95_lo": x["ci95"][0], "ci95_hi": x["ci95"][1], "T_lo_C": x["T_lo_C"],
                                "T_hi_C": x["T_hi_C"], "dT_K": x["dT_K"], "n_lo": x["n_lo"], "n_hi": x["n_hi"],
                                "s_setpoint_per_K": x["s_setpoint"], "law1_per_K": x["laws"]["law_1_kucinskis2022"],
                                "law2_per_K": x["laws"]["law_2_schimpe2018"]})
                rows.append(row)
    pd.DataFrame(rows).to_csv(os.path.join(out, "stage3_lever.csv"), index=False)
    fl = [{"source": f"Sandia {c}", "chemistry": c, **{k: v for k, v in flo[c].items()}} for c in CHEMS]
    fl += [{"source": "Geslin et al. 2025 (inputs/floor.json)", "chemistry": "NCA (SiOx-graphite)", "mean": floor["mean"],
            "median": floor["median"], "p95": floor["p95"], "pairs_with_difference": floor["n_pairs"]},
           {"source": "KIT LG HG2, Stage 2 (results/stage2_lever.json)", "chemistry": "NMC", "mean": kf["mean"],
            "median": kf["median"], "p95": kf["p95"], "pairs_with_difference": kf["pairs"]}]
    pd.DataFrame(fl).to_csv(os.path.join(out, "stage3_floor.csv"), index=False)
    pairs.to_csv(os.path.join(out, "stage3_pairs.csv"), index=False)
    o_a.to_csv(os.path.join(out, "stage3_pairs_oriented_measure_a.csv"), index=False)
    o_b.to_csv(os.path.join(out, "stage3_pairs_oriented_measure_b.csv"), index=False)
    res = {"thresholds": THRESH, "seed": SEED, "n_boot": N_BOOT, "registered_rate_C": REGISTERED,
           "registered_1C": lev, "secondary_2C": lev2, "kit_levers": kit_values(), "kit_floor": kf,
           "cells": int(len(cells)), "censored": int(cells.censored.astype(bool).sum()),
           "censored_cells": sorted(cells.loc[cells.censored.astype(bool), "cell_id"].tolist()),
           "checks_dropped": int(cells.checks_dropped.sum()),
           "cells_with_dropped_checks": cells.loc[cells.checks_dropped > 0, ["cell_id", "checks_dropped"]].to_dict("records"),
           "temperature_pass": tstats}
    json.dump(res, open(os.path.join(out, "stage3_lever.json"), "w"), indent=1, default=jsonable)
    json.dump({"floor": flo, "geslin_floor": {"mean": floor["mean"], "median": floor["median"], "p95": floor["p95"], "n_pairs": floor["n_pairs"]},
               "kit_floor": kf, "measure_a": share_a, "measure_b": share_b, "module": mod,
               "thermal_share_pairs": int(len(pts)), "pairs_total": int(len(pairs))},
              open(os.path.join(out, "stage3_thermal_share.json"), "w"), indent=1, default=jsonable)
    json.dump(sc, open(os.path.join(out, "stage3_scorecard.json"), "w"), indent=1, default=jsonable)
    return res, sc


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--out", default="results")
    a = ap.parse_args()
    res, sc = main(a.raw, a.work, a.out)
    print(json.dumps({"cells": res["cells"], "censored": res["censored"], "checks_dropped": res["checks_dropped"],
                      "scorecard": sc["summary"]}, default=jsonable))
