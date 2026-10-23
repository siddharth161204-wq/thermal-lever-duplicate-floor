#!/usr/bin/env python3
"""Stage 3 independent check, step 6 (after the freeze): compare a per-cell table of mine with the pipeline's
results/stage3_cells.csv and results/stage3_checks.csv, column by column.

Usage: python3 -I s6_compare_cells.py <my cells json> [<out json>]
"""
import csv, json, math, os, sys

REPO = "/home/claude/thermal-lever-duplicate-floor"

# (my key, pipeline column)
PAIRS = [("L_efc", "L_efc"), ("resolution_efc", "resolution_efc"), ("censored", "censored"),
         ("first_check_Ah", "q_first_Ah"), ("checks", "checks_found"), ("checks_dropped_half", "checks_dropped"),
         ("efc_total", "efc_total"), ("rows_ordered", None), ("measure_a_K", "swing_K"), ("measure_b_C", "rest_T_C"),
         ("T_cyc_C", "T_cyc_C"), ("cell_minus_env_K", "skin_minus_env_K"), ("units_used", "n_units"),
         ("units_in_life_skipped_missing_mean", "units_skipped_missing_T"), ("units_used_b", "units_with_rest_record"),
         ("L_sens_ref_first_cycling", "L_sens_efc"), ("L_sens_ref_first_check", "L_sens_efc")]


def num(x):
    if x in ("", None):
        return float("nan")
    if x in ("True", "False"):
        return x == "True"
    return float(x)


def compare(mine_path):
    mine = {c["cell_id"]: c for c in json.load(open(mine_path))}
    for c in mine.values():
        s = c.get("sensitivity", {})
        c.setdefault("L_sens_ref_first_check", s.get("L_sens_ref_first_check"))
        c.setdefault("L_sens_ref_first_cycling", s.get("L_sens_ref_first_cycling"))
    pipe = {r["cell_id"]: r for r in csv.DictReader(open(os.path.join(REPO, "results", "stage3_cells.csv")))}
    assert set(mine) == set(pipe), (set(mine) ^ set(pipe))
    res = {}
    for mk, pk in PAIRS:
        if pk is None:
            continue
        n = nboth = nmis = nexact = 0
        mx_abs = mx_rel = 0.0
        arg_abs = arg_rel = None
        mism = []
        for cid in sorted(mine):
            a = mine[cid].get(mk)
            b = num(pipe[cid][pk])
            if isinstance(b, bool) or isinstance(a, bool):
                n += 1
                if bool(a) != bool(b):
                    nmis += 1
                    mism.append(cid)
                continue
            a = float("nan") if a is None else float(a)
            if math.isnan(a) and math.isnan(b):
                nboth += 1
                continue
            if math.isnan(a) != math.isnan(b):
                nmis += 1
                mism.append(cid)
                continue
            n += 1
            d = abs(a - b)
            r = d / abs(b) if b != 0 else d
            nexact += d == 0
            if d > mx_abs:
                mx_abs, arg_abs = d, cid
            if r > mx_rel:
                mx_rel, arg_rel = r, cid
        res[f"{mk} vs {pk}"] = {"n_compared": n, "n_nan_both": nboth, "n_nan_mismatch": nmis, "nan_mismatch_cells": mism,
                                "n_exact": nexact, "max_abs": mx_abs, "argmax_abs": arg_abs, "max_rel": mx_rel, "argmax_rel": arg_rel}
    # checks: per cell, the list of (efc, capacity, n) against the pipeline's
    pc = {}
    for r in csv.DictReader(open(os.path.join(REPO, "results", "stage3_checks.csv"))):
        pc.setdefault(r["cell_id"], []).append(r)
    chk = {"cells_same_count": 0, "cells_diff_count": [], "max_abs_efc": 0.0, "max_abs_q": 0.0, "kept_flag_mismatch": 0}
    for cid, rows in pc.items():
        ml = mine[cid]["checks_list"]
        if len(ml) != len(rows):
            chk["cells_diff_count"].append([cid, len(ml), len(rows)])
            continue
        chk["cells_same_count"] += 1
        for a, b in zip(ml, rows):
            chk["max_abs_efc"] = max(chk["max_abs_efc"], abs(a["efc"] - float(b["efc"])))
            chk["max_abs_q"] = max(chk["max_abs_q"], abs(a["capacity_Ah"] - float(b["q_Ah"])))
            chk["kept_flag_mismatch"] += (not a.get("dropped_half", False)) != (b["kept"] == "True")
            if int(b["n_discharges"]) != a["n_discharges"]:
                chk.setdefault("n_discharges_mismatch", []).append([cid, a["efc"], a["n_discharges"], int(b["n_discharges"])])
    res["checks"] = chk
    return res


if __name__ == "__main__":
    r = compare(sys.argv[1])
    for k, v in r.items():
        if k == "checks":
            print("checks:", {kk: (vv if not isinstance(vv, list) or len(vv) < 8 else vv[:8] + ["..."]) for kk, vv in v.items()})
        else:
            print(f"{k:45s} n={v['n_compared']:3d} nan_both={v['n_nan_both']:2d} nan_mis={v['n_nan_mismatch']:2d} exact={v['n_exact']:3d} "
                  f"max_abs={v['max_abs']:.3e} ({(v['argmax_abs'] or '')[10:]}) max_rel={v['max_rel']:.3e} ({(v['argmax_rel'] or '')[10:]})"
                  + (f" mism={v['nan_mismatch_cells']}" if v['nan_mismatch_cells'] else ""))
    if len(sys.argv) > 2:
        json.dump(r, open(sys.argv[2], "w"), indent=1)
