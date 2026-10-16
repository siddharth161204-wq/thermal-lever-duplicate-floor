#!/usr/bin/env python3
"""Stage 2 checks.

V1. For the five cells fixed in DECISIONS.md, the check-based lifetime against the sensitivity lifetime from the
    per-cycle capacity at the cycling rate: they agree when the difference is within the check-based resolution.
V2. Measure (a) and measure (b) of three cells recomputed from the raw time series with separate code: the csv
    module, no pandas, none of steps.py, lifetimes.py or process_cells.py; labels, steps, counter ampere-hours,
    checks, the lifetime and the units are all rebuilt here. One worked unit per cell is printed for a hand check.
V3. Each law's secant at Stage 1's conditions by the Stage 2 code against Stage 1's registered values.

Usage:
  python3 src/stage2/verify.py v1 --cells results/stage2_cells.csv --ids ID1,ID2,... --out verification/stage2_v1.json
  python3 src/stage2/verify.py v2 --timeseries FILE --cell-id ID --params params.json --cells results/stage2_cells.csv
  python3 src/stage2/verify.py v3 --out verification/stage2_v3.json
"""
import argparse, csv, json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))


def v1(cells_csv, ids, out):
    import pandas as pd
    c = pd.read_csv(cells_csv).set_index("cell_id")
    rows = []
    for i in ids:
        r = c.loc[i]
        d = abs(r.L_efc - r.L_sens_efc)
        rows.append({"cell_id": i, "L_check_efc": float(r.L_efc), "L_sens_efc": float(r.L_sens_efc),
                     "resolution_efc": float(r.resolution_efc), "difference_efc": float(d),
                     "agree_within_resolution": bool(d <= r.resolution_efc)})
    res = {"check": "V1", "cells": rows, "all_agree": all(x["agree_within_resolution"] for x in rows)}
    json.dump(res, open(out, "w"), indent=1)
    return res


def _runs_from_csv(path, rated, i_rest, sign, cols):
    """Rebuild steps from the raw file with the csv module (independent of steps.py)."""
    runs, cur, prev = [], None, None
    with open(path, newline="") as f:
        rd = csv.reader(f)
        head = next(rd)
        ix = {k: head.index(v) for k, v in cols.items()}
        for row in rd:
            try:
                t = float(row[ix["time"]]); I = float(row[ix["current"]])
            except (ValueError, IndexError):
                continue
            V = float(row[ix["voltage"]]); qd = float(row[ix["q_dis"]]); T = float(row[ix["t_cell"]])
            Ic = -sign * I
            lab = 1 if Ic >= i_rest else (-1 if Ic <= -i_rest else 0)
            if cur is None or lab != cur["lab"]:
                if cur is not None:
                    runs.append(cur)
                cur = {"lab": lab, "t0": t, "t1": t, "intT": 0.0, "intI": 0.0, "T0": T, "Tlast": T, "Vlast": V,
                       "qd_before": prev["qd"] if prev else float("nan"), "qd_last": qd, "pT": T, "pI": abs(I), "pt": t}
            else:
                dt = t - cur["pt"]
                cur["intT"] += 0.5 * (cur["pT"] + T) * dt
                cur["intI"] += 0.5 * (cur["pI"] + abs(I)) * dt
                cur.update({"t1": t, "Tlast": T, "Vlast": V, "qd_last": qd, "pT": T, "pI": abs(I), "pt": t})
            prev = {"qd": qd}
        if cur is not None:
            runs.append(cur)
    for r in runs:
        dur = r["t1"] - r["t0"]
        r["dur"] = dur
        r["Tmean"] = r["intT"] / dur if dur > 0 else r["T0"]
        r["Imean"] = r["intI"] / dur if dur > 0 else r["pI"]
        if r["lab"] == -1:
            b = r["qd_before"]
            r["ah"] = r["qd_last"] if (not math.isfinite(b) or b > r["qd_last"] + 1e-9) else r["qd_last"] - b
    return runs


def v2(ts_path, cell_id, params, cells_csv):
    sys.path.insert(0, HERE)
    import sandia_io as IO               # only for the cell's conditions and the column names
    meta = IO.parse_cell_id(cell_id)
    rated, chem = meta["rated_Ah"], meta["chemistry"]
    cols = {k: v for k, v in IO.match_columns(IO.header(ts_path)).items() if k in ("time", "current", "voltage", "q_dis", "t_cell")}
    runs = _runs_from_csv(ts_path, rated, params["i_rest_frac"] * rated, params["discharge_sign"], cols)
    every = meta["c_discharge"] == 0.5 and meta["soc_min"] == 0 and meta["soc_max"] == 100
    cum, k, prev_c = 0.0, -1, False
    checks = []
    for r in runs:
        if r["lab"] != -1:
            continue
        r["efc_before"] = cum / rated
        cum += r["ah"]
        crate = r["Imean"] / rated
        r["check_like"] = (r["Vlast"] <= params["v_cut_low"][chem] + params["tol_v"]) and abs(crate - 0.5) <= params["tol_c"]
        if r["check_like"] and (every or not prev_c):
            k += 1
            checks.append({"efc": r["efc_before"], "q": []})
        if r["check_like"]:
            checks[-1]["q"].append(r["ah"])
        r["check"] = r["check_like"]
        prev_c = r["check_like"]
    q = [sum(c["q"]) / len(c["q"]) for c in checks]
    efc = [c["efc"] for c in checks]
    thr = params["thresholds"][chem] * q[0]
    L = float("nan")
    for i in range(1, len(q)):
        if q[i] < thr:
            if q[i - 1] >= thr:
                L = efc[i - 1] + (q[i - 1] - thr) / (q[i - 1] - q[i]) * (efc[i] - efc[i - 1])
            break
    limit = L if math.isfinite(L) else efc[-1]
    units = []
    for i, r in enumerate(runs):
        if r["lab"] != -1 or (r["check"] and not every):
            continue
        if r["efc_before"] > limit:
            break
        j = i - 1
        while j >= 0 and runs[j]["lab"] == 0:
            j -= 1
        if j < 0 or runs[j]["lab"] != 1:
            continue
        rest = list(range(j + 1, i))
        jj = j - 1
        while jj >= 0 and runs[jj]["lab"] == 0:
            rest.append(jj); jj -= 1
        jj = i + 1
        while jj < len(runs) and runs[jj]["lab"] == 0:
            rest.append(jj); jj += 1
        wc = runs[j]["dur"] + runs[i]["dur"]
        wr = sum(runs[x]["dur"] for x in rest)
        if not rest or wr <= 0 or wc <= 0:
            continue
        tc = (runs[j]["Tmean"] * runs[j]["dur"] + runs[i]["Tmean"] * runs[i]["dur"]) / wc
        if params.get("rest_baseline", "mean") == "end":
            tr = sum(runs[x]["Tlast"] for x in rest) / len(rest)
        else:
            tr = sum(runs[x]["Tmean"] * runs[x]["dur"] for x in rest) / wr
        units.append({"charge_T_mean": runs[j]["Tmean"], "charge_s": runs[j]["dur"], "discharge_T_mean": runs[i]["Tmean"],
                      "discharge_s": runs[i]["dur"], "rests": [(runs[x]["Tmean"], runs[x]["Tlast"], runs[x]["dur"]) for x in rest],
                      "T_cycling": tc, "T_rest": tr, "rise": tc - tr})
    rise = sum(u["rise"] for u in units) / len(units)
    rest_T = sum(u["T_rest"] for u in units) / len(units)
    import pandas as pd
    c = pd.read_csv(cells_csv).set_index("cell_id").loc[cell_id]
    return {"check": "V2", "cell_id": cell_id, "L_efc_by_hand": L, "L_efc_pipeline": float(c.L_efc),
            "rise_K_by_hand": rise, "rise_K_pipeline": float(c.rise_K), "rest_T_C_by_hand": rest_T,
            "rest_T_C_pipeline": float(c.rest_T_C), "n_units_by_hand": len(units), "n_units_pipeline": int(c.n_units),
            "agree": bool(abs(rise - c.rise_K) < 1e-6 and abs(rest_T - c.rest_T_C) < 1e-6 and len(units) == int(c.n_units)
                          and (abs(L - c.L_efc) < 1e-6 if math.isfinite(L) else bool(c.censored))),
            "worked_unit": units[0] if units else None}


def v3(out):
    sys.path.insert(0, HERE)
    import laws as LW
    res = LW.v3_against_stage1(json.load(open(os.path.join(ROOT, "results", "lever.json"))))
    res["check"] = "V3"
    res["pass"] = res["max_rel_diff_secant"] < 1e-9 and res["max_rel_diff_lifetime"] < 1e-12
    json.dump(res, open(out, "w"), indent=1)
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("which", choices=["v1", "v2", "v3"])
    ap.add_argument("--cells")
    ap.add_argument("--ids")
    ap.add_argument("--timeseries")
    ap.add_argument("--cell-id")
    ap.add_argument("--params")
    ap.add_argument("--out")
    a = ap.parse_args()
    if a.which == "v1":
        print(json.dumps(v1(a.cells, a.ids.split(","), a.out), indent=1))
    elif a.which == "v2":
        r = v2(a.timeseries, a.cell_id, json.load(open(a.params)), a.cells)
        print(json.dumps(r, indent=1))
        if a.out:
            json.dump(r, open(a.out, "w"), indent=1)
    else:
        print(json.dumps(v3(a.out), indent=1))
