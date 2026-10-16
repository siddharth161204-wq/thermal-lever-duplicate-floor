#!/usr/bin/env python3
"""Stage 2 inventory (D2 of the brief), from the files' structure only, before the predictions.

Per cell: chemistry, ambient, SOC window, charge and discharge rates, replicate, the files and their sizes, the
time series' columns, rows and cycle count, the temperature channels present (non-empty share of each column) with
the sampling interval, the step counts by kind, the capacity checks (count, EFC spacing, first check capacity) and the
lowest check capacity reached as a fraction of the first (used only for the threshold rule, DECISIONS.md), the
discharge-sign and counter checks, and the durations of the rests next to cycling steps. No temperature value is
averaged and no lifetime is computed here.

Usage:
  python3 src/stage2/inventory.py --root DATA_DIR --steps STEPS_DIR --params params.json --out inputs/sandia
"""
import argparse, json, os, sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sandia_io as IO      # noqa: E402
import steps as ST          # noqa: E402
import lifetimes as LT      # noqa: E402
import process_cells as PC  # noqa: E402


def channel_presence(path, chunksize=1_000_000):
    """Share of rows in which each temperature column holds a number. Presence only; no value is summarised."""
    m = IO.match_columns(IO.header(path))
    out = {}
    for key in ("t_cell", "t_env"):
        col = m.get(key)
        if not col:
            out[key] = {"column": None, "non_empty_share": 0.0}
            continue
        n = good = 0
        for ch in pd.read_csv(path, usecols=[col], chunksize=chunksize):
            v = pd.to_numeric(ch[col], errors="coerce")
            n += len(v)
            good += int(v.notna().sum())
        out[key] = {"column": col, "non_empty_share": good / n if n else 0.0}
    return out


def rest_durations(st):
    """Durations of the rest runs that touch a charge or a discharge step (structure)."""
    k = st["kind"].to_numpy()
    d = st["duration_s"].to_numpy()
    touch = np.zeros(len(k), bool)
    for i in range(len(k)):
        if k[i] == "rest" and ((i > 0 and k[i - 1] != "rest") or (i + 1 < len(k) and k[i + 1] != "rest")):
            touch[i] = True
    return d[touch]


def cell_row(c, steps_path, params):
    ts = c["files"].get("timeseries", [None])[0]
    cy = c["files"].get("cycle", [None])[0]
    row = {k: c[k] for k in ("cell_id", "chemistry", "ambient_C", "soc_min", "soc_max", "c_charge", "c_discharge", "replicate", "rated_Ah")}
    row["timeseries_file"] = os.path.basename(ts) if ts else None
    row["cycle_file"] = os.path.basename(cy) if cy else None
    if ts:
        cols = IO.header(ts)
        row["timeseries_columns"] = ";".join(cols)
        row.update({f"{k}_{x}": v for k, d in channel_presence(ts).items() for x, v in d.items()})
        row.update({f"sampling_{k}": v for k, v in ST.sampling_intervals(ts).items()})
    if cy:
        cd, m = IO.read_cycle_data(cy)
        row["cycle_file_rows"] = len(cd)
        row["cycle_file_max_cycle_index"] = float(cd["cycle"].max()) if "cycle" in cd else np.nan
    st = pd.read_csv(steps_path)
    row["timeseries_max_cycle_index"] = float(st["cycle_last"].max())
    for kind in ("charge", "discharge", "rest"):
        row[f"steps_{kind}"] = int((st["kind"] == kind).sum())
    dis = st[st["kind"] == "discharge"]
    row["discharge_counter_rises_share"] = float((dis["q_dis_max"] > dis["q_dis_first"]).mean()) if len(dis) else np.nan
    st["ah_ctr"] = PC.counter_ah(st)
    rel = ((st["ah_ctr"] - st["ah_int"]) / st["ah_ctr"]).where(st["kind"] == "discharge").dropna()
    row["ah_counter_vs_integral_median_rel"] = float(rel.median()) if len(rel) else np.nan
    row["ah_counter_vs_integral_p95_abs_rel"] = float(rel.abs().quantile(0.95)) if len(rel) else np.nan
    row["discharge_v_last_min"] = float(dis["v_last"].min()) if len(dis) else np.nan
    rd = rest_durations(st)
    row["adjacent_rest_s_median"] = float(np.median(rd)) if rd.size else np.nan
    row["adjacent_rest_s_p05"] = float(np.percentile(rd, 5)) if rd.size else np.nan
    st["ah"] = st[params["ah_column"]]
    every = c["c_discharge"] == 0.5 and c["soc_min"] == 0 and c["soc_max"] == 100
    d = LT.discharges(st, c["rated_Ah"], params["v_cut_low"][c["chemistry"]], params["tol_v"], params["tol_c"], every)
    ch = LT.checks_table(d)
    row["checks"] = len(ch)
    row["check_discharges"] = int((d["check_id"] >= 0).sum())
    row["check_efc_spacing_median"] = float(np.median(np.diff(ch["efc"]))) if len(ch) > 1 else np.nan
    row["efc_total"] = float((d["ah"].sum()) / c["rated_Ah"])
    row["first_check_Ah"] = float(ch["q_Ah"].iloc[0]) if len(ch) else np.nan
    row["lowest_check_frac_of_first"] = float(ch["q_Ah"].min() / ch["q_Ah"].iloc[0]) if len(ch) else np.nan
    return row


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--steps", required=True)
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    params = json.load(open(a.params))
    cells = IO.discover(a.root)
    rows = []
    for key, c in sorted(cells.items()):
        sp = os.path.join(a.steps, key + "_steps.csv")
        if os.path.exists(sp):
            rows.append(cell_row(c, sp, params))
    inv = pd.DataFrame(rows)
    os.makedirs(a.out, exist_ok=True)
    inv.to_csv(os.path.join(a.out, "inventory.csv"), index=False)
    print(f"{len(inv)} cells inventoried")
