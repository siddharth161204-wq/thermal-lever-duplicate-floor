#!/usr/bin/env python3
"""Stage 2: from each cell's step table to the per-cell table that analysis.py reads.

For every cell: its capacity checks, its lifetime at its chemistry's registered SoH threshold with the resolution
and censoring, the sensitivity lifetime for 0 to 100 % windows, and, when the step table carries the temperature
pass, measures (a) and (b). The parameters are the values fixed in DECISIONS.md (Stage 2), passed as JSON.

Usage:
  python3 src/stage2/process_cells.py --steps DIR --params params.json --out results/stage2_cells.csv \
      [--checks-out results/stage2_checks.csv]
"""
import argparse, glob, json, os, sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lifetimes as LT  # noqa: E402
import sandia_io as IO  # noqa: E402


def counter_ah(st):
    """Ampere-hours of each charge or discharge step from the cycler's own counters: the counter at the step's last
    sample minus its value at the sample before the step, or the last value alone when the counter was reset between
    them (a new cycle index starting with the step)."""
    out = np.full(len(st), np.nan)
    for kind, q in (("charge", "q_ch"), ("discharge", "q_dis")):
        m = (st["kind"] == kind).to_numpy()
        last, before = st[q + "_last"].to_numpy(), st[q + "_before"].to_numpy()
        reset = ~np.isfinite(before) | (before > last + 1e-9)
        out[m] = np.where(reset, last, last - before)[m]
    return out


def process(steps_dir, params, with_temperature):
    rows, all_checks = [], []
    suffix = "_steps_T.csv" if with_temperature else "_steps.csv"
    for path in sorted(glob.glob(os.path.join(steps_dir, "*" + suffix))):
        meta = IO.parse_cell_id(os.path.basename(path))
        if meta is None:
            continue
        st = pd.read_csv(path)
        chem = meta["chemistry"]
        st["ah_ctr"] = counter_ah(st)
        st["ah"] = st[params["ah_column"]]
        every = meta["c_discharge"] == 0.5 and meta["soc_min"] == 0 and meta["soc_max"] == 100
        d = LT.discharges(st, meta["rated_Ah"], params["v_cut_low"][chem], params["tol_v"], params["tol_c"], every)
        ch = LT.checks_table(d)
        frac = params["thresholds"][chem]
        life = LT.check_lifetime(ch, frac)
        row = dict(meta, soh_threshold=frac, **{k: v for k, v in life.items() if k != "bracket"})
        row["bracket_efc"] = json.dumps(life.get("bracket"))
        if meta["soc_min"] == 0 and meta["soc_max"] == 100:
            row.update(LT.sensitivity_lifetime(d, frac, meta["c_discharge"], params["tol_c"]))
        if with_temperature:
            prim = params.get("rest_baseline", "mean")
            row.update(LT.temperature_measures(st, d, life["L_efc"], life["efc_last_check"], prim, every))
            alt = LT.temperature_measures(st, d, life["L_efc"], life["efc_last_check"], "end" if prim == "mean" else "mean", every)
            row["rise_alt_K"], row["rest_alt_T_C"] = alt["rise_K"], alt["rest_T_C"]
        rows.append(row)
        ch.insert(0, "cell_id", meta["cell_id"])
        all_checks.append(ch)
    cells = pd.DataFrame(rows)
    checks = pd.concat(all_checks, ignore_index=True) if all_checks else pd.DataFrame()
    return cells, checks


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", required=True)
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--checks-out")
    ap.add_argument("--temperature", action="store_true")
    a = ap.parse_args()
    params = json.load(open(a.params))
    cells, checks = process(a.steps, params, a.temperature)
    cells.to_csv(a.out, index=False)
    if a.checks_out:
        checks.to_csv(a.checks_out, index=False)
    print(f"{len(cells)} cells; censored {int(cells.censored.sum()) if len(cells) else 0}")
