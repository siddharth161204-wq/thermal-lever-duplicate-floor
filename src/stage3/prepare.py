#!/usr/bin/env python3
"""Stage 3, structure only (DECISIONS.md D50 (1) and (2)): the rows of every time series in cycle order, the step pass
of src/stage2/steps.py repeated on the ordered rows of the cells whose rows moved or were dropped, and the ampere-hours
of every step from the counters with their restarts, beside src/stage2/process_cells.counter_ah. No temperature column
is averaged and no capacity check or lifetime is formed here.

Usage: python3 -I src/stage3/prepare.py --raw RAW_DIR --steps PASS1_STEPS_DIR --work WORK_DIR [--out results/stage3_prepare.json]
"""
import argparse, glob, json, os, sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "stage2"))
import rows as RW            # noqa: E402
import sandia_io as IO       # noqa: E402
import steps as ST           # noqa: E402
import process_cells as PC   # noqa: E402

I_REST = 0.005


def prepare(raw_dir, pass1_dir, work):
    os.makedirs(os.path.join(work, "clean"), exist_ok=True)
    os.makedirs(os.path.join(work, "steps3"), exist_ok=True)
    report = {}
    paths = sorted(glob.glob(os.path.join(raw_dir, "**", "*_timeseries.csv"), recursive=True))
    keys = [IO.parse_cell_id(os.path.basename(p))["key"] for p in paths]
    if len(set(keys)) != len(keys) or os.path.abspath(work).startswith(os.path.abspath(raw_dir) + os.sep):
        raise SystemExit("a cell appears twice under the raw folder, or the work folder lies inside it")
    for p in paths:
        meta = IO.parse_cell_id(os.path.basename(p))
        rated = meta["rated_Ah"]
        raw, arr, stats = RW.load_rows(p)
        changed = stats["rows_moved"] > 0 or stats["foreign_rows_dropped"] > 0
        if changed:
            cpath = os.path.join(work, "clean", meta["key"] + "_timeseries.csv")
            raw.to_csv(cpath, index=False)
            st, sstats = ST.segment(cpath, I_REST * rated, -1, False)
            stats["step_pass"] = "repeated on the ordered rows"
            stats["steps_pass1"] = int(len(pd.read_csv(os.path.join(pass1_dir, meta["key"] + "_steps.csv"), usecols=["step"])))
        else:
            st = pd.read_csv(os.path.join(pass1_dir, meta["key"] + "_steps.csv"))
            stats["step_pass"] = "pass 1"
        lab = RW.labels(arr["current"], rated, I_REST)
        s, e = RW.runs(lab)
        kinds = pd.Series(lab[s]).map({1: "charge", -1: "discharge", 0: "rest"}).to_numpy()
        if not (s.size == len(st) and np.array_equal(kinds, st.kind.to_numpy())
                and np.allclose(arr["time"][s], st.t_start_s.to_numpy(), rtol=0, atol=1e-6)):
            raise SystemExit(f"ordered rows do not reproduce the step table: {meta['key']}")
        st["ah_counter"] = RW.step_ah(arr, lab, s, e)
        st["ah_counter_ah_stage2"] = PC.counter_ah(st)
        st.to_csv(os.path.join(work, "steps3", meta["key"] + "_steps.csv"), index=False)
        dm = (st.kind == "discharge").to_numpy()
        diff = np.abs(st.ah_counter.to_numpy() - st.ah_counter_ah_stage2.to_numpy())
        dec = np.flatnonzero(np.diff(arr["q_dis"]) < -1e-9) + 1
        stats.update({"steps": int(len(st)), "discharges": int(dm.sum()),
                      "discharges_where_counter_ah_differs": int((diff[dm] > 1e-6).sum()),
                      "discharge_ah_total": float(st.ah_counter[dm].sum()),
                      "discharge_ah_total_counter_ah": float(np.nansum(st.ah_counter_ah_stage2[dm])),
                      "discharge_counter_decreases": int(dec.size),
                      "discharge_counter_decreases_without_cycle_change": int((arr["cycle"][dec] == arr["cycle"][dec - 1]).sum()),
                      "discharge_counter_decreases_without_cycle_change_in_discharge_steps":
                          int(((arr["cycle"][dec] == arr["cycle"][dec - 1]) & (lab[dec] == -1)).sum())})
        report[meta["key"]] = stats
    return report


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--steps", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--out")
    a = ap.parse_args()
    rep = prepare(a.raw, a.steps, a.work)
    tot = {k: sum(v.get(k, 0) for v in rep.values() if isinstance(v.get(k, 0), (int, float)))
           for k in ("late_rows", "rows_moved", "foreign_rows_dropped", "discharges", "discharges_where_counter_ah_differs",
                     "discharge_counter_decreases", "discharge_counter_decreases_without_cycle_change",
                     "discharge_counter_decreases_without_cycle_change_in_discharge_steps")}
    tot["cells_with_moved_or_dropped_rows"] = sum(1 for v in rep.values() if v["step_pass"] != "pass 1")
    tot["cells_where_counter_ah_differs"] = sum(1 for v in rep.values() if v["discharges_where_counter_ah_differs"])
    out = {"totals": tot, "cells": rep}
    if a.out:
        json.dump(out, open(a.out, "w"), indent=1)
    print(json.dumps(tot, indent=1))
