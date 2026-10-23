#!/usr/bin/env python3
"""Stage 3: the rows of a Sandia time series in cycle order, and the ampere-hours of each step from the cycler's
counters with their resets (DECISIONS.md D50).

Cycle order. In 20 files the cycler logs the last record of a cycle up to a second after the first record of the
next (one cycle below the running maximum of the cycle index); a stable sort by cycle index returns each such record
to the end of its own cycle. In one file (NMC, 25 C, 20 to 80 %, 0.5/3C, cell a) rows of another test, more than one
cycle below the running maximum, are interleaved with the cell's own; they are dropped.

Ampere-hours. The counters restart at most cycle-index changes. A step's ampere-hours are the sum, over the transitions
from the sample before the step to its last sample, of the counter's increase, a decrease being a restart after which
the new value counts in full. Where the counter restarts at most once, at the step's first sample, this equals
src/stage2/process_cells.counter_ah; it differs where a discharge's last sample already belongs to the next cycle.
"""
import os, sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "stage2"))
import sandia_io as IO  # noqa: E402

LOGICAL = ("time", "cycle", "current", "voltage", "q_ch", "q_dis", "t_cell", "t_env")


def load_rows(path):
    """The file's rows with time and current, in cycle order; returns the raw frame (original columns, cleaned order),
    the logical arrays and the counts of what was moved or dropped."""
    m = IO.match_columns(IO.header(path))
    raw = pd.read_csv(path, low_memory=False)
    a = {k: pd.to_numeric(raw[m[k]], errors="coerce").to_numpy(float) for k in LOGICAL}
    ok = np.isfinite(a["time"]) & np.isfinite(a["current"])
    if not np.all(np.isfinite(a["cycle"][ok])):
        raise SystemExit(f"cycle index missing in rows with time and current: {path}")
    idx = np.flatnonzero(ok)
    c = a["cycle"][idx]
    prevmax = np.r_[-np.inf, np.maximum.accumulate(c)[:-1]]
    foreign = (prevmax - c) > 1.0 + 1e-9                    # more than one cycle below the running maximum
    late = (c < prevmax) & ~foreign                         # the previous cycle's last record, logged late
    keep = idx[~foreign]
    order = np.argsort(a["cycle"][keep], kind="stable")
    sel = keep[order]
    moved = int(np.sum(sel != keep))
    t = a["time"][sel]
    back = np.diff(t)
    stats = {"rows": int(ok.sum()), "rows_without_time_or_current": int((~ok).sum()), "foreign_rows_dropped": int(foreign.sum()),
             "late_rows": int(late.sum()), "rows_moved": moved,
             "largest_time_step_back_s": float(-back.min()) if back.size and back.min() < 0 else 0.0,
             "time_steps_back": int((back < 0).sum())}
    if foreign.any():
        fi = idx[foreign]
        stats["foreign_row_range"] = [int(fi.min()), int(fi.max())]
        stats["foreign_cycles"] = [float(a["cycle"][fi].min()), float(a["cycle"][fi].max())]
    arr = {k: v[sel] for k, v in a.items()}
    return raw.iloc[sel].reset_index(drop=True), arr, stats


def labels(current, rated, i_rest_frac=0.005, discharge_sign=-1):
    """Sample labels exactly as src/stage2/steps.py: 1 charge, -1 discharge, 0 rest."""
    I = current * (-discharge_sign)
    i_rest = i_rest_frac * rated
    return np.where(I >= i_rest, 1, np.where(I <= -i_rest, -1, 0)).astype(np.int8)


def runs(lab):
    starts = np.flatnonzero(np.r_[True, lab[1:] != lab[:-1]])
    ends = np.r_[starts[1:], lab.size]
    return starts, ends


def increments(q):
    """Per sample, the counter's increase since the sample before; at a decrease (a restart) the new value."""
    inc = np.empty_like(q)
    if q.size:
        inc[0] = q[0]
        dq = np.diff(q)
        inc[1:] = np.where(dq >= -1e-9, dq, q[1:])
    return inc


def step_ah(arr, lab, starts, ends):
    """Ampere-hours of each step: discharge steps from the discharge counter, charge steps from the charge counter."""
    rid = np.repeat(np.arange(starts.size), ends - starts)
    out = np.full(starts.size, np.nan)
    for kind, q in ((1, "q_ch"), (-1, "q_dis")):
        s = np.bincount(rid, weights=np.nan_to_num(increments(arr[q])), minlength=starts.size)
        m = lab[starts] == kind
        out[m] = s[m]
    return out
