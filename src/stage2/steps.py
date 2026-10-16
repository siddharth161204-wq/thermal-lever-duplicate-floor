#!/usr/bin/env python3
"""Stage 2: one streaming pass over a Sandia time series that splits it into steps.

Each sample is labelled rest when |I| < i_rest_A, and charge or discharge by the sign of the current (the archive's
convention, discharge negative, is checked by the inventory against the capacity counters). Consecutive samples with
the same label form one step. For every step the pass records its time span, sample count, the archive's cycle index
at its first and last sample, the ampere-hours by trapezoidal integration of |I| over the step, the cycler's own
charge and discharge counters at the sample before the step, at its first and last sample and their maximum within
it, the mean |I| over the step's duration, and the voltage at its first and last sample with its minimum and maximum.

Without with_temperature no temperature column is read. With it, the time-weighted mean (trapezoidal), first, last,
minimum and maximum of the cell (skin) temperature and of the environment temperature are added per step. The
temperature pass runs only after the Stage 2 predictions are public (DECISIONS.md, Stage 2).
"""
import argparse, json, os, sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sandia_io as IO

BASE = ["time", "cycle", "current", "voltage", "q_ch", "q_dis"]
TEMP = ["t_cell", "t_env"]


def _runs(lab):
    starts = np.flatnonzero(np.r_[True, lab[1:] != lab[:-1]])
    ends = np.r_[starts[1:], lab.size]
    return starts, ends


def _aggregate(a, starts, ends, prev, with_temperature):
    """Aggregate complete runs [starts[k], ends[k]) of the arrays in a (dict of equal-length arrays)."""
    t, I = a["time"], a["current"]
    n = t.size
    nr = starts.size
    run_id = np.repeat(np.arange(nr), ends - starts)
    same = run_id[:-1] == run_id[1:]
    dt = np.diff(t)
    left = run_id[:-1][same]
    absI = np.abs(I)
    q_int = np.bincount(left, weights=(0.5 * (absI[:-1] + absI[1:]) * dt)[same], minlength=nr)
    t0, t1 = t[starts], t[ends - 1]
    dur = t1 - t0
    rows = {
        "t_start_s": t0, "t_end_s": t1, "duration_s": dur, "n_samples": ends - starts,
        "cycle_first": a["cycle"][starts], "cycle_last": a["cycle"][ends - 1],
        "ah_int": q_int / 3600.0,
        "mean_abs_I_A": np.where(dur > 0, q_int / np.where(dur > 0, dur, 1.0), absI[starts]),
        "v_first": a["voltage"][starts], "v_last": a["voltage"][ends - 1],
        "v_min": np.minimum.reduceat(a["voltage"], starts), "v_max": np.maximum.reduceat(a["voltage"], starts),
    }
    for q in ("q_ch", "q_dis"):
        v = a[q]
        rows[q + "_first"] = v[starts]
        rows[q + "_last"] = v[ends - 1]
        rows[q + "_max"] = np.fmax.reduceat(v, starts)
        before = np.r_[prev[q] if prev else np.nan, v[ends[:-1] - 1]] if nr else np.array([])
        rows[q + "_before"] = before
    rows["t_before_s"] = np.r_[prev["time"] if prev else np.nan, t[ends[:-1] - 1]] if nr else np.array([])
    if with_temperature:
        for key in TEMP:
            v = a.get(key)
            if v is None:
                for s in ("mean", "first", "last", "min", "max"):
                    rows[f"{key}_{s}"] = np.full(nr, np.nan)
                continue
            integ = np.bincount(left, weights=(0.5 * (v[:-1] + v[1:]) * dt)[same], minlength=nr)
            rows[f"{key}_mean"] = np.where(dur > 0, integ / np.where(dur > 0, dur, 1.0), v[starts])
            rows[f"{key}_first"] = v[starts]
            rows[f"{key}_last"] = v[ends - 1]
            rows[f"{key}_min"] = np.fmin.reduceat(v, starts)
            rows[f"{key}_max"] = np.fmax.reduceat(v, starts)
    last = {k: a[k][-1] for k in ("time", "q_ch", "q_dis")} if n else prev
    return rows, last


def segment(path, i_rest_A, discharge_sign=-1, with_temperature=False, chunksize=500_000):
    """Steps of one time-series file as a DataFrame; also returns the structural counts of the pass."""
    cols = BASE + (TEMP if with_temperature else [])
    buf, prev, out, last_t = None, None, [], None
    stats = {"rows": 0, "rows_dropped_nan_time_or_current": 0, "negative_time_steps": 0, "missing_columns": None}
    for chunk, missing in IO.iter_timeseries(path, cols, chunksize):
        stats["missing_columns"] = missing
        ok = np.isfinite(chunk["time"]) & np.isfinite(chunk["current"])
        stats["rows"] += int(ok.size)
        stats["rows_dropped_nan_time_or_current"] += int((~ok).sum())
        chunk = {k: v[ok] for k, v in chunk.items()}
        t = chunk["time"]
        if t.size:
            stats["negative_time_steps"] += int((np.diff(t) < 0).sum()) + int(last_t is not None and t[0] < last_t)
            last_t = t[-1]
        if buf is not None:
            chunk = {k: np.r_[buf[k], chunk[k]] for k in chunk}
        I = chunk["current"] * (-discharge_sign)          # positive = charge after this line
        lab = np.where(I >= i_rest_A, 1, np.where(I <= -i_rest_A, -1, 0)).astype(np.int8)
        starts, ends = _runs(lab)
        k = starts.size - 1                                # the last run may continue in the next chunk
        if k > 0:
            sub = {key: v[:starts[k]] for key, v in chunk.items()}
            rows, prev = _aggregate(sub, starts[:k], ends[:k], prev, with_temperature)
            rows["kind"] = lab[starts[:k]]
            out.append(pd.DataFrame(rows))
        buf = {key: v[starts[k]:] for key, v in chunk.items()}
    if buf is not None and buf["time"].size:
        I = buf["current"] * (-discharge_sign)
        lab = np.where(I >= i_rest_A, 1, np.where(I <= -i_rest_A, -1, 0)).astype(np.int8)
        starts, ends = _runs(lab)
        rows, prev = _aggregate(buf, starts, ends, prev, with_temperature)
        rows["kind"] = lab[starts]
        out.append(pd.DataFrame(rows))
    steps = pd.concat(out, ignore_index=True) if out else pd.DataFrame()
    if len(steps):
        steps["kind"] = steps["kind"].map({1: "charge", -1: "discharge", 0: "rest"})
        steps.insert(0, "step", np.arange(len(steps)))
        steps["gap_before_s"] = steps["t_start_s"] - steps["t_before_s"]
    return steps, stats


def sampling_intervals(path, chunksize=500_000, max_rows=None):
    """Structure only: the distribution of the time between consecutive samples, overall."""
    qs, n, prev = [], 0, None
    dts = []
    for chunk, _ in IO.iter_timeseries(path, ["time"], chunksize):
        t = chunk["time"][np.isfinite(chunk["time"])]
        if prev is not None and t.size:
            t = np.r_[prev, t]
        if t.size > 1:
            dts.append(np.diff(t))
        prev = t[-1] if t.size else prev
        n += t.size
        if max_rows and n >= max_rows:
            break
    d = np.concatenate(dts) if dts else np.array([])
    d = d[d > 0]
    if not d.size:
        return {}
    return {"dt_median_s": float(np.median(d)), "dt_p05_s": float(np.percentile(d, 5)),
            "dt_p95_s": float(np.percentile(d, 95)), "dt_max_s": float(d.max()), "n_intervals": int(d.size)}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("timeseries")
    ap.add_argument("--i-rest-A", type=float, required=True)
    ap.add_argument("--discharge-sign", type=int, default=-1)
    ap.add_argument("--temperature", action="store_true")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    st, stats = segment(a.timeseries, a.i_rest_A, a.discharge_sign, a.temperature)
    st.to_csv(a.out, index=False)
    print(json.dumps(stats))
