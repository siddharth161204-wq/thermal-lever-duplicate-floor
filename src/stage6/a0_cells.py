#!/usr/bin/env python3
"""Stage 6, Part A0 in the workspace: the cells and the sets, each cell's temperature measures, its contact flag and its
resistance (DECISIONS.md D86 to D88, D90), from the extracts that `a0_read.py` wrote where the files are. No lifetime or
capacity is read.

Usage: python3 -I src/stage6/a0_cells.py [extract directory] [output directory]
Defaults: inputs/matr/a0 and results. Writes stage6_cells.csv and stage6_sets.csv (the cycle life is added after
Part A0 is committed, D94).
"""
import csv, math, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s6_common as C                   # noqa: E402
import a0_traces as TR                  # noqa: E402

# D86: the four batch files, their papers, and the authors' own exclusions and carry-over: the Python loader,
# `Load Data.ipynb` of github.com/rdbraatz/data-driven-prediction-of-battery-cycle-life-before-capacity-degradation at
# 1ef13d2, by struct index from 0 as its keys b1c0 ... are; for the validation batch, the authors' validation code
# takes all 45 cells of the struct. The site's channel numbers for the same cells are checked in D86.
BATCHES = {
    "b1": {"date": "2017-05-12", "paper": "Severson", "authors_prefix": "b1c",
           "drop": {i: "excluded by the authors' loader: does not reach 80 % of nominal capacity"
                    for i in (8, 10, 12, 13, 22)}},
    "b2": {"date": "2017-06-30", "paper": "Severson", "authors_prefix": "b2c", "drop": {},
           "carried": {7: ("b1", 0), 8: ("b1", 1), 9: ("b1", 2), 15: ("b1", 3), 16: ("b1", 4)}},
    "b3": {"date": "2018-04-12", "paper": "Severson", "authors_prefix": "b3c",
           "drop": {i: "excluded by the authors' loader: noisy channel" for i in (37, 2, 23, 32, 42, 43)}},
    "b9": {"date": "2019-01-24", "paper": "Attia", "authors_prefix": None, "drop": {}},
}
WINDOW = (2, 100)
TEMPS = ("Tavg", "Tmax", "Tmin")


def load(extract_dir, label):
    cells = list(csv.DictReader(open(os.path.join(extract_dir, f"{label}_a0_cells.csv"))))
    cyc = {}
    for r in csv.DictReader(open(os.path.join(extract_dir, f"{label}_a0_cycles.csv"))):
        cyc.setdefault(int(r["index"]), []).append(r)
    return cells, cyc


def swap_temperatures(cells, cyc, date):
    """D87: the site's swapped thermocouples; each cell gets the other's temperature columns, cycle by cycle."""
    pair = TR.SITE_SWAPPED.get(date)
    if not pair:
        return []
    idx = {str(c["channel_id"]).strip(): int(c["index"]) for c in cells}
    if pair[0] not in idx or pair[1] not in idx:
        return []
    a, b = idx[pair[0]], idx[pair[1]]
    ra = {int(r["cycle"]): r for r in cyc.get(a, [])}
    rb = {int(r["cycle"]): r for r in cyc.get(b, [])}
    for n in set(ra) & set(rb):
        for k in TEMPS:
            ra[n][k], rb[n][k] = rb[n][k], ra[n][k]
    return [a, b]


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def series(rows, key, lo, hi, parity=None):
    """The finite values of one summary series over cycles lo to hi, optionally odd (1) or even (0) cycles only."""
    out = []
    for r in rows:
        c = int(r["cycle"])
        if lo <= c <= hi and (parity is None or c % 2 == parity):
            v = num(r[key])
            if math.isfinite(v):
                out.append(v)
    return np.array(out)


def mean(v):
    return float(v.mean()) if len(v) else float("nan")


def measures(rows):
    """D87 and D90: the registered temperature, its brackets and halves, and the resistance."""
    lo, hi = WINDOW
    ir = series(rows, "IR", lo, hi)
    irv = ir[ir > 0]
    ir_early = series(rows, "IR", 2, 10)
    ir_early = ir_early[ir_early > 0]
    m = {"T_reg": mean(series(rows, "Tavg", lo, hi)), "T_max": mean(series(rows, "Tmax", lo, hi)),
         "T_2_50": mean(series(rows, "Tavg", 2, 50)), "T_min": mean(series(rows, "Tmin", lo, hi)),
         "T_odd": mean(series(rows, "Tavg", lo, hi, 1)), "T_even": mean(series(rows, "Tavg", lo, hi, 0)),
         "T_51_100": mean(series(rows, "Tavg", 51, 100)), "T_cycles": int(len(series(rows, "Tavg", lo, hi))),
         "IR_ohm": mean(irv), "IR_valid_cycles": int(len(irv)), "IR_2_10_ohm": mean(ir_early)}
    m["T_rise"] = m["T_reg"] - m["T_min"]
    m["lnR"] = math.log(m["IR_ohm"]) if math.isfinite(m["IR_ohm"]) and m["IR_ohm"] > 0 else float("nan")
    m["lnR_2_10"] = math.log(m["IR_2_10_ohm"]) if math.isfinite(m["IR_2_10_ohm"]) and m["IR_2_10_ohm"] > 0 else float("nan")
    return m


def build(extract_dir, batches=BATCHES):
    """One row per struct entry of every batch: whether it enters and why not, its set, measures and contact flag."""
    rows, stats, swapped = [], {}, {}
    for label, cfg in batches.items():
        cells, cyc = load(extract_dir, label)
        swapped[cfg["date"]] = swap_temperatures(cells, cyc, cfg["date"])
        carried = cfg.get("carried", {})
        for c in cells:
            i = int(c["index"])
            r = {"cell": f"{cfg['date']}:{i}", "paper": cfg["paper"], "batch": cfg["date"], "label": label, "index": i,
                 "authors_key": f"{cfg['authors_prefix']}{i}" if cfg.get("authors_prefix") else "",
                 "policy": c["policy_readable"], "channel": str(c["channel_id"]).strip(), "barcode": c["barcode"],
                 "window_complete": c["window_complete"] == "True"}
            if i in cfg.get("drop", {}):
                r["enters"], r["why_not"] = False, cfg["drop"][i]
            elif i in carried:
                to = carried[i]
                r["enters"], r["why_not"] = False, (f"continuation of {batches[to[0]]['date']}:{to[1]}, appended to it "
                                                    "by the authors' loader")
            elif not r["window_complete"]:
                r["enters"], r["why_not"] = False, "cycles 2 to 100 not all in the file"
            else:
                r["enters"], r["why_not"] = True, ""
            r.update(measures(cyc.get(i, [])))
            stats[r["cell"]] = TR.trace_stats(cyc.get(i, []))
            rows.append(r)
    sets = {}
    for r in rows:
        if r["enters"]:
            sets.setdefault((r["batch"], r["policy"]), []).append(r)
    for k, members in sets.items():
        for r in members:
            r["set"], r["set_size"] = f"{k[0]}|{k[1]}", len(members)
    for r in rows:
        r.setdefault("set", "")
        r.setdefault("set_size", 0)
    TR.flag(rows, stats)
    return rows, sets


COLUMNS = ["cell", "paper", "batch", "index", "authors_key", "policy", "channel", "barcode", "enters", "why_not", "set",
           "set_size", "T_reg", "T_max", "T_2_50", "T_min", "T_rise", "T_odd", "T_even", "T_51_100", "T_cycles",
           "swing_med", "swing_early", "swing_late", "floor_med", "IR_ohm", "IR_valid_cycles", "lnR", "IR_2_10_ohm",
           "lnR_2_10", "contact_flag", "contact_reasons", "site_note"]


def fmt(v):
    return repr(v) if isinstance(v, float) else v


def write(rows, sets, out_dir, extra=()):
    os.makedirs(out_dir, exist_ok=True)
    cols = COLUMNS + list(extra)
    with open(os.path.join(out_dir, "stage6_cells.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: fmt(r.get(k, "")) for k in cols})
    with open(os.path.join(out_dir, "stage6_sets.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["set", "paper", "batch", "policy", "cells", "pairs", "channels", "T_reg_mean_C", "T_reg_sd_within_K",
                    "T_reg_range_K", "flagged_cells"])
        for (b, p), m in sorted(sets.items()):
            T = np.array([x["T_reg"] for x in m])
            w.writerow([f"{b}|{p}", m[0]["paper"], b, p, len(m), len(m) * (len(m) - 1) // 2,
                        " ".join(x["channel"] for x in m), repr(float(T.mean())),
                        repr(float(T.std(ddof=1))) if len(m) > 1 else "", repr(float(T.max() - T.min())),
                        sum(bool(x["contact_flag"]) for x in m)])


if __name__ == "__main__":
    ext = sys.argv[1] if len(sys.argv) > 1 else C.path("inputs/matr/a0")
    out = sys.argv[2] if len(sys.argv) > 2 else C.path("results")
    rows, sets = build(ext)
    write(rows, sets, out)
    print(len(rows), "entries,", sum(r["enters"] for r in rows), "enter,", len(sets), "sets,",
          sum(len(m) >= 2 for m in sets.values()), "of two or more,",
          sum(bool(r["contact_flag"]) for r in rows if r["enters"]), "flagged")
