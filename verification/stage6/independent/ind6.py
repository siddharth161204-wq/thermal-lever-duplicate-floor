#!/usr/bin/env python3
"""Independent check of Stage 6 (the lever within a fast-charging policy).

Written by the independent checker from the brief, predictions_stage6.md, DECISIONS.md D82 to D96 (and the earlier
entries they cite), and the committed inputs only: inputs/matr/a0/*, inputs/matr/life/*, inputs/matr/SOURCES.md and
MANIFEST.md, results/stage3_lever.csv, results/stage4_experiment.json, results/stage5_paper_numbers.json,
inputs/floor.json. Nothing of src/stage6, tests/stage6_*, results/stage6_* or verification/stage6 outside this
folder was opened, and inputs/matr/stage6_lifetimes.csv was not opened: the lifetimes are built here from
inputs/matr/life/*.csv and D86, D94.

Estimators are written here (fixed-effects slopes, two-regressor fixed effects, weighted between-set slope with batch
constants, Poisson-binomial exact test, percentile bootstrap, permutation test, the noncentral-t power by numerical
integration). scipy is used only for distribution functions (normal cdf, t quantile, log-gamma), quadrature and root
finding.

Usage: python3 -I verification/stage6/independent/ind6.py [--replicates 20]
Writes record.json, cells.csv and sets.csv beside this file.
"""
import argparse, csv, hashlib, json, math, os, sys, time

import numpy as np
from scipy import integrate, optimize, special, stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
MATR = os.path.join(ROOT, "inputs", "matr")

# ---------------------------------------------------------------------------------------------------------------
# Fixed facts from D84 to D88 (struct indices, channels, chamber)
# ---------------------------------------------------------------------------------------------------------------
BATCHES = [("b1", "2017-05-12", 1), ("b2", "2017-06-30", 1), ("b3", "2018-04-12", 1), ("b9", "2019-01-24", 2)]
BATCH_ORDER = {b[0]: i for i, b in enumerate(BATCHES)}
DATE = {b[0]: b[1] for b in BATCHES}
PAPER = {b[0]: b[2] for b in BATCHES}
EXCLUDE = {"b1": [8, 10, 12, 13, 22], "b2": [7, 8, 9, 15, 16], "b3": [2, 23, 32, 37, 42, 43], "b9": []}
CARRY = [(0, 7, 662), (1, 8, 981), (2, 9, 1060), (3, 15, 208), (4, 16, 482)]   # (b1 index, b2 continuation index, added length)
SWAP = {"b1": ("15", "16"), "b2": ("15", "16")}                                 # D87 (3)
DETACHED = {"b2": ("7", "21")}                                                    # D88 (f)
CHAMBER_C = 30.0

# Archived inputs (D82 (2)); read from the files below and checked against the full-precision values D82 and D96 quote
SANDIA_KEYS = {"lfp_1c_15_25": ("registered 1C", "LFP", "15 to 25 C"),
               "lfp_2c_15_25": ("secondary 2C", "LFP", "15 to 25 C"),
               "lfp_2c_25_35": ("secondary 2C", "LFP", "25 to 35 C")}

D_BOOT = 10000
ALPHA, POWER = 0.05, 0.80

# My own seeds (not the registered ones). Generator: numpy Philox counter-based bit generator, keyed by SeedSequence.
SEED_BASE = 5551117


def rng_for(*key):
    return np.random.Generator(np.random.Philox(np.random.SeedSequence([SEED_BASE, *key])))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def fnum(x):
    """Float or None for json (nan and inf become None)."""
    if x is None:
        return None
    x = float(x)
    return x if math.isfinite(x) else None


# ---------------------------------------------------------------------------------------------------------------
# 1. The entries, the traces, the cell table
# ---------------------------------------------------------------------------------------------------------------
def load_entries():
    entries = {}
    files = {}
    for lab, date, paper in BATCHES:
        p_cells = os.path.join(MATR, "a0", f"{lab}_a0_cells.csv")
        p_cyc = os.path.join(MATR, "a0", f"{lab}_a0_cycles.csv")
        p_struct = os.path.join(MATR, "a0", f"{lab}_a0_structure.json")
        p_life = os.path.join(MATR, "life", f"{lab}_life.csv")
        for p in (p_cells, p_cyc, p_struct, p_life):
            files[os.path.relpath(p, ROOT)] = {"bytes": os.path.getsize(p), "sha256": sha256(p)}
        desc = read_csv(p_cells)
        struct = json.load(open(p_struct))
        n = len(desc)
        if struct["cells"] != n or struct["batch_date"] != date:
            raise SystemExit(f"structure mismatch {lab}")
        tr = {k: np.full((n, 101), np.nan) for k in ("Tavg", "Tmax", "Tmin", "IR")}
        seen = np.zeros((n, 101), dtype=int)
        for r in read_csv(p_cyc):
            if r["label"] != lab:
                raise SystemExit("label")
            i, c = int(r["index"]), int(r["cycle"])
            if not (1 <= c <= 100):
                raise SystemExit("cycle out of range")
            seen[i, c] += 1
            for k in tr:
                tr[k][i, c] = float(r[k])
        if not np.all(seen[:, 1:] == 1):
            raise SystemExit(f"{lab}: cycles 1 to 100 not each present once")
        life = {int(r["index"]): r for r in read_csv(p_life)}
        if sorted(life) != list(range(n)):
            raise SystemExit(f"{lab}: life rows")
        for d in desc:
            i = int(d["index"])
            if d["label"] != lab or d["key"] != f"{lab}c{i}":
                raise SystemExit("key")
            lf = life[i]
            entries[(lab, i)] = {
                "label": lab, "index": i, "key": d["key"], "date": date, "paper": paper,
                "policy_readable": d["policy_readable"], "policy": d["policy"], "channel": d["channel_id"],
                "barcode": d["barcode"],
                "raw": {k: tr[k][i].copy() for k in tr},
                "life": {"cycle_life_stored": float(lf["cycle_life_stored"]) if lf["cycle_life_stored"] not in ("", None) else float("nan"),
                         "summary_entries": int(lf["summary_entries"]),
                         "first_below_0p88": int(lf["first_below_0p88"]) if lf["first_below_0p88"] not in ("", None) else None,
                         "first_below_0p88_confirmed": int(lf["first_below_0p88_confirmed"]) if lf["first_below_0p88_confirmed"] not in ("", None) else None,
                         "last_QDischarge": float(lf["last_QDischarge"])},
            }
    return entries, files


def apply_swaps(entries):
    """D87 (3): exchange Tavg, Tmax and Tmin cycle by cycle between channels 15 and 16 in 2017-05-12 and 2017-06-30."""
    done = []
    for lab, (ca, cb) in SWAP.items():
        ia = [k for k, e in entries.items() if e["label"] == lab and e["channel"] == ca]
        ib = [k for k, e in entries.items() if e["label"] == lab and e["channel"] == cb]
        if len(ia) != 1 or len(ib) != 1:
            raise SystemExit("swap channels not unique")
        ea, eb = entries[ia[0]], entries[ib[0]]
        tr_a = {k: ea["raw"][k].copy() for k in ("Tavg", "Tmax", "Tmin")}
        tr_b = {k: eb["raw"][k].copy() for k in ("Tavg", "Tmax", "Tmin")}
        ea["tr"] = {**tr_b, "IR": ea["raw"]["IR"].copy()}
        eb["tr"] = {**tr_a, "IR": eb["raw"]["IR"].copy()}
        ea["swapped_with"] = eb["key"]
        eb["swapped_with"] = ea["key"]
        done.append({"batch": DATE[lab], "keys": [ea["key"], eb["key"]], "channels": [ca, cb],
                     "policies": [ea["policy_readable"], eb["policy_readable"]]})
    for e in entries.values():
        if "tr" not in e:
            e["tr"] = {k: e["raw"][k].copy() for k in ("Tavg", "Tmax", "Tmin", "IR")}
            e["swapped_with"] = None
    return done


def entry_facts(entries):
    out = {}
    for lab, date, _ in BATCHES:
        es = [e for e in entries.values() if e["label"] == lab]
        chans = sorted(int(e["channel"]) for e in es)
        out[date] = {"entries": len(es), "channels_distinct": len(set(chans)) == len(chans),
                     "channels_missing_1_to_48": [c for c in range(1, 49) if c not in chans],
                     "excluded": [{"index": i, "channel": entries[(lab, i)]["channel"]} for i in EXCLUDE[lab]],
                     "stored_cycle_life_not_finite": [e["index"] for e in es if not math.isfinite(e["life"]["cycle_life_stored"])]}
    out["total_entries"] = len(entries)
    b2 = {e["channel"]: e["barcode"] for e in entries.values() if e["label"] == "b2"}
    out["b2_barcodes_channels_3_10"] = {"3": b2.get("3"), "10": b2.get("10")}
    return out


def carry_checks(entries):
    out = []
    for i1, i2, add in CARRY:
        a, b = entries[("b1", i1)], entries[("b2", i2)]
        out.append({"b1_index": i1, "b2_index": i2, "channel_b1": a["channel"], "channel_b2": b["channel"],
                    "barcode_b1": a["barcode"], "barcode_b2": b["barcode"], "b2_policy_readable": b["policy_readable"],
                    "added_length": add, "continuation_summary_entries": b["life"]["summary_entries"],
                    "continuation_stored_cycle_life": fnum(b["life"]["cycle_life_stored"]),
                    "length_equals_summary_entries": add == b["life"]["summary_entries"],
                    "same_channel": a["channel"] == b["channel"], "same_barcode": a["barcode"] == b["barcode"]})
    return out


def window(arr, lo, hi):
    """Cycles lo..hi inclusive (cycle n is column n)."""
    return arr[lo:hi + 1]


def median(a):
    return float(np.median(np.asarray(a, dtype=float)))


def cell_measures(e):
    tav, tmx, tmn, ir = (e["tr"][k] for k in ("Tavg", "Tmax", "Tmin", "IR"))
    w = slice(2, 101)
    m = {}
    m["T"] = float(np.mean(tav[w]))
    m["Tmax_m"] = float(np.mean(tmx[w]))
    m["T_2_50"] = float(np.mean(tav[2:51]))
    m["Tmin_m"] = float(np.mean(tmn[w]))
    m["T_minus_Tmin"] = m["T"] - m["Tmin_m"]
    m["T_odd"] = float(np.mean(tav[3:100:2]))         # cycles 3, 5, ..., 99
    m["T_even"] = float(np.mean(tav[2:101:2]))        # cycles 2, 4, ..., 100
    m["T_h1"] = float(np.mean(tav[2:51]))             # 2 to 50
    m["T_h2"] = float(np.mean(tav[51:101]))           # 51 to 100
    m["n_odd"], m["n_even"] = len(range(3, 100, 2)), len(range(2, 101, 2))
    irw = ir[w]
    ok = np.isfinite(irw) & (irw > 0)
    m["R_valid_cycles"] = int(ok.sum())
    m["R"] = float(np.mean(irw[ok])) if ok.any() else float("nan")
    ir10 = ir[2:11]
    ok10 = np.isfinite(ir10) & (ir10 > 0)
    m["R_2_10_valid_cycles"] = int(ok10.sum())
    m["R_2_10"] = float(np.mean(ir10[ok10])) if ok10.any() else float("nan")
    m["lnR"] = math.log(m["R"]) if math.isfinite(m["R"]) else float("nan")
    m["lnR_2_10"] = math.log(m["R_2_10"]) if math.isfinite(m["R_2_10"]) else float("nan")
    sw = tmx - tmn
    m["swing_med_2_100"] = median(sw[2:101])
    m["swing_mean_2_100"] = float(np.mean(sw[2:101]))
    m["swing_med_2_11"] = median(sw[2:12])
    m["swing_med_91_100"] = median(sw[91:101])
    steps = np.abs(np.diff(tav[2:101]))
    m["max_step_Tavg_2_100"] = float(np.max(steps))
    m["max_step_Tavg_1_100"] = float(np.max(np.abs(np.diff(tav[1:101]))))
    # longest run of identical consecutive Tavg values in cycles 2..100
    run, best = 1, 1
    for k in range(3, 101):
        if tav[k] == tav[k - 1]:
            run += 1
            best = max(best, run)
        else:
            run = 1
    m["longest_flat_run_2_100"] = best
    m["min_Tmin_2_100"] = float(np.min(tmn[w]))
    m["median_Tmin_2_100"] = median(tmn[w])
    m["median_Tmin_1_100"] = median(tmn[1:101])
    m["median_Tavg_2_100"] = median(tav[w])
    m["missing_in_window"] = int((~np.isfinite(tav[w])).sum() + (~np.isfinite(tmx[w])).sum() + (~np.isfinite(tmn[w])).sum())
    return m


def lifetime(e, entries):
    """D94 (2): the stored field; for a carried 2017-05-12 cell plus the loader's added length."""
    L = e["life"]["cycle_life_stored"]
    add = 0
    if e["label"] == "b1":
        for i1, i2, a in CARRY:
            if e["index"] == i1:
                add = a
    return L + add, add


def matlab_label(e, entries):
    """Beside, not used (D94 (2)): the first cycle below 0.88 Ah, or one more than the number of cycles, over the
    combined cycles of a carried cell. Two readings: literal (summary entry 1 counted), and skipping a first entry that
    is a zero placeholder (2017-05-12 holds no cycle 1, D84/SOURCES; its first summary entry reads below 0.88)."""
    lf = e["life"]
    n_own = lf["summary_entries"]
    fb = lf["first_below_0p88"]
    cont = None
    if e["label"] == "b1":
        for i1, i2, a in CARRY:
            if e["index"] == i1:
                cont = entries[("b2", i2)]["life"]
    def comb(first_own, skip_first):
        if first_own is not None and not (skip_first and first_own == 1):
            return first_own
        if cont is not None:
            if cont["first_below_0p88"] is not None:
                return n_own + cont["first_below_0p88"]
            return n_own + cont["summary_entries"] + 1
        return n_own + 1
    literal = comb(fb, False)
    # the skip reading only matters when the first below-0.88 entry is entry 1: then look for a later one, which the
    # extract does not give (first_below_0p88 is the first only); in 2017-05-12 the last capacity is above 0.88 in every
    # cell that has entry 1 below, except where noted, so the label is n + 1 (or the continuation's).
    if fb == 1:
        later_known = lf["last_QDischarge"] >= 0.88
        skip = comb(None, True) if later_known else None
    else:
        skip = literal
    return literal, skip


def build_cells(entries):
    excluded = {(lab, i) for lab, ix in EXCLUDE.items() for i in ix}
    cells = []
    for (lab, i), e in sorted(entries.items(), key=lambda kv: (BATCH_ORDER[kv[0][0]], kv[0][1])):
        if (lab, i) in excluded:
            continue
        m = cell_measures(e)
        if m["missing_in_window"]:
            raise SystemExit("a temperature of the window is missing")
        L, add = lifetime(e, entries)
        ml_lit, ml_skip = matlab_label(e, entries)
        c = {"key": e["key"], "label": lab, "batch": DATE[lab], "paper": PAPER[lab], "index": i,
             "policy_readable": e["policy_readable"], "channel": e["channel"], "barcode": e["barcode"],
             "set": f"{lab}|{e['policy_readable']}", "swapped_with": e["swapped_with"],
             "L_stored": e["life"]["cycle_life_stored"], "L_added": add, "L": L,
             "matlab_label_literal": ml_lit, "matlab_label_skip_zero_first": ml_skip, **m}
        c["lnL"] = math.log(L) if (math.isfinite(L) and L > 0) else float("nan")
        cells.append(c)
    return cells, sorted(excluded)


def sets_of(cells):
    sets = {}
    for c in cells:
        sets.setdefault(c["set"], []).append(c)
    return sets


def count_table(cells):
    def counts(sub):
        s = sets_of(sub)
        sizes = [len(v) for v in s.values()]
        multi = [n for n in sizes if n >= 2]
        hist = {}
        for n in multi:
            hist[n] = hist.get(n, 0) + 1
        return {"cells": len(sub), "sets": len(s), "sets_2plus": len(multi), "cells_in_sets_2plus": int(sum(multi)),
                "pairs": int(sum(n * (n - 1) // 2 for n in multi)), "single_cell_sets": int(sum(1 for n in sizes if n == 1)),
                "set_size_histogram_2plus": {str(k): hist[k] for k in sorted(hist)}}
    out = {"pooled": counts(cells), "paper1": counts([c for c in cells if c["paper"] == 1]),
           "paper2": counts([c for c in cells if c["paper"] == 2])}
    for lab, date, _ in BATCHES:
        out[date] = counts([c for c in cells if c["label"] == lab])
    return out


# ---------------------------------------------------------------------------------------------------------------
# 3. The contact rule (D88)
# ---------------------------------------------------------------------------------------------------------------
def contact_flags(cells):
    sets = sets_of(cells)
    detail = {}
    for sid, mem in sets.items():
        big = max(c["swing_med_2_100"] for c in mem) if len(mem) >= 2 else None
        for c in mem:
            r = {}
            r["a"] = c["swing_med_91_100"] < 0.5 * c["swing_med_2_11"]
            r["b_ratio"] = (c["swing_med_2_100"] / big) if big is not None else None
            r["b"] = (r["b_ratio"] is not None) and (c["swing_med_2_100"] < 0.5 * big)
            r["c"] = c["max_step_Tavg_2_100"] > 3.0
            r["d"] = c["longest_flat_run_2_100"] >= 5
            r["e"] = (c["min_Tmin_2_100"] < CHAMBER_C - 5.0) or (c["missing_in_window"] > 0)
            r["f"] = bool((c["label"] in SWAP and c["channel"] in SWAP[c["label"]]) or
                          (c["label"] in DETACHED and c["channel"] in DETACHED[c["label"]]))
            r["flagged"] = bool(r["a"] or r["b"] or r["c"] or r["d"] or r["e"] or r["f"])
            r["rules"] = [k for k in "abcdef" if r[k]]
            detail[c["key"]] = r
            c["flag"] = r["flagged"]
            c["flag_rules"] = "".join(r["rules"])
            c["b_ratio"] = r["b_ratio"]
    return detail


# ---------------------------------------------------------------------------------------------------------------
# Generic within-set and between-set machinery
# ---------------------------------------------------------------------------------------------------------------
def set_list(cells, xkey, zkey=None):
    """Sets (all sizes) with arrays of x, y = ln L and z; order: batch, then first appearance."""
    out = []
    for sid, mem in sets_of(cells).items():
        x = np.array([c[xkey] for c in mem], dtype=float)
        y = np.array([c["lnL"] for c in mem], dtype=float)
        z = np.array([c[zkey] for c in mem], dtype=float) if zkey else np.full(len(mem), np.nan)
        out.append({"id": sid, "label": mem[0]["label"], "batch": mem[0]["batch"], "paper": mem[0]["paper"],
                    "policy_readable": mem[0]["policy_readable"], "n": len(mem), "keys": [c["key"] for c in mem],
                    "x": x, "y": y, "z": z, "L": np.array([c["L"] for c in mem], dtype=float)})
    out.sort(key=lambda s: (BATCH_ORDER[s["label"]], s["keys"][0] if False else min(int(k.split("c")[-1]) for k in s["keys"])))
    return out


def within_sums(s):
    x, y = s["x"], s["y"]
    dx, dy = x - x.mean(), y - y.mean()
    r = {"Sxx": float(dx @ dx), "Sxy": float(dx @ dy), "Syy": float(dy @ dy)}
    # pairs form: oriented warmer first; d = ln(L_warm / L_cool)
    pdx, pdd, pxx, npair = 0.0, 0.0, 0.0, 0
    for i in range(s["n"]):
        for k in range(i + 1, s["n"]):
            if x[i] >= x[k]:
                w, c = i, k
            else:
                w, c = k, i
            dT = x[w] - x[c]
            d = y[w] - y[c]
            pdx += d * dT
            pxx += dT * dT
            pdd += d * d
            npair += 1
    r.update({"Pdx": pdx, "Pxx": pxx, "Pdd": pdd, "npairs": npair})
    # B2: cells with a valid z only
    ok = np.isfinite(s["z"])
    r["nR"] = int(ok.sum())
    if ok.sum() >= 2:
        xr, yr, zr = x[ok], y[ok], s["z"][ok]
        ax, ay, az = xr - xr.mean(), yr - yr.mean(), zr - zr.mean()
        r.update({"R_Sxx": float(ax @ ax), "R_Sxz": float(ax @ az), "R_Szz": float(az @ az), "R_Sxy": float(ax @ ay),
                  "R_Szy": float(az @ ay), "R_Syy": float(ay @ ay), "R_n": int(ok.sum())})
    else:
        r.update({"R_Sxx": 0.0, "R_Sxz": 0.0, "R_Szz": 0.0, "R_Sxy": 0.0, "R_Szy": 0.0, "R_Syy": 0.0, "R_n": int(ok.sum())})
    return r


def b1_point(sets2):
    """Fixed-effects slope with residuals computed cell by cell."""
    Sxx = sum(within_sums(s)["Sxx"] for s in sets2)
    Sxy = sum(within_sums(s)["Sxy"] for s in sets2)
    b = Sxy / Sxx
    rss = 0.0
    N = 0
    for s in sets2:
        dx, dy = s["x"] - s["x"].mean(), s["y"] - s["y"].mean()
        res = dy - b * dx
        rss += float(res @ res)
        N += s["n"]
    J = len(sets2)
    dof = N - J - 1
    sdc = math.sqrt(rss / dof)
    se = math.sqrt(rss / dof / Sxx)
    return {"slope": b, "s_w": -b, "N": N, "J": J, "dof": dof, "rss": rss, "sd_c_hat": sdc, "sqrt2_sd_c_hat": math.sqrt(2) * sdc,
            "se_slope": se, "t": b / se, "S_TT": Sxx, "dT_rms_cells_over_N": math.sqrt(Sxx / N),
            "dT_rms_cells_over_N_minus_J": math.sqrt(Sxx / (N - J))}


def b1_lstsq_check(sets2):
    """Same slope through an explicit dummy-variable least squares (numpy lstsq), as a check of the closed form."""
    rows = []
    J = len(sets2)
    for j, s in enumerate(sets2):
        for xi, yi in zip(s["x"], s["y"]):
            r = np.zeros(J + 1)
            r[j] = 1.0
            r[J] = xi
            rows.append((r, yi))
    A = np.array([r for r, _ in rows])
    yv = np.array([y for _, y in rows])
    coef, *_ = np.linalg.lstsq(A, yv, rcond=None)
    return float(coef[J])


def pairs_point(sets2):
    Pdx = sum(within_sums(s)["Pdx"] for s in sets2)
    Pxx = sum(within_sums(s)["Pxx"] for s in sets2)
    npairs = sum(within_sums(s)["npairs"] for s in sets2)
    b = Pdx / Pxx
    return {"slope": b, "s_w_pairs": -b, "pairs": npairs, "dT_rms_pairs": math.sqrt(Pxx / npairs)}


def b2_point(sets2):
    """ln L on T and ln R with one constant per set, over cells with a valid R; sets with >= 2 such cells."""
    use = [s for s in sets2 if within_sums(s)["R_n"] >= 2]
    if not use:
        return None
    W = [within_sums(s) for s in use]
    Sxx = sum(w["R_Sxx"] for w in W); Sxz = sum(w["R_Sxz"] for w in W); Szz = sum(w["R_Szz"] for w in W)
    Sxy = sum(w["R_Sxy"] for w in W); Szy = sum(w["R_Szy"] for w in W); Syy = sum(w["R_Syy"] for w in W)
    det = Sxx * Szz - Sxz * Sxz
    bT = (Szz * Sxy - Sxz * Szy) / det
    bR = (Sxx * Szy - Sxz * Sxy) / det
    N = sum(w["R_n"] for w in W)
    J = len(use)
    # s_w over the same cells
    bT_only = Sxy / Sxx
    # residual check computed cell by cell
    rss = 0.0
    for s in use:
        ok = np.isfinite(s["z"])
        x, y, z = s["x"][ok], s["y"][ok], s["z"][ok]
        res = (y - y.mean()) - bT * (x - x.mean()) - bR * (z - z.mean())
        rss += float(res @ res)
    # explicit lstsq check
    rows, yv = [], []
    for j, s in enumerate(use):
        ok = np.isfinite(s["z"])
        for xi, yi, zi in zip(s["x"][ok], s["y"][ok], s["z"][ok]):
            r = np.zeros(J + 2); r[j] = 1.0; r[J] = xi; r[J + 1] = zi
            rows.append(r); yv.append(yi)
    coef, *_ = np.linalg.lstsq(np.array(rows), np.array(yv), rcond=None)
    corr_T_lnR = Sxz / math.sqrt(Sxx * Szz)
    return {"lever_T_R_held": -bT, "slope_T_R_held": bT, "dlnL_dlnR": bR, "s_w_same_cells": -bT_only, "N": N, "J": J,
            "dof": N - J - 2, "rss": rss, "lstsq_slope_T": float(coef[J]), "lstsq_slope_lnR": float(coef[J + 1]),
            "within_corr_T_lnR": corr_T_lnR}


def poisson_binomial_pmf(ps):
    pmf = np.array([1.0])
    for p in ps:
        new = np.zeros(len(pmf) + 1)
        new[:-1] += pmf * (1 - p)
        new[1:] += pmf * p
        pmf = new
    return pmf


def exact_two_sided(pmf, k):
    """Outcomes no more probable than the observed one (scipy binomtest's rule, relative tolerance 1e-7)."""
    d = pmf[k]
    return float(min(1.0, pmf[pmf <= d * (1 + 1e-7)].sum()))


def b3_point(sets2):
    used, left_out, hits, ps, rows = 0, [], 0, [], []
    for s in sets2:
        x, L = s["x"], s["L"]
        xmax = x.max()
        Lmin = L.min()
        tie_T = int((x == xmax).sum()) > 1
        tie_L = int((L == Lmin).sum()) > 1
        if tie_T or tie_L:
            left_out.append({"set": s["id"], "tie_T": tie_T, "tie_L": tie_L})
            continue
        hit = int(np.argmax(x)) == int(np.argmin(L))
        used += 1
        hits += int(hit)
        ps.append(1.0 / s["n"])
        rows.append({"set": s["id"], "n": s["n"], "warmest_first_to_die": bool(hit)})
    pmf = poisson_binomial_pmf(ps)
    p = exact_two_sided(pmf, hits)
    # cross-check for equal set sizes with scipy's binomtest
    sizes = sorted(set(int(round(1 / q)) for q in ps))
    bt = None
    if len(sizes) == 1 and used > 0:
        bt = float(stats.binomtest(hits, used, 1.0 / sizes[0]).pvalue)
    return {"sets_used": used, "sets_left_out": len(left_out), "left_out": left_out, "hits": hits,
            "share": hits / used if used else None, "expected_hits": float(sum(ps)), "chance_share": float(sum(ps) / used) if used else None,
            "p_exact_two_sided": p, "scipy_binomtest_check": bt, "per_set": rows}


def c1_point(sets_all):
    """Set means of ln L on set means of T, weights n, one constant per batch."""
    by_b = {}
    for s in sets_all:
        by_b.setdefault(s["label"], []).append(s)
    Sxx = Sxy = Syy = 0.0
    xs, ys, ns, bs = [], [], [], []
    for lab, ss in by_b.items():
        n = np.array([s["n"] for s in ss], dtype=float)
        xb = np.array([s["x"].mean() for s in ss])
        yb = np.array([s["y"].mean() for s in ss])
        mx, my = (n @ xb) / n.sum(), (n @ yb) / n.sum()
        dx, dy = xb - mx, yb - my
        Sxx += float(n @ (dx * dx)); Sxy += float(n @ (dx * dy)); Syy += float(n @ (dy * dy))
        xs += list(dx); ys += list(dy); ns += list(n); bs += [lab] * len(ss)
    b = Sxy / Sxx
    # explicit weighted lstsq check with batch dummies
    labs = sorted(by_b, key=lambda l: BATCH_ORDER[l])
    rows, yv, wv = [], [], []
    for s in sets_all:
        r = np.zeros(len(labs) + 1); r[labs.index(s["label"])] = 1.0; r[-1] = s["x"].mean()
        rows.append(r); yv.append(s["y"].mean()); wv.append(s["n"])
    A, yv, wv = np.array(rows), np.array(yv), np.sqrt(np.array(wv, dtype=float))
    coef, *_ = np.linalg.lstsq(A * wv[:, None], yv * wv, rcond=None)
    return {"slope": b, "s_b": -b, "sets": len(sets_all), "batches": len(by_b), "cells": int(sum(s["n"] for s in sets_all)),
            "Vx": Sxx, "Sxy": Sxy, "Vy": Syy, "lstsq_slope": float(coef[-1])}


def c3_from(sw, Vx, Sxy, Vy):
    share = sw * sw * Vx / Vy
    resid = 1.0 - (Vy + 2.0 * sw * Sxy + sw * sw * Vx) / Vy
    return share, resid


# ---------------------------------------------------------------------------------------------------------------
# The bootstrap over sets (D93 (6)), vectorised through count matrices
# ---------------------------------------------------------------------------------------------------------------
def count_matrix(sets_all, rng, D=D_BOOT):
    """Strata: each batch's sets of two or more, and its single-cell sets, apart; within each stratum as many sets as
    it has, with replacement. Returns C (D x nsets) of how many times each set was drawn."""
    nset = len(sets_all)
    C = np.zeros((D, nset), dtype=np.int64)
    strata = []
    for lab, _, _ in BATCHES:
        for kind in ("2plus", "single"):
            idx = [j for j, s in enumerate(sets_all) if s["label"] == lab and ((s["n"] >= 2) if kind == "2plus" else (s["n"] == 1))]
            if idx:
                strata.append((lab, kind, idx))
    for lab, kind, idx in strata:
        m = len(idx)
        draw = rng.integers(0, m, size=(D, m))
        cnt = np.zeros((D, m), dtype=np.int64)
        rows = np.repeat(np.arange(D), m)
        np.add.at(cnt, (rows, draw.ravel()), 1)
        C[:, idx] = cnt
    return C, [(l, k, len(i)) for l, k, i in strata]


def boot_stats(sets_all, C, r_oe=None):
    """Every B and C quantity for each draw."""
    W = [within_sums(s) for s in sets_all]
    is2 = np.array([s["n"] >= 2 for s in sets_all], dtype=float)
    n = np.array([s["n"] for s in sets_all], dtype=float)
    g = lambda k: np.array([w[k] for w in W], dtype=float)
    Cf = C.astype(float)
    Sxx, Sxy, Syy = Cf @ (g("Sxx") * is2), Cf @ (g("Sxy") * is2), Cf @ (g("Syy") * is2)
    b1 = Sxy / Sxx
    N = Cf @ (n * is2)
    J = Cf @ is2
    rss = Syy - b1 * Sxy
    sdc = np.sqrt(np.maximum(rss, 0) / (N - J - 1))
    bp = (Cf @ (g("Pdx") * is2)) / (Cf @ (g("Pxx") * is2))
    # B2
    okR = np.array([w["R_n"] >= 2 for w in W], dtype=float) * is2
    rxx, rxz, rzz = Cf @ (g("R_Sxx") * okR), Cf @ (g("R_Sxz") * okR), Cf @ (g("R_Szz") * okR)
    rxy, rzy = Cf @ (g("R_Sxy") * okR), Cf @ (g("R_Szy") * okR)
    with np.errstate(invalid="ignore", divide="ignore"):
        det = rxx * rzz - rxz * rxz
        bT = (rzz * rxy - rxz * rzy) / det
        bR = (rxx * rzy - rxz * rxy) / det
        bT_only = rxy / rxx
    # C1 by batch, centred per batch for precision
    Vx = np.zeros(C.shape[0]); Sxyb = np.zeros(C.shape[0]); Vy = np.zeros(C.shape[0])
    for lab in sorted(set(s["label"] for s in sets_all), key=lambda l: BATCH_ORDER[l]):
        idx = [j for j, s in enumerate(sets_all) if s["label"] == lab]
        xb = np.array([sets_all[j]["x"].mean() for j in idx]); yb = np.array([sets_all[j]["y"].mean() for j in idx])
        xb = xb - xb.mean(); yb = yb - yb.mean()
        w = Cf[:, idx] * n[idx]
        Wt = w.sum(1)
        sx, sy = w @ xb, w @ yb
        Vx += w @ (xb * xb) - sx * sx / Wt
        Sxyb += w @ (xb * yb) - sx * sy / Wt
        Vy += w @ (yb * yb) - sy * sy / Wt
    bb = Sxyb / Vx
    sw = -b1
    sb = -bb
    share = sw * sw * Vx / Vy
    resid = 1.0 - (Vy + 2.0 * sw * Sxyb + sw * sw * Vx) / Vy
    out = {"s_w": sw, "s_w_pairs": -bp, "sd_c_hat": sdc, "sqrt2_sd_c_hat": math.sqrt(2) * sdc,
           "b2_lever_T_R_held": -bT, "b2_dlnL_dlnR": bR, "b2_s_w_same_cells": -bT_only,
           "s_b": sb, "s_b_minus_s_w": sb - sw, "c3_share": share, "c3_resid": resid}
    return out


def self_check_bootstrap(sets_all, C, bs, ndraw=25):
    """For the first draws, rebuild the resampled data set explicitly (every drawn copy its own set) and refit B1, B2 and
    C1 by dummy-variable least squares; returns the largest relative differences from the count-matrix results."""
    worst = {"s_w": 0.0, "s_b": 0.0, "b2_lever_T_R_held": 0.0}
    for d in range(min(ndraw, C.shape[0])):
        drawn = []
        for j, cnt in enumerate(C[d]):
            drawn += [sets_all[j]] * int(cnt)
        s2 = [s for s in drawn if s["n"] >= 2]
        b = b1_lstsq_check(s2)
        worst["s_w"] = max(worst["s_w"], abs(-b - bs["s_w"][d]) / max(1e-300, abs(bs["s_w"][d])))
        c1 = c1_point(drawn)
        worst["s_b"] = max(worst["s_b"], abs(c1["lstsq_slope"] * -1 - bs["s_b"][d]) / max(1e-300, abs(bs["s_b"][d])))
        b2 = b2_point(s2)
        if b2 is not None and np.isfinite(bs["b2_lever_T_R_held"][d]):
            worst["b2_lever_T_R_held"] = max(worst["b2_lever_T_R_held"],
                                            abs(-b2["lstsq_slope_T"] - bs["b2_lever_T_R_held"][d]) / abs(bs["b2_lever_T_R_held"][d]))
    return worst


def percentile_ci(v):
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return [None, None], 0
    lo, hi = np.percentile(v, [2.5, 97.5])
    return [float(lo), float(hi)], int(len(v))


def mc_sd_orderstat(v, q):
    """D77's measure: half the distance between the order statistics sqrt(D q (1 - q)) ranks either side of the
    percentile's rank (linear interpolation)."""
    v = np.sort(np.asarray(v, dtype=float)[np.isfinite(v)])
    D = len(v)
    if D == 0:
        return None
    r = q * (D - 1)
    h = math.sqrt(D * q * (1 - q))
    def at(rank):
        rank = min(max(rank, 0), D - 1)
        lo = int(math.floor(rank)); hi = min(lo + 1, D - 1)
        return v[lo] + (rank - lo) * (v[hi] - v[lo])
    return float(0.5 * (at(r + h) - at(r - h)))


# ---------------------------------------------------------------------------------------------------------------
# A0.5 reliability, A0.7 channel
# ---------------------------------------------------------------------------------------------------------------
def reliability(cells):
    out = {}
    for scope, sub in (("pooled", cells), ("paper1", [c for c in cells if c["paper"] == 1]),
                       ("paper2", [c for c in cells if c["paper"] == 2])):
        res = {}
        for name, (ka, kb) in (("r_oe", ("T_odd", "T_even")), ("r_h", ("T_h1", "T_h2"))):
            da, db = [], []
            for sid, mem in sets_of(sub).items():
                if len(mem) < 2:
                    continue
                a = np.array([c[ka] for c in mem]); b = np.array([c[kb] for c in mem])
                da += list(a - a.mean()); db += list(b - b.mean())
            da, db = np.array(da), np.array(db)
            r = float(da @ db / math.sqrt((da @ da) * (db @ db)))
            res[name] = r
            res[name + "_spearman_brown"] = 2 * r / (1 + r)
            res[name + "_cells"] = int(len(da))
            res[name + "_pearson_centered"] = float(np.corrcoef(da, db)[0, 1])
        out[scope] = res
    return out


def channel_test(cells, rng, nperm=10000):
    devs = {}
    for lab, _, _ in BATCHES:
        sub = [c for c in cells if c["label"] == lab]
        for sid, mem in sets_of(sub).items():
            if len(mem) < 2:
                continue
            T = np.array([c["T"] for c in mem])
            for c, d in zip(mem, T - T.mean()):
                devs.setdefault(lab, {})[c["channel"]] = float(d)
    labs = [b[0] for b in BATCHES if b[0] in devs]
    pos = {lab: {ch: i for i, ch in enumerate(sorted(devs[lab], key=lambda s: int(s)))} for lab in labs}
    vals = {lab: np.array([devs[lab][ch] for ch in sorted(devs[lab], key=lambda s: int(s))]) for lab in labs}
    pairs = []
    for ia in range(len(labs)):
        for ib in range(ia + 1, len(labs)):
            a, b = labs[ia], labs[ib]
            for ch in sorted(set(devs[a]) & set(devs[b]), key=lambda s: int(s)):
                pairs.append((a, pos[a][ch], b, pos[b][ch], ch))
    by_batchpair = {}
    for a, _, b, _, _ in pairs:
        by_batchpair[f"{DATE[a]} / {DATE[b]}"] = by_batchpair.get(f"{DATE[a]} / {DATE[b]}", 0) + 1

    def corr_of(v):
        x = np.array([v[a][i] for a, i, b, k, _ in pairs]); y = np.array([v[b][k] for a, i, b, k, _ in pairs])
        return float(np.corrcoef(x, y)[0, 1]), float(x @ y / math.sqrt((x @ x) * (y @ y)))

    r_obs, r_obs_unc = corr_of(vals)
    # vectorised permutations: index arrays per batch
    pa = np.array([[labs.index(a), i, labs.index(b), k] for a, i, b, k, _ in pairs])
    perms = {lab: np.argsort(rng.random((nperm, len(vals[lab]))), axis=1) for lab in labs}
    X = np.empty((nperm, len(pairs))); Y = np.empty((nperm, len(pairs)))
    for j, (ia_, i, ib_, k) in enumerate(pa):
        la, lb = labs[ia_], labs[ib_]
        X[:, j] = vals[la][perms[la][:, i]]
        Y[:, j] = vals[lb][perms[lb][:, k]]
    Xc = X - X.mean(1, keepdims=True); Yc = Y - Y.mean(1, keepdims=True)
    rp = (Xc * Yc).sum(1) / np.sqrt((Xc * Xc).sum(1) * (Yc * Yc).sum(1))
    ge = int((np.abs(rp) >= abs(r_obs) - 1e-15).sum())
    return {"pairs": len(pairs), "pairs_by_batch_pair": by_batchpair, "channels_with_pairs": len(set(p[4] for p in pairs)),
            "correlation": r_obs, "correlation_uncentred": r_obs_unc, "permutations": nperm,
            "p_two_sided": ge / nperm, "p_two_sided_plus_one": (ge + 1) / (nperm + 1),
            "carries_temperature": bool(r_obs > 0 and ge / nperm < 0.05), "perm_r_mean": float(rp.mean()), "perm_r_sd": float(rp.std())}


# ---------------------------------------------------------------------------------------------------------------
# A0.8 power: own noncentral-t evaluation, and simulation
# ---------------------------------------------------------------------------------------------------------------
def nct_two_sided_power(df, nc, tc):
    """P(|T'| > tc), T' = (Z + nc) / sqrt(V / df), V chi-square(df): the integral over v of
    [Phi(nc - tc sqrt(v/df)) + Phi(-nc - tc sqrt(v/df))] times the chi-square density, by adaptive quadrature."""
    lg = special.gammaln(df / 2.0)
    c0 = -(df / 2.0) * math.log(2.0) - lg

    def f(v):
        if v <= 0:
            return 0.0
        w = tc * math.sqrt(v / df)
        return (special.ndtr(nc - w) + special.ndtr(-nc - w)) * math.exp((df / 2.0 - 1.0) * math.log(v) - v / 2.0 + c0)

    sd = math.sqrt(2.0 * df)
    pts = sorted(set([0.0, max(df - 2.0, 0.0), max(df - 8 * sd, 0.0), df + 8 * sd, df + 80 * sd]))
    tot = 0.0
    for a, b in zip(pts[:-1], pts[1:]):
        if b > a:
            v, _ = integrate.quad(f, a, b, epsabs=1e-15, epsrel=1e-13, limit=400)
            tot += v
    return tot


def stage4_shape(m, q):
    """D61's registered shape: |e| log-normal with mean m and 95th percentile q (smaller root)."""
    z = float(stats.norm.ppf(0.95))
    L = math.log(q / m)
    disc = z * z - 2.0 * L
    sig = z - math.sqrt(disc)
    mu = math.log(m) - 0.5 * sig * sig
    return {"sigma": sig, "mu": mu, "sd_e": math.exp(mu + sig * sig), "z95": z,
            "mean_check": math.exp(mu + 0.5 * sig * sig), "p95_check": math.exp(mu + z * sig)}


def calc_detectable(df, S_TT, sd_c):
    tc = float(stats.t.ppf(1.0 - ALPHA / 2.0, df))
    f = lambda s: nct_two_sided_power(df, s * math.sqrt(S_TT) / sd_c, tc) - POWER
    hi = 1e-3
    while f(hi) < 0:
        hi *= 2
    s = optimize.brentq(f, 0.0, hi, xtol=1e-16, rtol=1e-12, maxiter=500)
    # cross-check of my integral against scipy's noncentral t at the root
    nc = s * math.sqrt(S_TT) / sd_c
    chk = float(stats.nct.sf(tc, df, nc) + stats.nct.cdf(-tc, df, nc))
    return s, tc, chk


def fe_design(sets2):
    """Within-set deviations of T and the set index of each cell."""
    dT, sidx = [], []
    for j, s in enumerate(sets2):
        dT += list(s["x"] - s["x"].mean())
        sidx += [j] * s["n"]
    return np.array(dT), np.array(sidx)


def sim_detectable(sets2, eps, df):
    """eps: (D, N) cell noises; FE test at each s; smallest s at which simulated power reaches 0.80, by bisection to
    1e-6 relative, and the exact smallest s by sweeping the draws' critical points."""
    dT, sidx = fe_design(sets2)
    J = len(sets2)
    S_TT = float(dT @ dT)
    # within-set demeaning of the noise
    D, N = eps.shape
    G = np.zeros((N, J))
    G[np.arange(N), sidx] = 1.0
    counts = G.sum(0)
    eps_dm = eps - ((eps @ G) / counts)[:, sidx]
    b0 = eps_dm @ dT / S_TT
    rss = (eps_dm * eps_dm).sum(1) - b0 * b0 * S_TT
    se = np.sqrt(rss / df / S_TT)
    tc = float(stats.t.ppf(1.0 - ALPHA / 2.0, df))
    c = tc * se

    def pw(s):
        return float(np.mean(np.abs(b0 - s) > c))

    lo, hi = 0.0, 0.01
    while pw(hi) < POWER:
        lo, hi = hi, hi * 2
    while (hi - lo) > 1e-6 * hi:
        mid = 0.5 * (lo + hi)
        if pw(mid) >= POWER:
            hi = mid
        else:
            lo = mid
    # exact sweep for s >= 0
    need = int(math.ceil(POWER * D - 1e-9))
    lo_end = b0 - c        # rejects for s < lo_end
    hi_start = b0 + c      # rejects for s > hi_start
    ev = []
    base = int(np.sum(lo_end > 0)) + int(np.sum(hi_start < 0))
    for le, hs in zip(lo_end, hi_start):
        if le > 0:
            ev.append((le, -1))
        if hs >= 0:
            ev.append((hs, +1))
    ev.sort()
    cnt = base
    exact = None
    if cnt >= need:
        exact = 0.0
    else:
        for x, dlt in ev:
            cnt += dlt
            if cnt >= need and dlt > 0:
                exact = x
                break
    return hi, exact, pw(hi)


def power_block(cells_scope, floors, shapes, rep_seeds, sim=True, scope_id=0):
    sets2 = [s for s in set_list(cells_scope, "T") if s["n"] >= 2]
    N = sum(s["n"] for s in sets2)
    J = len(sets2)
    df = N - J - 1
    dT, sidx = fe_design(sets2)
    S_TT = float(dT @ dT)
    # pairs form: within-set pairs
    Pxx, npairs = 0.0, 0
    Fnum, Fden = 0.0, 0.0
    for s in sets2:
        w = within_sums(s)
        Pxx += w["Pxx"]; npairs += w["npairs"]
        Fnum += s["n"] ** 2 * w["Sxx"]; Fden += 2 * s["n"] * w["Sxx"]
    F = Fnum / Fden
    dT_rms_pairs = math.sqrt(Pxx / npairs)
    out = {"N": N, "J": J, "df": df, "S_TT": S_TT, "dT_rms_over_N": math.sqrt(S_TT / N),
           "dT_rms_over_N_minus_J": math.sqrt(S_TT / (N - J)), "pairs": npairs, "dT_rms_pairs": dT_rms_pairs, "F_shared": F,
           "floors": {}}
    for fname, f in floors.items():
        sd_e = f["sd_e"]
        sd_c = sd_e / math.sqrt(2.0)
        s_calc, tc, chk = calc_detectable(df, S_TT, sd_c)
        r = {"sd_e": sd_e, "sd_c": sd_c, "tc": tc, "detectable_calc": s_calc,
             "power_at_root_scipy_nct_check": chk}
        # unscored: Stage 5's pairs form and shared-cell line; Stage 4's normal shape
        r["pairs_form_detectable"], _, _ = calc_detectable(npairs - 1, Pxx, sd_e)
        r["pairs_form_shared_cell_detectable"], _, _ = calc_detectable(npairs - 1, Pxx, sd_e * math.sqrt(F))
        sd_e_norm = f["mean"] * math.sqrt(math.pi / 2.0)
        r["normal_shape_sd_e"] = sd_e_norm
        r["normal_shape_detectable"], _, _ = calc_detectable(df, S_TT, sd_e_norm / math.sqrt(2.0))
        if sim:
            sims_ln, sims_nm, ex_ln, ex_nm = [], [], [], []
            for k, seed in enumerate(rep_seeds):
                rng = rng_for(80, scope_id, seed, 1 if fname == "Geslin" else 2)
                mag = np.exp(shapes[fname]["mu"] + shapes[fname]["sigma"] * rng.standard_normal((10000, N)))
                sign = rng.integers(0, 2, size=(10000, N)).astype(float) * 2.0 - 1.0
                eps = sign * mag / math.sqrt(2.0)
                v, ex, _ = sim_detectable(sets2, eps, df)
                sims_ln.append(v); ex_ln.append(ex)
                rng2 = rng_for(90, scope_id, seed, 1 if fname == "Geslin" else 2)
                eps2 = sd_c * rng2.standard_normal((10000, N))
                v2, ex2, _ = sim_detectable(sets2, eps2, df)
                sims_nm.append(v2); ex_nm.append(ex2)
            r["sim_lognormal_detectable"] = sims_ln[0]
            r["sim_lognormal_exact_sweep"] = ex_ln[0]
            r["sim_normal_detectable"] = sims_nm[0]
            r["sim_normal_exact_sweep"] = ex_nm[0]
            r["v5_rel_diff"] = (sims_ln[0] - s_calc) / s_calc
            r["v5_pass"] = abs(sims_ln[0] - s_calc) <= 0.05 * s_calc
            r["normal_check_rel_diff"] = (sims_nm[0] - s_calc) / s_calc
            r["replicates"] = len(rep_seeds)
            r["sim_lognormal_replicates"] = sims_ln
            r["sim_normal_replicates"] = sims_nm
            r["sim_lognormal_mc_sd"] = float(np.std(sims_ln, ddof=1)) if len(sims_ln) > 1 else None
            r["sim_normal_mc_sd"] = float(np.std(sims_nm, ddof=1)) if len(sims_nm) > 1 else None
            r["sim_lognormal_replicate_mean"] = float(np.mean(sims_ln))
            r["sim_normal_replicate_mean"] = float(np.mean(sims_nm))
        out["floors"][fname] = r
    return out


# ---------------------------------------------------------------------------------------------------------------
# One analysis variant: points and intervals for B and C
# ---------------------------------------------------------------------------------------------------------------
def scope_cells(cells, scope):
    if scope == "pooled":
        return list(cells)
    if scope == "paper1":
        return [c for c in cells if c["paper"] == 1]
    if scope == "paper2":
        return [c for c in cells if c["paper"] == 2]
    raise ValueError(scope)


def analysis(cells, scope, xkey="T", zkey="lnR", c1_2plus_only=False, seed_key=None, replicates=0, r_oe=None):
    sub = scope_cells(cells, scope)
    sets_all = set_list(sub, xkey, zkey)
    sets2 = [s for s in sets_all if s["n"] >= 2]
    pt = {}
    b1 = b1_point(sets2)
    b1["lstsq_check_slope"] = b1_lstsq_check(sets2)
    pt["B1"] = b1
    pt["B1_pairs"] = pairs_point(sets2)
    pt["B2"] = b2_point(sets2)
    pt["B3"] = b3_point(sets2)
    c_sets = sets2 if c1_2plus_only else sets_all
    c1 = c1_point(c_sets)
    pt["C1"] = c1
    pt["C2"] = {"s_b_minus_s_w": c1["s_b"] - b1["s_w"]}
    share, resid = c3_from(b1["s_w"], c1["Vx"], c1["Sxy"], c1["Vy"])
    pt["C3"] = {"share": share, "resid_form": resid, "between_policy_variance_sum_n_y2": c1["Vy"], "sum_n_x2": c1["Vx"],
                "sum_n_xy": c1["Sxy"]}
    # bootstrap
    boot = None
    if seed_key is not None:
        rng = rng_for(*seed_key)
        C, strata = count_matrix(c_sets, rng)
        bs = boot_stats(c_sets, C)
        selfcheck = self_check_bootstrap(c_sets, C, bs)
        ci = {}
        for k, v in bs.items():
            (lo, hi), nfin = percentile_ci(v)
            ci[k] = {"ci95": [lo, hi], "finite_draws": nfin,
                     "mc_sd_orderstat": [mc_sd_orderstat(v, 0.025), mc_sd_orderstat(v, 0.975)]}
        reps = {}
        if replicates:
            ends = {k: [] for k in bs}
            for rep in range(replicates):
                rr = rng_for(*seed_key, 1000 + rep)
                Cr, _ = count_matrix(c_sets, rr)
                br = boot_stats(c_sets, Cr)
                for k, v in br.items():
                    ends[k].append(percentile_ci(v)[0])
            for k, e in ends.items():
                if any(x[0] is None for x in e):
                    reps[k] = {"replicates": replicates, "sd_lo": None, "sd_hi": None, "mean_lo": None, "mean_hi": None}
                    continue
                e = np.array(e, dtype=float)
                reps[k] = {"replicates": replicates, "sd_lo": float(np.std(e[:, 0], ddof=1)), "sd_hi": float(np.std(e[:, 1], ddof=1)),
                           "mean_lo": float(np.mean(e[:, 0])), "mean_hi": float(np.mean(e[:, 1]))}
        boot = {"draws": D_BOOT, "seed_key": [SEED_BASE, *seed_key], "generator": "numpy Generator(Philox(SeedSequence([5551117, *key])))",
                "strata": [{"batch": DATE[l], "kind": k, "sets": m} for l, k, m in strata], "ci": ci, "replicate_runs": reps,
                "self_check_first_25_draws_max_rel_diff": selfcheck}
    res = {"scope": scope, "xkey": xkey, "zkey": zkey, "c1_2plus_only": c1_2plus_only,
           "counts": {"cells": len(sub), "sets": len(sets_all), "sets_2plus": len(sets2), "cells_2plus": int(sum(s["n"] for s in sets2))},
           "point": pt, "boot": boot}
    if r_oe is not None:
        res["divided_by_r_oe"] = divide_by_r(res, r_oe)
    return res


def divide_by_r(res, r):
    pt, bt = res["point"], res["boot"]
    out = {"r_oe": r}
    for name, val, ck in (("s_w", pt["B1"]["s_w"], "s_w"), ("s_w_pairs", pt["B1_pairs"]["s_w_pairs"], "s_w_pairs"),
                          ("s_b", pt["C1"]["s_b"], "s_b"), ("s_b_minus_s_w", pt["C2"]["s_b_minus_s_w"], "s_b_minus_s_w"),
                          ("b2_lever_T_R_held", pt["B2"]["lever_T_R_held"] if pt["B2"] else None, "b2_lever_T_R_held"),
                          ("b2_s_w_same_cells", pt["B2"]["s_w_same_cells"] if pt["B2"] else None, "b2_s_w_same_cells")):
        if val is None:
            out[name] = None
            continue
        e = {"point": val / r}
        if bt is not None:
            e["ci95"] = [(x / r if x is not None else None) for x in bt["ci"][ck]["ci95"]]
        out[name] = e
    return out


# ---------------------------------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------------------------------
def main(replicates):
    t0 = time.time()
    rec = {"what": "Independent check of Stage 6, written before any Stage 6 result was opened",
           "python": sys.version.split()[0], "numpy": np.__version__, "scipy": __import__("scipy").__version__}
    entries, files = load_entries()
    rec["inputs_read"] = files
    rec["entries"] = entry_facts(entries)
    rec["swaps"] = apply_swaps(entries)
    rec["carried"] = carry_checks(entries)
    cells, excluded = build_cells(entries)
    rec["excluded_entries"] = [{"batch": DATE[l], "index": i, "channel": entries[(l, i)]["channel"],
                                "policy_readable": entries[(l, i)]["policy_readable"],
                                "stored_cycle_life": fnum(entries[(l, i)]["life"]["cycle_life_stored"])} for l, i in excluded]
    rec["counts"] = count_table(cells)
    left_no_L = [c["key"] for c in cells if not math.isfinite(c["lnL"])]
    rec["cells_without_finite_L"] = left_no_L
    if left_no_L:
        raise SystemExit("a cell without a finite L; D94 (3) handling would be needed")
    # contact rule
    flags = contact_flags(cells)
    rec["contact"] = summarize_contact(cells, flags)
    # reliability
    rel = reliability(cells)
    rec["A0_5_reliability"] = rel
    rec["A0_5_reliability_without_flagged"] = reliability([c for c in cells if not c["flag"]])
    # resistance summary
    rec["A0_6_resistance"] = resistance_summary(cells)
    # channel
    rec["A0_7_channel"] = channel_test(cells, rng_for(70))
    rep_ch = []
    for rep in range(replicates):
        rep_ch.append(channel_test(cells, rng_for(70, 1000 + rep))["p_two_sided"])
    rec["A0_7_channel"]["p_replicates"] = rep_ch
    rec["A0_7_channel"]["p_mc_sd"] = float(np.std(rep_ch, ddof=1)) if len(rep_ch) > 1 else None
    # floors and power
    s4 = json.load(open(os.path.join(ROOT, "results", "stage4_experiment.json")))
    floors, shapes = {}, {}
    for fname in ("Geslin", "KIT"):
        fl = s4["shapes"][fname]["floor"]
        sh = stage4_shape(fl["mean"], fl["p95"])
        shapes[fname] = sh
        floors[fname] = {"mean": fl["mean"], "p95": fl["p95"], "sd_e": sh["sd_e"], "sd_e_archived": s4["shapes"][fname]["lognormal"]["sd_e"],
                         "sigma": sh["sigma"], "mu": sh["mu"], "sigma_archived": s4["shapes"][fname]["lognormal"]["sigma"],
                         "mu_archived": s4["shapes"][fname]["lognormal"]["mu"]}
    gfl = json.load(open(os.path.join(ROOT, "inputs", "floor.json")))
    floors["Geslin"]["floor_json_mean"] = gfl["mean"]; floors["Geslin"]["floor_json_p95"] = gfl["p95"]
    rec["floors"] = floors
    rec["A0_8_power"] = {}
    rep_seeds = list(range(1, max(replicates, 1) + 1))
    for si, scope in enumerate(("pooled", "paper1", "paper2")):
        rec["A0_8_power"][scope] = power_block(scope_cells(cells, scope), floors, shapes, rep_seeds, sim=True, scope_id=si + 1)
    rec["A0_8_power"]["pooled_without_flagged"] = power_block([c for c in cells if not c["flag"]], floors, shapes, rep_seeds[:1], sim=False)
    for scope in ("paper1", "paper2"):
        rec["A0_8_power"][scope + "_without_flagged"] = power_block([c for c in scope_cells(cells, scope) if not c["flag"]], floors, shapes, rep_seeds[:1], sim=False)
    print(f"A0 done {time.time() - t0:.0f} s", flush=True)
    # archived levers
    s3 = read_csv(os.path.join(ROOT, "results", "stage3_lever.csv"))
    sand = {}
    for name, (series, chem, sec) in SANDIA_KEYS.items():
        r = [x for x in s3 if x["series"] == series and x["chemistry"] == chem and x["secant"] == sec][0]
        sand[name] = {"series": series, "secant": sec, "s": float(r["s_per_K"]), "lo": float(r["ci95_lo"]), "hi": float(r["ci95_hi"])}
    rec["sandia_levers"] = sand
    s5 = json.load(open(os.path.join(ROOT, "results", "stage5_paper_numbers.json")))
    s5c1 = {(r["case"], r["floor"]): r["detectable_s_per_K"] for r in s5["C1"]}
    rec["stage5_detectable"] = {"KIT_KIT_floor": s5c1[("KIT", "KIT")], "KIT_Geslin_floor": s5c1[("KIT", "Geslin")]}
    # Part B and C, registered and brackets
    seeds = {}
    variants = {}
    k = 0
    for scope in ("pooled", "paper1", "paper2"):
        k += 1
        variants[f"registered|{scope}"] = analysis(cells, scope, "T", "lnR", seed_key=(1, k), replicates=replicates,
                                                   r_oe=rel[scope]["r_oe"])
    print(f"registered done {time.time() - t0:.0f} s", flush=True)
    br = [("Tmax", "Tmax_m", "lnR", False, False), ("Tavg_2_50", "T_2_50", "lnR", False, False), ("Tmin", "Tmin_m", "lnR", False, False),
          ("without_flagged", "T", "lnR", False, True),
          ("unscored_T_minus_Tmin", "T_minus_Tmin", "lnR", False, False), ("unscored_R_2_10", "T", "lnR_2_10", False, False),
          ("unscored_C1_sets_2plus_only", "T", "lnR", True, False)]
    for bi, (name, xkey, zkey, c2only, noflag) in enumerate(br):
        for si, scope in enumerate(("pooled", "paper1", "paper2")):
            use = [c for c in cells if not c["flag"]] if noflag else cells
            variants[f"{name}|{scope}"] = analysis(use, scope, xkey, zkey, c1_2plus_only=c2only, seed_key=(10 + bi, si + 1),
                                                   replicates=replicates)
    print(f"brackets done {time.time() - t0:.0f} s", flush=True)
    rec["parts_B_C"] = strip_variants(variants)
    rec["B5"] = {"runs": rec["A0_7_channel"]["carries_temperature"],
                 "note": "B5 runs only if A0.7 finds the channel carries a temperature (positive correlation, p < 0.05)."}
    rec["part_D"] = part_d(rec)
    rec["scorecard"] = scorecard(rec)
    rec["cells"] = [cell_row(c) for c in cells]
    rec["sets"] = sets_rows(cells)
    rec["runtime_s"] = time.time() - t0
    write_tables(rec)
    with open(os.path.join(HERE, "record.json"), "w") as f:
        json.dump(clean(rec), f, indent=1)
    print(json.dumps(clean(rec["scorecard"]), indent=1))
    print(f"done {time.time() - t0:.0f} s")


def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, (np.floating,)):
        return fnum(float(x))
    if isinstance(x, float):
        return fnum(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    return x


def strip_variants(variants):
    out = {}
    for k, v in variants.items():
        v = dict(v)
        pt = v["point"]
        if pt.get("B3"):
            pt["B3"] = dict(pt["B3"])
        out[k] = v
    return out


def summarize_contact(cells, flags):
    out = {"flagged_total": int(sum(c["flag"] for c in cells)), "by_batch": {}, "flagged_in_sets_2plus": 0, "rule_counts": {}}
    sets = sets_of(cells)
    for c in cells:
        if c["flag"] and len(sets[c["set"]]) >= 2:
            out["flagged_in_sets_2plus"] += 1
        for r in c["flag_rules"]:
            out["rule_counts"][r] = out["rule_counts"].get(r, 0) + 1
    for lab, date, _ in BATCHES:
        sub = [c for c in cells if c["label"] == lab]
        out["by_batch"][date] = {
            "flagged": [{"channel": c["channel"], "key": c["key"], "rules": c["flag_rules"], "b_ratio": c["b_ratio"]}
                        for c in sorted(sub, key=lambda c: int(c["channel"])) if c["flag"]],
            "b_ratios_sorted": sorted([c["b_ratio"] for c in sub if c["b_ratio"] is not None]),
            "max_step_Tavg_2_100": max(c["max_step_Tavg_2_100"] for c in sub),
            "max_step_Tavg_1_100": max(c["max_step_Tavg_1_100"] for c in sub),
            "longest_flat_run": max(c["longest_flat_run_2_100"] for c in sub),
            "min_Tmin_2_100": min(c["min_Tmin_2_100"] for c in sub),
            "median_Tmin_range": [min(c["median_Tmin_2_100"] for c in sub), max(c["median_Tmin_2_100"] for c in sub)],
            "median_Tmin_1_100_range": [min(c["median_Tmin_1_100"] for c in sub), max(c["median_Tmin_1_100"] for c in sub)],
            "mean_Tmin_2_100_range": [min(c["Tmin_m"] for c in sub), max(c["Tmin_m"] for c in sub)],
            "median_swing_2_100_range": [min(c["swing_med_2_100"] for c in sub), max(c["swing_med_2_100"] for c in sub)],
            "swing_collapse": [{"channel": c["channel"], "med_2_11": c["swing_med_2_11"], "med_91_100": c["swing_med_91_100"]}
                               for c in sub if c["swing_med_91_100"] < 0.5 * c["swing_med_2_11"]],
        }
    # D88's descriptive numbers: within-set correlation of T deviations with swing deviations
    for scope, sub in (("paper2", [c for c in cells if c["paper"] == 2]), ("paper1", [c for c in cells if c["paper"] == 1])):
        a, bm, bn = [], [], []
        for sid, mem in sets_of(sub).items():
            if len(mem) < 2:
                continue
            T = np.array([c["T"] for c in mem]); Sm = np.array([c["swing_med_2_100"] for c in mem])
            Sn = np.array([c["swing_mean_2_100"] for c in mem])
            a += list(T - T.mean()); bm += list(Sm - Sm.mean()); bn += list(Sn - Sn.mean())
        a, bm, bn = np.array(a), np.array(bm), np.array(bn)
        out[f"corr_dev_T_dev_swing_{scope}"] = {
            "median_swing_uncentred": float(a @ bm / math.sqrt((a @ a) * (bm @ bm))),
            "mean_swing_uncentred": float(a @ bn / math.sqrt((a @ a) * (bn @ bn))),
            "median_swing_pearson": float(np.corrcoef(a, bm)[0, 1]), "mean_swing_pearson": float(np.corrcoef(a, bn)[0, 1])}
    b9 = [c for c in cells if c["label"] == "b9"]
    rat = sorted(c["b_ratio"] for c in b9)
    out["b9_ratio_gap"] = {"below_half_max": max(r for r in rat if r < 0.5), "above_half_min": min(r for r in rat if r >= 0.5)}
    meds = sorted(c["swing_med_2_100"] for c in b9)
    out["b9_median_swings_sorted"] = meds
    out["lowest_Tmin_all_cells"] = min(c["min_Tmin_2_100"] for c in cells)
    out["largest_step_all_cells"] = max(c["max_step_Tavg_2_100"] for c in cells)
    return out


def resistance_summary(cells):
    out = {}
    for lab, date, _ in BATCHES:
        sub = [c for c in cells if c["label"] == lab]
        hist = {}
        for c in sub:
            hist[c["R_valid_cycles"]] = hist.get(c["R_valid_cycles"], 0) + 1
        Rs = [c["R"] for c in sub if math.isfinite(c["R"])]
        out[date] = {"valid_cycles_histogram": {str(k): hist[k] for k in sorted(hist)},
                     "R_min": min(Rs) if Rs else None, "R_max": max(Rs) if Rs else None, "cells_with_R": len(Rs)}
    Rs = [c["R"] for c in cells if math.isfinite(c["R"])]
    out["R_range_all"] = [min(Rs), max(Rs)]
    return out


def cell_row(c):
    keys = ["key", "label", "batch", "paper", "index", "policy_readable", "channel", "set", "swapped_with", "T", "Tmax_m", "T_2_50",
            "Tmin_m", "T_minus_Tmin", "T_odd", "T_even", "T_h1", "T_h2", "R", "R_valid_cycles", "R_2_10", "lnR", "lnR_2_10",
            "swing_med_2_100", "swing_med_2_11", "swing_med_91_100", "b_ratio", "max_step_Tavg_2_100", "longest_flat_run_2_100",
            "min_Tmin_2_100", "flag", "flag_rules", "L_stored", "L_added", "L", "lnL", "matlab_label_literal",
            "matlab_label_skip_zero_first"]
    return {k: c[k] for k in keys}


def sets_rows(cells):
    rows = []
    for s in set_list(cells, "T", "lnR"):
        mem = [c for c in cells if c["set"] == s["id"]]
        rows.append({"set": s["id"], "batch": s["batch"], "paper": s["paper"], "policy_readable": s["policy_readable"], "n": s["n"],
                     "keys": s["keys"], "channels": [c["channel"] for c in mem], "mean_T": float(s["x"].mean()),
                     "mean_lnL": float(s["y"].mean()), "n_flagged": int(sum(c["flag"] for c in mem)),
                     "L": [float(v) for v in s["L"]]})
    return rows


def write_tables(rec):
    cols = list(rec["cells"][0].keys())
    with open(os.path.join(HERE, "cells.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rec["cells"]:
            w.writerow({k: (repr(v) if isinstance(v, float) else v) for k, v in r.items()})
    cols = ["set", "batch", "paper", "policy_readable", "n", "n_flagged", "mean_T", "mean_lnL", "channels", "keys", "L"]
    with open(os.path.join(HERE, "sets.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rec["sets"]:
            w.writerow({k: (" ".join(map(str, r[k])) if isinstance(r[k], list) else (repr(r[k]) if isinstance(r[k], float) else r[k])) for k in cols})


def part_d(rec):
    pb = rec["parts_B_C"]
    row = {}
    for scope in ("pooled", "paper1", "paper2"):
        v = pb[f"registered|{scope}"]
        row[scope] = {"s_w": v["point"]["B1"]["s_w"], "s_w_ci95": v["boot"]["ci"]["s_w"]["ci95"],
                      "s_b": v["point"]["C1"]["s_b"], "s_b_ci95": v["boot"]["ci"]["s_b"]["ci95"],
                      "brackets_s_w": {b: pb[f"{b}|{scope}"]["point"]["B1"]["s_w"] for b in ("Tmax", "Tavg_2_50", "Tmin", "without_flagged")},
                      "brackets_s_b": {b: pb[f"{b}|{scope}"]["point"]["C1"]["s_b"] for b in ("Tmax", "Tavg_2_50", "Tmin", "without_flagged")},
                      "s_w_over_r_oe": v["divided_by_r_oe"]["s_w"]["point"], "s_b_over_r_oe": v["divided_by_r_oe"]["s_b"]["point"],
                      "detectable_KIT": rec["A0_8_power"][scope]["floors"]["KIT"]["detectable_calc"],
                      "detectable_Geslin": rec["A0_8_power"][scope]["floors"]["Geslin"]["detectable_calc"]}
    return {"this_stage": row, "sandia": rec["sandia_levers"], "stage5_detectable": rec["stage5_detectable"],
            "protocols": "Sandia: 0.5C charge, 0.5C to 3C discharge (levers at 1C and 2C), 15, 25 and 35 C air. These cells: "
                         "3.6C to 8C to 80 % then 1C CC-CV, 4C discharge to 2.0 V, 30 C chamber (SOURCES.md)."}


def scorecard(rec):
    pb = rec["parts_B_C"]
    reg = pb["registered|pooled"]
    sand = rec["sandia_levers"]
    pw = rec["A0_8_power"]["pooled"]["floors"]["KIT"]
    out = {}
    thr = sand["lfp_1c_15_25"]["s"]
    p1_calc = pw["detectable_calc"]
    p1_sim = pw.get("sim_lognormal_detectable")
    p1 = p1_calc < thr
    p1_scorable = True
    if not pw.get("v5_pass", True) and ((p1_sim < thr) != p1):
        p1_scorable = False
    out["P1"] = {"quantity": p1_calc, "threshold": thr, "holds": p1 if p1_scorable else None, "scorable": p1_scorable,
                 "v5_pass": pw.get("v5_pass"), "sim": p1_sim}
    sw = reg["point"]["B1"]["s_w"]; swci = reg["boot"]["ci"]["s_w"]["ci95"]
    out["P2"] = {"s_w": sw, "ci95": swci, "holds": bool(sw > 0 and swci[0] > 0)}
    b2 = reg["point"]["B2"]
    out["P3"] = {"lever_T_R_held": b2["lever_T_R_held"], "s_w_same_cells": b2["s_w_same_cells"],
                 "holds": bool(np.sign(b2["lever_T_R_held"]) == np.sign(b2["s_w_same_cells"]) and
                               abs(b2["lever_T_R_held"]) >= 0.5 * abs(b2["s_w_same_cells"])),
                 "ratio": b2["lever_T_R_held"] / b2["s_w_same_cells"]}
    sb = reg["point"]["C1"]["s_b"]; sbci = reg["boot"]["ci"]["s_b"]["ci95"]
    out["P4"] = {"s_b": sb, "ci95": sbci, "holds": bool(sb > 0 and sbci[0] > 0)}
    dci = reg["boot"]["ci"]["s_b_minus_s_w"]["ci95"]
    out["P5"] = {"s_b_minus_s_w": reg["point"]["C2"]["s_b_minus_s_w"], "ci95": dci, "holds": bool(dci[0] > 0)}
    roe = rec["A0_5_reliability"]["pooled"]["r_oe"]
    out["P6"] = {"r_oe": roe, "holds": bool(roe > 0.8)}
    lo7, hi7 = sand["lfp_2c_25_35"]["lo"], sand["lfp_2c_25_35"]["hi"]
    out["P7"] = {"s_w": sw, "interval": [lo7, hi7], "holds": bool(lo7 <= sw <= hi7)}
    q8 = reg["point"]["B1"]["sqrt2_sd_c_hat"]
    g, kf = rec["floors"]["Geslin"]["sd_e"], rec["floors"]["KIT"]["sd_e"]
    out["P8"] = {"sqrt2_sd_c_hat": q8, "interval": [g, kf], "interval_archived": [rec["floors"]["Geslin"]["sd_e_archived"], rec["floors"]["KIT"]["sd_e_archived"]],
                 "holds": bool(g <= q8 <= kf)}
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--replicates", type=int, default=20)
    a = ap.parse_args()
    main(a.replicates)
