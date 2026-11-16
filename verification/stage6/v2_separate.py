#!/usr/bin/env python3
"""V2 of Stage 6 (DECISIONS.md D96 (1)): Parts B and C by separate code written from D86 to D95.

Imports nothing of src/stage6. Reads the committed extracts (`inputs/matr/a0/`) and lifetimes
(`inputs/matr/stage6_lifetimes.csv`) itself; applies the authors' exclusions and carry-over (D86) and the site's swapped
thermocouples (D87 (3)); forms each cell's T and ln R (D87, D90); fits every slope by numpy.linalg.lstsq on explicit
design matrices with one dummy column per set (or per batch for C1); and draws its own bootstrap with Python's
random.Random (seed 20261117), refitting each draw, a set drawn twice entering as two sets. Passes where every point
value agrees with results/stage6_within.json and results/stage6_between.json to 1e-10 relative and every interval end
lies within 4 sqrt(2) Monte Carlo standard deviations of the pipeline's end.

Usage: python3 -I verification/stage6/v2_separate.py [draws]
"""
import csv, json, math, os, random, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
EXT = os.path.join(ROOT, "inputs", "matr", "a0")
FILES = {"b1": ("2017-05-12", "Severson"), "b2": ("2017-06-30", "Severson"), "b3": ("2018-04-12", "Severson"),
         "b9": ("2019-01-24", "Attia")}
OUT = {"b1": {8, 10, 12, 13, 22}, "b2": {7, 8, 9, 15, 16}, "b3": {2, 23, 32, 37, 42, 43}, "b9": set()}
SWAP = {"2017-05-12": ("15", "16"), "2017-06-30": ("15", "16")}


def read_cells():
    cells = []
    for lab, (date, paper) in FILES.items():
        desc = list(csv.DictReader(open(os.path.join(EXT, f"{lab}_a0_cells.csv"))))
        series = {}
        for r in csv.DictReader(open(os.path.join(EXT, f"{lab}_a0_cycles.csv"))):
            series.setdefault(int(r["index"]), {})[int(r["cycle"])] = r
        chan = {int(d["index"]): d["channel_id"].strip() for d in desc}
        if date in SWAP:
            a = [i for i, c in chan.items() if c == SWAP[date][0]]
            b = [i for i, c in chan.items() if c == SWAP[date][1]]
            if a and b:
                sa, sb = series[a[0]], series[b[0]]
                for n in set(sa) & set(sb):
                    for k in ("Tavg", "Tmax", "Tmin"):
                        sa[n][k], sb[n][k] = sb[n][k], sa[n][k]
        for d in desc:
            i = int(d["index"])
            if i in OUT[lab]:
                continue
            s = series[i]
            T = [float(s[n]["Tavg"]) for n in range(2, 101)]
            ir = [float(s[n]["IR"]) for n in range(2, 101)]
            ir = [v for v in ir if math.isfinite(v) and v > 0]
            cells.append({"cell": f"{date}:{i}", "batch": date, "paper": paper, "policy": d["policy_readable"],
                          "T": sum(T) / len(T), "lnR": math.log(sum(ir) / len(ir)) if ir else None})
    return cells


def lifetimes():
    out = {}
    for r in csv.DictReader(open(os.path.join(ROOT, "inputs", "matr", "stage6_lifetimes.csv"))):
        if r["cycle_life"] not in ("", "nan"):
            out[r["cell"]] = float(r["cycle_life"])
    return out


def fit_within(groups, cols):
    """groups: list of lists of cell dicts (a set drawn twice appears twice). cols: names of the regressors.
    Returns the coefficients of y on cols with one dummy per group, the residual sd and its df."""
    rows, y, dummies = [], [], []
    for g, members in enumerate(groups):
        for c in members:
            rows.append([c[k] for k in cols])
            y.append(c["y"])
            dummies.append(g)
    X = np.zeros((len(rows), len(cols) + len(groups)))
    X[:, :len(cols)] = np.array(rows)
    X[np.arange(len(rows)), len(cols) + np.array(dummies)] = 1.0
    coef, *_ = np.linalg.lstsq(X, np.array(y), rcond=None)
    res = np.array(y) - X @ coef
    df = len(rows) - len(groups) - len(cols)
    return coef[:len(cols)], math.sqrt(float(res @ res) / df), df


def fit_between(groups):
    """Set means of y on set means of T, weighted by cells, one dummy per batch."""
    batches = sorted(set(m[0]["batch"] for m in groups))
    X, y, w = [], [], []
    for m in groups:
        row = [sum(c["T"] for c in m) / len(m)] + [1.0 if m[0]["batch"] == b else 0.0 for b in batches]
        X.append(row)
        y.append(sum(c["y"] for c in m) / len(m))
        w.append(len(m))
    X, y, sw = np.array(X), np.array(y), np.sqrt(np.array(w, dtype=float))
    coef, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
    xt = X[:, 0] - (X[:, 1:] @ np.linalg.lstsq(X[:, 1:] * sw[:, None], X[:, 0] * sw, rcond=None)[0])
    yt = y - (X[:, 1:] @ np.linalg.lstsq(X[:, 1:] * sw[:, None], y * sw, rcond=None)[0])
    return coef[0], float((np.array(w) * xt * xt).sum()), float((np.array(w) * yt * yt).sum()), float((np.array(w) * xt * yt).sum())


def pairs_slope(groups):
    num = den = 0.0
    for m in groups:
        for i in range(len(m)):
            for j in range(i + 1, len(m)):
                a, b = (m[i], m[j]) if m[i]["T"] > m[j]["T"] else (m[j], m[i])
                dT, d = a["T"] - b["T"], a["y"] - b["y"]
                num += dT * d
                den += dT * dT
    return -num / den


def sign_test(groups):
    hits, ps, ties = 0, [], 0
    for m in groups:
        if len(m) < 2:
            continue
        Ts, ys = [c["T"] for c in m], [c["y"] for c in m]
        if Ts.count(max(Ts)) > 1 or ys.count(min(ys)) > 1:
            ties += 1
            continue
        hits += int(Ts.index(max(Ts)) == ys.index(min(ys)))
        ps.append(1.0 / len(m))
    dist = [1.0]
    for p in ps:
        new = [0.0] * (len(dist) + 1)
        for k, v in enumerate(dist):
            new[k] += v * (1 - p)
            new[k + 1] += v * p
        dist = new
    pval = sum(v for v in dist if v <= dist[hits] * (1 + 1e-7))
    return {"sets": len(ps), "warmest_first": hits, "p_two_sided": min(1.0, pval), "left_out_ties": ties}


def resistance_stats(groups):
    rm = [[c for c in m if c["lnR"] is not None] for m in groups if len(m) >= 2]
    rm = [m for m in rm if len(m) >= 2]
    if not rm:
        return {}
    (bT, bR), _, _ = fit_within(rm, ["T", "r"])
    (bT1,), _, _ = fit_within(rm, ["T"])
    return {"s_w_R_held": -bT, "lnL_per_lnR": bR, "s_w_same_cells": -bT1}


def statistics(groups, with_R):
    multi = [m for m in groups if len(m) >= 2]
    (b,), sd, df = fit_within(multi, ["T"])
    out = {"s_w": -b, "sd_cell": sd, "s_pairs": pairs_slope(multi)}
    bb, CTT, Cyy, CTy = fit_between(groups)
    out["s_b"] = -bb
    out["s_b_minus_s_w"] = out["s_b"] - out["s_w"]
    out["share"] = out["s_w"] ** 2 * CTT / Cyy
    out["share_residual_form"] = 1.0 - (Cyy + 2 * out["s_w"] * CTy + out["s_w"] ** 2 * CTT) / Cyy
    if with_R:
        out.update(resistance_stats(groups))
    return out


def draw(groups, rnd):
    strata = {}
    for m in groups:
        strata.setdefault((m[0]["batch"], len(m) >= 2), []).append(m)
    out = []
    for key in sorted(strata):
        g = strata[key]
        out.extend(g[rnd.randrange(len(g))] for _ in range(len(g)))
    return out


def main(B=10000):
    cells = read_cells()
    life = lifetimes()
    for c in cells:
        c["y"] = math.log(life[c["cell"]])
        c["r"] = c["lnR"]
    W = json.load(open(os.path.join(ROOT, "results", "stage6_within.json")))["registered"]
    Bt = json.load(open(os.path.join(ROOT, "results", "stage6_between.json")))["registered"]
    rnd = random.Random(20261117)
    report, ok = {"check": "V2 (D96 (1))", "rows": []}, True
    for scope, papers in (("pooled", ("Severson", "Attia")), ("Severson", ("Severson",)), ("Attia", ("Attia",))):
        sets = {}
        for c in cells:
            if c["paper"] in papers:
                sets.setdefault((c["batch"], c["policy"]), []).append(c)
        groups = [sets[k] for k in sorted(sets)]
        pt = statistics(groups, True)
        draws = [statistics(draw(groups, rnd), False) for _ in range(B)]
        rgroups = [[c for c in m if c["lnR"] is not None] for m in groups if len(m) >= 2]
        rgroups = [m for m in rgroups if len(m) >= 2]
        if rgroups:
            rdraws = [resistance_stats(draw(rgroups, rnd)) for _ in range(B)]
            for d, r in zip(draws, rdraws):
                d.update(r)
        pipe = {"s_w": W[scope]["B1"]["s_w_per_K"], "s_pairs": W[scope]["B1"]["pairs_form_s_per_K"],
                "sd_cell": W[scope]["B4"]["sd_cell"], "s_b": Bt[scope]["C1"]["s_b_per_K"],
                "s_b_minus_s_w": Bt[scope]["C2"]["s_b_minus_s_w_per_K"], "share": Bt[scope]["C3"]["share"],
                "share_residual_form": Bt[scope]["C3"]["share_residual_form"]}
        if "s_w_R_held_per_K" in W[scope].get("B2", {}):
            pipe.update({"s_w_R_held": W[scope]["B2"]["s_w_R_held_per_K"], "lnL_per_lnR": W[scope]["B2"]["lnL_per_lnR"],
                         "s_w_same_cells": W[scope]["B2"]["s_w_same_cells_per_K"]})
        for k, p in pipe.items():
            d = abs(pt[k] - p["value"])
            good = d <= 1e-10 * max(1.0, abs(p["value"]))
            row = {"scope": scope, "quantity": k, "separate": pt[k], "pipeline": p["value"], "abs_diff": d}
            if k in draws[0]:
                lo, hi = np.percentile([x[k] for x in draws], [2.5, 97.5])
                tol = [4 * math.sqrt(2) * m for m in p["mc_sd_of_ends"]]
                ends = [abs(lo - p["ci95"][0]), abs(hi - p["ci95"][1])]
                good = good and ends[0] <= tol[0] and ends[1] <= tol[1]
                row.update({"separate_ci95": [float(lo), float(hi)], "pipeline_ci95": p["ci95"],
                            "ends_in_mc_sd": [e / m if m > 0 else None for e, m in zip(ends, p["mc_sd_of_ends"])]})
            row["pass"] = bool(good)
            ok &= good
            report["rows"].append(row)
            print(f"{scope:9s} {k:20s} {pt[k]:+.12f} {p['value']:+.12f} diff {d:.1e} {row.get('separate_ci95', '')} {good}")
        st = sign_test([m for m in groups if len(m) >= 2])
        pst = W[scope]["B3"]
        same = all(st[k] == pst[k] for k in ("sets", "warmest_first", "left_out_ties")) and \
            abs(st["p_two_sided"] - pst["p_two_sided"]) <= 1e-10
        ok &= same
        report["rows"].append({"scope": scope, "quantity": "B3", "separate": st, "pipeline": pst, "pass": bool(same)})
    report["pass"] = bool(ok)
    json.dump(report, open(os.path.join(HERE, "v2.json"), "w"), indent=1)
    print("V2 pass:", ok)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 10000)
