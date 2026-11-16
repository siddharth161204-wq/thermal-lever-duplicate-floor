#!/usr/bin/env python3
"""Stage 6, Part A0 from the cell table: the counts (A0.3), the reliability of the registered temperature (A0.5), the
resistance's validity (A0.6), the channel (A0.7) and the power (A0.8), with V5 and the scores of P1 and P6
(DECISIONS.md D85 to D92, D95). No lifetime or capacity is read.

Usage: python3 -I src/stage6/a0_analysis.py [extract directory] [output directory] [--quick]
Writes stage6_a0.json and stage6_power.json in the output directory (results/ by default), verification/stage6/v5.json
and the A0 part of stage6_scorecard.json.
"""
import json, math, os, sys, time

import numpy as np
from scipy import optimize, stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s6_common as C                   # noqa: E402
import a0_cells as A                    # noqa: E402

SCOPES = {"Severson": ("Severson",), "Attia": ("Attia",), "pooled": ("Severson", "Attia")}
N_SIM = 10000
N_PERM = 10000


def within_cells(rows, papers, flagged=True):
    out = [r for r in rows if r["enters"] and r["paper"] in papers and r["set_size"] >= 2]
    if not flagged:
        keep = [r for r in out if not r.get("contact_flag")]
        sizes = {}
        for r in keep:
            sizes[r["set"]] = sizes.get(r["set"], 0) + 1
        out = [r for r in keep if sizes[r["set"]] >= 2]
    return out


def deviations(cells, key):
    g = {}
    for r in cells:
        g.setdefault(r["set"], []).append(r[key])
    mean = {k: float(np.mean(v)) for k, v in g.items()}
    return np.array([r[key] - mean[r["set"]] for r in cells])


def counts(rows):
    out = {}
    for name, papers in SCOPES.items():
        ent = [r for r in rows if r["enters"] and r["paper"] in papers]
        sets = {}
        for r in ent:
            sets.setdefault(r["set"], 0)
            sets[r["set"]] += 1
        multi = {k: v for k, v in sets.items() if v >= 2}
        out[name] = {"cells_entering": len(ent), "sets": len(sets), "sets_of_two_or_more": len(multi),
                     "cells_in_them": int(sum(multi.values())), "pairs": int(sum(v * (v - 1) // 2 for v in multi.values())),
                     "set_sizes": {str(k): int(sum(1 for v in sets.values() if v == k)) for k in sorted(set(sets.values()))}}
    by_batch = {}
    for r in rows:
        b = by_batch.setdefault(r["batch"], {"entries": 0, "enter": 0, "why_not": {}})
        b["entries"] += 1
        if r["enters"]:
            b["enter"] += 1
        else:
            b["why_not"][r["why_not"]] = b["why_not"].get(r["why_not"], 0) + 1
    out["by_batch"] = by_batch
    return out


def reliability(rows):
    out = {}
    for name, papers in SCOPES.items():
        cells = within_cells(rows, papers)
        res = {"cells": len(cells)}
        for lab, (a, b) in {"odd_even": ("T_odd", "T_even"), "halves": ("T_2_50", "T_51_100")}.items():
            da, db = deviations(cells, a), deviations(cells, b)
            r = float((da * db).sum() / math.sqrt((da * da).sum() * (db * db).sum()))
            res[lab] = {"r": r, "spearman_brown": 2 * r / (1 + r)}
        out[name] = res
    return out


def resistance(rows):
    out = {}
    for b in sorted(set(r["batch"] for r in rows)):
        ent = [r for r in rows if r["enters"] and r["batch"] == b]
        out[b] = {"cells": len(ent), "all_99_valid": sum(r["IR_valid_cycles"] == 99 for r in ent),
                  "some_invalid": sum(0 < r["IR_valid_cycles"] < 99 for r in ent),
                  "none_valid": sum(r["IR_valid_cycles"] == 0 for r in ent),
                  "valid_cycles_min": min((r["IR_valid_cycles"] for r in ent), default=None),
                  "IR_ohm_range": [min((r["IR_ohm"] for r in ent if math.isfinite(r["IR_ohm"])), default=None),
                                   max((r["IR_ohm"] for r in ent if math.isfinite(r["IR_ohm"])), default=None)]}
    return out


def channel(rows, seed_k=70):
    """A0.7: the correlation, over every channel and every two batches in which it held a cell of a set of two or more,
    between that cell's within-set deviation of the registered temperature in the earlier batch and in the later one;
    a permutation test, the channel labels shuffled within each batch."""
    cells = within_cells(rows, ("Severson", "Attia"))
    dev = deviations(cells, "T_reg")
    recorded = {b: all(str(r["channel"]).strip() not in ("", "None") for r in rows if r["batch"] == b)
                for b in sorted(set(r["batch"] for r in rows))}
    by = {}
    for r, d in zip(cells, dev):
        by.setdefault(r["batch"], {})[str(r["channel"]).strip()] = d
    batches = sorted(by)

    def corr(maps):
        x, y = [], []
        for i, a in enumerate(batches):
            for b in batches[i + 1:]:
                for ch in maps[a]:
                    if ch in maps[b]:
                        x.append(maps[a][ch])
                        y.append(maps[b][ch])
        x, y = np.array(x), np.array(y)
        if len(x) < 3:
            return float("nan"), len(x)
        return float(np.corrcoef(x, y)[0, 1]), len(x)

    r_obs, n_pairs = corr(by)
    rng = np.random.default_rng([C.SEED, seed_k])
    perm = []
    for _ in range(N_PERM):
        pm = {}
        for b in batches:
            chs = list(by[b].keys())
            vals = np.array([by[b][c] for c in chs])
            pm[b] = dict(zip(chs, vals[rng.permutation(len(vals))]))
        perm.append(corr(pm)[0])
    perm = np.array(perm)
    p = float((1 + np.sum(np.abs(perm) >= abs(r_obs))) / (1 + N_PERM)) if math.isfinite(r_obs) else float("nan")
    # the instrument for B5: each cell's channel's mean deviation in the other batches
    inst = {}
    for r in cells:
        ch = str(r["channel"]).strip()
        other = [by[b][ch] for b in batches if b != r["batch"] and ch in by[b]]
        inst[r["cell"]] = float(np.mean(other)) if other else float("nan")
    carries = bool(math.isfinite(r_obs) and r_obs > 0 and p < 0.05)
    return {"recorded_by_batch": recorded, "pairs": n_pairs, "r": r_obs, "p_permutation_two_sided": p,
            "permutations": N_PERM, "seed": [C.SEED, seed_k], "channel_carries_temperature": carries,
            "cells_with_instrument": int(sum(math.isfinite(v) for v in inst.values())), "instrument": inst}


def sim_detectable(Tdev, groups, shape, sd_cell, seed, n_sim=N_SIM, normal=False):
    """V5: the detectable |s| by simulation, cell noise e / sqrt(2) with e from the shape (or normal of sd_cell)."""
    rng = np.random.default_rng(seed)
    N = len(Tdev)
    J = len(np.unique(groups))
    df = N - J - 1
    tc = stats.t.ppf(1 - C.ALPHA / 2, df)
    Sxx = float((Tdev * Tdev).sum())
    u, se = np.empty(n_sim), np.empty(n_sim)
    gi = np.unique(groups, return_inverse=True)[1]
    for k in range(n_sim):
        e = (sd_cell * rng.standard_normal(N)) if normal else C.shape_draws(shape, N, rng) / math.sqrt(2.0)
        ed = e - (np.bincount(gi, e) / np.bincount(gi))[gi]
        uk = float((Tdev * ed).sum() / Sxx)
        res = ed - uk * Tdev
        u[k] = uk
        se[k] = math.sqrt(float(res @ res) / df / Sxx)
    power = lambda s: float(np.mean(np.abs(u - s) > tc * se))
    lo, hi = 0.0, 1e-3
    while power(hi) < C.POWER:
        hi *= 2.0
    while hi - lo > 1e-6 * hi:
        mid = 0.5 * (lo + hi)
        if power(mid) >= C.POWER:
            hi = mid
        else:
            lo = mid
    return hi


def power(rows, archive, quick=False):
    out, v5 = [], []
    k = 0
    for name, papers in SCOPES.items():
        cells = within_cells(rows, papers)
        Tdev = deviations(cells, "T_reg")
        groups = np.array([r["set"] for r in cells])
        N, J = len(cells), len(np.unique(groups))
        Sxx = float((Tdev * Tdev).sum())
        sizes = {}
        for r in cells:
            sizes[r["set"]] = sizes.get(r["set"], 0) + 1
        n_pairs = int(sum(v * (v - 1) // 2 for v in sizes.values()))
        sxx_set = {}
        for r, d in zip(cells, Tdev):
            sxx_set[r["set"]] = sxx_set.get(r["set"], 0.0) + d * d
        pair_ss = sum(sizes[s] * v for s, v in sxx_set.items())          # sum over pairs of dT^2 = n_j Sxx_j
        F = sum(sizes[s] ** 2 * v for s, v in sxx_set.items()) / (2.0 * pair_ss)
        dT_rms = math.sqrt(pair_ss / n_pairs)
        for floor in ("KIT", "Geslin"):
            sh = archive["shapes"][floor]
            sd_c = sh["sd_cell"]
            calc = C.fe_detectable(Sxx, N, J, sd_c)
            k += 1
            sim = sim_detectable(Tdev, groups, sh, sd_c, [C.SEED, 80 + k], 2000 if quick else N_SIM)
            simn = sim_detectable(Tdev, groups, sh, sd_c, [C.SEED, 90 + k], 2000 if quick else N_SIM, normal=True)
            row = {"scope": name, "floor": floor, "registered": True, "cells": N, "sets": J, "df": N - J - 1,
                   "Sxx_K2": Sxx, "T_dev_rms_K": math.sqrt(Sxx / N), "sd_cell": sd_c, "sd_e_pair": sh["sd_e"],
                   "detectable_s_per_K": calc, "sim_registered_shape_detectable_s_per_K": sim,
                   "sim_normal_detectable_s_per_K": simn, "v5_rel_diff": sim / calc - 1.0,
                   "v5_pass": abs(sim / calc - 1.0) <= 0.05, "seeds": [[C.SEED, 80 + k], [C.SEED, 90 + k]],
                   "pairs": n_pairs, "pair_dT_rms_K": dT_rms, "shared_cell_F": F,
                   "pairs_form_detectable_s_per_K": C.pairs_detectable(n_pairs, dT_rms, sh["sd_e"]),
                   "pairs_form_shared_cells_detectable_s_per_K": C.pairs_detectable(n_pairs, dT_rms, sh["sd_e"] * math.sqrt(F)),
                   "normal_shape_detectable_s_per_K": C.fe_detectable(Sxx, N, J, sh["sd_e_normal_shape"] / math.sqrt(2.0))}
            cu = within_cells(rows, papers, flagged=False)           # unscored: without the flagged cells (D92 (5))
            if cu:
                du = deviations(cu, "T_reg")
                Ju = len(set(r["set"] for r in cu))
                row["without_flagged"] = {"cells": len(cu), "sets": Ju, "T_dev_rms_K": math.sqrt(float((du * du).mean())),
                                          "detectable_s_per_K": C.fe_detectable(float((du * du).sum()), len(cu), Ju, sd_c)}
            out.append(row)
            v5.append({k2: row[k2] for k2 in ("scope", "floor", "detectable_s_per_K", "sim_registered_shape_detectable_s_per_K",
                                              "sim_normal_detectable_s_per_K", "v5_rel_diff", "v5_pass")})
    return out, v5


def main(ext, out_dir, quick=False, batches=None):
    t0 = time.time()
    arch = C.archive()
    rows, sets = A.build(ext, batches or A.BATCHES)
    res = {"counts": counts(rows), "reliability": reliability(rows), "resistance": resistance(rows),
           "channel": channel(rows)}
    pw, v5 = power(rows, arch, quick)
    pooled_kit = next(p for p in pw if p["scope"] == "pooled" and p["floor"] == "KIT")
    lever = C.LEVERS_D82[("registered 1C", "15 to 25 C")][0]
    p1 = pooled_kit["detectable_s_per_K"] < lever
    p1_sim = pooled_kit["sim_registered_shape_detectable_s_per_K"] < lever
    r6 = res["reliability"]["pooled"]["odd_even"]["r"]
    score = {"P1": {"quantity": "pooled detectable |s| under the KIT floor, % per K", "value": 100 * pooled_kit["detectable_s_per_K"],
                    "threshold": 100 * lever, "holds": bool(p1) if (pooled_kit["v5_pass"] or p1 == p1_sim) else None,
                    "note": "" if (pooled_kit["v5_pass"] or p1 == p1_sim) else "not scored: V5 fails and the simulation gives the other verdict"},
             "P6": {"quantity": "pooled within-set correlation of the registered temperature on odd and on even cycles",
                    "value": r6, "threshold": 0.8, "holds": bool(r6 > 0.8)}}
    os.makedirs(out_dir, exist_ok=True)
    json.dump(res, open(os.path.join(out_dir, "stage6_a0.json"), "w"), indent=1)
    json.dump({"method": "D92", "cases": pw}, open(os.path.join(out_dir, "stage6_power.json"), "w"), indent=1)
    vdir = C.path("verification/stage6") if batches is None else out_dir
    os.makedirs(vdir, exist_ok=True)
    json.dump({"check": "V5 (D95)", "cases": v5, "pass": all(x["v5_pass"] for x in v5)},
              open(os.path.join(vdir, "v5.json"), "w"), indent=1)
    json.dump({"scored_at": "end of Part A0", "P1": score["P1"], "P6": score["P6"]},
              open(os.path.join(out_dir, "stage6_scorecard.json"), "w"), indent=1)
    print(json.dumps(res["counts"], indent=1))
    print(json.dumps(res["reliability"], indent=1))
    print(json.dumps({k: v for k, v in res["channel"].items() if k != "instrument"}, indent=1))
    for p in pw:
        print(p["scope"], p["floor"], f"N {p['cells']} J {p['sets']} Trms {p['T_dev_rms_K']:.3f} detectable {100*p['detectable_s_per_K']:.3f} "
              f"sim {100*p['sim_registered_shape_detectable_s_per_K']:.3f} ({100*p['v5_rel_diff']:+.1f} %) normal sim "
              f"{100*p['sim_normal_detectable_s_per_K']:.3f}; pairs form {100*p['pairs_form_detectable_s_per_K']:.3f}")
    print(json.dumps(score, indent=1), f"{time.time() - t0:.0f} s")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    main(args[0] if args else C.path("inputs/matr/a0"), args[1] if len(args) > 1 else C.path("results"),
         "--quick" in sys.argv)
