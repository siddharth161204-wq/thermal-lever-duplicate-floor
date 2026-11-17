#!/usr/bin/env python3
"""Stage 6, post hoc (DECISIONS.md D97 (4) and D98), after every registered result was written. Nothing here is scored.

(1) The channel test of D91 (`a0_analysis.channel`) with each cell's within-set deviation of the Tmin measure, and of
its median swing, in place of T: does a channel's floor or swing repeat from batch to batch? 5,000 permutations of the
channel labels within each batch, each from a generator seeded [20261117, 1].
(2) Leave one set out: B1's s_w with each set of two or more left out in turn, for T and for the Tmin bracket, and C1's
s_b with each set left out in turn, pooled; the range and the set at each end.
(3) B1's s_w within each batch, with the percentile interval of 10,000 draws of that batch's sets (strata as in
D93 (6)), seeds [20261117, 10 + i] for the batches in date order.

Usage: python3 -I src/stage6/s6_posthoc.py [extract directory] [lifetimes CSV] [results directory]
Writes stage6_posthoc.json in the results directory.
"""
import json, math, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s6_common as C                   # noqa: E402
import s6_lever as Lv                   # noqa: E402
import s6_parts as P                    # noqa: E402
import a0_cells as A                    # noqa: E402
import a0_analysis as AN                # noqa: E402

SEED_PH = 20261117
N_PERM_PH = 5000
BOTH = ("Severson", "Attia")


def channel_key(rows, key, seed):
    """`a0_analysis.channel` with `key` in place of T_reg and N_PERM_PH permutations."""
    cells = AN.within_cells(rows, BOTH)
    dev = AN.deviations(cells, key)
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
        return float(np.corrcoef(x, y)[0, 1]), len(x)

    r_obs, n_pairs = corr(by)
    rng = np.random.default_rng(seed)
    perm = []
    for _ in range(N_PERM_PH):
        pm = {}
        for b in batches:
            chs = list(by[b].keys())
            vals = np.array([by[b][c] for c in chs])
            pm[b] = dict(zip(chs, vals[rng.permutation(len(vals))]))
        perm.append(corr(pm)[0])
    perm = np.array(perm)
    p = float((1 + np.sum(np.abs(perm) >= abs(r_obs))) / (1 + N_PERM_PH))
    return {"measure": key, "pairs": n_pairs, "r": r_obs, "p_permutation_two_sided": p,
            "permutations": N_PERM_PH, "seed": list(seed)}


def point(cells):
    sets = Lv.set_table(cells, "T", "y")
    o = Lv.ones(sets)
    w = Lv.within(sets, o)
    b = Lv.between(sets, o, w["s_w"])
    return float(w["s_w"][0]), float(b["s_b"][0]), sets


def leave_one_out(cells, what):
    full_w, full_b, sets = point(cells)
    ids = [s["set"] for s in sets if (s["n"] >= 2 if what == "s_w" else True)]
    vals = []
    for sid in ids:
        sw, sb, _ = point([c for c in cells if c["set"] != sid])
        vals.append((sw if what == "s_w" else sb, sid))
    vals.sort()
    return {"full": full_w if what == "s_w" else full_b, "sets_left_out_in_turn": len(ids),
            "min": vals[0][0], "set_at_min": vals[0][1], "max": vals[-1][0], "set_at_max": vals[-1][1]}


def main(ext, lif, out_dir):
    rows, _ = A.build(ext)
    life = P.read_life(lif)
    out = {"note": "post hoc, after every registered result; nothing here is scored (DECISIONS.md D97 (4), D98)",
           "channel": {"T_min": channel_key(rows, "T_min", (SEED_PH, 1)),
                       "swing_med": channel_key(rows, "swing_med", (SEED_PH, 1))}}
    loo = {}
    for tkey in ("T_reg", "T_min"):
        cells = P.cells_for(rows, life, BOTH, tkey)
        loo[tkey] = {"s_w": leave_one_out(cells, "s_w")}
        if tkey == "T_reg":
            loo[tkey]["s_b"] = leave_one_out(cells, "s_b")
    out["leave_one_set_out_pooled"] = loo
    per_batch = {}
    cells = P.cells_for(rows, life, BOTH, "T_reg")
    for i, b in enumerate(sorted(set(c["batch"] for c in cells))):
        cb = [c for c in cells if c["batch"] == b]
        sets = Lv.set_table(cb, "T", "y")
        multi = [s for s in sets if s["n"] >= 2]
        if multi:
            per_batch[b] = {"cells": int(sum(s["n"] for s in multi)), "sets": len(multi),
                            "s_w": float(Lv.within(sets, Lv.ones(sets))["s_w"][0])}
            if len(multi) >= 2:
                W = Lv.draw_weights(sets, C.N_BOOT, np.random.default_rng([SEED_PH, 10 + i]))
                ci, mc, used = Lv.percentile_ci(Lv.within(sets, W)["s_w"])
                per_batch[b].update({"ci95": ci, "draws_used": used, "seed": [SEED_PH, 10 + i]})
            else:
                per_batch[b]["ci95"] = None
                per_batch[b]["note"] = "one set of two or more: every draw takes it, so there is no interval"
    out["s_w_by_batch"] = per_batch
    json.dump(out, open(os.path.join(out_dir, "stage6_posthoc.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else C.path("inputs/matr/a0"),
         sys.argv[2] if len(sys.argv) > 2 else C.path("inputs/matr/stage6_lifetimes.csv"),
         sys.argv[3] if len(sys.argv) > 3 else C.path("results"))
