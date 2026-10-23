#!/usr/bin/env python3
"""Stage 3 independent check, step 4: every reading choice of s1_cells.py run end to end, to show which scored
numbers and verdicts depend on a reading. Also an unregistered floor sensitivity without the NMC 20-80 % checks that
are coulomb-limited cycling discharges ending near the cut-off.

Usage: python3 -I s4_variants.py <out_json>
"""
import copy, json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib
import s1_cells
import s3_groups

VARIANTS = [
    {},
    {"check_v": "one_sided"},
    {"grouping": "G2"},
    {"units": "U2"},
    {"units": "U3"},
    {"units": "U5"},
    {"zero_dt": "nan"},
    {"rests": "B2"},
]


def summarise(cells, g):
    lev = {k: {kk: v.get(kk) for kk in ("s_per_K", "T_lo_C", "T_hi_C", "s_setpoint_per_K", "law2_s_per_K", "ci95")}
           for k, v in g["levers"].items()}
    return {"verdicts": g["verdicts"], "levers": lev,
            "floor_mean": {c: g["floor"][c]["mean"] for c in g["floor"]},
            "P6": g["P6"], "share_a_all": {k: g["thermal_share"]["a"]["all"][k] for k in ("n_entering", "k_warmer_shorter", "n_ties", "n_pairs")},
            "P4_alternatives": g["P4_alternatives"], "module_divergence": g["module"]["divergence"]}


def main(out_path):
    idx = lib.index()
    res = {}
    base = None
    for opt in VARIANTS:
        tag = "_".join(f"{k}={v}" for k, v in sorted(opt.items())) or "primary"
        cells = [s1_cells.process(cid, idx, opt) for cid in sorted(idx)]
        g = s3_groups.analyse(copy.deepcopy(cells))
        rec = summarise(cells, g)
        if base is None:
            base = {c["cell_id"]: c for c in cells}
        diffs = {}
        for k in ("L_efc", "resolution_efc", "T_cyc_C", "measure_a_K", "measure_b_C"):
            mx, arg = 0.0, None
            for c in cells:
                a, b = c[k], base[c["cell_id"]][k]
                if (isinstance(a, float) and math.isnan(a)) and (isinstance(b, float) and math.isnan(b)):
                    continue
                d = abs(a - b) if not (math.isnan(a) or math.isnan(b)) else float("inf")
                if d > mx:
                    mx, arg = d, c["cell_id"]
            diffs[k] = {"max_abs_vs_primary": mx, "cell": arg}
        rec["per_cell_max_diff_vs_primary"] = diffs
        res[tag] = rec
        print(tag, rec["verdicts"], {k: round(v["s_per_K"], 6) if v["s_per_K"] is not None and not math.isnan(v["s_per_K"]) else None
                                     for k, v in rec["levers"].items() if "1C" in k}, flush=True)

    # unregistered floor sensitivity: drop checks of 20-80 % cells whose every discharge is the 60 % coulomb limit
    cells = [s1_cells.process(cid, idx, {}) for cid in sorted(idx)]
    changed = {}
    for c in cells:
        if c["window"] != "20-80":
            continue
        lim = 0.6 * c["rated_Ah"]
        bad = [k for k, ch in enumerate(c["checks_list"]) if all(abs(x - lim) <= 0.005 * lim for x in ch["discharge_ah"])]
        if bad:
            keep = [ch for k, ch in enumerate(c["checks_list"]) if k not in bad]
            for ch in keep:
                ch.pop("dropped_half", None)
            lt = s1_cells.lifetime(keep, lib.THRESHOLD[c["chemistry"]])
            changed[c["cell_id"]] = {"L_registered_rule": c["L_efc"], "L_without_coulomb_limited_checks": lt["L_efc"],
                                     "resolution_without": lt["resolution_efc"], "checks_removed": len(bad),
                                     "removed_check_efc": [c["checks_list"][k]["efc"] for k in bad]}
            c["L_efc"], c["resolution_efc"] = lt["L_efc"], lt["resolution_efc"]
    g = s3_groups.analyse(copy.deepcopy(cells))
    res["floor_without_coulomb_limited_checks"] = {"cells_changed": changed,
                                                   "floor": g["floor"], "P5": g["verdicts"]["P5"]}
    print("floor without coulomb-limited checks:", {k: round(v["mean"], 5) for k, v in g["floor"].items()}, g["verdicts"]["P5"])
    with open(out_path, "w") as fh:
        json.dump(res, fh, indent=1, default=float)


if __name__ == "__main__":
    main(sys.argv[1])
