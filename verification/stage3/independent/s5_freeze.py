#!/usr/bin/env python3
"""Stage 3 independent check, step 5: run steps 1 to 4 and write my outputs to frozen/, then record the UTC time and
the SHA-256 of every output, script and input in frozen/freeze_manifest.json. Nothing of the pipeline's is read here.

Usage: python3 -I s5_freeze.py
"""
import csv, copy, datetime, hashlib, json, math, os, platform, sys
import numpy as np, pandas, scipy

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib, s1_cells, s3_groups, s4_variants
import s2_v2_purepython as s2

FROZEN = os.path.join(HERE, "frozen")
REPO = "/home/claude/thermal-lever-duplicate-floor"
V2_CELLS = ["SNL_18650_LFP_15C_0-100_0.5/1C_a", "SNL_18650_LFP_15C_0-100_0.5/2C_a", "SNL_18650_LFP_25C_0-100_0.5/3C_a"]
EXTRA_CELLS = ["SNL_18650_NMC_35C_0-100_0.5/1C_a", "SNL_18650_LFP_35C_0-100_0.5/1C_c",
               "SNL_18650_NMC_25C_20-80_0.5/3C_a", "SNL_18650_NCA_25C_0-100_0.5/0.5C_b"]

CELL_COLS = ["cell_id", "chemistry", "ambient_C", "window", "c_charge", "c_discharge", "replicate", "rated_Ah",
             "rows_logged", "late_records", "rows_other_record_dropped", "rows_ordered", "steps_total", "steps_charge",
             "steps_discharge", "steps_rest", "discharge_ah_total", "charge_ah_total", "efc_total",
             "registered_counter_ah_total", "registered_counter_share", "discharges_registered_rule_differs",
             "qualifying_discharges", "checks", "checks_dropped_half", "first_check_Ah", "threshold", "L_efc",
             "resolution_efc", "censored", "censor_efc", "L_sens_ref_first_check", "L_sens_ref_first_cycling",
             "units_total", "units_in_life", "units_used", "units_used_b", "units_in_life_skipped_missing_mean",
             "measure_a_K", "measure_b_C", "T_cyc_C", "cell_minus_env_K", "first_unit_cycle"]


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for ch in iter(lambda: fh.read(1 << 20), b""):
            h.update(ch)
    return h.hexdigest()


def first_cells_by_rate(cells):
    """V2 selection (D49): the first 1C, 2C and 3C cells in cell-ID order with an uncensored lifetime."""
    sel = {}
    for c in sorted(cells, key=lambda r: r["cell_id"]):
        if c["c_discharge"] in (1.0, 2.0, 3.0) and not c["censored"] and c["c_discharge"] not in sel:
            sel[c["c_discharge"]] = c["cell_id"]
    return sel


def logged_counter_share(cid, idx):
    """The registered counter rule on the rows as logged (as the inventory was computed), for comparison."""
    a = lib.load(cid, idx)
    m = lib.parse_id(cid)
    lab = lib.labels(a["I"], m["rated_Ah"])
    st, en = lib.steps(lab)
    dis = lab[st] == lib.DISCHARGE
    reg = lib.counter_ah_registered(a["Qd"], st, en)[dis]
    ps = lib.step_ah_per_sample(a["Qd"], st, en)[dis]
    return float(reg.sum()), float(ps.sum())


def main():
    os.makedirs(FROZEN, exist_ok=True)
    idx = lib.index()
    cells = [s1_cells.process(cid, idx, {}) for cid in sorted(idx)]
    for c in cells:
        c["L_sens_ref_first_check"] = c["sensitivity"].get("L_sens_ref_first_check")
        c["L_sens_ref_first_cycling"] = c["sensitivity"].get("L_sens_ref_first_cycling")
    with open(os.path.join(FROZEN, "cells_full.json"), "w") as fh:
        json.dump(cells, fh, indent=0, default=float)
    with open(os.path.join(FROZEN, "cells.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(CELL_COLS)
        for c in cells:
            w.writerow([repr(c.get(k)) if isinstance(c.get(k), float) else c.get(k) for k in CELL_COLS])
    with open(os.path.join(FROZEN, "checks.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cell_id", "check", "first_discharge_step", "n_discharges", "capacity_Ah", "fraction_of_first", "efc",
                    "cycle_index", "t_start_s", "dropped_half"])
        for c in cells:
            for k, ch in enumerate(c["checks_list"]):
                w.writerow([c["cell_id"], k, ch["first_discharge_step"], ch["n_discharges"], repr(ch["capacity_Ah"]),
                            repr(ch["capacity_Ah"] / c["first_check_Ah"]), repr(ch["efc"]), ch["cycle_index"], ch["t_start_s"],
                            ch["dropped_half"]])

    g = s3_groups.analyse(copy.deepcopy(cells))
    with open(os.path.join(FROZEN, "groups.json"), "w") as fh:
        json.dump(g, fh, indent=1, default=float)
    with open(os.path.join(FROZEN, "pairs.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        cols = ["group", "chemistry", "c_discharge", "a", "b", "censored_any", "rel_diff", "abs_diff_efc", "resolution_efc", "below_resolution"]
        w.writerow(cols)
        for p in g["pairs_all"]:
            w.writerow([repr(p.get(k)) if isinstance(p.get(k), float) else p.get(k) for k in cols])
    for meas in ("a", "b"):
        with open(os.path.join(FROZEN, f"oriented_{meas}.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            cols = ["group", "chemistry", "c_discharge", "a", "b", "status", "warmer", "cooler", "dT", "predicted", "observed",
                    "resolution", "tie", "enters_sign_test", "warmer_shorter"]
            w.writerow(cols)
            for q in g["oriented_pairs"][meas]:
                w.writerow([repr(q.get(k)) if isinstance(q.get(k), float) else q.get(k) for k in cols])

    # V2 (registered) and the extra cells, by the pure-Python implementation, and the agreement with s1
    sel = first_cells_by_rate(cells)
    v2_ids = [sel[1.0], sel[2.0], sel[3.0]]
    pure = [s2.analyse(cid) for cid in v2_ids + EXTRA_CELLS]
    byid = {c["cell_id"]: c for c in cells}
    agree = {}
    for r in pure:
        m = byid[r["cell_id"]]
        rec = {}
        for k in ("L_efc", "resolution_efc", "measure_a_K", "measure_b_C", "T_cyc_C", "discharge_ah_total",
                  "registered_counter_ah_total", "first_check_Ah"):
            a, b = r[k], m[k]
            if a is None:
                rec[k] = "both missing" if (isinstance(b, float) and math.isnan(b)) else "MISMATCH"
            else:
                rec[k] = {"pure_python": a, "vectorised": b, "abs_diff": abs(a - b), "rel_diff": abs(a - b) / abs(b) if b else abs(a - b)}
        for k in ("units_used", "checks", "checks_dropped_half", "late_records", "rows_other_record_dropped", "steps_charge",
                  "steps_discharge", "steps_rest"):
            rec[k] = {"pure_python": r[k], "vectorised": m[k], "equal": r[k] == m[k]}
        agree[r["cell_id"]] = rec
    with open(os.path.join(FROZEN, "v2_cells.json"), "w") as fh:
        json.dump({"v2_selection_rule": "first 1C, 2C and 3C cells in cell-ID order with an uncensored lifetime (D49)",
                   "v2_cells": v2_ids, "extra_cells": EXTRA_CELLS, "pure_python": pure,
                   "agreement_with_vectorised": agree}, fh, indent=1, default=float)

    # D50 counts
    d50 = {"late_records_total": sum(c["late_records"] for c in cells),
           "files_with_late_records": sum(1 for c in cells if c["late_records"]),
           "late_records_per_file": {c["cell_id"]: c["late_records"] for c in cells if c["late_records"]},
           "late_record_max_lag_s": {c["cell_id"]: c["late_detail"]["max_lag_behind_first_record_of_next_cycle_s"] for c in cells if c.get("late_detail")},
           "other_record_rows": {c["cell_id"]: c["other_detail"] for c in cells if c.get("other_detail")},
           "discharge_steps_total": sum(c["steps_discharge"] for c in cells),
           "discharge_counter_decreases_ordered_rows": sum(c["discharge_counter_decreases"] for c in cells),
           "decreases_without_cycle_index_change": {c["cell_id"]: {"n": c["discharge_counter_decreases_same_cycle"],
                                                                   "step_kinds (0 rest, 1 charge, -1 discharge)": c["discharge_counter_decreases_same_cycle_step_kinds"]}
                                                    for c in cells if c["discharge_counter_decreases_same_cycle"]},
           "discharges_where_registered_rule_differs": sum(c["discharges_registered_rule_differs"] for c in cells),
           "cells_where_registered_rule_differs": [c["cell_id"] for c in cells if c["discharges_registered_rule_differs"]],
           "registered_rule_share_of_discharge_ah": {c["cell_id"]: c["registered_counter_share"] for c in cells if c["discharges_registered_rule_differs"]},
           "per_sample_over_trapezoid_discharges_5_samples_or_more": [
               min(c["per_sample_over_trapezoid_min_max"][0] for c in cells if c["per_sample_over_trapezoid_min_max"][0] is not None),
               max(c["per_sample_over_trapezoid_min_max"][1] for c in cells if c["per_sample_over_trapezoid_min_max"][1] is not None)],
           "registered_rule_on_rows_as_logged": {}}
    for cid in ("SNL_18650_LFP_35C_0-100_0.5/1C_a", "SNL_18650_NMC_35C_0-100_0.5/1C_a", "SNL_18650_NMC_35C_0-100_0.5/2C_a",
                "SNL_18650_NCA_35C_0-100_0.5/1C_c", "SNL_18650_NCA_35C_0-100_0.5/1C_d"):
        reg, ps = logged_counter_share(cid, idx)
        d50["registered_rule_on_rows_as_logged"][cid] = {"registered_ah": reg, "per_sample_ah_logged_rows": ps,
                                                         "per_sample_ah_ordered_rows": byid[cid]["discharge_ah_total"],
                                                         "registered_share_of_ordered_per_sample": reg / byid[cid]["discharge_ah_total"]}
    with open(os.path.join(FROZEN, "d50_counts.json"), "w") as fh:
        json.dump(d50, fh, indent=1, default=float)

    # every reading variant end to end
    s4_variants.main(os.path.join(FROZEN, "variants.json"))

    # freeze
    t = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    outs = sorted(f for f in os.listdir(FROZEN) if f != "freeze_manifest.json")
    man = {"frozen_utc": t, "outputs_sha256": {f: sha(os.path.join(FROZEN, f)) for f in outs},
           "scripts_sha256": {f: sha(os.path.join(HERE, f)) for f in sorted(os.listdir(HERE)) if f.endswith(".py")},
           "inputs_sha256": {"inputs/sandia/inventory.csv": sha(os.path.join(REPO, "inputs/sandia/inventory.csv")),
                             "inputs/floor.json": sha(os.path.join(REPO, "inputs/floor.json")),
                             "results/module_spread.csv": sha(os.path.join(REPO, "results/module_spread.csv")),
                             "results/stage2_lever.csv": sha(os.path.join(REPO, "results/stage2_lever.csv")),
                             "results/stage2_lever.json": sha(os.path.join(REPO, "results/stage2_lever.json")),
                             "src/thermal_lever.py": sha(os.path.join(REPO, "src/thermal_lever.py"))},
           "raw_files_sha256": {cid: v["sha256"] for cid, v in idx.items()},
           "software": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pandas.__version__,
                        "scipy": scipy.__version__},
           "pipeline_files_opened_before_freeze": []}
    with open(os.path.join(FROZEN, "freeze_manifest.json"), "w") as fh:
        json.dump(man, fh, indent=1)
    print("frozen at", t)
    for f, h in man["outputs_sha256"].items():
        print(f, h)


if __name__ == "__main__":
    main()
