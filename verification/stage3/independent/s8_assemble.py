#!/usr/bin/env python3
"""Stage 3 independent check, step 8 (after the freeze): assemble verification/stage3/v2_independent.json from the
frozen outputs, the comparisons of s6 and s7, the D50 counts against results/stage3_prepare.json, the V2 cells and
the three judged points.

Usage: python3 -I s8_assemble.py
"""
import copy, csv, json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib, s1_cells, s3_groups, s6_compare_cells, s7_compare_groups

REPO = "/home/claude/thermal-lever-duplicate-floor"
FROZEN = os.path.join(HERE, "frozen")
OUT = os.path.join(REPO, "verification", "stage3", "v2_independent.json")
SCRATCH = os.environ.get("S3V2_SCRATCH", "/tmp/claude-0/-home-claude/6d63e4c8-32e9-55cd-b414-99eb5512ec94/scratchpad/work")

PIPELINE_FILES_BEFORE_OPENING = {
    # SHA-256 and modification time, recorded at 20:58:59 UTC, after the freeze and before any was opened
    "results/stage3_cells.csv": ["f646126e3155558b2008365098c5f3a155d3372e2d59d2ab7ce1b2e222859506", "2026-10-23 20:13:12.125859807 +0000"],
    "results/stage3_checks.csv": ["92f8cf8f6506dcc52e60a31a7886053d509f8b8589e940c3744faf715d7e3187", "2026-10-23 20:13:12.177859809 +0000"],
    "results/stage3_lever.csv": ["2cadaad5babc3f14b19fa6d7642fa3c8fea474512b3854425181262979ff1eb7", "2026-10-23 20:13:12.184625740 +0000"],
    "results/stage3_lever.json": ["54581949864772721c86419788993cb888a6e87505b58aa07dde02a72c69685d", "2026-10-23 20:13:12.209859810 +0000"],
    "results/stage3_pairs.csv": ["672476bc058fe323a069d6b9341162d0d8c523a323a1218764aa8915f29bc4c8", "2026-10-23 20:13:12.193859809 +0000"],
    "results/stage3_pairs_oriented_measure_a.csv": ["45720d4201b02bcf12382e687fb75f4e4f47ff44f2c6c1175363f7f83c121bac", "2026-10-23 20:13:12.197085078 +0000"],
    "results/stage3_pairs_oriented_measure_b.csv": ["24b31be1d1167f14b843e20ffd0ac1d8c144985a8a32d80c983fe6111f69f8fd", "2026-10-23 20:13:12.198410099 +0000"],
    "results/stage3_floor.csv": ["59b35df48ed3606d8ce333766853af8d6775e3a9b6aa07f19134f2543b5bb535", "2026-10-23 20:13:12.188232240 +0000"],
    "results/stage3_thermal_share.json": ["d09319c23d600b228ed40212636f37c3b4cc6e665fa01a407a63c16d3b5adf5b", "2026-10-23 20:13:12.211811949 +0000"],
    "results/stage3_scorecard.json": ["769f2454b8577e67744c0fcbe4abe5f08d3c07bf9037647d54803792fa79d125", "2026-10-23 20:13:12.212233137 +0000"],
    "results/stage3_prepare.json": ["ee11f28dd7abb8b2ae36c81282b04c4a2158be573fe4b438639647b49899f863", "2026-10-23 20:02:42.586013822 +0000"],
}


def clean(x):
    if isinstance(x, float):
        return None if math.isnan(x) else (("inf" if x > 0 else "-inf") if math.isinf(x) else x)
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, (np.floating,)):
        return clean(float(x))
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    return x


def d50_vs_prepare(cells):
    p = json.load(open(os.path.join(REPO, "results", "stage3_prepare.json")))
    by = {c["cell_id"].replace("/", "-"): c for c in cells}
    pairs = [("late_rows", "late_records"), ("foreign_rows_dropped", "rows_other_record_dropped"), ("steps", "steps_total"),
             ("discharges", "steps_discharge"), ("discharges_where_counter_ah_differs", "discharges_registered_rule_differs"),
             ("discharge_counter_decreases", "discharge_counter_decreases"),
             ("discharge_counter_decreases_without_cycle_change", "discharge_counter_decreases_same_cycle")]
    mism = []
    mx_ah = mx_reg = 0.0
    for k, v in p["cells"].items():
        c = by[k]
        for a, b in pairs:
            if v[a] != c[b]:
                mism.append([k, a, v[a], c[b]])
        mx_ah = max(mx_ah, abs(v["discharge_ah_total"] - c["discharge_ah_total"]))
        mx_reg = max(mx_reg, abs(v["discharge_ah_total_counter_ah"] - c["registered_counter_ah_total"]))
    mine_tot = {"late_rows": sum(c["late_records"] for c in cells), "foreign_rows_dropped": sum(c["rows_other_record_dropped"] for c in cells),
                "discharges": sum(c["steps_discharge"] for c in cells),
                "discharges_where_counter_ah_differs": sum(c["discharges_registered_rule_differs"] for c in cells),
                "discharge_counter_decreases": sum(c["discharge_counter_decreases"] for c in cells),
                "discharge_counter_decreases_without_cycle_change": sum(c["discharge_counter_decreases_same_cycle"] for c in cells),
                "cells_where_counter_ah_differs": sum(1 for c in cells if c["discharges_registered_rule_differs"])}
    return {"per_cell_count_mismatches": mism, "cells_compared": len(p["cells"]),
            "max_abs_diff_discharge_ah_total_Ah": mx_ah, "max_abs_diff_counter_ah_total_Ah": mx_reg,
            "totals_mine": mine_tot, "totals_pipeline": {k: p["totals"][k] for k in mine_tot}}


def frozen_share_vs_pipeline(g):
    """My own thermal-share layout (frozen) against the pipeline, sign-test fields and NCA/NMC slopes."""
    P = json.load(open(os.path.join(REPO, "results", "stage3_thermal_share.json")))
    names = {"all": "all", "rate_1C_or_below": "1C or below", "rate_2C_and_3C": "2C and 3C", "chem_LFP": "chemistry LFP",
             "chem_NCA": "chemistry NCA", "chem_NMC": "chemistry NMC", "one_pair_per_group": "one pair per group (warmest and coolest)"}
    out = {}
    for meas in ("a", "b"):
        for mk, pk in names.items():
            m, p = g["thermal_share"][meas][mk], P["measure_" + meas][pk]
            rec = {"n_entering": [m["n_entering"], p["warmer_shorter"] + p["warmer_longer"]],
                   "warmer_shorter": [m["k_warmer_shorter"], p["warmer_shorter"]], "ties": [m["n_ties"], p["ties"]],
                   "share": [clean(m["share"]), clean(p["share"])], "p_one_sided": [clean(m["p_one_sided"]), clean(p["sign_test_p_one_sided"])],
                   "p_two_sided": [clean(m["p_two_sided"]), clean(p["sign_test_p_two_sided"])]}
            if "slope" in p:
                rec["slope"] = [m.get("slope"), p["slope"]]
                rec["slope_ci95"] = [m.get("slope_ci95"), p["slope_ci95"]]
                rec["rms_ratio"] = [m.get("rms_pred_over_obs"), p["sd_ratio_pred_over_obs"]]
            out[f"{meas}/{pk}"] = rec
    return out


def main():
    man = json.load(open(os.path.join(FROZEN, "freeze_manifest.json")))
    frozen_cells = json.load(open(os.path.join(FROZEN, "cells_full.json")))
    frozen_groups = json.load(open(os.path.join(FROZEN, "groups.json")))
    v2 = json.load(open(os.path.join(FROZEN, "v2_cells.json")))
    d50 = json.load(open(os.path.join(FROZEN, "d50_counts.json")))
    variants = json.load(open(os.path.join(FROZEN, "variants.json")))

    # per-cell comparisons
    cmp_frozen = s6_compare_cells.compare(os.path.join(FROZEN, "cells_full.json"))
    idx = lib.index()
    pr_cells = [s1_cells.process(cid, idx, {"check_v": "one_sided", "units": "U3"}) for cid in sorted(idx)]
    pr_path = os.path.join(SCRATCH, "cells_pipeline_readings.json")
    json.dump(pr_cells, open(pr_path, "w"), default=float)
    cmp_pr = s6_compare_cells.compare(pr_path)

    # group comparisons (both runs)
    grp_path = os.path.join(SCRATCH, "compare_groups.json")
    s7_compare_groups.main(grp_path)
    grp = json.load(open(grp_path))

    # V2 cells against the pipeline
    pipe_cells = {r["cell_id"]: r for r in csv.DictReader(open(os.path.join(REPO, "results", "stage3_cells.csv")))}
    v2cmp = {}
    for r in v2["pure_python"]:
        p = pipe_cells[r["cell_id"]]
        rec = {}
        for mk, pk in (("L_efc", "L_efc"), ("resolution_efc", "resolution_efc"), ("measure_a_K", "swing_K"),
                       ("measure_b_C", "rest_T_C"), ("T_cyc_C", "T_cyc_C"), ("units_used", "n_units"),
                       ("first_check_Ah", "q_first_Ah"), ("checks", "checks_found"), ("checks_dropped_half", "checks_dropped")):
            a = r[mk]
            b = float(p[pk]) if p[pk] != "" else float("nan")
            if a is None:
                rec[mk] = {"mine": None, "pipeline": clean(b), "agree": math.isnan(b)}
            else:
                rec[mk] = {"mine": a, "pipeline": b, "abs_diff": abs(a - b)}
        rec["discharge_ah_total_Ah"] = {"mine": r["discharge_ah_total"],
                                        "pipeline_stage3_prepare": json.load(open(os.path.join(REPO, "results", "stage3_prepare.json")))["cells"][r["cell_id"].replace("/", "-")]["discharge_ah_total"]}
        v2cmp[r["cell_id"]] = rec

    L = frozen_groups["levers"]
    lf = L["LFP_1C_25-35"]
    point_i = {
        "question": "P4 not scored because one LFP 35 C cell has no T_cyc: faithful to D43, D44 and D48?",
        "facts": {"cell": "SNL_18650_LFP_35C_0-100_0.5/1C_c", "L_efc": frozen_cells[[c["cell_id"] for c in frozen_cells].index("SNL_18650_LFP_35C_0-100_0.5/1C_c")]["L_efc"],
                  "units_in_lifetime_window": 237, "units_with_temperature_in_window": 0,
                  "first_unit_with_temperature": {"cycle_index": 262, "efc_at_discharge_start": 241.41},
                  "censored_cells_in_secant": 0,
                  "D44_determinability": "D44 defines a secant as not determinable only when a cell in either group is censored; none is",
                  "pipeline_reason": "a 25 to 35 C lever not determinable"},
        "reading_strict_D43": {"T_35_C": None, "s_per_K": None, "P4": "not scored (pipeline)"},
        "reading_T_cyc_over_cells_with_one": {"T_25_C": lf["T_lo_C"], "T_35_C": lf["alt_T_over_cells_with_one"]["T_hi_C"],
                                              "s_per_K": lf["alt_T_over_cells_with_one"]["s_per_K"],
                                              "law2_same_secant_per_K": lf["alt_T_over_cells_with_one"]["law2_s_per_K"],
                                              "bootstrap_ci95_own_seed_unregistered": lf["alt_T_over_cells_with_one"].get("ci95_own_seed"),
                                              **{k: frozen_groups["P4_alternatives"]["T_cyc_over_cells_with_one"][k] for k in ("ratio_to_law2", "s_NCA", "s_NMC", "exceeds_both", "same_sign_as_law2", "within_factor_3", "verdict")}},
        "reading_setpoint_10K": {"s_per_K": lf["s_setpoint_per_K"], "law2_setpoint_25_to_35_per_K": lf["law2_s_setpoint_per_K"],
                                 "law2_logged_cells_with_one_per_K": lf["alt_T_over_cells_with_one"]["law2_s_per_K"],
                                 "ratio_to_law2_setpoint": lf["s_setpoint_per_K"] / lf["law2_s_setpoint_per_K"],
                                 "ratio_to_law2_logged": lf["s_setpoint_per_K"] / lf["alt_T_over_cells_with_one"]["law2_s_per_K"],
                                 **{k: frozen_groups["P4_alternatives"]["setpoint_10K"][k] for k in ("s_NCA", "s_NMC", "exceeds_both", "same_sign_as_law2", "within_factor_3", "verdict")}},
        "reading_drop_the_cell": frozen_groups["P4_alternatives"]["drop_cell_without_T_cyc"],
        "same_three_readings_with_the_pipeline_cell_table": s3_groups.analyse(copy.deepcopy(pr_cells))["P4_alternatives"],
        "pipeline_T_group_35C_in_stage3_lever_json": json.load(open(os.path.join(REPO, "results", "stage3_lever.json")))["registered_1C"]["LFP"]["T_group_C"]["35C"],
        "judgement": ("Defensible as a literal reading but not forced. D43 makes the group temperature the mean T_cyc of the "
                      "uncensored cells, which is undefined here, so the registered logged lever cannot be formed; D48's "
                      "'not determinable' then applies in its ordinary sense. D44's own determinability rule is about "
                      "censoring only, so the scorecard's reason misnames the cause. The case was not foreseen, and the two "
                      "natural repairs give opposite verdicts: mean T_cyc over the three cells that have one gives "
                      "s = 0.1287 per K and P4 holds (ratio to Law 2 2.92, inside [1/3, 3] by 0.08); the 10 K setpoint gives "
                      "s = 0.1382 per K and P4 does not hold (ratio 3.09 to Law 2 at the setpoints, 3.14 to Law 2 at the "
                      "logged temperatures). The note should report all three and say the verdict depends on this reading."),
    }
    sa = frozen_groups["thermal_share"]["a"]["all"]
    sb = frozen_groups["thermal_share"]["b"]["all"]
    or_a = [q for q in frozen_groups["oriented_pairs"]["a"] if q.get("enters_sign_test")]
    point_ii = {
        "question": "P7 scored 'did not hold' on one untied pair out of 60, with 59 ties: is that what D48 says?",
        "D48_text": "P7 holds if the share of D46 (a) over all pairs lies in [0.35, 0.65]",
        "facts": {"pairs_oriented": 60, "ties": sa["n_ties"], "untied": sa["n_entering"], "warmer_shorter": sa["k_warmer_shorter"],
                  "share_measure_a": sa["share"], "share_measure_b_same_pair": sb["share"],
                  "p_one_sided": sa["p_one_sided"], "p_two_sided": sa["p_two_sided"],
                  "the_untied_pair": {k: or_a[0].get(k) for k in ("warmer", "cooler", "dT", "observed", "resolution")} if or_a else None,
                  "untied_pair_lifetimes_efc": {c["cell_id"]: c["L_efc"] for c in frozen_cells if or_a and c["cell_id"] in (or_a[0]["warmer"], or_a[0]["cooler"])},
                  "untied_margin_efc": (lambda Ls: abs(Ls[0] - Ls[1]) - or_a[0]["resolution"])(
                      [c["L_efc"] for c in frozen_cells if c["cell_id"] in (or_a[0]["warmer"], or_a[0]["cooler"])]) if or_a else None},
        "judgement": ("Yes, literally. D48 sets no minimum number of pairs for P7 (it sets five per rate group for P8 only), "
                      "and its 'not scored' clause covers only a lever that is not determinable, which P7 does not use. "
                      "With one untied pair the share can only be 0 or 1, so P7 could not have held whatever the data; the "
                      "sign test gives p = 1.0 one-sided and two-sided, and measure (b) gives the opposite share (1/1) on "
                      "the same pair. The verdict is faithful to D48 but carries no information; the note should say it "
                      "rests on one pair untied by 13 EFC (248.7 EFC difference against a 235.2 EFC resolution)."),
    }
    point_iii = {
        "question": "Measures for 0.5C-discharge groups excluded from the thermal share (D43, D50 (5)): applied consistently?",
        "facts": {"oriented_pairs_a_and_b_contain_0p5C": False, "P6_pairs": 36, "P6_pairs_all_1C": True,
                  "D47_pairs_all": 60, "D47_pairs_NCA": 16, "floor_pairs_total": 86,
                  "cells_with_measures_used_false": "the 24 cells discharged at 0.5C plus LFP 35 C 1C c (no measures)",
                  "P6_if_the_21_partial_window_0p5C_pairs_with_measure_a_were_included": {"pairs": 57, "fraction_below_1K": 0.8947368421052632, "verdict": "holds"}},
        "judgement": ("Yes. No 0.5C-discharge cell appears in either oriented-pair file, in P6, in the D47 pair numbers or in "
                      "any slope; the floor keeps all 86 pairs. Including the 0.5C partial-window pairs, as D48's P6 wording "
                      "alone would, raises P6's fraction to 0.895 and changes no verdict. In stage3_pairs.csv the 26 pairs at "
                      "0.5C carry the rate-group label '1C or below', which no thermal number uses."),
    }

    flags = json.load(open(os.path.join(HERE, "flags.json")))
    out = {
        "meta": {"check": "V2 of DECISIONS.md D49 as registered, plus an independent read of every registered Stage 3 number "
                          "and verdict (D38 to D50, predictions_stage3.md)",
                 "date": "2026-10-23 (UTC)", "repository": REPO, "branch": "stage3", "results_commit_checked": "9c6bb84",
                 "independence": ("Code written from DECISIONS.md, predictions_stage3.md, inputs/sandia/inventory.csv, inputs/floor.json, "
                                  "results/module_spread.csv, results/stage2_lever.csv and .json, and src/thermal_lever.py (Stage 1 laws, "
                                  "imported as allowed). Nothing under src/stage2, src/stage3 or tests was opened, nor /tmp/sandia/steps, "
                                  "verification/stage3/*.json or *.log, results/stage3_* or results/stage3_posthoc.json before the freeze. "
                                  "The raw files were read only, with python3 -I, from scripts outside the data folders."),
                 "freeze": {"frozen_utc": man["frozen_utc"], "outputs_sha256": man["outputs_sha256"],
                            "scripts_sha256_at_freeze": man["scripts_sha256"], "inputs_sha256": man["inputs_sha256"],
                            "software": man["software"]},
                 "pipeline_files_hashes_before_opening": PIPELINE_FILES_BEFORE_OPENING,
                 "pipeline_files_opened_after_freeze_at_utc": "2026-10-23T20:58:59Z",
                 "note_on_commit_767bb9c": ("At 20:42:06 UTC, during this check, commit 767bb9c (not made by me) added early versions of "
                                            "lib.py, s0_cache.py, s1_cells.py and s3_groups.py from this folder; the files here were "
                                            "changed after that and the frozen versions are the ones hashed above. That commit also "
                                            "added results/stage3_posthoc.json and src/stage3/posthoc.py, neither of which I opened, "
                                            "and changed CITATION.cff; it changed none of the inputs or results compared here."),
                 "raw_files_sha256": man["raw_files_sha256"]},
        "summary": {
            "verdict": ("Every registered number and all ten verdicts reproduce. Run with the two readings the pipeline uses "
                        "where D41 and D43 leave a choice, my separate code reproduces all 1,306 group-level fields compared "
                        "(levers, bootstrap intervals, laws, floor, pairs, oriented pairs, thermal share, module, scorecard) and "
                        "the per-cell table to floating-point rounding. Under my own primary readings, frozen before any "
                        "pipeline file was opened, all ten verdicts agree, every number a verdict rests on agrees to 6.5e-5 "
                        "relative or better, the reported intervals to 2.3e-4 and P6's reported median difference to 7.3e-3; "
                        "the larger differences are confined to unscored 2C and D47 numbers and are explained by two readings "
                        "(F1, F2). V2 as registered passes: the three V2 cells agree with the pipeline to 7.1e-15 under either "
                        "reading. Three flags concern verdicts: P4 depends on an unforeseen reading (F5), P7 rests on one pair "
                        "(F6), and the NMC floor that P5 scores is doubled by spurious checks without changing P5 (F3)."),
            "verdicts_mine_frozen": frozen_groups["verdicts"],
            "verdicts_pipeline": {k: ("not scored" if not v["scored"] else ("holds" if v["held"] else "did not hold"))
                                  for k, v in json.load(open(os.path.join(REPO, "results", "stage3_scorecard.json"))).items() if k.startswith("P")},
            "largest_differences": {
                "pipeline_readings_run_groups": grp["pipeline_readings"]["overall"],
                "pipeline_readings_run_cells": {k: {kk: v[kk] for kk in ("max_abs", "argmax_abs", "max_rel", "argmax_rel", "n_compared", "n_nan_mismatch")} for k, v in cmp_pr.items() if k != "checks" and "L_sens" not in k},
                "frozen_primary_groups": grp["frozen_primary"]["summary"],
                "frozen_primary_cells": {k: {kk: v[kk] for kk in ("max_abs", "argmax_abs", "max_rel", "argmax_rel", "n_compared", "n_nan_mismatch")} for k, v in cmp_frozen.items() if k != "checks"},
                "v2_cells_max_abs": max(x.get("abs_diff", 0.0) for r in v2cmp.values() for x in r.values() if isinstance(x, dict) and "abs_diff" in x)},
            "flags_changing_or_possibly_changing_a_scored_number_or_verdict": [f["id"] for f in flags if f.get("changes_scored")]},
        "readings": json.load(open(os.path.join(HERE, "readings.json"))),
        "v2_registered": {"selection": v2["v2_selection_rule"], "cells": v2["v2_cells"], "extra_cells_unregistered": v2["extra_cells"],
                          "pure_python_vs_vectorised": v2["agreement_with_vectorised"], "against_pipeline": v2cmp,
                          "pure_python_records": v2["pure_python"]},
        "d50": {"mine": d50, "against_stage3_prepare": d50_vs_prepare(frozen_cells)},
        "per_cell": {"frozen_primary_vs_pipeline": cmp_frozen, "pipeline_readings_vs_pipeline": cmp_pr},
        "groups": {"frozen_primary_vs_pipeline": {"summary": grp["frozen_primary"]["summary"], "issues": grp["frozen_primary"]["issues"],
                                                  "rows": grp["frozen_primary"]["rows"]},
                   "frozen_primary_own_share_layout_vs_pipeline": frozen_share_vs_pipeline(frozen_groups),
                   "pipeline_readings_vs_pipeline": {"summary": grp["pipeline_readings"]["summary"], "issues": grp["pipeline_readings"]["issues"],
                                                     "rows": grp["pipeline_readings"]["rows"],
                                                     "share_cis_if_every_subset_drew": grp["pipeline_readings"]["share_ci_if_every_subset_draws"]["overall"]}},
        "my_frozen_results": {"levers": frozen_groups["levers"], "rng_stream": frozen_groups["rng_stream"], "floor": frozen_groups["floor"],
                              "thermal_share": frozen_groups["thermal_share"], "P6": frozen_groups["P6"], "module": frozen_groups["module"],
                              "verdicts": frozen_groups["verdicts"], "P4_alternatives": frozen_groups["P4_alternatives"],
                              "kit_beside_nmc": frozen_groups["kit_beside_nmc"], "kit_floor": frozen_groups["kit_floor"]},
        "reading_variants_end_to_end": variants,
        "points": {"i_P4": point_i, "ii_P7": point_ii, "iii_0p5C_exclusion": point_iii},
        "flags": flags,
    }
    with open(OUT, "w") as fh:
        json.dump(clean(out), fh, indent=1)
    print("written", OUT, os.path.getsize(OUT))


if __name__ == "__main__":
    main()
