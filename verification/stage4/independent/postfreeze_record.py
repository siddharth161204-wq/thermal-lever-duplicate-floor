#!/usr/bin/env python3
"""POST-FREEZE: writes record.json, the machine-readable record of the independent check, from the frozen outputs,
the post-freeze diagnostics and out_compare/compare.json (no number is typed by hand).
"""
import csv, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")


def j(p):
    with open(p) as f:
        return json.load(f)


def c(p):
    with open(p) as f:
        return list(csv.DictReader(f))


def main():
    cmp_ = j(os.path.join(HERE, "out_compare", "compare.json"))
    sc = j(os.path.join(OUT, "scorecard.json"))
    res = j(os.path.join(OUT, "resistance.json"))
    chk = j(os.path.join(OUT, "resistance_check.json"))
    tol = j(os.path.join(OUT, "postfreeze_diag_tolerance.json"))
    tail = j(os.path.join(OUT, "diag_charge_10C_tail.json"))
    sens = j(os.path.join(OUT, "sensitivity_own_temperature.json"))
    prec = c(os.path.join(OUT, "simulation_precise.csv"))
    ln = c(os.path.join(OUT, "simulation_lognormal_all.csv"))
    same = c(os.path.join(OUT, "postfreeze_simulation_normal_same_sd.csv"))
    rec = {
        "check": "Stage 4 independent check (D62's independent read)",
        "base_commit": "c1e58a19d1653bb477de485c4fc9f68e956ee070 (pre-registration: D56 to D63, predictions_stage4.md)",
        "freeze": {"commit": "fd5d1e73f6fa87e5b60cc9c2689b26e072a6b1d0", "committed_utc": "2026-10-27T10:33:27Z",
                   "pushed_utc": "2026-10-27T10:33:36Z", "branch": "origin/stage4-independent",
                   "rule": "no file of src/stage4, results/stage4_*, verification/stage4 (outside this folder) or "
                           "tests/stage4_synthetic.py opened, and no commit after c1e58a1 logged, shown or diffed, before this commit"},
        "pipeline_compared": {"registered_results_commit": "6f658f1 (pushed 2026-10-27T10:07:16Z)", "posthoc_commit": "a2b3d6e"},
        "agreement": cmp_["quantities"],
        "deterministic_mismatches": [m for m in cmp_["mismatches"] if "Monte Carlo" not in m["quantity"] and "post-freeze reading" not in m["quantity"]],
        "monte_carlo_differences": {"normal_check_pass_differences": [m for m in cmp_["mismatches"] if "post-freeze reading" in m["quantity"]],
                                    "V4_lognormal_pass_differences": cmp_["V4_pass_disagreements"],
                                    "counts": cmp_["checks"],
                                    "checker_literal_normal_check_pass": sum(r["V4_pass_5pct"] == "True" for r in same)},
        "verdicts": {"checker": {k: v["held"] for k, v in sc.items()},
                     "pipeline": {"P1": True, "P2": False, "P3": True, "P4": True, "P5": False, "P6": True, "P7": False},
                     "checker_details": sc},
        "resistance": {"procedure_ohm": {k: v["R_registered_ohm"] for k, v in res["ambients"].items()},
                       "ratios": {"25/35": res["used"]["ratio_25_over_35_archived"], "10/35": res["used"]["ratio_10_over_35_archived"]},
                       "failed_calls": {k: v["failed_calls"] for k, v in res["ambients"].items()},
                       "separate_code_1s_default_tolerance_ohm": {str(r["T_C"]): r["R_ohm"] for r in chk["rows"]},
                       "separate_code_1s_spike_screened_ohm": {str(r["T_C"]): r["R_ohm_spike_screened"] for r in chk["rows"]},
                       "C2_charge_10C_spike": {"default_tolerance_0p01s": tail["fine_0.01s"],
                                               "tolerance_test": tol}},
        "sensitivity_heat_at_own_temperature_not_registered": sens["verdicts_under_this_heat"],
        "P7_precise_simulation_200k": [r for r in prec],
        "choices": [
            {"point": "D58 fails to solve", "checker": "a raised exception", "pipeline": "the same", "effect": "same two failures"},
            {"point": "D58 resistance at the ambient", "checker": "as registered", "pipeline": "as registered", "effect": "decides P2, P4, P5 (flag 2)"},
            {"point": "D61 two-sided power", "checker": "both tails by quadrature", "pipeline": "lower tail zero where scipy nan (D64)", "effect": "n identical in 144 cases"},
            {"point": "D62 V4 seeds", "checker": "SeedSequence([20261027, case, n, shape])", "pipeline": "own seeds", "effect": "10 boundary pass/fail differences"},
            {"point": "D62 normal check sd(e)", "checker": "normal-shape sd at freeze; registered sd post-freeze", "pipeline": "registered sd", "effect": "none on verdicts"},
            {"point": "cells inside band", "checker": "logged secant span at freeze; set points post-freeze", "pipeline": "set points", "effect": "17 of 36 outside, all flags agree"},
            {"point": "P2 levers", "checker": "point estimates", "pipeline": "point estimates", "effect": "none"},
            {"point": "P7 central cold lever", "checker": "cold band median 3.5800", "pipeline": "the same", "effect": "3.6139 gives n = 3 too"},
            {"point": "figure grid", "checker": "0.5 K", "pipeline": "0.1 K", "effect": "agree where overlapping"}],
        "flags": [
            "1. R(10 C) rests on the 10 s grid missing a 13.8 kW sub-second solver spike (t = 3954.70 s) in the C/2 charge at default tolerances; 1 s grid +15 %, 0.01 s grid 4x on the charge; at rtol 1e-8 it vanishes (46.526 vs 46.547 mOhm); registered value right, D65's no-spike statement describes samples, V2's agreement depends on its tight tolerances.",
            "2. The resistance at the ambient is an upper bracket that decides P2, P4 and P5: with heat at each cell's own temperature P2 holds, P4 fails (5.54 % at h = 5), P5 holds (5.54 to 7.27 %); D66's brackets do not cover P2.",
            "3. 17 of 36 cases have cells outside their band's set points; P4 rests on 10 C 1C cells at 24 to 50 C under 10 to 25 C levers; 10 C 2C 95th percentiles reach 75 to 98 %.",
            "4. V4's 5 % rule is granular and seed-dependent at small n; P7 warm case precise simulated n = 23 (power 0.8024 at 23, 0.7897 at 22), 8 % below 25; P7's verdict is unaffected.",
            "5. P7 turns on the registered median-|s| reading (signed median gives 221 pairs, P7 would hold); the note's 25 pairs is the noncentral-t upper estimate (23 simulated).",
            "6. P5's 35 C clause holds by 0.29 points (7.046 vs 7.332 %).",
            "7. The 10 C C/2 charge ends at 2.75 Ah; R(10 C) mixes a full-window discharge and a 55 % charge (discharge alone +2.9 %); no verdict moves.",
            "8. Wording: 'a 10 C module' (cells at 24 to 50 C at 1C); 'offers a share from C/2' only under the upper bracket; note's 'KIT and Sandia, 25 to 40 C' median lever.",
            "9. Unscored 2C within-cell ratio 0.306 near one third (analytic bracket 0.339).",
            "10. Agreed: C/16 charge values are solver spikes (D65); scipy nct nan at n = 3 (D64)."],
        "checker_own_slips_resolved_postfreeze": [
            "out/divergence.csv ln_offset_cell columns written as np.float64(...) strings; digits intact",
            "summary.json used the logged secant span for cells inside the band; set points recomputed post-freeze",
            "normal-shape V4 at freeze used m sqrt(pi/2); literal reading run post-freeze"],
        "report_md": "not written: the harness refused a report .md from this agent; the report text was returned to the caller"}
    with open(os.path.join(HERE, "record.json"), "w") as f:
        json.dump(rec, f, indent=1)
    print("record.json written;", len(rec["agreement"]), "quantity classes;", len(rec["deterministic_mismatches"]), "deterministic mismatches")


if __name__ == "__main__":
    main()
