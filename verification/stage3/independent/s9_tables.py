#!/usr/bin/env python3
"""Stage 3 independent check, step 9: print the report's comparison tables (markdown) from v2_independent.json and the
pipeline files, so that no number in the report is typed by hand.

Usage: python3 -I s9_tables.py > tables.md
"""
import csv, json, math, os

REPO = "/home/claude/thermal-lever-duplicate-floor"
D = json.load(open(os.path.join(REPO, "verification", "stage3", "v2_independent.json")))
PL = json.load(open(os.path.join(REPO, "results", "stage3_lever.json")))
PS = json.load(open(os.path.join(REPO, "results", "stage3_scorecard.json")))
PT = json.load(open(os.path.join(REPO, "results", "stage3_thermal_share.json")))


def e(x, n=3):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "n/a"
    if isinstance(x, str):
        return x
    if x == 0:
        return "0"
    return f"{x:.{n}e}"


def g(x, n=5):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "n/a"
    return f"{x:.{n}g}"


def rows_of(run):
    return {r[0]: r for r in D["groups"][run]["rows"]}


print("### Per cell, 86 cells\n")
print("| Column (mine vs pipeline) | Frozen primary: max abs (cell) | max rel | Pipeline readings: max abs | max rel |")
print("|---|---|---|---|---|")
pf, pr = D["per_cell"]["frozen_primary_vs_pipeline"], D["per_cell"]["pipeline_readings_vs_pipeline"]
for k in pf:
    if k == "checks" or "L_sens" in k:
        continue
    a, b = pf[k], pr[k]
    cell = (a["argmax_abs"] or "")[10:]
    print(f"| {k} | {e(a['max_abs'])} ({cell}) | {e(a['max_rel'])} | {e(b['max_abs'])} | {e(b['max_rel'])} |")
print()
print(f"Checks: frozen primary has the pipeline's check list in {pf['checks']['cells_same_count']} of 86 cells "
      f"(different count: {pf['checks']['cells_diff_count']}); the pipeline-readings run has it in {pr['checks']['cells_same_count']} of 86, "
      f"with largest differences {e(pr['checks']['max_abs_efc'])} EFC and {e(pr['checks']['max_abs_q'])} Ah.\n")

print("### Levers (D44), registered 1C and secondary 2C\n")
print("| Secant | s mine (frozen) | s pipeline | rel diff | 95 % interval mine | interval pipeline | Pipeline-readings run max rel over all lever fields |")
print("|---|---|---|---|---|---|---|")
lev = D["my_frozen_results"]["levers"]
R1, R2 = rows_of("frozen_primary_vs_pipeline"), rows_of("pipeline_readings_vs_pipeline")
for ser in ("1C", "2C"):
    for chem in ("LFP", "NCA", "NMC"):
        for sec in ("25-35", "15-25"):
            m = lev[f"{chem}_{ser}_{sec}"]
            p = PL["registered_1C" if ser == "1C" else "secondary_2C"][chem]["s_25_35" if sec == "25-35" else "s_15_25"]
            if p is None:
                print(f"| {ser} {chem} {sec} C | not formed (T_35 undefined) | not determinable | | | | |")
                continue
            pre = f"lever/{ser}/{chem}/{sec}/"
            mx = max((r[4] for k, r in R2.items() if k.startswith(pre)), default=0)
            rel = abs(m["s_per_K"] - p["s"]) / abs(p["s"])
            print(f"| {ser} {chem} {sec} C | {m['s_per_K']:.7f} | {p['s']:.7f} | {e(rel)} | {m['ci95'][0]:.5f} to {m['ci95'][1]:.5f} | "
                  f"{p['ci95'][0]:.5f} to {p['ci95'][1]:.5f} | {e(mx)} |")
print()
print("| Secant | T_lo C | T_hi C | setpoint s | Law 1 | Law 2 | Law 1 setpoint | Law 2 setpoint |")
print("|---|---|---|---|---|---|---|---|")
for ser in ("1C", "2C"):
    for chem in ("LFP", "NCA", "NMC"):
        for sec in ("25-35", "15-25"):
            p = PL["registered_1C" if ser == "1C" else "secondary_2C"][chem]["s_25_35" if sec == "25-35" else "s_15_25"]
            if p is None:
                pl = PL["registered_1C"][chem]
                print(f"| {ser} {chem} {sec} C (pipeline) | {pl['T_group_C']['25C']:.4f} | undefined (35.0863 over 3 cells) | {lev['LFP_1C_25-35']['s_setpoint_per_K']:.5f} (mine) | | | | |")
                continue
            print(f"| {ser} {chem} {sec} C (pipeline) | {p['T_lo_C']:.4f} | {p['T_hi_C']:.4f} | {p['s_setpoint']:.5f} | "
                  f"{p['laws']['law_1_kucinskis2022']:.5f} | {p['laws']['law_2_schimpe2018']:.5f} | "
                  f"{p['laws_setpoint']['law_1_kucinskis2022']:.5f} | {p['laws_setpoint']['law_2_schimpe2018']:.5f} |")
print()

print("### Floor (D45)\n")
print("| Chemistry | pairs (groups) | mean mine | mean pipeline | median | 95th pct mine | 95th pct pipeline | below resolution | max rel, pipeline-readings run |")
print("|---|---|---|---|---|---|---|---|---|")
fl = D["my_frozen_results"]["floor"]
for chem in ("LFP", "NCA", "NMC"):
    m = fl[chem]
    p = PT["floor"][chem]
    mx = max((r[4] for k, r in R2.items() if k.startswith(f"floor/{chem}/")), default=0)
    print(f"| {chem} | {m['pairs_observed']} ({m['groups_with_pairs']}) | {m['mean']:.6f} | {p['mean']:.6f} | {p['median']:.6f} | "
          f"{m['p95']:.5f} | {p['p95']:.5f} | {m['pairs_below_resolution']} of {m['pairs_observed']} | {e(mx)} |")
print()

print("### Thermal share (D46), sign test\n")
print("| Subset | measure | pairs | ties | untied | warmer shorter | share | p one-sided | p two-sided | same in both runs |")
print("|---|---|---|---|---|---|---|---|---|---|")
for meas in ("a", "b"):
    ps = PT["measure_" + meas]
    for name in ("all", "1C or below", "2C and 3C", "chemistry LFP", "chemistry NCA", "chemistry NMC", "one pair per group (warmest and coolest)"):
        p = ps[name]
        k1 = f"share_{meas}/{name}/"
        same = all(r[3] == 0 or r[5] == "both missing" for k, r in R2.items() if k.startswith(k1) and any(x in k for x in ("/pairs", "/ties", "/warmer", "/share", "/p1", "/p2")))
        same1 = all(r[3] == 0 or r[5] == "both missing" for k, r in R1.items() if k.startswith(k1) and any(x in k for x in ("/pairs", "/ties", "/warmer", "/share", "/p1", "/p2")))
        n = p["warmer_shorter"] + p["warmer_longer"]
        print(f"| {name} | ({meas}) | {p['pairs']} | {p['ties']} | {n} | {p['warmer_shorter']} | {g(p['share'], 3)} | {g(p['sign_test_p_one_sided'], 3)} | "
              f"{g(p['sign_test_p_two_sided'], 3)} | {'yes' if (same and same1) else 'NO'} |")
print()
print("| Slope subset | measure | slope pipeline | interval pipeline | slope, pipeline-readings run | interval, pipeline-readings run | ratio (c) pipeline | my frozen slope (own layout) |")
print("|---|---|---|---|---|---|---|---|")
own = D["my_frozen_results"]["thermal_share"]
for meas in ("a", "b"):
    for name, mk in (("chemistry NCA", "chem_NCA"), ("chemistry NMC", "chem_NMC")):
        p = PT["measure_" + meas][name]
        r_s = R2[f"share_{meas}/{name}/slope"][1]
        r_lo = R2[f"share_{meas}/{name}/ci_lo"][1]
        r_hi = R2[f"share_{meas}/{name}/ci_hi"][1]
        print(f"| {name} | ({meas}) | {p['slope']:.6f} | {p['slope_ci95'][0]:.5f} to {p['slope_ci95'][1]:.5f} | {r_s:.6f} | {r_lo:.5f} to {r_hi:.5f} | "
              f"{p['sd_ratio_pred_over_obs']:.5f} | {own[meas][mk]['slope']:.5f} |")
print()

print("### Module (D47) and P9\n")
pm = PT["module"]
mm = D["my_frozen_results"]["module"]
print("| Quantity | mine (frozen) | pipeline | rel diff |")
print("|---|---|---|---|")
for lab, a, b in (("abs(s_NCA) x 0.391279 K", mm["divergence"], pm["implied_divergence"]),
                  ("interval, low", mm["divergence_ci95_from_s_ci"][0], pm["implied_divergence_ci95"][0]),
                  ("interval, high", mm["divergence_ci95_from_s_ci"][1], pm["implied_divergence_ci95"][1]),
                  ("NCA pairs (16), mean", mm["NCA_pairs"]["mean"], pm["pairs_nca"]["implied_mean"]),
                  ("NCA pairs, max", mm["NCA_pairs"]["max"], pm["pairs_nca"]["implied_max"]),
                  ("all pairs (60), mean", mm["all_pairs"]["mean"], pm["pairs_all"]["implied_mean"]),
                  ("all pairs, max", mm["all_pairs"]["max"], pm["pairs_all"]["implied_max"])):
    print(f"| {lab} | {a:.6f} | {b:.6f} | {e(abs(a - b) / abs(b))} |")
print()

print("### Scorecard (D48)\n")
print("| Prediction | number it rests on, mine (frozen) | pipeline | my verdict | pipeline verdict |")
print("|---|---|---|---|---|")
L = D["my_frozen_results"]["levers"]
V = D["summary"]["verdicts_mine_frozen"]
VP = D["summary"]["verdicts_pipeline"]
fr = D["my_frozen_results"]
items = [("P1", f"NMC 25-35 C s = {L['NMC_1C_25-35']['s_per_K']:.6f}", f"{PS['P1']['s_per_K']:.6f}"),
         ("P2", f"NMC 15-25 C s = {L['NMC_1C_15-25']['s_per_K']:.6f}", f"{PS['P2']['s_15_25_per_K']:.6f}"),
         ("P3", f"NCA 25-35 C s = {L['NCA_1C_25-35']['s_per_K']:.6f}", f"{PS['P3']['s_per_K']:.6f}"),
         ("P4", "LFP 25-35 C s undefined", PS["P4"]["reason"]),
         ("P5", "means " + ", ".join(f"{c} {fr['floor'][c]['mean']:.5f}" for c in ("LFP", "NCA", "NMC")),
          ", ".join(f"{c} {PS['P5']['chemistries'][c]['mean']:.5f}" for c in ("LFP", "NCA", "NMC"))),
         ("P6", f"{fr['P6']['n_below_1K']} of {fr['P6']['n_pairs']} = {fr['P6']['share_below_1K']:.4f}; median {fr['P6']['median_dT_K']:.4f} K",
          f"{PS['P6']['fraction_below_1K']:.4f} of {PS['P6']['pairs']}; median {PS['P6']['median_dT_K']:.4f} K"),
         ("P7", f"{fr['thermal_share']['a']['all']['k_warmer_shorter']} of {fr['thermal_share']['a']['all']['n_entering']} untied, {fr['thermal_share']['a']['all']['n_ties']} ties",
          f"{PS['P7']['warmer_shorter']} of {PS['P7']['warmer_shorter'] + PS['P7']['warmer_longer']}, {PS['P7']['ties']} ties"),
         ("P8", f"untied pairs {fr['thermal_share']['a']['rate_1C_or_below']['n_entering']} (1C or below) and {fr['thermal_share']['a']['rate_2C_and_3C']['n_entering']} (2C and 3C)",
          f"{PS['P8']['pairs_entering_sign_test']['1C or below']} and {PS['P8']['pairs_entering_sign_test']['2C and 3C']}"),
         ("P9", f"{fr['module']['divergence']:.6f}", f"{PS['P9']['implied_divergence']:.6f}"),
         ("P10", ", ".join(f"{c} {L[c + '_1C_15-25']['s_per_K']:.5f}" for c in ("LFP", "NCA", "NMC")),
          ", ".join(f"{c} {PS['P10']['s_15_25_per_K'][c]:.5f}" for c in ("LFP", "NCA", "NMC")))]
for k, a, b in items:
    print(f"| {k} | {a} | {b} | {V[k]} | {VP[k]} |")
print()

print("### D50 counts\n")
d = D["d50"]
print("| Count | mine | pipeline (stage3_prepare.json) | D50 text |")
print("|---|---|---|---|")
tm, tp = d["against_stage3_prepare"]["totals_mine"], d["against_stage3_prepare"]["totals_pipeline"]
txt = {"late_rows": "9,776", "foreign_rows_dropped": "4,313", "discharges": "279,973", "discharges_where_counter_ah_differs": "25,412",
       "discharge_counter_decreases": "275,861", "discharge_counter_decreases_without_cycle_change": "2", "cells_where_counter_ah_differs": "18"}
for k in tm:
    print(f"| {k} | {tm[k]:,} | {tp[k]:,} | {txt.get(k, '')} |")
print()
print(f"Per-cell count mismatches against stage3_prepare.json: {len(d['against_stage3_prepare']['per_cell_count_mismatches'])}; "
      f"largest per-cell difference in discharge ampere-hours {e(d['against_stage3_prepare']['max_abs_diff_discharge_ah_total_Ah'])} Ah, "
      f"in the registered-rule total {e(d['against_stage3_prepare']['max_abs_diff_counter_ah_total_Ah'])} Ah.\n")
print("| Cell | registered rule share of discharge ampere-hours |")
print("|---|---|")
for k, v in d["mine"]["registered_rule_share_of_discharge_ah"].items():
    print(f"| {k[10:]} | {100 * v:.2f} % |")
