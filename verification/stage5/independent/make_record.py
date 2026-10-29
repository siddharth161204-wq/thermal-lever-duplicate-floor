#!/usr/bin/env python3
"""Independent Stage 5 check: record.json (base commit, freeze rule, readings, seeds, headline numbers, flags and the
SHA-256 of every output), built from out/*.json so that the record cannot differ from the outputs.

Usage: PYTHONDONTWRITEBYTECODE=1 python3 verification/stage5/independent/make_record.py
"""
import json, os, platform, sys

import numpy, scipy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ind_common as C   # noqa: E402


def j(name):
    return json.load(open(os.path.join(C.OUT, name)))


def r6(x):
    return None if x is None else float(f"{x:.6g}")


def main():
    P = j("pairs.json")
    B = j("B.json")["analyses"]
    Cj = j("C.json")
    S = j("scorecard.json")
    X = j("B_crosscheck.json")
    V3 = j("v3_hand.json")
    S4 = j("crosscheck_stage4_power.json")
    MC = j("C_mc_spread.json")
    INF = j("influence.json")
    c1 = {(r["case"], r["floor"]): r for r in Cj["C1"]}
    c2 = {(r["dataset"], r["band"], r["floor"]): r for r in Cj["C2"]}
    mc1 = {(r["case"], r["floor"]): r for r in MC["C1"]}

    def b_line(name):
        a = B[name]
        out = {"pairs": a["n_pairs"], "sets": a["n_sets"]}
        if a.get("beta") is not None:
            out["beta"] = a["beta"]
            out["beta_ci95"] = a["beta_ci95"]
            out["resid_sd"] = a["resid_sd"]
            out["resid_sd_ci95"] = a["resid_sd_ci95"]
        if a.get("s_w_pct_per_K") is not None:
            out["s_w_pct_per_K"] = a["s_w_pct_per_K"]
            out["s_w_ci95_pct_per_K"] = a["s_w_ci95_pct_per_K"]
        if "sign" in a:
            out["sign"] = {k: a["sign"][k] for k in ("agree", "pairs_tested", "left_out", "share", "p_two_sided")}
        if "floor_verdicts" in a:
            out["floor_verdicts"] = {k: v["inside_interval"] for k, v in a["floor_verdicts"].items()}
        return out

    readings = [
        "Sets for every bootstrap are KIT parameter sets (one protocol at one temperature) and Sandia groups (chemistry, "
        "ambient, window, rates), as D74 reads the brief's 'resampling protocols, not pairs'; the pooled draws are "
        "stratified (KIT sets and Sandia groups drawn separately). Resampling protocols instead is reported as a "
        "sensitivity (out/B.json, 'sens ... protocol clusters').",
        "The 25 C band is the one whose archived interval (hi - lo) is narrower, protocol by protocol for KIT and chemistry "
        "by chemistry for Sandia; a lever that is not determinable has no interval and cannot be the narrower (LFP 25 C "
        "goes to 15 to 25 C). Sandia 15 C sets take the 15 to 25 C lever and 35 C sets the 25 to 35 C lever.",
        "Every Sandia pair takes its chemistry's registered 1C lever for its band, the 2C and 3C groups included "
        "(D72's literal wording); the secondary 2C series is not used.",
        "KIT tie: |L_warm - L_cool| strictly below the larger res90_efc of the two cells. Sandia tie: below_resolution as "
        "archived.",
        "Orientation: KIT as archived (warm = higher T_C, checked against stage2_cells.csv); Sandia by the larger rise "
        "(measure (a)) and, for the measure (b) bracket, by the larger rest temperature; no exact tie occurs in either.",
        "d = ln(L_warm / L_cool) as the log of the ratio; p = -s dT with s in natural-log units per K.",
        "Untied bracket: the dataset restricted to its untied pairs, the bootstrap drawing the sets that hold one (2 KIT "
        "sets); drawing all 36 KIT sets (sets without an untied pair adding nothing, empty draws left out and counted) "
        "gives the same interval and is reported as a sensitivity. Pooled untied equals KIT untied (Sandia has no "
        "untied scored pair).",
        "B4's interval uses B1's draws: in each draw beta is refitted over the drawn pairs, the residuals of those pairs "
        "give sd* = sqrt(sum r^2 / (n* - 1)) with n* the draw's pair count. 'Lies inside that interval' is read as the "
        "closed interval.",
        "B2 for the 0 C KIT sets uses the same bootstrap over its 12 sets.",
        "B3 is computed for KIT, Sandia and pooled, each over all scored pairs and per band; a pair with d = 0 or p = 0 "
        "would be left out (none is).",
        "Extras not named in D74, computed from the same draws and labelled so: beta and B4 per band (D76 mentions 'B1 to "
        "B4 per dataset and band'), s_w over a whole dataset.",
        "C1 calculation: noncentral t with df n - 1 and nc = |s| sqrt(sum dT^2) / sd(e), both tails, from scipy.stats.nct "
        "directly; where scipy returns nan for the lower tail the quadrature value is used (this happened only at "
        "Brent's trial points nc = 11.8 and 23.7, never at a reported value). Brent's method on nc, then |s| = nc sd(e) / "
        "sqrt(sum dT^2).",
        "C1 simulation: one set of 10,000 datasets per case and shape; rejection |t| > t_0.975(n - 1); 'reaches 0.80' as "
        "simulated power >= 0.80; bisection on |s| from [0, 2 s_calc] (doubling if needed) until hi - lo <= 1e-6 hi, "
        "the reported value being hi; the exact step of the simulated power is reported beside. The random sign is "
        "+1 or -1 with probability 1/2 each.",
        "The normal-e simulation uses the registered shape's sd(e); the normal-shape calculation uses m sqrt(pi / 2).",
        "Shared-cell F = (sum over cells of w^2) / (2 sum over pairs of dT^2) per selection; the pooled F sums over the "
        "cells of both datasets; the shared-cell line keeps df = n - 1.",
        "C2's simulated power is C1's per-band simulation (the same 10,000 datasets) evaluated at the band's median |s|.",
        "C2's simulated n follows Stage 4's bisection (2 first; then lo = 2, hi = max(4 n_calc, 3), doubling while hi "
        "fails), fresh 10,000 datasets at each n tried, each drawing n dT values with replacement from the band's "
        "observed dT and e from the registered shape.",
        "C2 is computed for every dataset, band and floor; D75 names KIT under the KIT floor and Sandia and pooled under "
        "both floors; KIT under the Geslin floor is an extra.",
        "P6 compares the detectable |s| with the full-precision band medians (1.302942364773262 and 3.5800471015633453 % "
        "per K), 'exceeds' strict. P1 and P2 read the interval as closed; P7 as width ratio >= 2; P8 as upper end <= 3.",
        "V3: D77's 'first Sandia NMC pair at 25 C in the order of results/stage3_pairs.csv' is literally the 0-100 %, "
        "0.5/0.5C pair, which has no measure (a) and is not scored; its parenthesis names the 0-100 %, 0.5/1C pair of "
        "cells a and b, which is the one computed.",
    ]
    seeds = {
        "base_seed": C.BASE_SEED,
        "generator": "numpy.random.Generator(numpy.random.SFC64(numpy.random.SeedSequence([5550129, *key])))",
        "B (partb.py)": {"KIT": [1], "Sandia": [2], "pooled": [3], "KIT warm": [4], "KIT cold": [5], "Sandia warm": [6],
                         "Sandia cold": [7], "pooled warm": [8], "pooled cold": [9], "KIT 0 C": [10], "untied KIT": [11],
                         "untied KIT warm": [12], "untied pooled": [13], "untied pooled warm": [14],
                         "Sandia measure b": [15], "Sandia measure b warm": [151], "Sandia measure b cold": [152],
                         "Sandia without NCA 15C 1C": [16], "sens untied all sets": [17],
                         "sens protocol clusters": [20, 21, 22, 23, 24]},
        "B replicates (crosscheck_b.py)": "[900 + i, r] for 200 replicates, [950 + i] for 1,000,000 draws",
        "C1 (partc.py)": "[100 + case index, floor index (0 Geslin, 1 KIT), 1 log-normal or 2 normal]",
        "C2 simulated n (partc.py)": "[200 + 10 band index + dataset index, floor index, 1, n]",
        "C replicates (mc_spread_c.py)": "C1 [700 + case, floor, r]; C2 [800 + (band == cold), r, n]",
        "not used": "the pipeline's default_rng([20261029, k]) and V2's random.Random(20261030)",
    }
    headline = {
        "counts": {"KIT": {"in_file": 142, "scored": P["counts"]["KIT"]["scored_pairs"], "sets": P["counts"]["KIT"]["scored_sets"],
                           "warm": P["counts"]["KIT"]["by_band"]["warm"], "cold": P["counts"]["KIT"]["by_band"]["cold"],
                           "not_scored": P["counts"]["KIT"]["not_scored_by_reason"]},
                   "Sandia": {"in_file": 86, "with_measure_a": P["counts"]["Sandia"]["pairs_with_measure_a_both"],
                              "scored": P["counts"]["Sandia"]["scored_pairs"], "sets": P["counts"]["Sandia"]["scored_sets"],
                              "warm": P["counts"]["Sandia"]["by_band"]["warm"], "cold": P["counts"]["Sandia"]["by_band"]["cold"],
                              "not_scored": P["counts"]["Sandia"]["not_scored_by_reason"]},
                   "pooled": P["counts"]["pooled"]},
        "B": {k: b_line(k) for k in ("KIT", "Sandia", "pooled", "KIT warm", "KIT cold", "Sandia warm", "Sandia cold",
                                     "pooled warm", "pooled cold", "KIT 0 C", "untied KIT", "untied pooled",
                                     "Sandia measure b", "Sandia without NCA 15C 1C")},
        "C1": {f"{k[0]} | {k[1]}": {"pairs": v["n_pairs"], "calc_pct_per_K": v["detectable_s_calc_pct_per_K"],
                                    "sim_lognormal_pct_per_K": v["detectable_s_sim_lognormal_pct_per_K"],
                                    "v4_rel_diff": v["v4_rel_diff"], "v4_pass": v["v4_pass"],
                                    "sim_normal_same_sd_pct_per_K": v["detectable_s_sim_normal_same_sd_pct_per_K"],
                                    "calc_normal_shape_pct_per_K": v["detectable_s_calc_normal_shape_pct_per_K"],
                                    "F": v["F_shared_cell"], "calc_shared_cell_pct_per_K": v["detectable_s_calc_shared_cell_pct_per_K"]}
               for k, v in c1.items()},
        "C2": {f"{k[0]} {k[1]} | {k[2]}": {"named": v["named_in_D75"], "power_calc": v["power_calc"],
                                           "power_sim_lognormal": v["power_sim_lognormal"], "n_calc": v["n_required_calc"],
                                           "n_sim_lognormal": v["n_required_sim_lognormal"], "n_shared_cell": v["n_required_shared_cell"]}
               for k, v in c2.items()},
        "verdicts": S["verdicts"],
    }
    kb, pb = X["mc_spread"]["KIT"]["beta"], X["mc_spread"]["pooled"]["beta"]
    pw = X["mc_spread"]["pooled warm"]["s_w"]
    sn = X["mc_spread"]["Sandia without NCA 15C 1C"]["beta"]
    v4fail = [f"{k[0]} | {k[1]} ({v['v4_rel_diff']:+.4f}; passes in {mc1[k]['v4_pass_share']:.0%} of 20 replicate "
              f"simulations, mean {mc1[k]['rel_diff_mean']:+.4f})" for k, v in c1.items() if not v["v4_pass"]]
    near5 = [f"{k[0]} | {k[1]} (replicate pass share {mc1[k]['v4_pass_share']:.2f}, mean {mc1[k]['rel_diff_mean']:+.4f})"
             for k in c1 if 0.0 < mc1[k]["v4_pass_share"] < 1.0]
    mcc2 = {x["band"]: x for x in MC["C2"]}
    flags = [
        f"Ties: by D72's rule 140 of the 142 KIT pairs (104 of the 106 scored) are ties, so the KIT untied bracket holds "
        f"two pairs (P049-3/P049-2 and P052-3/P052-1, both 25 C, warm band) in two sets and no cold pair; its numbers "
        f"(beta {B['untied KIT']['beta']:.1f}, s_w {B['untied KIT']['s_w_pct_per_K']:.0f} % per K, residual sd "
        f"{B['untied KIT']['resid_sd']:.3f}) carry no information, and the pooled untied bracket is the same two pairs. "
        f"D72 reports the empty Sandia bracket but not the near-empty KIT one.",
        f"Monte Carlo noise in the intervals: this run's KIT beta interval [{kb['registered_run_ci95'][0]:.3f}, "
        f"{kb['registered_run_ci95'][1]:.3f}] lies {abs(kb['registered_run_ci95'][0]-kb['lo_mean'])/kb['lo_sd']:.1f} and "
        f"{abs(kb['registered_run_ci95'][1]-kb['hi_mean'])/kb['hi_sd']:.1f} Monte Carlo sd from the mean of 200 replicate "
        f"10,000-draw bootstraps [{kb['lo_mean']:.3f}, {kb['hi_mean']:.3f}] (end sd {kb['lo_sd']:.3f}, {kb['hi_sd']:.3f}); "
        f"1,000,000 draws give [{kb['big_ci95'][0]:.3f}, {kb['big_ci95'][1]:.3f}]. The pooled lower end "
        f"{pb['registered_run_ci95'][0]:.3f} lies {abs(pb['registered_run_ci95'][0]-pb['lo_mean'])/pb['lo_sd']:.1f} sd "
        f"from its replicate mean {pb['lo_mean']:.3f}. The seed was fixed before the run and is kept; the pipeline's "
        f"intervals should be compared within several Monte Carlo sd, not to rounding. P1, P2 and P8 hold in all 200 "
        f"replicates and at 1,000,000 draws.",
        f"Pooled warm s_w: this run's interval [{100*pw['registered_run_ci95'][0]:.3f}, {100*pw['registered_run_ci95'][1]:.3f}] "
        f"% per K excludes zero, but {pw['contains_share']['0.0']:.1%} of the replicates and the 1,000,000-draw interval "
        f"[{100*pw['big_ci95'][0]:.3f}, {100*pw['big_ci95'][1]:.3f}] include it: any statement that the pooled warm "
        f"within-set lever is negative rests on Monte Carlo noise (unscored).",
        f"B4, pooled against the KIT floor: the interval's lower end {B['pooled']['resid_sd_ci95'][0]:.4f} lies "
        f"{B['pooled']['resid_sd_ci95'][0]-C.SD_E_TASK['KIT']:.4f} above sd(e) = {C.SD_E_TASK['KIT']:.5f}, so the pooled "
        f"noise is recorded as not the KIT floor; stable over Monte Carlo (no replicate includes it; 1,000,000 draws give "
        f"{X['mc_spread']['pooled']['sd']['big_ci95'][0]:.4f}) but a margin of 1.4 %. KIT and Sandia include the KIT "
        f"floor's sd(e) and exclude Geslin's in every replicate. For KIT the comparison is partly circular (D74): the KIT "
        f"floor comes from cells 1 and 2 of these same sets.",
        f"P3 fails on the point estimate: the KIT warm s_w is {B['KIT warm']['s_w_pct_per_K']:.2f} % per K "
        f"[{B['KIT warm']['s_w_ci95_pct_per_K'][0]:.1f}, {B['KIT warm']['s_w_ci95_pct_per_K'][1]:.1f}]. In the same band the "
        f"registered per-band sign test (unscored) gives {B['KIT warm']['sign']['agree']} of {B['KIT warm']['sign']['pairs_tested']} "
        f"agreeing, two-sided p = {B['KIT warm']['sign']['p_two_sided']:.4f}, and the per-band beta (an extra, not a B1 "
        f"quantity) is {B['KIT warm']['beta']:.2f} [{B['KIT warm']['beta_ci95'][0]:.2f}, {B['KIT warm']['beta_ci95'][1]:.2f}], "
        f"excluding both 0 and 1: the KIT beta near zero ({B['KIT']['beta']:.3f}) mixes a negative warm band with a cold band "
        f"at {B['KIT cold']['beta']:.2f} [{B['KIT cold']['beta_ci95'][0]:.2f}, {B['KIT cold']['beta_ci95'][1]:.2f}].",
        f"P4 holds on a point estimate of {B['KIT cold']['s_w_pct_per_K']:.2f} % per K whose interval "
        f"[{B['KIT cold']['s_w_ci95_pct_per_K'][0]:.1f}, {B['KIT cold']['s_w_ci95_pct_per_K'][1]:.1f}] spans zero widely: "
        f"deterministic, but the sign carries little information.",
        f"P7 fails by a wide margin (width ratio {S['details']['P7']['ratio']:.3f}; replicate means "
        f"{X['mc_spread']['Sandia']['beta']['width_mean']:.3f} and {X['mc_spread']['KIT']['beta']['width_mean']:.3f}; "
        f"1,000,000 draws {S['details']['P7']['big_draw_ratio']:.3f}): the Sandia swing differences (dT rms "
        f"{Cj['design']['Sandia']['dT_rms_K']:.2f} K against {Cj['design']['KIT']['dT_rms_K']:.2f} K for KIT) give its 56 "
        f"pairs more leverage than KIT's 106; the ties the prediction names do not decide the width.",
        f"Sandia without the NCA 15 C 1C pair (bracket): beta {B['Sandia without NCA 15C 1C']['beta']:.3f} "
        f"[{B['Sandia without NCA 15C 1C']['beta_ci95'][0]:.3f}, {B['Sandia without NCA 15C 1C']['beta_ci95'][1]:.3f}] "
        f"excludes zero in every replicate (lower end sd {sn['lo_sd']:.3f}; 1,000,000 draws {sn['big_ci95'][0]:.3f}), "
        f"against {B['Sandia']['beta']:.3f} [{B['Sandia']['beta_ci95'][0]:.3f}, {B['Sandia']['beta_ci95'][1]:.3f}] with it: "
        f"one pair whose cell b sensor did not follow the cell (dT 3.98 K, the largest in its band) moves the Sandia "
        f"interval off zero. Unscored.",
        f"Sandia measure (b) orients {P['counts']['Sandia']['scored_orientation_b_differs']} of the 56 pairs the other way "
        f"from measure (a); beta {B['Sandia measure b']['beta']:.3f} [{B['Sandia measure b']['beta_ci95'][0]:.3f}, "
        f"{B['Sandia measure b']['beta_ci95'][1]:.3f}].",
        "Pooling mixes two temperature measures: KIT dT is the difference of mean cell temperatures, Sandia dT the "
        "difference of within-cycle swings (D43, D72); and the Sandia 2C and 3C groups carry the 1C lever (D72), though "
        "Stage 3's secondary 2C levers differ (NCA 25 to 35 C +1.30 against -2.04 % per K at 1C; NMC -7.14 against -2.44).",
        f"V4 (C1): within 5 % in {Cj['v4']['pass']} of {Cj['v4']['cases']} cases with this run's seeds; it fails in "
        + ("; ".join(v4fail) if v4fail else "none") + ". Seed-dependent V4 verdicts (replicate pass share strictly "
        "between 0 and 1): " + ("; ".join(near5) if near5 else "none") + ". P6's case passes in every replicate.",
        f"P6 holds with margin: KIT calculated detectable |s| {c1[('KIT', 'KIT')]['detectable_s_calc_pct_per_K']:.3f} % per K "
        f"(simulated {c1[('KIT', 'KIT')]['detectable_s_sim_lognormal_pct_per_K']:.3f}) against 3.580; under the Geslin "
        f"floor, beside it, the value is {c1[('KIT', 'Geslin')]['detectable_s_calc_pct_per_K']:.3f} % per K, below the cold "
        f"median, so the verdict rests on the KIT floor being the one scored.",
        f"C1's calculation keeps df = n - 1 and independent pairs, while the pairs of a set share cells (F = "
        f"{Cj['design']['KIT']['F']:.4f} KIT, {Cj['design']['Sandia']['F']:.4f} Sandia, {Cj['design']['pooled']['F']:.4f} "
        f"pooled); the shared-cell line raises KIT's detectable |s| to "
        f"{c1[('KIT', 'KIT')]['detectable_s_calc_shared_cell_pct_per_K']:.3f} % per K, still at df = n - 1.",
        f"C2: with dT resampled from a skewed band the simulated n exceeds the calculation at small n: Sandia cold under "
        f"the Geslin floor {c2[('Sandia', 'cold', 'Geslin')]['n_required_sim_lognormal']} against "
        f"{c2[('Sandia', 'cold', 'Geslin')]['n_required_calc']}, pooled cold "
        f"{c2[('pooled', 'cold', 'Geslin')]['n_required_sim_lognormal']} against {c2[('pooled', 'cold', 'Geslin')]['n_required_calc']}. "
        f"The named cases agree closely: KIT warm {c2[('KIT', 'warm', 'KIT')]['n_required_sim_lognormal']} against "
        f"{c2[('KIT', 'warm', 'KIT')]['n_required_calc']} (5 replicates {mcc2['warm']['replicate_n_sim']}), KIT cold "
        f"{c2[('KIT', 'cold', 'KIT')]['n_required_sim_lognormal']} against {c2[('KIT', 'cold', 'KIT')]['n_required_calc']} "
        f"(5 replicates {mcc2['cold']['replicate_n_sim']}).",
        f"Influence (checker's leave-one-set-out sensitivity, out/influence.json): one set carries most of each fit. KIT "
        f"P036 (10 C, 10-100 1.67/1.00C, lever -11.59 % per K) holds {INF['KIT']['top_sets_by_p2_share'][0]['p2_share']:.0%} "
        f"of KIT's sum of p^2; without it KIT beta is {INF['KIT']['without_top_set']['beta']:.3f} "
        f"[{INF['KIT']['without_top_set']['beta_ci95'][0]:.2f}, {INF['KIT']['without_top_set']['beta_ci95'][1]:.2f}]. "
        f"Sandia NMC 25 C 0-100 0.5/3C (a 3C group carrying the 1C lever; its cell c has the smallest swing, 6.0 K against "
        f"10.1 to 11.3 K, and the shortest life) holds {INF['Sandia']['top_sets_by_p2_share'][0]['p2_share']:.0%} of "
        f"Sandia's and {INF['pooled']['top_sets_by_p2_share'][0]['p2_share']:.0%} of the pooled sum of p^2 (Sandia as a "
        f"whole 72 %); without it Sandia beta is {INF['Sandia']['without_top_set']['beta']:.3f} "
        f"[{INF['Sandia']['without_top_set']['beta_ci95'][0]:.2f}, {INF['Sandia']['without_top_set']['beta_ci95'][1]:.2f}] "
        f"and pooled {INF['pooled']['without_top_set']['beta']:.3f} [{INF['pooled']['without_top_set']['beta_ci95'][0]:.2f}, "
        f"{INF['pooled']['without_top_set']['beta_ci95'][1]:.2f}], against {B['Sandia']['beta']:.3f} and "
        f"{B['pooled']['beta']:.3f} with it. No verdict changes (P7's ratio would be about 1.3), but the Sandia and "
        f"pooled point values describe one group more than the dataset.",
        f"C2 simulated n, Monte Carlo: this run's KIT cold value {c2[('KIT', 'cold', 'KIT')]['n_required_sim_lognormal']} "
        f"is below all five replicates ({mcc2['cold']['replicate_n_sim']}): at n = 198 this run's simulated power is "
        f"{c2[('KIT', 'cold', 'KIT')]['sim_power_at_n_sim']:.4f}, while other seeds average about 0.793 there, so the "
        f"bisection's first crossing landed a few pairs low; the simulated n carries a spread of a few pairs near 200 and "
        f"about 1 % near 4,200.",
        "D77's V3 wording: the literal first NMC pair at 25 C in stage3_pairs.csv is the 0.5/0.5C pair without measure (a); "
        "the parenthesis names the 0.5/1C pair, which this check computed.",
        "Brief against D72: the brief's A1 gives every one of the 142 KIT pairs a lever and P1's expectation cites 142 "
        "pairs at a median 0.6 K; D72 leaves the 36 pairs at 0 C unscored (no secant reaches 0 C), so the scored test has "
        "106, and 0.59 K is the median within-set temperature range, while the median pair dT is 0.33 K over the 142 "
        "pairs and 0.26 K over the 106 scored. The predictions file states both.",
        "Brief against D74: the brief resamples 'protocols, not pairs'; D74 resamples sets. With protocols as clusters "
        f"(the 12 KIT protocols; the 11 Sandia chemistry, window and rate combinations) the KIT interval is [{B['sens KIT, protocol clusters']['beta_ci95'][0]:.3f}, "
        f"{B['sens KIT, protocol clusters']['beta_ci95'][1]:.3f}], Sandia [{B['sens Sandia, protocol clusters']['beta_ci95'][0]:.3f}, "
        f"{B['sens Sandia, protocol clusters']['beta_ci95'][1]:.3f}], pooled [{B['sens pooled, protocol clusters']['beta_ci95'][0]:.3f}, "
        f"{B['sens pooled, protocol clusters']['beta_ci95'][1]:.3f}]: no verdict would change.",
        f"scipy's noncentral t returned nan for the lower tail {Cj['nct_lower_tail_nan_count']} times, all at Brent's trial "
        "points nc = 11.8 and 23.7 far from the root (about 2.8), filled by quadrature; no reported value is affected "
        "(D64's issue, harmless here).",
    ]
    outputs = {}
    for fn in sorted(os.listdir(C.OUT)):
        outputs["out/" + fn] = C.sha256_file(os.path.join(C.OUT, fn))
    for fn in sorted(os.listdir(C.HERE)):
        if fn.endswith((".py", ".sh")):
            outputs[fn] = C.sha256_file(os.path.join(C.HERE, fn))
    rec = {
        "check": "Stage 5 independent check (D77's independent read)",
        "base_commit": C.BASE_COMMIT + " (pre-registration: DECISIONS.md D71 to D78 and predictions_stage5.md)",
        "branch": "stage5-independent",
        "freeze": {"rule": "Until the commit that adds this record (message starting 'Stage 5 independent check: outputs "
                           "frozen') was committed and pushed to origin/stage5-independent: nothing under src/stage5/, "
                           "results/stage5_*, verification/stage5/ outside this folder, or tests/stage5_synthetic.py was "
                           "opened, listed or imported, and no commit after da613a1 was logged, shown, diffed or fetched "
                           "(no git fetch, no git log of origin/stage5). Stage 1 to 4 code and results were read; "
                           "src/stage4/experiment.py was imported only after the checker's own power code had run, as a "
                           "cross-check (crosscheck_stage4_power.py).",
                   "commit_and_push_time": "given in the hand-back message; a commit cannot hold its own hash"},
        "environment": {"python": platform.python_version(), "numpy": numpy.__version__, "scipy": scipy.__version__},
        "inputs": P["input_hashes"],
        "checks": {"v1_oriented_measure_a": P["v1_oriented_measure_a"],
                   "kit_dT_and_lifetimes_against_stage2_cells": "exact (stop otherwise)",
                   "sd_e_rebuilt": {k: v["sd_e"] for k, v in P["shapes_rebuilt"].items()},
                   "sd_e_equal_to_task_values": {k: v["sd_e"] == C.SD_E_TASK[k] for k, v in P["shapes_rebuilt"].items()},
                   "v3_by_hand": V3["against_own_pairs_csv"],
                   "part_b_second_method_max_diff": X["point_second_method_max_diff"],
                   "part_c_against_stage4_power_t": {"C1_max_rel_diff": S4["C1_max_rel_diff"],
                                                     "C2_all_n_equal": S4["C2_all_n_equal"],
                                                     "C2_max_power_abs_diff": S4["C2_max_power_abs_diff"]}},
        "readings": readings,
        "seeds": seeds,
        "headline": headline,
        "flags": flags,
        "outputs_sha256": outputs,
    }
    with open(os.path.join(C.HERE, "record.json"), "w") as f:
        json.dump(rec, f, indent=1)
    print(json.dumps(headline["verdicts"]))
    for fl in flags:
        print("-", fl)


if __name__ == "__main__":
    main()
