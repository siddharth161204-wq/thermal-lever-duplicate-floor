#!/usr/bin/env python3
"""Stage 5 independent check, comparison phase: D79 and D80 (DECISIONS.md at 67685ae), the README's Stage 5 section,
the note's Stage 5 paragraph and the Zenodo description at 67685ae, read against the registered results (101ca54), the
post-hoc numbers (157827f) and the checker's own numbers (frozen at ea3a5b4, with the post-freeze Monte Carlo study).

Every quoted number listed below is checked by code: its rounding against the pipeline value, and the checker's value
and 1,000,000-draw value beside it. Statements judged wrong, overstated or missing are recorded with a suggested wording.

Output: out_compare/text_review.json.
Usage: PYTHONDONTWRITEBYTECODE=1 python3 verification/stage5/independent/compare_text.py
"""
import json, math, os, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ind_common as C   # noqa: E402
import partc as P        # noqa: E402

PIPE = "101ca54a199def64ab7227a41277e2408b7f5e91"
DRAFT = "67685ae92dd66fd10359456d21fde8e66a87bc1c"
OUTC = os.path.join(C.HERE, "out_compare")


def gshow(rev, rel):
    return subprocess.run(["git", "-C", C.ROOT, "show", f"{rev}:{rel}"], capture_output=True, check=True).stdout


def rounds_to(value, quoted, decimals):
    return round(value, decimals) == round(quoted, decimals) or abs(round(value, decimals) - quoted) < 0.5 * 10 ** -decimals - 1e-15


def main():
    pb = json.loads(gshow(PIPE, "results/stage5_beta.json"))
    pp = json.loads(gshow(PIPE, "results/stage5_power.json"))
    v2 = json.loads(gshow(PIPE, "verification/stage5/v2.json"))
    v3 = json.loads(gshow(PIPE, "verification/stage5/v3.json"))
    texts = {"DECISIONS.md": gshow(DRAFT, "DECISIONS.md").decode(), "README.md": gshow(DRAFT, "README.md").decode(),
             "note/thermal_lever_note.md": gshow(DRAFT, "note/thermal_lever_note.md").decode(), ".zenodo.json": gshow(DRAFT, ".zenodo.json").decode()}
    X = json.load(open(os.path.join(C.OUT, "B_crosscheck.json")))["mc_spread"]
    MP = json.load(open(os.path.join(C.OUT, "pairs.json")))
    MCS = json.load(open(os.path.join(C.OUT, "C_mc_spread.json")))
    D80 = json.load(open(os.path.join(OUTC, "postfreeze_d80_checks.json")))
    FN = json.load(open(os.path.join(OUTC, "postfreeze_flag_numbers.json")))
    pc1 = {(r["case"], r["floor"]): r for r in pp["C1"]}
    pc2 = {(r["dataset"], r["band"], r["floor"]): r for r in pp["C2"]}
    sh = C.shapes()

    checks = []

    def num(where, quote, quoted, pipeline_value, decimals, checker_value=None, big=None):
        ok = rounds_to(pipeline_value, quoted, decimals)
        present = all(quote in t for t in texts.values() if False) or any(quote in t for t in texts.values())
        checks.append({"where": where, "quote": quote, "quoted": quoted, "pipeline": pipeline_value, "checker": checker_value,
                       "big_draw": big, "rounds_correctly_from_pipeline": ok, "quote_found_in_text": present})

    kw = pc2[("KIT", "warm", "KIT")]
    # the spreads and ratios
    num("D80 (5), README, note, Zenodo", "0.35 K", 0.35, kw["dT_rms_K"], 2, D80["dT_rms_K"]["KIT warm"])
    num("D80 (5), README, note", "0.44 K", 0.44, pc1[("KIT", "KIT")]["dT_rms_K"], 2)
    num("D80 (6), README", "1.54", 1.54, pc1[("Sandia", "KIT")]["dT_rms_K"], 2)
    num("D80 (6), README", "3.5 times", 3.5, pc1[("Sandia", "KIT")]["dT_rms_K"] / pc1[("KIT", "KIT")]["dT_rms_K"], 1)
    # beta and intervals
    for lab, path, q, dec in (("KIT beta", ("dataset", "KIT"), (-0.09, -2.05, 1.37), 2), ("Sandia beta", ("dataset", "Sandia"), (1.55, -0.74, 2.50), 2),
                              ("pooled beta", ("dataset", "pooled"), (1.10, -0.49, 1.96), 2),
                              ("measure (b) beta", ("brackets", "Sandia measure (b)"), (0.70, -2.03, 2.64), 2),
                              ("without NCA pair beta", ("brackets", "Sandia without the NCA 15 C 1C pair"), (1.81, 0.11, 2.90), 2)):
        d = pb[path[0]][path[1]]
        key = {"KIT beta": "KIT", "Sandia beta": "Sandia", "pooled beta": "pooled", "without NCA pair beta": "Sandia without NCA 15C 1C"}.get(lab)
        bigci = X[key]["beta"]["big_ci95"] if key else json.load(open(os.path.join(OUTC, "postfreeze_mc_extra.json")))["Sandia measure b beta"]["big_ci95"]
        num(f"{lab} (D80, README, note)", f"{q[0]:+.2f}".replace("+", ""), q[0], d["beta"], dec)
        num(f"{lab} interval lower end", f"{q[1]:.2f}", q[1], d["beta_boot"]["ci95"][0], dec, None, bigci[0])
        num(f"{lab} interval upper end", f"{q[2]:.2f}", q[2], d["beta_boot"]["ci95"][1], dec, None, bigci[1])
    # s_w
    for lab, path, q, dec, key in (("KIT warm s_w", ("band", "KIT", "warm"), (-5.19, -16.35, 8.18), 2, "KIT warm"),
                                   ("KIT cold s_w", ("band", "KIT", "cold"), (-0.79, -9.44, 4.06), 2, "KIT cold"),
                                   ("Sandia warm s_w", ("band", "Sandia", "warm"), (-4.35, -7.12, 0.46), 2, "Sandia warm"),
                                   ("Sandia cold s_w", ("band", "Sandia", "cold"), (7.10, -4.08, 18.48), 2, "Sandia cold"),
                                   ("pooled warm s_w", ("band", "pooled", "warm"), (-4.41, -7.14, -0.10), 2, "pooled warm"),
                                   ("pooled cold s_w", ("band", "pooled", "cold"), (4.42, -3.26, 8.10), 2, "pooled cold"),
                                   ("KIT 0 C s_w", ("kit_0C_unscored",), (-1.28, -3.50, 0.65), 2, "KIT 0 C")):
        d = pb
        for q_ in path:
            d = d[q_]
        bigci = [100 * v for v in X[key]["s_w"]["big_ci95"]]
        num(f"{lab} (README)", f"{q[0]}", q[0], d["s_w_pct_per_K"], dec)
        num(f"{lab} lower end", f"{q[1]}", q[1], d["s_w_boot_pct"]["ci95"][0], dec, None, bigci[0])
        num(f"{lab} upper end", f"{q[2]}", q[2], d["s_w_boot_pct"]["ci95"][1], dec, None, bigci[1])
    # sign tests and noise
    num("KIT sign p", "p = 0.14", 0.14, pb["dataset"]["KIT"]["sign"]["p_two_sided"], 2)
    num("KIT warm sign p", "p = 0.012", 0.012, pb["band"]["KIT"]["warm"]["sign"]["p_two_sided"], 3)
    num("KIT cold sign p", "p = 0.41", 0.41, pb["band"]["KIT"]["cold"]["sign"]["p_two_sided"], 2)
    num("Sandia sign p", "p = 0.35", 0.35, pb["dataset"]["Sandia"]["sign"]["p_two_sided"], 2)
    num("KIT share", "42 %", 42, 100 * pb["dataset"]["KIT"]["sign"]["share"], 0)
    for lab, path, q in (("KIT resid sd", ("dataset", "KIT"), (0.131, 0.092, 0.168)), ("Sandia resid sd", ("dataset", "Sandia"), (0.158, 0.091, 0.222)),
                         ("pooled resid sd", ("dataset", "pooled"), (0.142, 0.105, 0.175))):
        d = pb[path[0]][path[1]]
        key = path[1]
        num(lab, f"{q[0]}", q[0], d["resid_sd"], 3)
        num(lab + " lower end", f"{q[1]}", q[1], d["resid_sd_boot"]["ci95"][0], 3, None, X[key]["sd"]["big_ci95"][0])
        num(lab + " upper end", f"{q[2]}", q[2], d["resid_sd_boot"]["ci95"][1], 3, None, X[key]["sd"]["big_ci95"][1])
    num("resid sd over KIT sd(e)", "1.26 times", 1.26, pb["dataset"]["KIT"]["resid_sd"] / sh["KIT"]["sd_e"], 2)
    num("resid sd over Geslin sd(e)", "2.9 times", 2.9, pb["dataset"]["KIT"]["resid_sd"] / sh["Geslin"]["sd_e"], 1)
    # C1 and C2
    c = pc1[("KIT", "KIT")]
    num("KIT detectable", "6.56", 6.56, 100 * c["detectable_s_per_K"], 2)
    num("KIT detectable (note, Zenodo)", "6.6", 6.6, 100 * c["detectable_s_per_K"], 1)
    num("KIT simulated (Monte Carlo)", "6.47", 6.47, 100 * c["sim_lognormal_detectable_s_per_K"], 2,
        [r for r in MCS["C1"] if r["case"] == "KIT" and r["floor"] == "KIT"][0]["sim_mean_pct_per_K"])
    num("KIT shared cells", "8.03", 8.03, 100 * c["detectable_s_shared_cells_per_K"], 2)
    num("KIT against Geslin", "2.81", 2.81, 100 * pc1[("KIT", "Geslin")]["detectable_s_per_K"], 2)
    num("over warm median", "5.0 times", 5.0, c["detectable_s_per_K"] / abs(C.S_MED_TASK["warm"]), 1)
    num("over cold median", "1.8 times", 1.8, c["detectable_s_per_K"] / abs(C.S_MED_TASK["cold"]), 1)
    num("all 142 KIT pairs", "3.85", 3.85, 100 * pc1[("KIT all 142 pairs", "KIT")]["detectable_s_per_K"], 2)
    num("Sandia detectable, KIT floor", "2.58", 2.58, 100 * pc1[("Sandia", "KIT")]["detectable_s_per_K"], 2)
    num("warm power", "0.065", 0.065, kw["power_calc"], 3)
    kc = pc2[("KIT", "cold", "KIT")]
    num("cold power", "0.21", 0.21, kc["power_calc"], 2)
    num("warm pairs", "4,221", 4221, kw["pairs_needed_calc"], 0)
    num("warm pairs simulated (Monte Carlo)", "4,259", 4259, kw["pairs_needed_sim_lognormal_resampled"], 0,
        [r for r in MCS["C2"] if r["band"] == "warm"][0]["mean"])
    num("warm shared", "6,330", 6330, kw["pairs_needed_shared_cells"], 0)
    num("warm Geslin", "774", 774, pc2[("KIT", "warm", "Geslin")]["pairs_needed_calc"], 0)
    num("cold pairs", "205", 205, kc["pairs_needed_calc"], 0)
    num("cold simulated (Monte Carlo)", "207", 207, kc["pairs_needed_sim_lognormal_resampled"], 0,
        [r for r in MCS["C2"] if r["band"] == "cold"][0]["mean"])
    num("cold shared", "306", 306, kc["pairs_needed_shared_cells"], 0)
    num("cold Geslin", "39", 39, pc2[("KIT", "cold", "Geslin")]["pairs_needed_calc"], 0)
    num("Stage 4 warm chamber, KIT floor, 2 K", "128", 128, D80["stage4_warm_chamber_central_2K"]["KIT"]["n_pairs"], 0)
    num("Stage 4 warm chamber, Geslin floor, 2 K", "25", 25, D80["stage4_warm_chamber_central_2K"]["Geslin"]["n_pairs"], 0)
    # P7 widths
    wS = pb["dataset"]["Sandia"]["beta_boot"]["width"]
    wK = pb["dataset"]["KIT"]["beta_boot"]["width"]
    num("Sandia width", "3.24", 3.24, wS, 2, None, X["Sandia"]["beta"]["big_ci95"][1] - X["Sandia"]["beta"]["big_ci95"][0])
    num("KIT width", "3.42", 3.42, wK, 2, None, X["KIT"]["beta"]["big_ci95"][1] - X["KIT"]["beta"]["big_ci95"][0])
    num("width ratio", "0.95 times", 0.95, wS / wK, 2)
    num("P3 interval width", "about 25 % per K wide", 25, pb["band"]["KIT"]["warm"]["s_w_boot_pct"]["ci95"][1] - pb["band"]["KIT"]["warm"]["s_w_boot_pct"]["ci95"][0], 0)
    # ties and one pair per set
    num("median interval over lifetime", "50 %", 50, 100 * D80["pipeline_posthoc_157827f"]["ties"]["KIT"]["median_resolution_over_lifetime"], 0,
        100 * D80["ties"]["KIT"]["median_res_over_L_distinct_cells"])
    num("median difference over interval", "9 %", 9, 100 * D80["ties"]["KIT"]["median_dL_over_res"], 0)
    num("Sandia median difference over interval", "10 %", 10, 100 * D80["ties"]["Sandia"]["median_dL_over_res"], 0)
    for b, q in (("all", (15, 36, 0.41)), ("warm", (8, 24, 0.15)), ("cold", (7, 12, 0.77))):
        r = D80["one_pair_per_set"]["KIT"][b]
        checks.append({"where": f"D80 (2) KIT {b}, one pair per set", "quoted": q, "checker": [r["agree"], r["pairs"], r["p_two_sided"]],
                       "rounds_correctly_from_pipeline": (r["agree"], r["pairs"]) == q[:2] and rounds_to(r["p_two_sided"], q[2], 2)})
    r = D80["one_pair_per_set"]["Sandia"]["all"]
    checks.append({"where": "D80 (2) Sandia, one pair per set", "quoted": (10, 21, 1.0), "checker": [r["agree"], r["pairs"], r["p_two_sided"]],
                   "rounds_correctly_from_pipeline": (r["agree"], r["pairs"]) == (10, 21) and rounds_to(r["p_two_sided"], 1.0, 1)})
    # V checks quoted in the README and D79
    reg = [r for r in pp["C1"] if r["registered"]]
    rels = [-r["v4_rel_diff"] * 100 for r in reg]
    num("V4 registered range low", "1.2", 1.2, min(rels), 1)
    num("V4 registered range high", "3.7", 3.7, max(rels), 1)
    num("V4 Sandia cold KIT", "5.2", 5.2, -100 * pc1[("Sandia cold", "KIT")]["v4_rel_diff"], 1,
        -100 * [r for r in MCS["C1"] if r["case"] == "Sandia cold" and r["floor"] == "KIT"][0]["rel_diff_mean"])
    num("V4 Sandia cold Geslin", "6.9", 6.9, -100 * pc1[("Sandia cold", "Geslin")]["v4_rel_diff"], 1,
        -100 * [r for r in MCS["C1"] if r["case"] == "Sandia cold" and r["floor"] == "Geslin"][0]["rel_diff_mean"])
    num("normal check", "1.2 %", 1.2, 100 * max(abs(r["normal_check_rel_diff"]) for r in pp["C1"]), 1)
    num("V2 point", "2e-16", 2e-16, max(r["abs_diff"] for r in v2["rows"]), 16)
    num("V2 ends", "2.2 Monte Carlo standard deviations", 2.2, max(max(r["ends_in_mc_sd"]) for r in v2["rows"]), 1)
    num("V3", "6e-15", 6e-15, max(max(c_["rel_diff"].values()) for c_ in v3["cases"]), 15)
    # V1's claim of identity with the first release, from the checker's own hashes
    v1_first = all(h["first_archived"]["sha256"] == h["sha256_at_v1.3.0"] for h in MP["input_hashes"])
    # the figure: where the KIT curve of panel (b) meets the warm median
    n_fig = P.n_required(abs(C.S_MED_TASK["warm"]) * pc1[("KIT", "KIT")]["dT_rms_K"], sh["KIT"]["sd_e"])
    # the F of the Sandia selections
    F = {r["case"]: r["shared_cell_F"] for r in pp["C1"] if r["floor"] == "KIT"}

    pw = X["pooled warm"]["s_w"]
    findings = [
        {"id": 1, "kind": "mis-rounded", "where": "D80 (5); README 'What the data could detect' and 'What it means for the thesis'; note; Zenodo description",
         "quote": "at its 0.35 K spread / at the warm band's 0.35 K",
         "fact": f"the KIT warm band's dT root mean square is {kw['dT_rms_K']:.4f} K (mean {D80['dT_mean_K']['KIT warm']:.3f} K)",
         "suggest": "0.34 K"},
        {"id": 2, "kind": "Monte Carlo, wording", "where": "README 'The within-set lever' (pooled warm), and any later use",
         "quote": "pooled warm -4.41 (-7.14 to -0.10)",
         "fact": (f"the registered draw's upper end, {100 * pb['band']['pooled']['warm']['s_w_boot']['ci95'][1]:.2f} % per K, lies "
                  f"{abs(pb['band']['pooled']['warm']['s_w_boot']['ci95'][1] - pw['hi_mean']) / pw['hi_sd']:.1f} Monte Carlo sd below the mean "
                  f"of 200 independent 10,000-draw runs ({100 * pw['hi_mean']:+.2f}, sd {100 * pw['hi_sd']:.2f}); {pw['contains_share']['0.0']:.1%} "
                  f"of those runs include zero; 1,000,000 draws give {100 * pw['big_ci95'][0]:.2f} to {100 * pw['big_ci95'][1]:+.2f} % per K"),
         "suggest": ("pooled warm -4.41 (-7.14 to -0.10 in the registered draw; with 1,000,000 draws -7.11 to +0.03, so the interval's "
                     "upper end lies at zero within its Monte Carlo error and the pooled warm lever is not resolved from zero)")},
        {"id": 3, "kind": "overstated", "where": "README 'The noise' and 'What it means for the thesis'; note; Zenodo; D80 (3)",
         "quote": "which is the floor / their spread is the duplicate floor / the residual spread is the KIT floor's / the within-set noise is the floor of D31",
         "fact": (f"by D74's rule the KIT noise is recorded as the KIT floor because sd(e) = 0.104 lies inside 0.092 to 0.168; the point "
                  f"value is 1.26 times sd(e), the KIT floor is computed from cells 1 and 2 of these same sets (D74), and the pooled "
                  f"interval (0.105 to 0.175) excludes it"),
         "suggest": ("'the residual spread is consistent with the KIT floor (0.131 against its sd(e) of 0.104, from cells of these same "
                     "sets); the pooled interval lies just above it'")},
        {"id": 4, "kind": "overstated", "where": "note (Sandia sentence)",
         "quote": "and show no thermal lever",
         "fact": "beta 1.55 (-0.74 to +2.50) neither shows nor excludes a lever; the brief's pitfall is that a beta at or above one would not show one",
         "suggest": "'and cannot show a thermal lever' (as the README and D80 (4) say: 'does not show')"},
        {"id": 5, "kind": "missing", "where": "D80 (4) and (6); README 'Sandia, and what beta does not show' and the P7 sentence; note",
         "quote": "Sandia beta 1.55; P7 '... because its pairs span 3.5 times KIT's dT'",
         "fact": (f"the 23 Sandia pairs at 2C or 3C discharge (13 of 21 groups), which D72 gives the 1C lever, hold "
                  f"{FN['flag9']['share of Sandia sum p^2 at 2C or 3C']:.0%} of Sandia's sum of p^2, and one group (NMC 25 C 0-100 0.5/3C) "
                  f"{FN['flag7']['Sandia NMC 25 C 0-100 0.5/3C']['p2_share_of_Sandia']:.0%}; over the 1C groups alone beta is "
                  f"{FN['flag9']['sensitivity (post-freeze, not registered): Sandia beta over the 1C groups only']['beta']:.2f}, without that group "
                  f"{FN['flag7']['Sandia NMC 25 C 0-100 0.5/3C']['Sandia beta without it']:.2f}; the interval width is set by sum p^2 (Sandia "
                  f"{D80['sum_p2']['Sandia']:.4f} against KIT {D80['sum_p2']['KIT']:.4f}), not by dT alone, and the two dT are different "
                  f"measures (swing against mean cell temperature)"),
         "suggest": "record in D81 (flags 7 and 9) and add one clause where Sandia's beta and P7 are explained"},
        {"id": 6, "kind": "missing", "where": "D80 (6) P3; README P3 sentence",
         "quote": "so its sign carries little",
         "fact": (f"in the same band the registered per-band sign test is 24 of 70 (p = 0.012) and the per-band beta (not registered) is "
                  f"{FN['flag5']['KIT warm beta (checker extra)']:.2f} ({FN['flag5']['KIT warm beta ci95 (1,000,000)'][0]:.1f} to "
                  f"{FN['flag5']['KIT warm beta ci95 (1,000,000)'][1]:.2f} with 1,000,000 draws; with protocol clusters "
                  f"{FN['flag5']['KIT warm beta ci95, 12 protocol clusters (1,000,000)'][0]:.1f} to "
                  f"{FN['flag5']['KIT warm beta ci95, 12 protocol clusters (1,000,000)'][1]:.2f}; leaving out any one set "
                  f"{FN['flag5']['KIT warm beta, leaving out one set: range'][0]:.1f} to {FN['flag5']['KIT warm beta, leaving out one set: range'][1]:.1f}), "
                  f"excluding both 0 and 1, while the one-pair-per-set sign test gives 8 of 24 (p = 0.15)"),
         "suggest": "record in D81 (flag 5); the P3 sentence can stay, with 'the warm band's within-set differences run against its levers more often than not (D81)'"},
        {"id": 7, "kind": "seed-dependent", "where": "D79 (3); README checks paragraph",
         "quote": "fails in two unscored ones, the 20 Sandia cold-band pairs under either floor (5.2 and 6.9 % below)",
         "fact": ("under the KIT floor the case sits at the 5 % line: the checker's run gives 4.7 % below (passes) and 15 of its 20 replicate "
                  "simulations pass (mean 4.7 %, sd 0.6 %); under the Geslin floor it fails in 19 of 20 (mean 5.9 %)"),
         "suggest": "'fails in this run in the 20 Sandia cold-band pairs under either floor (5.2 and 6.9 % below; under the KIT floor the verdict depends on the seed)'"},
        {"id": 8, "kind": "imprecise", "where": "README Limits",
         "quote": "1.1 to 1.8 for the Sandia groups",
         "fact": (f"F is 1 for a two-cell group and 2 for a four-cell group; 1.13, 1.72 and 1.85 are the F of the Sandia cold band, all 56 "
                  f"pairs and the warm band ({F['Sandia cold']:.4f}, {F['Sandia']:.4f}, {F['Sandia warm']:.4f})"),
         "suggest": "'1.5 for KIT (three-cell sets) and 1.72 over the Sandia pairs (1.13 cold band, 1.85 warm; 1 for a two-cell group, 2 for a four-cell one)'"},
        {"id": 9, "kind": "missing", "where": "note; Zenodo description",
         "quote": "the smallest |s| they could detect is 6.6 % per K ... the warm median needs 4,221 pairs, against Stage 4's 128",
         "fact": "all three numbers are against the KIT floor; against the Geslin floor they are 2.81 % per K, 774 and 25",
         "suggest": "add 'against the KIT floor' once"},
        {"id": 10, "kind": "clarity", "where": "note Figure 4 caption (b)",
         "quote": "against the number of pairs at the observed temperature spread of the KIT and Sandia pairs",
         "fact": (f"panel (b)'s KIT curve is at the 106 pairs' 0.44 K spread, where the warm median is reached at {n_fig} pairs "
                  f"(KIT floor, calculation); the text's 4,221 is at the warm band's own 0.34 K"),
         "suggest": "add 'all 106 KIT pairs (0.44 K rms)' to the caption, or say in the text which spread 4,221 belongs to"},
        {"id": 11, "kind": "precision", "where": "README bullets (intervals at two decimals)",
         "quote": "-16.35 to +8.18; -9.44 to +4.06; -7.12 to +0.46; -4.08 to +18.48 ...",
         "fact": (f"the second decimal of a 10,000-draw interval end is Monte Carlo noise: the end sd is about "
                  f"{100 * X['KIT warm']['s_w']['lo_sd']:.2f} and {100 * X['KIT warm']['s_w']['hi_sd']:.2f} % per K for the KIT warm band "
                  f"(1,000,000 draws: {100 * X['KIT warm']['s_w']['big_ci95'][0]:.2f} to {100 * X['KIT warm']['s_w']['big_ci95'][1]:+.2f}) and 0.04 "
                  f"for the KIT beta ends (1,000,000 draws: {X['KIT']['beta']['big_ci95'][0]:.2f} to {X['KIT']['beta']['big_ci95'][1]:+.2f})"),
         "suggest": "quote the s_w interval ends to 0.1 % per K, as D80 and the note do, or say the registered draw's ends carry Monte Carlo error of a few hundredths"},
        {"id": 12, "kind": "pending", "where": "README file list", "quote": "D79 to D81 after its results",
         "fact": "D81 is not yet written at 67685ae", "suggest": "keep once D81 is committed"},
        {"id": 13, "kind": "minor omission", "where": "D79 (1)", "quote": "analyse was made to return the counts of an empty subset",
         "fact": "the change in 101ca54 also drops the residual for a subset of one pair; no registered subset has one pair, so nothing moves",
         "suggest": "'... to return the counts of an empty subset and to skip the residual of a one-pair subset'"},
        {"id": 14, "kind": "verified", "where": "README checks paragraph (V1)",
         "quote": "identical (SHA-256) to its version at v1.3.0 and at the release that first archived it",
         "fact": f"the checker's hashes agree: every D71 input equals its first-archived version: {v1_first}", "suggest": "none"},
        {"id": 15, "kind": "style", "where": "README How to run (Stage 5 lines)", "quote": "python3 tests/stage5_synthetic.py ...",
         "fact": "the Stage 4 lines run with python3 -I, the Stage 5 lines without it", "suggest": "use -I as for Stage 4, if the modules allow it"},
    ]
    bad = [c_ for c_ in checks if not c_["rounds_correctly_from_pipeline"]]
    out = {"texts_read_at": DRAFT, "results_read_at": PIPE, "numbers_checked": len(checks), "numbers_not_matching": bad,
           "checks": checks, "findings": findings, "figure_panel_b_warm_median_pairs_at_KIT_spread": n_fig}
    os.makedirs(OUTC, exist_ok=True)
    with open(os.path.join(OUTC, "text_review.json"), "w") as f:
        json.dump(out, f, indent=1)
    print("numbers checked:", len(checks), "not matching:", [(b_["where"], b_["quoted"], b_["pipeline"]) for b_ in bad])
    print("figure panel (b), KIT spread, warm median:", n_fig)
    for fd in findings:
        print(fd["id"], fd["kind"], "|", fd["where"], "|", fd["fact"][:300])


if __name__ == "__main__":
    main()
