#!/usr/bin/env python3
"""Stage 5 independent check, comparison phase (written after the freeze at ea3a5b4): the checker's frozen outputs
against the registered pipeline results of 101ca54 (pushed 2026-10-29T04:22:03Z), read with git show.

Deterministic classes (largest absolute and relative difference, and where): the per-pair table; B1 to B4 point values;
counts and sign tests; C1 calculated values (20 cases, with the normal shape and the shared cells); C2 calculated power
and pairs needed; the verdicts. Monte Carlo quantities (bootstrap interval ends, simulated detectable |s|, simulated
power and pairs) are set against the checker's estimate of their Monte Carlo spread (replicate runs with independent
seeds: out/B_crosscheck.json, out/C_mc_spread.json, out_compare/postfreeze_mc_extra.json), as z = (value - replicate
mean) / replicate sd for each run and as the difference over sqrt(2) sd, the spread of a difference of two runs.

Output: out_compare/compare.json, out_compare/compare_table.csv.
Usage: PYTHONDONTWRITEBYTECODE=1 python3 verification/stage5/independent/compare.py
"""
import csv, io, json, math, os, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ind_common as C   # noqa: E402

PIPE = "101ca54a199def64ab7227a41277e2408b7f5e91"
FREEZE = "ea3a5b46a891e74d504c818d30b35d34b4351ff5"
OUTC = os.path.join(C.HERE, "out_compare")


def gshow(rev, rel):
    return subprocess.run(["git", "-C", C.ROOT, "show", f"{rev}:{rel}"], capture_output=True, check=True).stdout


def fnum(x):
    if x in (None, ""):
        return None
    return float(x)


def fbool(x):
    if isinstance(x, bool):
        return x
    return str(x) == "True"


class Num:
    """Largest absolute and relative difference of matched numbers, with where they occur."""

    def __init__(self, name):
        self.name = name
        self.n = 0
        self.max_abs = 0.0
        self.max_rel = 0.0
        self.where_abs = ""
        self.where_rel = ""
        self.one_sided = []

    def add(self, where, a, b):
        if a is None and b is None:
            return
        if (a is None) != (b is None):
            self.one_sided.append({"where": where, "checker": a, "pipeline": b})
            return
        self.n += 1
        d = abs(a - b)
        m = max(abs(a), abs(b))
        r = d / m if m > 0 else 0.0
        if d > self.max_abs:
            self.max_abs, self.where_abs = d, where
        if r > self.max_rel:
            self.max_rel, self.where_rel = r, where

    def out(self):
        return {"compared": self.n, "max_abs_diff": self.max_abs, "where_max_abs": self.where_abs,
                "max_rel_diff": self.max_rel, "where_max_rel": self.where_rel, "present_on_one_side_only": self.one_sided}


class Cat:
    def __init__(self, name):
        self.name = name
        self.n = 0
        self.mismatches = []

    def add(self, where, a, b):
        self.n += 1
        if a != b:
            self.mismatches.append({"where": where, "checker": a, "pipeline": b})

    def out(self):
        return {"compared": self.n, "mismatches": len(self.mismatches), "list": self.mismatches[:40]}


def mc_entry(name, checker, pipeline, mean, sd, big=None):
    z_c = (checker - mean) / sd if sd else None
    z_p = (pipeline - mean) / sd if sd else None
    d = pipeline - checker
    return {"quantity": name, "checker": checker, "pipeline": pipeline, "replicate_mean": mean, "replicate_sd": sd,
            "z_checker": z_c, "z_pipeline": z_p, "diff": d, "diff_over_sd": d / sd if sd else None,
            "diff_over_sqrt2_sd": d / (math.sqrt(2) * sd) if sd else None, "big_draw_value": big}


def main():
    mine_pairs = list(csv.DictReader(open(os.path.join(C.OUT, "pairs.csv"))))
    pipe_pairs = list(csv.DictReader(io.StringIO(gshow(PIPE, "results/stage5_pairs.csv").decode())))
    pb = json.loads(gshow(PIPE, "results/stage5_beta.json"))
    pp = json.loads(gshow(PIPE, "results/stage5_power.json"))
    ps = json.loads(gshow(PIPE, "results/stage5_scorecard.json"))
    MB = json.load(open(os.path.join(C.OUT, "B.json")))["analyses"]
    MCj = json.load(open(os.path.join(C.OUT, "C.json")))
    MS = json.load(open(os.path.join(C.OUT, "scorecard.json")))
    MP = json.load(open(os.path.join(C.OUT, "pairs.json")))
    X = json.load(open(os.path.join(C.OUT, "B_crosscheck.json")))["mc_spread"]
    MCS = json.load(open(os.path.join(C.OUT, "C_mc_spread.json")))
    EX = json.load(open(os.path.join(OUTC, "postfreeze_mc_extra.json")))
    res = {"pipeline_commit": PIPE, "checker_freeze_commit": FREEZE,
           "pipeline_files_sha256": {rel: __import__("hashlib").sha256(gshow(PIPE, rel)).hexdigest() for rel in
                                     ("results/stage5_pairs.csv", "results/stage5_beta.json", "results/stage5_power.json",
                                      "results/stage5_scorecard.json")}}

    # ---------------- per-pair table
    def key_mine(x):
        cells = {x["warm"], x["cool"]} if x["dataset"] == "KIT" else {x["cell_a"], x["cell_b"]}
        return (x["dataset"], x["set"], frozenset(cells))

    def key_pipe(x):
        return (x["dataset"], x["set"], frozenset({x["warm"], x["cool"]}))

    km = {key_mine(x): x for x in mine_pairs}
    kp = {key_pipe(x): x for x in pipe_pairs}
    if len(km) != len(mine_pairs) or len(kp) != len(pipe_pairs):
        C.stop("duplicate pair keys")
    pr = {"rows_checker": len(mine_pairs), "rows_pipeline": len(pipe_pairs),
          "only_checker": [str(sorted(k[2])) for k in km if k not in kp], "only_pipeline": [str(sorted(k[2])) for k in kp if k not in km]}
    cats = {f: Cat(f) for f in ("scored", "band", "tie", "orientation (warm, cool)", "orientation by measure (b)")}
    nums = {f: Num(f) for f in ("dT_K", "d", "p", "s_per_K", "s_lo_per_K", "s_hi_per_K", "L_warm", "L_cool", "resolution_efc",
                                "dT_b_K", "d_b", "p_b")}
    band_label_notes = []
    for k, m in km.items():
        p = kp[k]
        where = f"{k[0]} {k[1]} {'/'.join(sorted(k[2]))}"
        sm, sp = fbool(m["scored"]), fbool(p["scored"])
        cats["scored"].add(where, sm, sp)
        cats["band"].add(where, m["band"], p["band"])
        if m["band"] != p["band"]:
            band_label_notes.append(where)
        cats["tie"].add(where, fbool(m["tie"]), fbool(p["tie"]))
        oriented = m["dT_K"] not in ("", None)
        if oriented:
            cats["orientation (warm, cool)"].add(where, (m["warm"], m["cool"]), (p["warm"], p["cool"]))
            nums["dT_K"].add(where, fnum(m["dT_K"]), fnum(p["dT_K"]))
            nums["d"].add(where, fnum(m["d"]), fnum(p["d"]))
            nums["L_warm"].add(where, fnum(m["L_warm"]), fnum(p["L_warm_efc"]))
            nums["L_cool"].add(where, fnum(m["L_cool"]), fnum(p["L_cool_efc"]))
        nums["resolution_efc"].add(where, fnum(m["resolution_efc"]), fnum(p["resolution_efc"]))
        if sm and sp:
            nums["p"].add(where, fnum(m["p"]), fnum(p["p"]))
            nums["s_per_K"].add(where, fnum(m["s_per_K"]), fnum(p["s_per_K"]))
            nums["s_lo_per_K"].add(where, fnum(m["s_ci95_lo"]), fnum(p["s_lo_per_K"]))
            nums["s_hi_per_K"].add(where, fnum(m["s_ci95_hi"]), fnum(p["s_hi_per_K"]))
            if k[0] == "Sandia":
                cats["orientation by measure (b)"].add(where, (m["warm_b"], m["cool_b"]), (p["warm_b"], p["cool_b"]))
                nums["dT_b_K"].add(where, fnum(m["dT_b_K"]), fnum(p["dT_b_K"]))
                nums["d_b"].add(where, fnum(m["d_b"]), fnum(p["d_b"]))
                nums["p_b"].add(where, fnum(m["p_b"]), fnum(p["p_b"]))
        else:
            # the lever beside an unscored pair is a label, not an input of any number; reported where both give one
            if fnum(m["s_per_K"]) is not None and fnum(p["s_per_K"]) is not None:
                nums["s_per_K"].add(where + " (unscored)", fnum(m["s_per_K"]), fnum(p["s_per_K"]))
    pr["categorical"] = {k: v.out() for k, v in cats.items()}
    pr["numeric"] = {k: v.out() for k, v in nums.items()}
    res["pairs"] = pr

    # ---------------- B1 to B4 point values
    def g(d, *path):
        for q in path:
            if d is None or q not in d:
                return None
            d = d[q]
        return d

    beta_map = {"KIT": ("dataset", "KIT"), "Sandia": ("dataset", "Sandia"), "pooled": ("dataset", "pooled"),
                "untied KIT": ("brackets", "KIT untied"), "Sandia measure b": ("brackets", "Sandia measure (b)"),
                "Sandia without NCA 15C 1C": ("brackets", "Sandia without the NCA 15 C 1C pair")}
    sw_map = {"KIT warm": ("band", "KIT", "warm"), "KIT cold": ("band", "KIT", "cold"), "Sandia warm": ("band", "Sandia", "warm"),
              "Sandia cold": ("band", "Sandia", "cold"), "pooled warm": ("band", "pooled", "warm"),
              "pooled cold": ("band", "pooled", "cold"), "KIT 0 C": ("kit_0C_unscored",),
              "untied KIT warm": ("brackets", "KIT untied, by band", "warm")}
    sd_map = {k: beta_map[k] for k in ("KIT", "Sandia", "pooled", "untied KIT")}
    b_num = {q: Num(q) for q in ("beta", "s_w_per_K", "resid_sd")}
    for k, path in beta_map.items():
        b_num["beta"].add(k, MB[k]["beta"], g(pb, *path, "beta"))
    for k, path in sw_map.items():
        b_num["s_w_per_K"].add(k, MB[k]["s_w_per_K"], g(pb, *path, "s_w_per_K"))
    for k, path in sd_map.items():
        b_num["resid_sd"].add(k, MB[k]["resid_sd"], g(pb, *path, "resid_sd"))
    res["B_point"] = {k: v.out() for k, v in b_num.items()}
    res["B_point_note"] = ("the pipeline gives no per-band beta or residual sd and no s_w over a whole dataset; those are "
                           "extras of the checker and are not compared. Its pooled untied bracket points to the KIT untied "
                           "numbers ('same_pairs_as'), which the checker's pooled untied bracket equals exactly.")

    # ---------------- counts and sign tests
    cnt = Cat("counts")
    pc = pb["counts"]
    mk = MP["counts"]
    for label, mine_v in (("KIT pairs in file", mk["KIT"]["pairs_in_file"]), ("KIT scored", mk["KIT"]["scored_pairs"]),
                          ("KIT scored sets", mk["KIT"]["scored_sets"]),
                          ("KIT 0 C pairs, not scored", sum(v for r, v in mk["KIT"]["not_scored_by_reason"].items())),
                          ("KIT warm", mk["KIT"]["by_band"]["warm"]["pairs"]), ("KIT cold", mk["KIT"]["by_band"]["cold"]["pairs"]),
                          ("KIT ties among scored", mk["KIT"]["scored_ties"]), ("KIT ties among all 142", mk["KIT"]["ties_all_142"]),
                          ("Sandia pairs in file", mk["Sandia"]["pairs_in_file"]),
                          ("Sandia with measure (a)", mk["Sandia"]["pairs_with_measure_a_both"]),
                          ("Sandia scored", mk["Sandia"]["scored_pairs"]), ("Sandia scored groups", mk["Sandia"]["scored_sets"]),
                          ("Sandia warm", mk["Sandia"]["by_band"]["warm"]["pairs"]), ("Sandia cold", mk["Sandia"]["by_band"]["cold"]["pairs"]),
                          ("Sandia ties among scored", mk["Sandia"]["scored_ties"]),
                          ("Sandia ties among the 60 with measure (a)", mk["Sandia"]["ties_among_measure_a_pairs"]),
                          ("Sandia with measure (a) but no lever", mk["Sandia"]["not_scored_by_reason"].get("LFP 25 to 35 C lever not determinable (D52)")),
                          ("pooled scored", mk["pooled"]["scored_pairs"])):
        cnt.add(label, mine_v, pc.get(label))
    for k, path in list(beta_map.items()) + list(sw_map.items()):
        cnt.add(k + " pairs", MB[k]["n_pairs"], g(pb, *path, "pairs"))
        cnt.add(k + " sets", MB[k]["n_sets"], g(pb, *path, "sets"))
    sign_map = dict(beta_map)
    sign_map.update({k: v for k, v in sw_map.items() if k != "KIT 0 C"})
    sg_cat = Cat("sign counts")
    sg_num = Num("sign p two-sided")
    for k, path in sign_map.items():
        s_m = MB[k].get("sign")
        s_p = g(pb, *path, "sign")
        if s_m is None or s_p is None:
            continue
        sg_cat.add(k, (s_m["agree"], s_m["pairs_tested"], s_m["left_out"]), (s_p["agree"], s_p["pairs"], s_p["left_out_zero"]))
        sg_num.add(k, s_m["p_two_sided"], s_p["p_two_sided"])
        sg_num.add(k + " share", s_m["share"], s_p["share"])
    res["counts"] = cnt.out()
    res["sign"] = {"counts": sg_cat.out(), "p_and_share": sg_num.out()}

    # ---------------- C1 calculated
    pc1 = {(r["case"].replace("KIT all 142 pairs", "KIT all 142"), r["floor"]): r for r in pp["C1"]}
    c1n = {q: Num(q) for q in ("detectable_s_per_K", "normal_shape_per_K", "F", "shared_cells_per_K", "dT_rms_K", "sum_dT2_K2",
                                "dT_median_K", "sd_e", "power_at_detectable")}
    c1c = Cat("C1 counts")
    for r in MCj["C1"]:
        k = (r["case"], r["floor"])
        q = pc1[k]
        w = f"{k[0]} | {k[1]}"
        c1c.add(w, (r["n_pairs"], r["n_sets"]), (q["pairs"], q["sets"]))
        c1n["detectable_s_per_K"].add(w, r["detectable_s_calc_pct_per_K"] / 100, q["detectable_s_per_K"])
        c1n["normal_shape_per_K"].add(w, r["detectable_s_calc_normal_shape_pct_per_K"] / 100, q["detectable_s_normal_shape_per_K"])
        c1n["F"].add(w, r["F_shared_cell"], q["shared_cell_F"])
        c1n["shared_cells_per_K"].add(w, r["detectable_s_calc_shared_cell_pct_per_K"] / 100, q["detectable_s_shared_cells_per_K"])
        c1n["dT_rms_K"].add(w, r["dT_rms_K"], q["dT_rms_K"])
        c1n["sum_dT2_K2"].add(w, r["sum_dT2_K2"], q["sum_dT2_K2"])
        c1n["dT_median_K"].add(w, r["dT_median_K"], q["dT_median_K"])
        c1n["sd_e"].add(w, r["sd_e"], q["sd_e"])
        c1n["power_at_detectable"].add(w, r["power_quadrature_at_calc"], q["power_check_at_detectable"])
    res["C1_calculated"] = {"counts": c1c.out(), **{k: v.out() for k, v in c1n.items()}}

    # ---------------- C2 calculated
    pc2 = {(r["dataset"], r["band"], r["floor"]): r for r in pp["C2"]}
    c2n = {q: Num(q) for q in ("power_calc", "dT_rms_K", "sd_e", "F", "abs_s_per_K")}
    c2c = Cat("C2 pairs needed (calculated, shared-cell, normal shape) and pairs available")
    for r in MCj["C2"]:
        k = (r["dataset"], r["band"], r["floor"])
        q = pc2[k]
        w = " | ".join(k)
        c2n["power_calc"].add(w, r["power_calc"], q["power_calc"])
        c2n["dT_rms_K"].add(w, r["dT_rms_K"], q["dT_rms_K"])
        c2n["sd_e"].add(w, r["sd_e"], q["sd_e"])
        c2n["F"].add(w, r["F_shared_cell"], q["shared_cell_F"])
        c2n["abs_s_per_K"].add(w, r["band_median_abs_s_per_K"], q["abs_s_per_K"])
        c2c.add(w, (r["n_required_calc"], r["n_required_shared_cell"], r["n_required_calc_normal_shape"], r["n_pairs_available"]),
                (q["pairs_needed_calc"], q["pairs_needed_shared_cells"], q["pairs_needed_normal_shape"], q["pairs"]))
    res["C2_calculated"] = {"pairs_needed": c2c.out(), **{k: v.out() for k, v in c2n.items()}}

    # ---------------- verdicts
    vc = Cat("verdicts")
    for k in ("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"):
        vc.add(k, MS["verdicts"][k], ps["scorecard"][k]["held"])
    vc.add("P6 scored", MS["details"]["P6"]["scored"], ps["scorecard"]["P6"]["scored"])
    fl = Cat("B4 floor verdicts (rest on bootstrap intervals)")
    for k in ("KIT", "Sandia", "pooled", "untied KIT"):
        for f in ("KIT", "Geslin"):
            fl.add(f"{k} against the {f} floor", MB[k]["floor_verdicts"][f]["inside_interval"],
                   g(pb, *beta_map[k], "floor_inside_interval", f, "inside"))
    v4 = Cat("V4 pass per C1 case (Monte Carlo)")
    mcs1 = {(r["case"], r["floor"]): r for r in MCS["C1"]}
    v4_notes = []
    for r in MCj["C1"]:
        k = (r["case"], r["floor"])
        q = pc1[k]
        v4.add(f"{k[0]} | {k[1]}", r["v4_pass"], q["v4_pass"])
        if r["v4_pass"] != q["v4_pass"]:
            v4_notes.append({"case": f"{k[0]} | {k[1]}", "checker_rel": r["v4_rel_diff"], "pipeline_rel": q["v4_rel_diff"],
                             "replicate_pass_share": mcs1[k]["v4_pass_share"], "replicate_rel_mean": mcs1[k]["rel_diff_mean"],
                             "replicate_rel_sd": mcs1[k]["rel_diff_sd"]})
    res["verdicts"] = {"P1_to_P8": vc.out(), "floor_verdicts": fl.out(), "v4_pass": v4.out(), "v4_differences": v4_notes}

    # ---------------- Monte Carlo quantities, in units of the checker's replicate spread
    mc = []
    ends = (("lo", 0), ("hi", 1))
    for k, path in beta_map.items():
        ci_p = g(pb, *path, "beta_boot", "ci95")
        ci_m = MB[k]["bootstrap"]["beta_ci95"]
        if k == "Sandia measure b":
            spread = EX["Sandia measure b beta"]
        elif k in X:
            spread = X[k]["beta"]
        else:
            spread = None
        for e, i in ends:
            if spread:
                mc.append(mc_entry(f"beta {e} {k}", ci_m[i], ci_p[i], spread[e + "_mean"], spread[e + "_sd"], spread["big_ci95"][i]))
            else:
                mc.append({"quantity": f"beta {e} {k}", "checker": ci_m[i], "pipeline": ci_p[i], "diff": ci_p[i] - ci_m[i],
                           "note": "two pairs in two sets: the ends are the two single-pair values in almost every run"})
    for k, path in sw_map.items():
        ci_p = g(pb, *path, "s_w_boot", "ci95")
        ci_m = MB[k]["bootstrap"]["s_w_ci95"]
        spread = X.get(k, {}).get("s_w")
        for e, i in ends:
            if spread:
                mc.append(mc_entry(f"s_w {e} {k} (per K)", ci_m[i], ci_p[i], spread[e + "_mean"], spread[e + "_sd"], spread["big_ci95"][i]))
            else:
                mc.append({"quantity": f"s_w {e} {k} (per K)", "checker": ci_m[i], "pipeline": ci_p[i], "diff": ci_p[i] - ci_m[i],
                           "note": "two pairs in two sets: degenerate"})
    for k, path in sd_map.items():
        ci_p = g(pb, *path, "resid_sd_boot", "ci95")
        ci_m = MB[k]["bootstrap"]["resid_sd_ci95"]
        spread = X.get(k, {}).get("sd")
        for e, i in ends:
            if spread:
                mc.append(mc_entry(f"resid_sd {e} {k}", ci_m[i], ci_p[i], spread[e + "_mean"], spread[e + "_sd"], spread["big_ci95"][i]))
            else:
                mc.append({"quantity": f"resid_sd {e} {k}", "checker": ci_m[i], "pipeline": ci_p[i], "diff": ci_p[i] - ci_m[i],
                           "note": "two pairs in two sets: degenerate"})
    # P7's width ratio
    wr_m = MS["details"]["P7"]["ratio"]
    wr_p = ps["scorecard"]["P7"]["ratio"]
    mc.append({"quantity": "P7 width ratio Sandia/KIT", "checker": wr_m, "pipeline": wr_p, "diff": wr_p - wr_m,
               "replicate_means_ratio": X["Sandia"]["beta"]["width_mean"] / X["KIT"]["beta"]["width_mean"],
               "big_draw_ratio": MS["details"]["P7"]["big_draw_ratio"]})
    ex_n = {(r["case"], r["floor"]): r for r in EX["C1 normal same sd"]}
    for r in MCj["C1"]:
        k = (r["case"], r["floor"])
        q = pc1[k]
        s = mcs1[k]
        mc.append(mc_entry(f"C1 simulated log-normal {k[0]} | {k[1]} (% per K)", r["detectable_s_sim_lognormal_pct_per_K"],
                           100 * q["sim_lognormal_detectable_s_per_K"], s["sim_mean_pct_per_K"], s["sim_sd_pct_per_K"]))
        sn = ex_n[k]
        mc.append(mc_entry(f"C1 simulated normal {k[0]} | {k[1]} (% per K)", r["detectable_s_sim_normal_same_sd_pct_per_K"],
                           100 * q["sim_normal_same_sd_detectable_s_per_K"], sn["mean_pct_per_K"], sn["sd_pct_per_K"]))
    mcs2 = {(r["dataset"], r["band"], r["floor"]): r for r in MCS["C2"]}
    ex2 = {(r["dataset"], r["band"], r["floor"]): r for r in EX["C2 n_sim"]}
    for r in MCj["C2"]:
        k = (r["dataset"], r["band"], r["floor"])
        q = pc2[k]
        pbar = 0.5 * (r["power_sim_lognormal"] + q["power_sim_lognormal"])
        sd = math.sqrt(pbar * (1 - pbar) / 10000)
        e = mc_entry(f"C2 simulated power {' | '.join(k)}", r["power_sim_lognormal"], q["power_sim_lognormal"], pbar, sd)
        e["note"] = "spread: the binomial sd of a 10,000-dataset estimate; replicate_mean is the mean of the two runs"
        mc.append(e)
        s = mcs2.get(k) or ex2.get(k)
        mc.append(mc_entry(f"C2 simulated pairs {' | '.join(k)}", r["n_required_sim_lognormal"], q["pairs_needed_sim_lognormal_resampled"],
                           s["mean"], s["sd"]) | {"replicates": s["replicate_n_sim"]})
    zero_spread = []
    for e in mc:
        if e.get("diff_over_sqrt2_sd") is not None:
            e["within_3_sqrt2_sd"] = abs(e["diff_over_sqrt2_sd"]) <= 3.0
        elif e.get("replicate_sd") == 0:
            zero_spread.append({"quantity": e["quantity"], "checker": e["checker"], "pipeline": e["pipeline"],
                                "replicates": e.get("replicates"), "note": "all replicates equal; the spread is below one unit"})
    res["monte_carlo"] = mc
    sp = [e for e in mc if e.get("diff_over_sqrt2_sd") is not None]
    worst = max(sp, key=lambda e: abs(e["diff_over_sqrt2_sd"]))
    wp = max(sp, key=lambda e: abs(e["z_pipeline"]))
    wc = max(sp, key=lambda e: abs(e["z_checker"]))
    res["monte_carlo_summary"] = {
        "quantities": len(mc), "with_spread": len(sp),
        "max_abs_diff_over_sqrt2_sd": abs(worst["diff_over_sqrt2_sd"]), "where": worst["quantity"],
        "over_3_sqrt2_sd": [e["quantity"] for e in sp if abs(e["diff_over_sqrt2_sd"]) > 3],
        "over_2_sqrt2_sd": [e["quantity"] for e in sp if abs(e["diff_over_sqrt2_sd"]) > 2],
        "pipeline_max_abs_z": abs(wp["z_pipeline"]), "pipeline_where": wp["quantity"],
        "pipeline_count_abs_z_over_2": sum(abs(e["z_pipeline"]) > 2 for e in sp),
        "checker_max_abs_z": abs(wc["z_checker"]), "checker_where": wc["quantity"],
        "checker_count_abs_z_over_2": sum(abs(e["z_checker"]) > 2 for e in sp),
        "expected_count_over_2_by_chance": round(0.0455 * len(sp), 1),
        "zero_spread_cases": zero_spread,
        "reading": ("z is a run's distance from the mean of the checker's independent replicate runs, in replicate sd; the pipeline's "
                    "values are the test of consistency, the checker's own registered run is one draw of the same distribution")}

    # ---------------- summary table
    det_classes = [("pairs: categorical", sum(v["mismatches"] for v in pr["categorical"].values()), None, None, None),
                   ("pairs: numeric", None, max(v["max_abs_diff"] for v in pr["numeric"].values()),
                    max(v["max_rel_diff"] for v in pr["numeric"].values()),
                    max(pr["numeric"].items(), key=lambda kv: kv[1]["max_rel_diff"])[0])]
    for k, v in res["B_point"].items():
        det_classes.append((f"B point {k}", None, v["max_abs_diff"], v["max_rel_diff"], v["where_max_rel"]))
    det_classes.append(("counts", res["counts"]["mismatches"], None, None, None))
    det_classes.append(("sign counts", res["sign"]["counts"]["mismatches"], None, None, None))
    det_classes.append(("sign p and share", None, res["sign"]["p_and_share"]["max_abs_diff"], res["sign"]["p_and_share"]["max_rel_diff"],
                        res["sign"]["p_and_share"]["where_max_rel"]))
    for k, v in res["C1_calculated"].items():
        if k == "counts":
            det_classes.append(("C1 counts", v["mismatches"], None, None, None))
        else:
            det_classes.append((f"C1 {k}", None, v["max_abs_diff"], v["max_rel_diff"], v["where_max_rel"]))
    for k, v in res["C2_calculated"].items():
        if k == "pairs_needed":
            det_classes.append(("C2 pairs needed (calc, shared, normal shape)", v["mismatches"], None, None, None))
        else:
            det_classes.append((f"C2 {k}", None, v["max_abs_diff"], v["max_rel_diff"], v["where_max_rel"]))
    det_classes.append(("verdicts P1 to P8", res["verdicts"]["P1_to_P8"]["mismatches"], None, None, None))
    det_classes.append(("B4 floor verdicts", res["verdicts"]["floor_verdicts"]["mismatches"], None, None, None))
    res["deterministic_all_agree"] = (all((c[1] in (None, 0)) for c in det_classes) and
                                      all((c[3] is None or c[3] <= 1e-12) for c in det_classes))
    res["table"] = [{"class": c[0], "mismatches": c[1], "max_abs_diff": c[2], "max_rel_diff": c[3], "where_max_rel": c[4]}
                    for c in det_classes]
    os.makedirs(OUTC, exist_ok=True)
    with open(os.path.join(OUTC, "compare.json"), "w") as f:
        json.dump(res, f, indent=1)
    with open(os.path.join(OUTC, "compare_table.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["class", "mismatches", "max_abs_diff", "max_rel_diff", "where_max_rel"])
        w.writeheader()
        for row in res["table"]:
            w.writerow(row)
    for row in res["table"]:
        print(row)
    print("deterministic all agree:", res["deterministic_all_agree"])
    print("MC:", json.dumps(res["monte_carlo_summary"]))
    print("V4 differences:", json.dumps(v4_notes))
    print("pairs one-sided:", {k: len(v["present_on_one_side_only"]) for k, v in pr["numeric"].items()})
    print("band mismatches:", pr["categorical"]["band"]["mismatches"], pr["categorical"]["band"]["list"][:6])


if __name__ == "__main__":
    main()
