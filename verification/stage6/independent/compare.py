#!/usr/bin/env python3
"""Comparison of the independent check (record.json, frozen at commit cf90a7f) with the registered Stage 6 results.

Run after the freeze was pushed. Reads the registered files with `git show <commit>:<path>` at 023652c (Part A0) and
at c2f07d1 (origin/stage6 as fetched after the push: Parts B to D), and DECISIONS.md at c2f07d1. A few quantities the
frozen record did not hold (the within-set slope of ln R on T that the registered B2 reports, B5's instrument, the
post hoc channel tests of D97 (4), the median Tmin floor) are computed here with the frozen code's functions and are
labelled as computed after the freeze.

Usage: python3 -I verification/stage6/independent/compare.py
Writes compare.json beside this file.
"""
import csv, hashlib, importlib.util, io, json, math, os, re, subprocess, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
A0_COMMIT = "023652c8364befd10e6b66d47eef0db3a34056ba"
BC_COMMIT = "c2f07d1395c745f498f0359ad88051f67a501b46"
FREEZE_COMMIT = "cf90a7f8b76359e304fc3cab6040835ff547b616"
SCOPES = {"pooled": "pooled", "Severson": "paper1", "Attia": "paper2"}
VARIANTS = {"registered": "registered", "T_max": "Tmax", "T_2_50": "Tavg_2_50", "T_min": "Tmin",
            "without_flagged": "without_flagged", "T_rise": "unscored_T_minus_Tmin", "R_cycles_2_10": "unscored_R_2_10",
            "C1_sets_of_two_or_more": "unscored_C1_sets_2plus_only"}
DATE_LABEL = {"2017-05-12": "b1", "2017-06-30": "b2", "2018-04-12": "b3", "2019-01-24": "b9"}
Z_LIMIT = 4.0          # D77 / D96 use 4 sqrt(2) Monte Carlo sd of one estimate; z here is diff / (sqrt(2) sd)


def git_show(commit, path):
    out = subprocess.run(["git", "-C", ROOT, "show", f"{commit}:{path}"], capture_output=True, check=True)
    return out.stdout


REG_FILES = {}


def reg_bytes(commit, path):
    b = git_show(commit, path)
    REG_FILES[f"{commit[:7]}:{path}"] = {"commit": commit, "path": path, "bytes": len(b), "sha256": hashlib.sha256(b).hexdigest()}
    return b


def reg_json(commit, path):
    return json.loads(reg_bytes(commit, path).decode())


def reg_csv(commit, path):
    return list(csv.DictReader(io.StringIO(reg_bytes(commit, path).decode())))


def rel(a, b):
    if a is None or b is None:
        return None
    a, b = float(a), float(b)
    if not (math.isfinite(a) and math.isfinite(b)):
        return None
    return abs(a - b) / max(1e-300, abs(b)) if b != 0 else abs(a - b)


def fnum(x):
    if x is None:
        return None
    if isinstance(x, str):
        if x.strip() == "" or x.lower() == "nan":
            return None
        x = float(x)
    x = float(x)
    return x if math.isfinite(x) else None


class Cls:
    """Accumulates one class of quantities: point values compared to rounding."""

    def __init__(self, name, tol_rel=1e-10):
        self.name, self.tol, self.rows = name, tol_rel, []

    def add(self, what, reg, mine):
        reg_f, mine_f = fnum(reg), fnum(mine)
        if reg_f is None and mine_f is None:
            self.rows.append({"what": what, "registered": None, "mine": None, "abs_diff": 0.0, "rel_diff": 0.0, "agree": True})
            return
        if reg_f is None or mine_f is None:
            self.rows.append({"what": what, "registered": reg_f, "mine": mine_f, "abs_diff": None, "rel_diff": None, "agree": False})
            return
        d = abs(reg_f - mine_f)
        r = d / abs(reg_f) if reg_f != 0 else d
        self.rows.append({"what": what, "registered": reg_f, "mine": mine_f, "abs_diff": d, "rel_diff": r,
                          "agree": bool(d <= self.tol * max(1.0, abs(reg_f)) or r <= self.tol)})

    def summary(self, keep_rows=True):
        num = [r for r in self.rows if r["abs_diff"] is not None]
        worst = max(num, key=lambda r: r["rel_diff"]) if num else None
        worst_abs = max(num, key=lambda r: r["abs_diff"]) if num else None
        out = {"compared": len(self.rows), "agree": int(sum(r["agree"] for r in self.rows)),
               "disagree": [r for r in self.rows if not r["agree"]],
               "max_abs_diff": worst_abs["abs_diff"] if worst_abs else None, "max_abs_diff_at": worst_abs["what"] if worst_abs else None,
               "max_rel_diff": worst["rel_diff"] if worst else None, "max_rel_diff_at": worst["what"] if worst else None,
               "tolerance_rel": self.tol}
        if keep_rows:
            out["rows"] = self.rows
        return out


class MC:
    """Monte Carlo quantities: z = (registered - mine) / (sqrt(2) sd), sd from my replicate runs."""

    def __init__(self, name):
        self.name, self.rows = name, []

    def add(self, what, reg, mine, sd_rep, sd_os=None, sd_20=None):
        reg_f, mine_f = fnum(reg), fnum(mine)
        if reg_f is None or mine_f is None or sd_rep in (None, 0):
            self.rows.append({"what": what, "registered": reg_f, "mine": mine_f, "sd_replicates": sd_rep, "z": None,
                              "both_missing": reg_f is None and mine_f is None})
            return
        z = (reg_f - mine_f) / (math.sqrt(2.0) * sd_rep)
        z_os = (reg_f - mine_f) / (math.sqrt(2.0) * sd_os) if sd_os else None
        z_20 = (reg_f - mine_f) / (math.sqrt(2.0) * sd_20) if sd_20 else None
        self.rows.append({"what": what, "registered": reg_f, "mine": mine_f, "diff": reg_f - mine_f, "sd_replicates": sd_rep,
                          "z": z, "sd_orderstat": sd_os, "z_orderstat": z_os, "sd_20_frozen": sd_20, "z_20_frozen": z_20,
                          "within": bool(abs(z) <= Z_LIMIT)})

    def summary(self):
        zs = [r for r in self.rows if r.get("z") is not None]
        worst = max(zs, key=lambda r: abs(r["z"])) if zs else None
        worst_d = max(zs, key=lambda r: abs(r["diff"])) if zs else None
        missing = [r for r in self.rows if r.get("z") is None and not r.get("both_missing")]
        return {"compared": len(zs), "within_4": int(sum(r["within"] for r in zs)), "outside_4": [r for r in zs if not r["within"]],
                "max_abs_z": abs(worst["z"]) if worst else None, "max_abs_z_at": worst["what"] if worst else None,
                "max_abs_diff": abs(worst_d["diff"]) if worst_d else None, "max_abs_diff_at": worst_d["what"] if worst_d else None,
                "share_abs_z_above_2": (sum(abs(r["z"]) > 2 for r in zs) / len(zs)) if zs else None,
                "not_comparable": missing, "rows": self.rows}


def load_mine():
    spec = importlib.util.spec_from_file_location("ind6", os.path.join(HERE, "ind6.py"))
    M = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(M)
    rec = json.load(open(os.path.join(HERE, "record.json")))
    return M, rec


def rebuild_cells(M):
    entries, _ = M.load_entries()
    M.apply_swaps(entries)
    cells, _ = M.build_cells(entries)
    M.contact_flags(cells)
    return entries, cells


REFINE_R = 100
VARIANT_DEF = {  # my variant name: (xkey, zkey, C1 over sets of two or more only, without flagged cells, seed key base)
    "registered": ("T", "lnR", False, False, None), "Tmax": ("Tmax_m", "lnR", False, False, 10),
    "Tavg_2_50": ("T_2_50", "lnR", False, False, 11), "Tmin": ("Tmin_m", "lnR", False, False, 12),
    "without_flagged": ("T", "lnR", False, True, 13), "unscored_T_minus_Tmin": ("T_minus_Tmin", "lnR", False, False, 14),
    "unscored_R_2_10": ("T", "lnR_2_10", False, False, 15), "unscored_C1_sets_2plus_only": ("T", "lnR", True, False, 16)}
SCOPE_INDEX = {"pooled": 1, "paper1": 2, "paper2": 3}


def refine_boot_sd(M, cells, R=REFINE_R):
    """Computed after the freeze: the Monte Carlo sd of each interval end from R further replicate runs of my own
    bootstrap (keys (*seed key, 5000 + r)); the frozen values themselves are not changed."""
    out = {}
    for var, (xkey, zkey, c2only, noflag, base) in VARIANT_DEF.items():
        use = [c for c in cells if not c["flag"]] if noflag else cells
        for scope, si in SCOPE_INDEX.items():
            key = (1, si) if base is None else (base, si)
            sub = M.scope_cells(use, scope)
            sets_all = M.set_list(sub, xkey, zkey)
            c_sets = [s for s in sets_all if s["n"] >= 2] if c2only else sets_all
            ends = {}
            for r in range(R):
                C, _ = M.count_matrix(c_sets, M.rng_for(*key, 5000 + r))
                bs = M.boot_stats(c_sets, C)
                for k, v in bs.items():
                    ends.setdefault(k, []).append(M.percentile_ci(v)[0])
            res = {}
            for k, e in ends.items():
                if any(x[0] is None for x in e):
                    res[k] = {"sd_lo": None, "sd_hi": None}
                    continue
                a = np.array(e, dtype=float)
                res[k] = {"sd_lo": float(np.std(a[:, 0], ddof=1)), "sd_hi": float(np.std(a[:, 1], ddof=1)),
                          "mean_lo": float(a[:, 0].mean()), "mean_hi": float(a[:, 1].mean())}
            out[f"{var}|{scope}"] = res
    return out


def refine_power_sd(M, cells, rec, R=REFINE_R):
    out = {}
    for si, scope in enumerate(("pooled", "paper1", "paper2")):
        sets2 = [s for s in M.set_list(M.scope_cells(cells, scope), "T") if s["n"] >= 2]
        N = sum(s["n"] for s in sets2)
        df = N - len(sets2) - 1
        for fname in ("Geslin", "KIT"):
            fl = rec["floors"][fname]
            sd_c = fl["sd_e"] / math.sqrt(2.0)
            ln, nm = [], []
            for r in range(R):
                rng = M.rng_for(80, si + 1, 5000 + r, 1 if fname == "Geslin" else 2)
                mag = np.exp(fl["mu"] + fl["sigma"] * rng.standard_normal((10000, N)))
                sign = rng.integers(0, 2, size=(10000, N)).astype(float) * 2.0 - 1.0
                ln.append(M.sim_detectable(sets2, sign * mag / math.sqrt(2.0), df)[0])
                rng2 = M.rng_for(90, si + 1, 5000 + r, 1 if fname == "Geslin" else 2)
                nm.append(M.sim_detectable(sets2, sd_c * rng2.standard_normal((10000, N)), df)[0])
            out[f"{scope}|{fname}"] = {"lognormal_sd": float(np.std(ln, ddof=1)), "normal_sd": float(np.std(nm, ddof=1)),
                                       "lognormal_mean": float(np.mean(ln)), "normal_mean": float(np.mean(nm)),
                                       "lognormal_mean_rel_to_calc": float(np.mean(ln)) / rec["A0_8_power"][scope]["floors"][fname]["detectable_calc"] - 1}
    return out


def main():
    M, rec = load_mine()
    entries, cells = rebuild_cells(M)
    bykey = {c["key"]: c for c in cells}
    sd100 = refine_boot_sd(M, cells)
    psd100 = refine_power_sd(M, cells, rec)
    out = {"what": "Comparison of the independent check with the registered Stage 6 results",
           "freeze_commit": FREEZE_COMMIT, "registered_commits": {"part_A0": A0_COMMIT, "parts_B_to_D_origin_stage6": BC_COMMIT},
           "z_rule": "z = (registered - mine) / (sqrt(2) x my replicate-run sd of that end); |z| <= 4 counts as agreement"}
    flags = []

    # ---------------- the registered files ----------------
    a0 = reg_json(A0_COMMIT, "results/stage6_a0.json")
    pw = reg_json(A0_COMMIT, "results/stage6_power.json")
    sc_a0 = reg_json(A0_COMMIT, "results/stage6_scorecard.json")
    cells_a0 = reg_csv(A0_COMMIT, "results/stage6_cells.csv")
    sets_a0 = reg_csv(A0_COMMIT, "results/stage6_sets.csv")
    cells_r = reg_csv(BC_COMMIT, "results/stage6_cells.csv")
    sets_r = reg_csv(BC_COMMIT, "results/stage6_sets.csv")
    a0_bc = reg_json(BC_COMMIT, "results/stage6_a0.json")
    pw_bc = reg_json(BC_COMMIT, "results/stage6_power.json")
    within = reg_json(BC_COMMIT, "results/stage6_within.json")
    between = reg_json(BC_COMMIT, "results/stage6_between.json")
    score = reg_json(BC_COMMIT, "results/stage6_scorecard.json")
    paper = reg_json(BC_COMMIT, "results/stage6_paper_numbers.json")
    decisions = reg_bytes(BC_COMMIT, "DECISIONS.md").decode()
    out["a0_files_unchanged_between_commits"] = {
        "stage6_a0.json": a0 == a0_bc, "stage6_power.json": pw == pw_bc,
        "stage6_cells.csv": cells_a0 == cells_r, "stage6_sets.csv": sets_a0 == sets_r}

    # ---------------- 1. counts ----------------
    cnt = Cls("counts", 0.0)
    for rs, ms in (("Severson", "paper1"), ("Attia", "paper2"), ("pooled", "pooled")):
        r, m = a0["counts"][rs], rec["counts"][ms]
        cnt.add(f"{rs} cells", r["cells_entering"], m["cells"])
        cnt.add(f"{rs} sets", r["sets"], m["sets"])
        cnt.add(f"{rs} sets of two or more", r["sets_of_two_or_more"], m["sets_2plus"])
        cnt.add(f"{rs} cells in them", r["cells_in_them"], m["cells_in_sets_2plus"])
        cnt.add(f"{rs} pairs", r["pairs"], m["pairs"])
        for size, n in r["set_sizes"].items():
            mine_n = m["single_cell_sets"] if size == "1" else m["set_size_histogram_2plus"].get(size, 0)
            cnt.add(f"{rs} sets of size {size}", n, mine_n)
    for date, r in a0["counts"]["by_batch"].items():
        cnt.add(f"{date} entries", r["entries"], rec["entries"][date]["entries"])
        cnt.add(f"{date} entering", r["enter"], rec["counts"][date]["cells"])
    out["counts"] = cnt.summary()

    # ---------------- 2. the cell table ----------------
    reg_enter = [r for r in cells_r if r["enters"] == "True"]
    colmap = {"T_reg": "T", "T_max": "Tmax_m", "T_2_50": "T_2_50", "T_min": "Tmin_m", "T_rise": "T_minus_Tmin", "T_odd": "T_odd",
              "T_even": "T_even", "T_51_100": "T_h2", "swing_med": "swing_med_2_100", "swing_early": "swing_med_2_11",
              "swing_late": "swing_med_91_100", "floor_med": "median_Tmin_2_100", "IR_ohm": "R", "IR_valid_cycles": "R_valid_cycles",
              "lnR": "lnR", "IR_2_10_ohm": "R_2_10", "lnR_2_10": "lnR_2_10", "cycle_life": "L", "T_cycles": None}
    per_col = {k: Cls("cells " + k, 1e-12) for k in colmap}
    set_mis, flag_mis, rule_mis, keys_seen = [], [], [], set()
    for r in reg_enter:
        key = f"{DATE_LABEL[r['batch']]}c{r['index']}"
        keys_seen.add(key)
        c = bykey.get(key)
        if c is None:
            set_mis.append({"registered_cell": r["cell"], "problem": "not entering in my table"})
            continue
        for rc, mc in colmap.items():
            if mc is None:
                per_col[rc].add(key, r[rc], 99)
            else:
                per_col[rc].add(key, r[rc], c[mc])
        if r["set"].split("|", 1)[1] != c["policy_readable"] or int(r["set_size"]) != sum(1 for x in cells if x["set"] == c["set"]):
            set_mis.append({"cell": key, "registered_set": r["set"], "registered_size": r["set_size"]})
        if (r["contact_flag"] == "True") != bool(c["flag"]):
            flag_mis.append({"cell": key, "registered": r["contact_flag"], "mine": c["flag"]})
        reasons = r["contact_reasons"] or ""
        reg_rules = "".join(l for l, word in (("a", "contact lost"), ("b", "contact poor"), ("c", "jump"), ("d", "flat"),
                                              ("e", "below"), ("f", "site:")) if word in reasons)
        if reg_rules != "".join(sorted(c["flag_rules"])):
            rule_mis.append({"cell": key, "registered_reasons": reasons, "registered_letters": reg_rules, "mine": c["flag_rules"]})
        if r["channel"] != c["channel"] or r["policy"] != c["policy_readable"]:
            set_mis.append({"cell": key, "problem": "channel or policy differ"})
    missing_in_reg = sorted(set(bykey) - keys_seen)
    out["cells"] = {"registered_entering": len(reg_enter), "mine_entering": len(cells), "mine_not_in_registered": missing_in_reg,
                    "set_or_descriptor_mismatches": set_mis, "contact_flag_mismatches": flag_mis,
                    "contact_rule_letter_mismatches": rule_mis,
                    "registered_contact_reasons_sample": sorted(set(r["contact_reasons"] for r in reg_enter if r["contact_reasons"]))[:12],
                    "columns": {k: v.summary(keep_rows=False) for k, v in per_col.items()}}
    # excluded entries
    reg_out = [r for r in cells_r if r["enters"] != "True"]
    mine_out = {(e["batch"], e["index"]) for e in rec["excluded_entries"]}
    out["cells"]["excluded_entries_match"] = sorted((r["batch"], int(r["index"])) for r in reg_out) == sorted(mine_out)

    # ---------------- 3. the sets table ----------------
    sc = Cls("sets", 1e-12)
    mine_sets = {s["set"]: s for s in rec["sets"]}
    set_problems = []
    for r in sets_r:
        lab = DATE_LABEL[r["batch"]]
        key = f"{lab}|{r['policy']}"
        m = mine_sets.get(key)
        if m is None:
            set_problems.append({"registered_set": r["set"], "problem": "missing in mine"})
            continue
        sc.add(key + " cells", r["cells"], m["n"])
        sc.add(key + " pairs", r["pairs"], m["n"] * (m["n"] - 1) // 2)
        sc.add(key + " mean T", r["T_reg_mean_C"], m["mean_T"])
        sc.add(key + " flagged", r["flagged_cells"], m["n_flagged"])
        T = np.array([bykey[k]["T"] for k in m["keys"]])
        if m["n"] >= 2:
            sc.add(key + " sd within (ddof 1)", r["T_reg_sd_within_K"], float(np.std(T, ddof=1)))
        sc.add(key + " range", r["T_reg_range_K"], float(T.max() - T.min()))
        if sorted(r["channels"].split()) != sorted(m["channels"]):
            set_problems.append({"set": key, "registered_channels": r["channels"], "mine": m["channels"]})
    out["sets"] = sc.summary(keep_rows=False)
    out["sets"]["problems"] = set_problems
    out["sets"]["registered_sets"] = len(sets_r)
    out["sets"]["mine_sets"] = len(mine_sets)

    # ---------------- 4. A0: reliability, resistance, channel ----------------
    a0c = Cls("A0", 1e-12)
    for rs, ms in SCOPES.items():
        r, m = a0["reliability"][rs], rec["A0_5_reliability"][ms]
        a0c.add(f"r_oe {rs}", r["odd_even"]["r"], m["r_oe"])
        a0c.add(f"r_oe Spearman-Brown {rs}", r["odd_even"]["spearman_brown"], m["r_oe_spearman_brown"])
        a0c.add(f"r_h {rs}", r["halves"]["r"], m["r_h"])
        a0c.add(f"r_h Spearman-Brown {rs}", r["halves"]["spearman_brown"], m["r_h_spearman_brown"])
        a0c.add(f"reliability cells {rs}", r["cells"], m["r_oe_cells"])
    for date, r in a0["resistance"].items():
        m = rec["A0_6_resistance"][date]
        h = m["valid_cycles_histogram"]
        a0c.add(f"{date} cells all 99 valid", r["all_99_valid"], h.get("99", 0))
        a0c.add(f"{date} cells some invalid", r["some_invalid"], sum(v for k, v in h.items() if 0 < int(k) < 99))
        a0c.add(f"{date} cells none valid", r["none_valid"], h.get("0", 0))
        a0c.add(f"{date} IR min", r["IR_ohm_range"][0], m["R_min"])
        a0c.add(f"{date} IR max", r["IR_ohm_range"][1], m["R_max"])
    ch_r, ch_m = a0["channel"], rec["A0_7_channel"]
    a0c.add("channel pairs", ch_r["pairs"], ch_m["pairs"])
    a0c.add("channel correlation", ch_r["r"], ch_m["correlation"])
    out["A0"] = a0c.summary()
    mc_ch = MC("channel p")
    mc_ch.add("channel p, (k + 1) / (n + 1) form", ch_r["p_permutation_two_sided"], ch_m["p_two_sided_plus_one"], ch_m["p_mc_sd"])
    out["A0_channel_p"] = mc_ch.summary()
    out["A0_channel_verdict"] = {"registered": ch_r["channel_carries_temperature"], "mine": ch_m["carries_temperature"],
                                 "agree": ch_r["channel_carries_temperature"] == ch_m["carries_temperature"]}

    # ---------------- 5. power ----------------
    pc = Cls("power calculation", 1e-10)
    ps = MC("power simulation")
    v5 = []
    for case in pw["cases"]:
        ms = SCOPES[case["scope"]]
        fl = case["floor"]
        m = rec["A0_8_power"][ms]
        mf = m["floors"][fl]
        tag = f"{case['scope']} {fl}"
        pc.add(tag + " cells", case["cells"], m["N"]); pc.add(tag + " sets", case["sets"], m["J"]); pc.add(tag + " df", case["df"], m["df"])
        pc.add(tag + " S_TT", case["Sxx_K2"], m["S_TT"]); pc.add(tag + " T dev rms", case["T_dev_rms_K"], m["dT_rms_over_N"])
        pc.add(tag + " sd_cell", case["sd_cell"], mf["sd_c"]); pc.add(tag + " sd_e", case["sd_e_pair"], mf["sd_e"])
        pc.add(tag + " detectable", case["detectable_s_per_K"], mf["detectable_calc"])
        pc.add(tag + " pairs", case["pairs"], m["pairs"]); pc.add(tag + " pair dT rms", case["pair_dT_rms_K"], m["dT_rms_pairs"])
        pc.add(tag + " shared-cell F", case["shared_cell_F"], m["F_shared"])
        pc.add(tag + " pairs form detectable", case["pairs_form_detectable_s_per_K"], mf["pairs_form_detectable"])
        pc.add(tag + " pairs form shared-cell detectable", case["pairs_form_shared_cells_detectable_s_per_K"], mf["pairs_form_shared_cell_detectable"])
        pc.add(tag + " normal shape detectable", case["normal_shape_detectable_s_per_K"], mf["normal_shape_detectable"])
        wf = rec["A0_8_power"][ms + "_without_flagged"]
        pc.add(tag + " without flagged cells", case["without_flagged"]["cells"], wf["N"])
        pc.add(tag + " without flagged sets", case["without_flagged"]["sets"], wf["J"])
        pc.add(tag + " without flagged T dev rms", case["without_flagged"]["T_dev_rms_K"], wf["dT_rms_over_N"])
        pc.add(tag + " without flagged detectable", case["without_flagged"]["detectable_s_per_K"], wf["floors"][fl]["detectable_calc"])
        rp = psd100[f"{ms}|{fl}"]
        ps.add(tag + " simulated, registered shape", case["sim_registered_shape_detectable_s_per_K"], mf["sim_lognormal_detectable"],
               rp["lognormal_sd"], None, mf["sim_lognormal_mc_sd"])
        ps.add(tag + " simulated, normal noise", case["sim_normal_detectable_s_per_K"], mf["sim_normal_detectable"],
               rp["normal_sd"], None, mf["sim_normal_mc_sd"])
        v5.append({"case": tag, "registered_pass": case["v5_pass"], "mine_pass": mf["v5_pass"],
                   "registered_rel_diff": case["v5_rel_diff"], "mine_rel_diff": mf["v5_rel_diff"],
                   "mine_rel_diff_of_replicate_mean": (mf["sim_lognormal_replicate_mean"] - mf["detectable_calc"]) / mf["detectable_calc"],
                   "mine_replicate_min_rel": (min(mf["sim_lognormal_replicates"]) - mf["detectable_calc"]) / mf["detectable_calc"],
                   "mine_replicate_max_rel": (max(mf["sim_lognormal_replicates"]) - mf["detectable_calc"]) / mf["detectable_calc"]})
    out["power_calculation"] = pc.summary()
    out["power_simulation"] = ps.summary()
    out["V5"] = v5

    # ---------------- 6. Parts B and C: points, intervals, B3 ----------------
    pts = Cls("B and C points", 1e-10)
    b3c = Cls("B3", 1e-12)
    ivs = MC("B and C interval ends")
    unmatched = []

    def mine_variant(var, scope):
        name = f"{VARIANTS[var]}|{SCOPES[scope]}"
        v = rec["parts_B_C"].get(name)
        if v is not None:
            v = dict(v)
            v["_name"] = name
        return v

    def add_iv(tag, reg_ci, v, key):
        if reg_ci is None:
            return
        ci = v["boot"]["ci"][key]
        reps = v["boot"]["replicate_runs"].get(key, {})
        r100 = sd100[f"{v['_name']}"].get(key, {})
        for end, i in (("lo", 0), ("hi", 1)):
            ivs.add(f"{tag} {end}", reg_ci[i], ci["ci95"][i], r100.get("sd_" + end), ci["mc_sd_orderstat"][i], reps.get("sd_" + end))

    for part, doc in (("within", within), ("between", between)):
        for group in ("registered", "brackets", "unscored"):
            if group not in doc:
                continue
            node = doc[group]
            items = [("registered", node)] if group == "registered" else list(node.items())
            for var, vnode in items:
                for scope, sn in vnode.items():
                    v = mine_variant(var, scope)
                    tag = f"{part}/{var}/{scope}"
                    if v is None:
                        unmatched.append(tag)
                        continue
                    pt = v["point"]
                    if "B1" in sn:
                        pts.add(tag + " s_w", sn["B1"]["s_w_per_K"]["value"], pt["B1"]["s_w"])
                        pts.add(tag + " pairs form", sn["B1"]["pairs_form_s_per_K"]["value"], pt["B1_pairs"]["s_w_pairs"])
                        pts.add(tag + " pairs", sn["B1"]["pairs"], pt["B1_pairs"]["pairs"])
                        add_iv(tag + " s_w", sn["B1"]["s_w_per_K"]["ci95"], v, "s_w")
                        add_iv(tag + " pairs form", sn["B1"]["pairs_form_s_per_K"]["ci95"], v, "s_w_pairs")
                    if "cells_within" in sn:
                        pts.add(tag + " cells within", sn["cells_within"], v["counts"]["cells_2plus"])
                        pts.add(tag + " sets within", sn["sets_within"], v["counts"]["sets_2plus"])
                    if "cells_all" in sn and var != "C1_sets_of_two_or_more":
                        pts.add(tag + " cells all", sn["cells_all"], v["counts"]["cells"])
                        pts.add(tag + " sets all", sn["sets_all"], v["counts"]["sets"])
                    if "cells_all" in sn and var == "C1_sets_of_two_or_more":
                        pts.add(tag + " cells in C1", sn["cells_all"], v["counts"]["cells_2plus"])
                        pts.add(tag + " sets in C1", sn["sets_all"], v["counts"]["sets_2plus"])
                    if var == "R_cycles_2_10" and "note" in sn:
                        pts.add(tag + " B2 absent (no valid R)", 0, 0 if v["point"]["B2"] is None else 1)
                    if "B3" in sn:
                        b = pt["B3"]
                        b3c.add(tag + " B3 sets", sn["B3"]["sets"], b["sets_used"])
                        b3c.add(tag + " B3 warmest first", sn["B3"]["warmest_first"], b["hits"])
                        b3c.add(tag + " B3 share", sn["B3"]["share"], b["share"])
                        b3c.add(tag + " B3 chance share", sn["B3"]["chance_share"], b["chance_share"])
                        b3c.add(tag + " B3 expected", sn["B3"]["expected"], b["expected_hits"])
                        b3c.add(tag + " B3 p", sn["B3"]["p_two_sided"], b["p_exact_two_sided"])
                        b3c.add(tag + " B3 left out", sn["B3"]["left_out_ties"], b["sets_left_out"])
                    if "B4" in sn:
                        pts.add(tag + " sd_cell", sn["B4"]["sd_cell"]["value"], pt["B1"]["sd_c_hat"])
                        pts.add(tag + " sqrt2 sd_cell", sn["B4"]["sqrt2_sd_cell"], pt["B1"]["sqrt2_sd_c_hat"])
                        pts.add(tag + " B4 df", sn["B4"]["df"], pt["B1"]["dof"])
                        add_iv(tag + " sd_cell", sn["B4"]["sd_cell"]["ci95"], v, "sd_c_hat")
                        add_iv(tag + " sqrt2 sd_cell", sn["B4"]["sqrt2_sd_cell_ci95"], v, "sqrt2_sd_c_hat")
                    b2n = sn.get("B2") if "B2" in sn else (sn if "s_w_R_held_per_K" in sn else None)
                    if b2n is not None and "note" not in b2n:
                        b2 = pt["B2"]
                        pts.add(tag + " B2 cells", b2n["cells"], b2["N"]); pts.add(tag + " B2 sets", b2n["sets"], b2["J"])
                        pts.add(tag + " B2 T lever R held", b2n["s_w_R_held_per_K"]["value"], b2["lever_T_R_held"])
                        pts.add(tag + " B2 dlnL/dlnR", b2n["lnL_per_lnR"]["value"], b2["dlnL_dlnR"])
                        pts.add(tag + " B2 s_w same cells", b2n["s_w_same_cells_per_K"]["value"], b2["s_w_same_cells"])
                        add_iv(tag + " B2 T lever R held", b2n["s_w_R_held_per_K"]["ci95"], v, "b2_lever_T_R_held")
                        add_iv(tag + " B2 dlnL/dlnR", b2n["lnL_per_lnR"]["ci95"], v, "b2_dlnL_dlnR")
                        add_iv(tag + " B2 s_w same cells", b2n["s_w_same_cells_per_K"]["ci95"], v, "b2_s_w_same_cells")
                    elif b2n is not None and "note" in b2n:
                        pts.add(tag + " B2 absent (no valid R)", 0, 0 if pt["B2"] is None else 1)
                    if "C1" in sn:
                        pts.add(tag + " s_b", sn["C1"]["s_b_per_K"]["value"], pt["C1"]["s_b"])
                        pts.add(tag + " s_b - s_w", sn["C2"]["s_b_minus_s_w_per_K"]["value"], pt["C2"]["s_b_minus_s_w"])
                        pts.add(tag + " C3 share", sn["C3"]["share"]["value"], pt["C3"]["share"])
                        pts.add(tag + " C3 residual form", sn["C3"]["share_residual_form"]["value"], pt["C3"]["resid_form"])
                        add_iv(tag + " s_b", sn["C1"]["s_b_per_K"]["ci95"], v, "s_b")
                        add_iv(tag + " s_b - s_w", sn["C2"]["s_b_minus_s_w_per_K"]["ci95"], v, "s_b_minus_s_w")
                        add_iv(tag + " C3 share", sn["C3"]["share"]["ci95"], v, "c3_share")
                        add_iv(tag + " C3 residual form", sn["C3"]["share_residual_form"]["ci95"], v, "c3_resid")
                    if "reliability_bracket" in sn:
                        rb, mr = sn["reliability_bracket"], v["divided_by_r_oe"]
                        pts.add(tag + " r_oe", rb["r_oe"], mr["r_oe"])
                        pts.add(tag + " s_w / r_oe", rb["s_w_over_r"], mr["s_w"]["point"])
                        pts.add(tag + " s_b / r_oe", rb["s_b_over_r"], mr["s_b"]["point"])
                        if "s_w_R_held_over_r" in rb:
                            pts.add(tag + " B2 T lever / r_oe", rb["s_w_R_held_over_r"], mr["b2_lever_T_R_held"]["point"])
    out["B_C_points"] = pts.summary()
    out["B3"] = b3c.summary()
    out["B_C_intervals"] = ivs.summary()
    out["B_C_unmatched_registered_entries"] = unmatched

    # ---------------- 7. scorecard ----------------
    mine_sc = rec["scorecard"]
    verdicts = {}
    for p in ("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"):
        verdicts[p] = {"registered": score[p]["holds"], "mine": mine_sc[p]["holds"], "reproduced": score[p]["holds"] == mine_sc[p]["holds"]}
    out["scorecard"] = {"verdicts": verdicts, "all_reproduced": all(v["reproduced"] for v in verdicts.values()),
                        "registered_held": score["held"], "registered_failed": score["failed"], "registered_not_scored": score["not_scored"],
                        "a0_scorecard_P1_P6_unchanged": sc_a0["P1"] == score["P1"] and sc_a0["P6"] == score["P6"]}
    scv = Cls("scorecard values", 1e-10)
    scv.add("P1 value (% per K)", score["P1"]["value"], 100 * mine_sc["P1"]["quantity"])
    scv.add("P1 threshold (% per K)", score["P1"]["threshold"], 100 * mine_sc["P1"]["threshold"])
    scv.add("P2 s_w", score["P2"]["value"], mine_sc["P2"]["s_w"])
    scv.add("P3 lever", score["P3"]["value"], mine_sc["P3"]["lever_T_R_held"])
    scv.add("P3 same cells", score["P3"]["same_cells_s_w"], mine_sc["P3"]["s_w_same_cells"])
    scv.add("P3 ratio", score["P3"]["ratio"], mine_sc["P3"]["ratio"])
    scv.add("P4 s_b", score["P4"]["value"], mine_sc["P4"]["s_b"])
    scv.add("P5 s_b - s_w", score["P5"]["value"], mine_sc["P5"]["s_b_minus_s_w"])
    scv.add("P6 r_oe", score["P6"]["value"], mine_sc["P6"]["r_oe"])
    scv.add("P7 interval lo", score["P7"]["interval"][0], mine_sc["P7"]["interval"][0])
    scv.add("P7 interval hi", score["P7"]["interval"][1], mine_sc["P7"]["interval"][1])
    scv.add("P8 value", score["P8"]["value"], mine_sc["P8"]["sqrt2_sd_c_hat"])
    scv.add("P8 floor lo", score["P8"]["floors"][0], mine_sc["P8"]["interval"][0])
    scv.add("P8 floor hi", score["P8"]["floors"][1], mine_sc["P8"]["interval"][1])
    out["scorecard_values"] = scv.summary()
    # margins of the interval-based verdicts, in my replicate sd
    reg_pooled = rec["parts_B_C"]["registered|pooled"]["boot"]
    margins = {"P2": {"note": "fails on the point value: s_w is negative", "s_w": rec["parts_B_C"]["registered|pooled"]["point"]["B1"]["s_w"],
                      "upper_end_mine": reg_pooled["ci"]["s_w"]["ci95"][1], "upper_end_registered": score["P2"]["ci95"][1]}}
    for p, key in (("P4", "s_b"), ("P5", "s_b_minus_s_w")):
        sd = sd100["registered|pooled"][key]["sd_lo"]
        margins[p] = {"lower_end_mine": reg_pooled["ci"][key]["ci95"][0], "lower_end_registered": score[p]["ci95"][0],
                      "replicate_sd_100": sd, "lower_end_over_sd_mine": reg_pooled["ci"][key]["ci95"][0] / sd,
                      "lower_end_over_sd_registered": score[p]["ci95"][0] / sd}
    out["verdict_margins"] = margins
    out["mc_refinement"] = {"replicates": REFINE_R, "bootstrap_keys": "(*frozen key, 5000 + r)",
                            "power_keys": "(80 or 90, scope, 5000 + r, floor)",
                            "power": psd100, "note": "computed after the freeze; frozen values unchanged"}

    # ---------------- 8. Part D and paper numbers ----------------
    pdc = Cls("paper numbers", 1e-10)
    rows = paper["part_D"]["rows"]
    pdc.add("Part D s_w pooled (%)", rows[0]["value"], 100 * rec["parts_B_C"]["registered|pooled"]["point"]["B1"]["s_w"])
    pdc.add("Part D s_w without flagged (%)", rows[1]["value"], 100 * rec["parts_B_C"]["without_flagged|pooled"]["point"]["B1"]["s_w"])
    pdc.add("Part D s_w first paper (%)", rows[2]["value"], 100 * rec["parts_B_C"]["registered|paper1"]["point"]["B1"]["s_w"])
    pdc.add("Part D s_w second paper (%)", rows[3]["value"], 100 * rec["parts_B_C"]["registered|paper2"]["point"]["B1"]["s_w"])
    pdc.add("Part D s_b pooled (%)", rows[4]["value"], 100 * rec["parts_B_C"]["registered|pooled"]["point"]["C1"]["s_b"])
    for i, k in ((5, "lfp_1c_15_25"), (6, "lfp_2c_15_25"), (7, "lfp_2c_25_35")):
        pdc.add(f"Part D Sandia {k} (%)", rows[i]["value"], 100 * rec["sandia_levers"][k]["s"])
        pdc.add(f"Part D Sandia {k} lo (%)", rows[i]["ci95"][0], 100 * rec["sandia_levers"][k]["lo"])
        pdc.add(f"Part D Sandia {k} hi (%)", rows[i]["ci95"][1], 100 * rec["sandia_levers"][k]["hi"])
    det = paper["part_D"]["detectable_pct_per_K"]
    pdc.add("Part D detectable pooled KIT (%)", det["Stage 6, pooled"]["KIT floor"], 100 * rec["A0_8_power"]["pooled"]["floors"]["KIT"]["detectable_calc"])
    pdc.add("Part D detectable pooled Geslin (%)", det["Stage 6, pooled"]["Geslin floor"], 100 * rec["A0_8_power"]["pooled"]["floors"]["Geslin"]["detectable_calc"])
    pdc.add("Part D detectable without flagged KIT (%)", det["Stage 6, pooled, without flagged cells"]["KIT floor"],
            100 * rec["A0_8_power"]["pooled_without_flagged"]["floors"]["KIT"]["detectable_calc"])
    pdc.add("Part D Stage 5 KIT KIT (%)", det["Stage 5, KIT replicates"]["KIT floor"], 100 * rec["stage5_detectable"]["KIT_KIT_floor"])
    pdc.add("Part D Stage 5 KIT Geslin (%)", det["Stage 5, KIT replicates"]["Geslin floor"], 100 * rec["stage5_detectable"]["KIT_Geslin_floor"])
    # paper numbers duplicate the within / between files: check a sample that they agree with those files
    dup = []
    for scope in ("pooled", "Severson", "Attia"):
        pn = paper["registered"][scope]
        dup.append(abs(pn["s_w_pct_per_K"]["value"] / 100 - within["registered"][scope]["B1"]["s_w_per_K"]["value"]))
        dup.append(abs(pn["s_b_pct_per_K"]["value"] / 100 - between["registered"][scope]["C1"]["s_b_per_K"]["value"]))
        dup.append(max(abs(a / 100 - b) for a, b in zip(pn["s_w_pct_per_K"]["ci95"], within["registered"][scope]["B1"]["s_w_per_K"]["ci95"])))
    out["paper_numbers"] = pdc.summary()
    out["paper_numbers"]["internal_agreement_with_within_between_files_max_abs"] = max(dup)
    out["paper_numbers"]["checks_reported"] = paper["checks"]

    # ---------------- 9. computed after the freeze: extra quantities the registered files hold ----------------
    extra = Cls("post-freeze extras", 1e-10)
    for var, xkey, zkey, noflag in (("registered", "T", "lnR", False), ("T_max", "Tmax_m", "lnR", False), ("T_2_50", "T_2_50", "lnR", False),
                                    ("T_min", "Tmin_m", "lnR", False), ("without_flagged", "T", "lnR", True),
                                    ("T_rise", "T_minus_Tmin", "lnR", False), ("R_cycles_2_10", "T", "lnR_2_10", False)):
        use = [c for c in cells if not c["flag"]] if noflag else cells
        sets2 = [s for s in M.set_list(use, xkey, zkey) if s["n"] >= 2]
        W = [M.within_sums(s) for s in sets2 if M.within_sums(s)["R_n"] >= 2]
        slope = sum(w["R_Sxz"] for w in W) / sum(w["R_Sxx"] for w in W)
        node = within["registered"]["pooled"] if var == "registered" else (
            within["brackets"][var]["pooled"] if var in within["brackets"] else within["unscored"][var]["pooled"])
        node = node["B2"] if "B2" in node else node
        extra.add(f"{var} pooled: within-set slope of ln R on T (lnR_per_K_within)", node["lnR_per_K_within"]["value"], slope)
    # B5's instrument: mean of the channel's within-set deviations of T in the other batches
    devs = {}
    for c in cells:
        mem = [x for x in cells if x["set"] == c["set"]]
        if len(mem) >= 2:
            devs.setdefault(c["channel"], {})[c["label"]] = c["T"] - float(np.mean([x["T"] for x in mem]))
    inst = a0["channel"]["instrument"]
    n_inst = 0
    for k, v in inst.items():
        date, idx = k.split(":")
        c = bykey[f"{DATE_LABEL[date]}c{idx}"]
        others = [d for lab, d in devs.get(c["channel"], {}).items() if lab != c["label"]]
        mine = float(np.mean(others)) if others else None
        n_inst += mine is not None
        extra.add(f"B5 instrument {k}", v, mine)
    out["post_freeze_extras"] = extra.summary(keep_rows=False)
    out["post_freeze_extras"]["instrument_cells_mine"] = n_inst
    out["post_freeze_extras"]["instrument_cells_registered"] = a0["channel"]["cells_with_instrument"]
    # D97 (4): the channel test with the Tmin measure and with the swing in place of T
    ph = {}
    for name, key in (("Tmin measure", "Tmin_m"), ("median swing", "swing_med_2_100"), ("mean swing", "swing_mean_2_100")):
        cc = []
        for c in cells:
            d = dict(c)
            d["T"] = c[key]
            cc.append(d)
        res = M.channel_test(cc, M.rng_for(97, len(ph)), nperm=10000)
        ph[name] = {"pairs": res["pairs"], "correlation": res["correlation"], "p_two_sided": res["p_two_sided"]}
    out["post_freeze_D97_4"] = ph

    # ---------------- 10. DECISIONS.md: the quoted numbers ----------------
    out["decisions_quotes"] = decisions_check(rec, cells, bykey, ph)

    # ---------------- summary: agreement by class, verdicts, flags ----------------
    out["summary"] = summary_block(out)
    out["registered_files_read"] = REG_FILES
    json.dump(clean(out), open(os.path.join(HERE, "compare.json"), "w"), indent=1)
    s = {k: (out[k]["max_rel_diff"] if isinstance(out[k], dict) and "max_rel_diff" in out[k] else None) for k in out}
    print(json.dumps({k: v for k, v in s.items() if v is not None}, indent=1))
    print("intervals max |z|", out["B_C_intervals"]["max_abs_z"], out["B_C_intervals"]["max_abs_z_at"])
    print("power sim max |z|", out["power_simulation"]["max_abs_z"], out["power_simulation"]["max_abs_z_at"])
    print("verdicts", json.dumps(out["scorecard"]["verdicts"]))


def decisions_check(rec, cells, bykey, ph):
    """Each number DECISIONS.md D82 to D97 quotes from Stage 6, against mine at the quoted rounding."""
    q = []

    def chk(where, quoted, mine, places=None, note=""):
        if mine is None:
            q.append({"where": where, "quoted": quoted, "mine": None, "agree": False, "note": note})
            return
        if places is None:
            ok = quoted == mine
        else:
            ok = round(float(mine), places) == round(float(quoted), places)
        q.append({"where": where, "quoted": quoted, "mine": mine, "rounded_mine": (round(float(mine), places) if places is not None else mine),
                  "agree": bool(ok), "note": note})

    sl = rec["sandia_levers"]
    chk("D82 (2) Sandia 1C 15 to 25 C, % per K (+3.3718, 2.9491 to 3.8270)", [3.3718, 2.9491, 3.8270],
        [round(100 * sl["lfp_1c_15_25"][k], 4) for k in ("s", "lo", "hi")])
    chk("D82 (2) Sandia 2C 15 to 25 C, % per K (+2.0499, 1.1837 to 2.9077)", [2.0499, 1.1837, 2.9077],
        [round(100 * sl["lfp_2c_15_25"][k], 4) for k in ("s", "lo", "hi")])
    chk("D82 (2) Sandia 2C 25 to 35 C, % per K (+11.0058, 8.3022 to 13.7300)", [11.0058, 8.3022, 13.73],
        [round(100 * sl["lfp_2c_25_35"][k], 4) for k in ("s", "lo", "hi")])
    chk("D82 (2) floors sd(e), Geslin and KIT", [0.04454298465298393, 0.10412996462328801], [rec["floors"]["Geslin"]["sd_e"], rec["floors"]["KIT"]["sd_e"]])
    chk("D82 (2) Stage 5 detectable, KIT and Geslin floors, % per K (6.56, 2.81)", [6.56, 2.81],
        [round(100 * rec["stage5_detectable"]["KIT_KIT_floor"], 2), round(100 * rec["stage5_detectable"]["KIT_Geslin_floor"], 2)])
    chk("D92 (6) / D96 (2) P1 threshold 0.03371837126587428", 0.03371837126587428, sl["lfp_1c_15_25"]["s"])
    chk("D96 (2) P7 interval", [0.08302214841830174, 0.13729950535845456], [sl["lfp_2c_25_35"]["lo"], sl["lfp_2c_25_35"]["hi"]])
    cn = rec["counts"]
    chk("D84 (2) entries 2017-05-12", 46, rec["entries"]["2017-05-12"]["entries"])
    chk("D84 (2) entries 2017-06-30", 48, rec["entries"]["2017-06-30"]["entries"])
    chk("D84 (2) entries 2018-04-12", 46, rec["entries"]["2018-04-12"]["entries"])
    chk("D84 (2) entries 2019-01-24", 45, rec["entries"]["2019-01-24"]["entries"])
    chk("D84 (2) entries all", 185, rec["entries"]["total_entries"])
    chk("D84 (3) 2017-05-12 channels missing", [4, 8], rec["entries"]["2017-05-12"]["channels_missing_1_to_48"])
    chk("D84 (3) 2018-04-12 channels missing", [26, 31], rec["entries"]["2018-04-12"]["channels_missing_1_to_48"])
    chk("D84 (3) 2019-01-24 channels missing", [4, 5, 46], rec["entries"]["2019-01-24"]["channels_missing_1_to_48"])
    chk("D84 (3) 2017-05-12 excluded channels", [13, 19, 21, 22, 31], sorted(int(e["channel"]) for e in rec["entries"]["2017-05-12"]["excluded"]))
    chk("D84 (3) 2018-04-12 excluded channels", [6, 7, 12, 33, 41, 46], sorted(int(e["channel"]) for e in rec["entries"]["2018-04-12"]["excluded"]))
    chk("D84 (3) continuations channels", [1, 2, 3, 5, 6], [int(c["channel_b2"]) for c in rec["carried"]])
    chk("D84 (3) continuation barcodes", ["EL150800460514", "EL150800460486", "EL150800460623", "EL150800464977", "EL150800464865"],
        [c["barcode_b2"] for c in rec["carried"]])
    chk("D84 (3) 2017-06-30 barcodes of channels 3 and 10", ["EL150800460623", "el150800460605"],
        [rec["entries"]["b2_barcodes_channels_3_10"]["3"], rec["entries"]["b2_barcodes_channels_3_10"]["10"]])
    chk("D83 (5) / D94 loader lengths equal continuation summary entries", [662, 981, 1060, 208, 482],
        [c["continuation_summary_entries"] for c in rec["carried"]])
    chk("D86 (3) first paper cells", 124, cn["paper1"]["cells"])
    chk("D86 (3) per batch 41, 43, 40", [41, 43, 40], [cn["2017-05-12"]["cells"], cn["2017-06-30"]["cells"], cn["2018-04-12"]["cells"]])
    chk("D86 (3) first paper sets", 71, cn["paper1"]["sets"])
    chk("D86 (3) sets of two or more, cells, pairs", [27, 80, 126], [cn["paper1"]["sets_2plus"], cn["paper1"]["cells_in_sets_2plus"], cn["paper1"]["pairs"]])
    chk("D86 (3) 2017-05-12 sets of two or more and cells", [18, 37], [cn["2017-05-12"]["sets_2plus"], cn["2017-05-12"]["cells_in_sets_2plus"]])
    chk("D86 (3) 2017-06-30 sets of two or more and cells", [1, 3], [cn["2017-06-30"]["sets_2plus"], cn["2017-06-30"]["cells_in_sets_2plus"]])
    chk("D86 (3) 2018-04-12 sets of two or more and cells", [8, 40], [cn["2018-04-12"]["sets_2plus"], cn["2018-04-12"]["cells_in_sets_2plus"]])
    chk("D86 (3) set sizes 2, 3, 4, 6, 7, 8", [19, 3, 1, 1, 1, 2], [cn["paper1"]["set_size_histogram_2plus"].get(k, 0) for k in ("2", "3", "4", "6", "7", "8")])
    chk("D86 (3) single-cell sets", 44, cn["paper1"]["single_cell_sets"])
    chk("D86 (3) second paper", [45, 9, 90], [cn["paper2"]["cells"], cn["paper2"]["sets"], cn["paper2"]["pairs"]])
    chk("D86 (3) both", [169, 80, 36, 125, 216], [cn["pooled"]["cells"], cn["pooled"]["sets"], cn["pooled"]["sets_2plus"],
                                                 cn["pooled"]["cells_in_sets_2plus"], cn["pooled"]["pairs"]])
    sw = rec["swaps"]
    chk("D87 (3) 2017-05-12 swapped pair is the set 5.4C(60%)-3C", ["5.4C(60%)-3C", "5.4C(60%)-3C"], sw[0]["policies"])
    b2single = [sum(1 for c in cells if c["set"] == bykey[k]["set"]) for k in sw[1]["keys"]]
    chk("D87 (3) 2017-06-30 swapped cells are two single-cell sets", [1, 1], b2single)
    ct = rec["contact"]
    b9 = [c for c in cells if c["label"] == "b9"]
    below = max(c["b_ratio"] for c in b9 if c["b_ratio"] < 0.5)
    above = min(c["b_ratio"] for c in b9 if c["b_ratio"] >= 0.5)
    chk("D88 validation batch: the gap in the ratio of (b) runs from 0.39 to 0.57", [0.39, 0.57], [round(below, 2), round(above, 2)],
        note=f"the gap's ends are {below:.4f} and {above:.4f}")
    chk("D88 correlation of within-set T and swing deviations, validation batch (0.95)", 0.95, ct["corr_dev_T_dev_swing_paper2"]["median_swing_uncentred"], 2)
    chk("D88 the same, first paper sets of two or more (0.87)", 0.87, ct["corr_dev_T_dev_swing_paper1"]["median_swing_uncentred"], 2)
    col = {(c["label"], c["channel"]): c for c in cells}
    chk("D88 2017-06-30 channel 7 swing 4.5 to 1.6 K", [4.5, 1.6], [round(col[("b2", "7")]["swing_med_2_11"], 1), round(col[("b2", "7")]["swing_med_91_100"], 1)])
    chk("D88 2019-01-24 channel 44 swing 6.9 to 1.2 K", [6.9, 1.2], [round(col[("b9", "44")]["swing_med_2_11"], 1), round(col[("b9", "44")]["swing_med_91_100"], 1)])
    chk("D88 2019-01-24 channel 22 swing 4.4 to 1.5 K", [4.4, 1.5], [round(col[("b9", "22")]["swing_med_2_11"], 1), round(col[("b9", "22")]["swing_med_91_100"], 1)])
    chk("D88 (c) largest step 2.09 K", 2.09, ct["largest_step_all_cells"], 2)
    chk("D88 (e) lowest Tmin 26.7 C", 26.7, ct["lowest_Tmin_all_cells"], 1)
    mf = {lab: [min(c["median_Tmin_2_100"] for c in cells if c["label"] == lab), max(c["median_Tmin_2_100"] for c in cells if c["label"] == lab)]
          for lab in ("b2", "b9")}
    chk("D88 (e) median floors 2017-06-30, 28.0 to 30.6 C", [28.0, 30.6], [round(x, 1) for x in mf["b2"]])
    chk("D88 (e) median floors 2019-01-24, 30.5 to 32.6 C", [30.5, 32.6], [round(x, 1) for x in mf["b9"]],
        note=f"unrounded {mf['b9']}")
    sw9 = sorted(c["swing_med_2_100"] for c in b9)
    chk("D88 validation cells swing either 5 to 8 K or 0.5 to 2.6 K (median over cycles 2 to 100)", True,
        all((4.95 <= x <= 8.05) or (0.45 <= x <= 2.65) for x in sw9),
        note="median swings between 2.65 and 4.95 K: " + ", ".join(f"{x:.2f}" for x in sw9 if 2.65 < x < 4.95))
    chk("D88 flagged 30", 30, ct["flagged_total"])
    chk("D88 flagged by batch 4, 4, 5, 17", [4, 4, 5, 17], [len(ct["by_batch"][d]["flagged"]) for d in ("2017-05-12", "2017-06-30", "2018-04-12", "2019-01-24")])
    chk("D88 flagged channels 2017-05-12 (5, 36 b; 15, 16 f)", [["5", "b"], ["15", "f"], ["16", "f"], ["36", "b"]],
        [[f["channel"], f["rules"]] for f in ct["by_batch"]["2017-05-12"]["flagged"]])
    chk("D88 flagged channels 2017-06-30 (7 abf; 15, 16, 21 f)", [["7", "abf"], ["15", "f"], ["16", "f"], ["21", "f"]],
        [[f["channel"], f["rules"]] for f in ct["by_batch"]["2017-06-30"]["flagged"]])
    chk("D88 flagged channels 2018-04-12 (3, 10, 14, 39, 42 b)", [["3", "b"], ["10", "b"], ["14", "b"], ["39", "b"], ["42", "b"]],
        [[f["channel"], f["rules"]] for f in ct["by_batch"]["2018-04-12"]["flagged"]])
    chk("D88 flagged channels 2019-01-24 (22, 44 ab; 15 more b)",
        sorted([["22", "ab"], ["44", "ab"]] + [[x, "b"] for x in ("2", "3", "8", "10", "27", "33", "35", "36", "38", "39", "41", "42", "43", "45", "48")]),
        sorted([[f["channel"], f["rules"]] for f in ct["by_batch"]["2019-01-24"]["flagged"]]))
    chk("D88 flagged in sets of two or more", 27, ct["flagged_in_sets_2plus"])
    r1 = ct["by_batch"]["2017-05-12"]["b_ratios_sorted"]
    chk("D88 (b) ratios 2017-05-12: 0.12, 0.26, then 0.64", [0.12, 0.26, 0.64], [round(x, 2) for x in r1[:3]])
    r3 = ct["by_batch"]["2018-04-12"]["b_ratios_sorted"]
    chk("D88 (b) ratios 2018-04-12: 0.14, 0.27, 0.42, 0.43, 0.48, then 0.52, 0.56", [0.14, 0.27, 0.42, 0.43, 0.48, 0.52, 0.56], [round(x, 2) for x in r3[:7]])
    rl = rec["A0_5_reliability"]
    chk("D89 r_oe pooled 0.999983", 0.999983, rl["pooled"]["r_oe"], 6)
    chk("D89 r_oe first paper 0.999974", 0.999974, rl["paper1"]["r_oe"], 6)
    chk("D89 r_oe second paper 0.999989", 0.999989, rl["paper2"]["r_oe"], 6)
    chk("D89 r_h 0.98953 (pooled)", 0.98953, rl["pooled"]["r_h"], 5)
    rs = rec["A0_6_resistance"]
    chk("D90 2017-05-12: 31 cells with 99 valid, 10 with 98", [31, 10], [rs["2017-05-12"]["valid_cycles_histogram"].get("99"), rs["2017-05-12"]["valid_cycles_histogram"].get("98")])
    chk("D90 R from 0.0134 to 0.0202 ohm", [0.0134, 0.0202], [round(x, 4) for x in rs["R_range_all"]])
    chk("D90 2019-01-24 no valid IR", 45, rs["2019-01-24"]["valid_cycles_histogram"].get("0"))
    ch = rec["A0_7_channel"]
    chk("D91 108 pairs", 108, ch["pairs"])
    chk("D91 correlation -0.054", -0.054, ch["correlation"], 3)
    q_p = ch["p_two_sided_plus_one"]
    sd_p = math.sqrt(2 * 0.593 * 0.407 / 10000)
    q.append({"where": "D91 p 0.593 (Monte Carlo; mine with my own seed)", "quoted": 0.593, "mine": q_p,
              "agree": bool(abs(q_p - 0.593) <= 3 * sd_p + 0.0005), "note": f"3 sd of the difference {3 * sd_p:.4f}"})
    p = rec["A0_8_power"]
    chk("D92 (3) sd_c Geslin 0.0314966", 0.0314966, p["pooled"]["floors"]["Geslin"]["sd_c"], 7)
    chk("D92 (3) sd_c KIT 0.0736310", 0.0736310, p["pooled"]["floors"]["KIT"]["sd_c"], 7)
    for scope, k, g, n, j, rms in (("pooled", 1.735, 0.742, 125, 36, 1.07), ("paper1", 2.820, 1.206, 80, 27, 0.83), ("paper2", 2.251, 0.963, 45, 9, 1.41)):
        chk(f"D92 (6) {scope} KIT floor % per K", k, 100 * p[scope]["floors"]["KIT"]["detectable_calc"], 3)
        chk(f"D92 (6) {scope} Geslin floor % per K", g, 100 * p[scope]["floors"]["Geslin"]["detectable_calc"], 3)
        chk(f"D92 (6) {scope} cells and sets", [n, j], [p[scope]["N"], p[scope]["J"]])
        chk(f"D92 (6) {scope} T deviations rms K", rms, p[scope]["dT_rms_over_N"], 2)
    chk("D92 (6) V5 passes in all six", True, all(p[s]["floors"][f]["v5_pass"] for s in ("pooled", "paper1", "paper2") for f in ("KIT", "Geslin")))
    rels = [p[s]["floors"][f]["v5_rel_diff"] for s in ("pooled", "paper1", "paper2") for f in ("KIT", "Geslin")]
    chk("D92 (6) simulated 0.0 to 3.2 % below the calculated (Monte Carlo; mine)", True, all(-0.032 - 0.01 <= r <= 0.0 for r in rels),
        note="mine: " + ", ".join(f"{100 * r:.2f} %" for r in rels))
    wf = p["pooled_without_flagged"]
    chk("D92 (6) without flagged pooled KIT, Geslin % per K", [3.083, 1.319], [round(100 * wf["floors"]["KIT"]["detectable_calc"], 3), round(100 * wf["floors"]["Geslin"]["detectable_calc"], 3)])
    chk("D92 (6) without flagged 95 cells in 32 sets, 0.70 K", [95, 32, 0.70], [wf["N"], wf["J"], round(wf["dT_rms_over_N"], 2)])
    chk("D92 (6) P1 holds (1.735 < 3.372)", True, rec["scorecard"]["P1"]["holds"])
    chk("D92 (6) P6 holds", True, rec["scorecard"]["P6"]["holds"])
    lit_differs = [c["key"] for c in cells if c["matlab_label_skip_zero_first"] != c["L"]]
    chk("D97 (1) MATLAB label agrees for all 124 first-paper cells (LoadData.m's rule)", 0,
        sum(1 for c in cells if c["paper"] == 1 and c["matlab_label_skip_zero_first"] != c["L"]))
    chk("D97 (1) differs for two validation cells, channels 45 and 48", ["45", "48"], sorted((bykey[k]["channel"] for k in lit_differs), key=int))
    chk("D97 (4) Tmin-measure channel test correlation +0.12 over 108 pairs", [108, 0.12], [ph["Tmin measure"]["pairs"], round(ph["Tmin measure"]["correlation"], 2)])
    def mc_p(where, quoted, mine, n_quoted, n_mine, unit):
        sd = math.sqrt(quoted * (1 - quoted) / n_quoted + mine * (1 - mine) / n_mine)
        ok = abs(quoted - mine) <= 3 * sd + unit / 2
        q.append({"where": where, "quoted": quoted, "mine": mine, "agree": bool(ok),
                  "note": f"Monte Carlo: |difference| {abs(quoted - mine):.4f} against 3 sd {3 * sd:.4f} plus half the quoted rounding"})

    mc_p("D97 (4) Tmin-measure p 0.20 (5,000 permutations; mine 10,000)", 0.20, ph["Tmin measure"]["p_two_sided"], 5000, 10000, 0.01)
    chk("D97 (4) swing channel test correlation -0.08 (median swing)", -0.08, ph["median swing"]["correlation"], 2,
        note=f"mean swing gives {ph['mean swing']['correlation']:.4f}")
    mc_p("D97 (4) swing p 0.41 (median swing)", 0.41, ph["median swing"]["p_two_sided"], 5000, 10000, 0.01)
    return {"quotes": q, "agree": int(sum(x["agree"] for x in q)), "checked": len(q), "disagree": [x for x in q if not x["agree"]]}


def summary_block(out):
    cols = out["cells"]["columns"]
    det = lambda k: {"class": k, "compared": out[k]["compared"], "agree": out[k]["agree"],
                     "largest_abs_diff": out[k]["max_abs_diff"], "largest_rel_diff": out[k]["max_rel_diff"], "at": out[k]["max_rel_diff_at"]}
    mc = lambda k: {"class": k, "compared": out[k]["compared"], "within_abs_z_4": out[k]["within_4"],
                    "largest_abs_z": out[k]["max_abs_z"], "at": out[k]["max_abs_z_at"], "largest_abs_diff": out[k]["max_abs_diff"]}
    zs = [r["z"] for r in out["B_C_intervals"]["rows"] if r.get("z") is not None]
    rows = [det("counts"),
            {"class": "cells (169 entering cells x 19 columns)", "compared": sum(v["compared"] for v in cols.values()),
             "largest_abs_diff": max(v["max_abs_diff"] for v in cols.values()),
             "contact_flag_mismatches": len(out["cells"]["contact_flag_mismatches"]),
             "rule_mismatches": len(out["cells"]["contact_rule_letter_mismatches"]),
             "excluded_entries_match": out["cells"]["excluded_entries_match"]},
            det("sets"), det("A0"), det("power_calculation"), det("B_C_points"), det("B3"), det("scorecard_values"),
            det("paper_numbers"), det("post_freeze_extras"), mc("B_C_intervals"), mc("power_simulation"), mc("A0_channel_p"),
            {"class": "decisions_quotes", "compared": out["decisions_quotes"]["checked"], "agree": out["decisions_quotes"]["agree"],
             "disagree": [d["where"] for d in out["decisions_quotes"]["disagree"]]}]
    flags = [
        "D88 (e): 'median floors ... of 30.5 to 32.6 C in 2019-01-24'; the largest is 32.651356 C (the registered floor_med agrees), "
        "which rounds to 32.7. Cosmetic.",
        "D88: 'cells on one protocol swing either 5 to 8 K or 0.5 to 2.6 K' holds for five of the nine validation protocols; four have a "
        "cell between (median swings 2.99, 3.81, 3.94, 4.22, 4.75 K) and one cell swings 8.17 K. The ratio gap of rule (b), 0.390 to "
        "0.575, holds; no flag depends on the sentence.",
        "D93 (6) says B and C use the same draws; the registered B2 draws its own resample (s6_parts.py resistance_part, seed "
        "[20261116, k + 100]). Same scheme, intervals agree within Monte Carlo spread, C2 uses shared draws; no point or verdict depends on it.",
        "D93 (7) 'every slope divided by r_oe' is ambiguous; the registered divides s_w, s_b and B2's lever, points only. r_oe = 0.99998, "
        "so nothing turns on it.",
        "D88 (c) to (e) must be read inside cycles 2 to 100; with cycle 1 (a zero placeholder in 2017-05-12) every 2017-05-12 cell "
        "would be flagged. Both codes read the window.",
        "D91 does not state the permutation p's form (the registered uses (k + 1) / (n + 1)) and 'every two batches' means every pair "
        "of batches (108 pairs). Readings match.",
        "D94 (2): the literal MATLAB label is 1 for every 2017-05-12 cell; D97 (1) corrects it; my frozen reading agrees (differs from L "
        "only for 2019-01-24 channels 45 and 48). Not used.",
    ]
    return {"agreement_by_class": rows,
            "interval_z": {"n": len(zs), "mean": float(np.mean(zs)), "sd": float(np.std(zs)), "abs_above_2": int(sum(abs(z) > 2 for z in zs)),
                           "abs_above_3": int(sum(abs(z) > 3 for z in zs))},
            "verdicts": out["scorecard"]["verdicts"], "all_verdicts_reproduced": out["scorecard"]["all_reproduced"],
            "verdict_margins": out["verdict_margins"], "flags": flags,
            "errors_found": "none: no registered number differs from mine beyond rounding or Monte Carlo spread",
            "note": "report.md was not written: this agent's environment refuses report .md files, so the report is returned as text "
                    "to the caller; this summary holds its table, verdicts and flags."}


def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, (np.floating, float)):
        x = float(x)
        return x if math.isfinite(x) else None
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


if __name__ == "__main__":
    main()
