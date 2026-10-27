#!/usr/bin/env python3
"""POST-FREEZE: the checker's frozen outputs (commit fd5d1e7) against the pipeline's registered Stage 4 results.

The pipeline's files are read from git: the registered results and checks from commit 6f658f1 (unchanged at the
stage4 tip), the post-hoc files from a2b3d6e. For every quantity: the number of values compared, the largest
absolute and relative difference, and every disagreement in a verdict, flag, integer or label.
Output: out_compare/compare.json (and a short printout)
"""
import csv, io, json, math, os, subprocess, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = os.path.join(HERE, "out")
CMP = os.path.join(HERE, "out_compare")
REG, POST = "6f658f1", "a2b3d6e"


def git(rev, path):
    return subprocess.run(["git", "show", f"{rev}:{path}"], cwd=ROOT, check=True, capture_output=True).stdout.decode()


def pcsv(path, rev=REG):
    return list(csv.DictReader(io.StringIO(git(rev, path))))


def pjson(path, rev=REG):
    return json.loads(git(rev, path))


def mcsv(name):
    with open(os.path.join(OUT, name)) as f:
        return list(csv.DictReader(f))


def mjson(name):
    with open(os.path.join(OUT, name)) as f:
        return json.load(f)


def fl(x):
    """float, accepting the 'np.float64(...)' wrapper that the checker's frozen divergence.csv wrote for the per-cell
    ln-offsets (a formatting slip in that file; the digits inside are the full-precision values)."""
    if isinstance(x, str) and x.startswith("np.float64(") and x.endswith(")"):
        x = x[len("np.float64("):-1]
    return float(x)


class Acc:
    def __init__(self):
        self.q = {}
        self.mismatch = []

    def num(self, qty, a, b, where=""):
        a, b = fl(a), fl(b)
        d = self.q.setdefault(qty, {"n": 0, "max_abs": 0.0, "max_rel": 0.0, "where_max_abs": "", "nan_both": 0})
        if math.isnan(a) and math.isnan(b):
            d["nan_both"] += 1
            return
        d["n"] += 1
        ad = abs(a - b)
        rel = ad / abs(b) if b != 0 else (0.0 if a == 0 else float("inf"))
        if ad > d["max_abs"]:
            d["max_abs"], d["where_max_abs"] = ad, where
        d["max_rel"] = max(d["max_rel"], rel)

    def same(self, qty, a, b, where=""):
        d = self.q.setdefault(qty, {"n": 0, "disagree": 0})
        d["n"] += 1
        if a != b:
            d["disagree"] += 1
            self.mismatch.append({"quantity": qty, "where": where, "checker": a, "pipeline": b})


def b(x):
    return str(x).strip() in ("True", "true", "1")


def main():
    os.makedirs(CMP, exist_ok=True)
    A = Acc()
    # ---------------- resistance ----------------
    pr = pjson("results/stage4_resistance.json")
    mr = mjson("resistance.json")
    for T in ("35", "25", "10"):
        A.num("R_used_ohm", mr["ambients"][T]["R_registered_ohm"] if T != "35" else 0.01871279570954531, pr["R_used_ohm"][T], T)
        A.num("R_rerun_ohm (procedure at each ambient)", mr["ambients"][T]["R_registered_ohm"], pr["runs"][T]["registered"]["R_DC_ohm"], T)
        for x, y in zip(mr["ambients"][T]["pulses"], pr["runs"][T]["pulses"]):
            A.num("pulse R_ohm", x["R_ohm"], y["R_ohm"], f"{T} {x['c_rate']} {x['pulse_s']}")
        for x, y in zip(mr["ambients"][T]["cycle_heat"], pr["runs"][T]["cycle_heat"]):
            assert x["c_rate"] == y["c_rate"] and x["direction"] == y["direction"]
            A.same("cycle_heat failed", x["status"] != "ok", "failed" in y, f"{T} {x['c_rate']} {x['direction']}")
            if x["status"] == "ok" and "R_eff_irreversible_ohm" in y:
                A.num("cycle_heat R_eff_ohm", x["R_eff_irreversible_ohm"], y["R_eff_irreversible_ohm"], f"{T} {x['c_rate']} {x['direction']}")
                A.num("cycle_heat charge_Ah", x["charge_Ah"], y["charge_Ah"], f"{T} {x['c_rate']} {x['direction']}")
    A.num("ratio 25/35", mr["used"]["ratio_25_over_35_archived"], pr["ratio_to_35C"]["25"])
    A.num("ratio 10/35", mr["used"]["ratio_10_over_35_archived"], pr["ratio_to_35C"]["10"])
    A.same("35 C rerun identical to archived", mr["used"]["rerun_minus_archived_ohm"] == 0.0, pr["rerun_35C_identical_to_archived"])
    # ---------------- module ----------------
    pm = {(float(r["ambient_C"]), float(r["crate"]), float(r["h"])): r for r in pcsv("results/stage4_module.csv")}
    mm = {(float(r["T_amb_C"]), float(r["crate"]), float(r["h"])): r for r in mcsv("module.csv")}
    assert set(pm) == set(mm) and len(mm) == 36
    for k in sorted(mm):
        x, y = mm[k], pm[k]
        w = f"{k[0]:g} C {k[1]:g}C h={k[2]:g}"
        for qa, qb in (("R_ohm", "R_ohm"), ("q_cell_W", "q_cell_W"), ("T_max_C", "T_max_C"), ("T_min_C", "T_min_C"),
                       ("T_mean_C", "T_mean_C"), ("spread_K", "spread_K"), ("hot_minus_mean_K", "hot_minus_mean_K")):
            A.num("module " + qb + " (own solver)", x[qa], y[qb], w)
        for i in range(8):
            A.num("module T_cell_C (own solver)", x[f"T_cell{i}_C"], y[f"T_cell{i}_C"], f"{w} cell {i}")
            A.num("module T_cell_C (module_steady cross-check)", x[f"ms_T_cell{i}_C"], y[f"T_cell{i}_C"], f"{w} cell {i}")
            A.num("module offset_K (own solver)", x[f"offset{i}_K"], y[f"offset_cell{i}_K"], f"{w} cell {i}")
        A.num("module spread_K (module_steady cross-check)", x["ms_spread_K"], y["spread_K"], w)
        A.same("module i_hot", int(x["i_hot"]), int(y["i_hot"]), w)
        A.same("module i_cold", int(x["i_cold"]), int(y["i_cold"]), w)
        A.same("module newton_iterations (module_steady)", int(x["ms_newton_iterations"]), int(y["newton_iterations"]), w)
        A.num("module balance_W (module_steady)", x["ms_balance_W"], y["balance_W"], w)
        A.same("module converged", b(x["ms_converged_D58"]), b(y["converged"]), w)
    # cells inside the band: the pipeline reads the set points (D66); the checker had written the logged secant span
    sp = {"warm": (25.0, 40.0), "cold": (10.0, 25.0)}
    inside = {}
    for k in sorted(mm):
        band = "warm" if k[0] in (35.0, 25.0) else "cold"
        ins = float(mm[k]["T_min_C"]) >= sp[band][0] and float(mm[k]["T_max_C"]) <= sp[band][1]
        inside[k] = ins
        A.same("cells inside band set points (checker recomputed post-freeze with KIT set points)", ins, b(pm[k]["cells_inside_band_setpoints"]),
               f"{k[0]:g} C {k[1]:g}C h={k[2]:g}")
    # ---------------- resolved cell ----------------
    prc = {(float(r["ambient_C"]), float(r["crate"])): r for r in pcsv("results/stage4_resolved_cell.csv")}
    for r in mcsv("resolved_cell.csv"):
        k = (float(r["T_amb_C"]), float(r["crate"]))
        y = prc[k]
        w = f"{k[0]:g} C {k[1]:g}C"
        for qa in ("q_W", "T_max_C", "T_min_C", "T_vol_mean_C", "T_surface_mean_C", "dT_cell_K", "dT_hot_minus_mean_K", "module_spread_K"):
            A.num("resolved " + qa, r[qa], y[qa], w)
        A.num("resolved within_over_module", r["ratio_cell_over_module"], y["within_over_module"], w)
        A.same("resolved newton_iterations", int(r["newton_iterations"]), int(y["newton_iterations"]), w)
        A.same("resolved cg_all_converged", b(r["cg_all_converged"]), b(y["cg_all_converged"]), w)
    # ---------------- divergence ----------------
    pd_ = pcsv("results/stage4_divergence.csv")
    plev = {(float(r["ambient_C"]), float(r["crate"]), float(r["h"]), r["lever"]): r for r in pd_ if r["row_kind"] == "lever"}
    pmed = {(float(r["ambient_C"]), float(r["crate"]), float(r["h"])): r for r in pd_ if r["row_kind"] == "band median"}
    pp95 = {(float(r["ambient_C"]), float(r["crate"]), float(r["h"])): r for r in pd_ if r["row_kind"] == "band 95th percentile"}
    mlev = {(float(r["T_amb_C"]), float(r["crate"]), float(r["h"]), r["lever"]): r for r in mcsv("divergence.csv")}
    assert set(plev) == set(mlev), (len(plev), len(mlev))
    for k in sorted(mlev):
        x, y = mlev[k], plev[k]
        w = f"{k[0]:g} C {k[1]:g}C h={k[2]:g} {k[3]}"
        A.num("divergence |s| x spread", x["implied_divergence"], y["divergence"], w)
        A.num("divergence interval low", x["implied_divergence_lo"], y["divergence_lo"], w)
        A.num("divergence interval high", x["implied_divergence_hi"], y["divergence_hi"], w)
        A.num("lever s_per_K", x["s_per_K"], y["s_per_K"], w)
        for i in range(8):
            A.num("divergence ln-offset per cell", x[f"ln_offset_cell{i}"], y[f"offset_lnL_cell{i}"], w)
    mc = {(float(r["T_amb_C"]), float(r["crate"]), float(r["h"])): r for r in mcsv("divergence_cases.csv")}
    pmap = pjson("results/stage4_map.json")
    pcase = {(float(r["ambient_C"]), float(r["crate"]), float(r["h"])): r for r in pmap["cases"]}
    for k in sorted(mc):
        x = mc[k]
        w = f"{k[0]:g} C {k[1]:g}C h={k[2]:g}"
        A.num("band median", x["median"], pmed[k]["divergence"], w)
        A.num("band 95th percentile", x["p95"], pp95[k]["divergence"], w)
        A.num("band median (map.json)", x["median"], pcase[k]["median"], w)
        A.num("band max", x["max"], pcase[k]["max"], w)
        A.same("median lever", x["median_lever"], pcase[k]["median_lever"], w)
        for fk in ("Geslin", "KIT"):
            A.same(f"offers {fk}", b(x[f"offers_{fk}"]), bool(pcase[k][f"offers_{fk}"]), w)
            A.same(f"detectable {fk}", b(x[f"detectable_{fk}"]), bool(pcase[k][f"detectable_{fk}"]), w)
            A.same(f"levers above {fk} mean", int(x[f"levers_above_{fk}_mean"]), int(pcase[k][f"levers_above_mean_{fk}"]), w)
            A.same(f"levers above {fk} p95", int(x[f"levers_above_{fk}_p95"]), int(pcase[k][f"levers_above_p95_{fk}"]), w)
    mmap = mjson("map.json")["map"]
    for key, v in mmap.items():
        pv = pmap["map"][v["floor"]][f"{v['ambient_C']:g} C"][f"h = {v['h']:g}"]
        A.same("map smallest rate offering", v["smallest_rate_offers"], pv["smallest_rate_offering_a_share"], key)
        A.same("map smallest rate detectable", v["smallest_rate_detectable"], pv["smallest_rate_offering_a_detectable_share"], key)
    # ---------------- experiment ----------------
    pe = pcsv("results/stage4_experiment.csv")
    me = mcsv("experiment.csv")
    sim_ln = {int(r["case_index"]): r for r in mcsv("simulation_lognormal_all.csv")}
    sim_no = {int(r["case_index"]): r for r in mcsv("simulation_normal_all.csv")}

    def mkey_m(r):
        ch = "" if r["group"] == "C1" else ("warm" if "warm" in r["group"] else "cold")
        return (r["group"][:2], ch, r["lever"].split(" (")[0] if r["lever"].startswith("KIT 10 to 25 C median") else r["lever"],
                float(r["dT_K"]), r["floor"])

    mreg = {}
    for i, r in enumerate(me):
        if r["shape"] == "lognormal":
            mreg.setdefault(mkey_m(r), {})["ln"] = (i, r)
        else:
            mreg.setdefault(mkey_m(r), {})["no"] = (i, r)
    v4_disagree = []
    unmatched = []
    sim_same = {}
    p_same = os.path.join(OUT, "postfreeze_simulation_normal_same_sd.csv")
    if os.path.exists(p_same):
        sim_same = {int(r["case_index"]): r for r in mcsv("postfreeze_simulation_normal_same_sd.csv")}
    for r in pe:
        lev = r["lever"]
        pch = "" if not r["chamber"] else ("warm" if r["chamber"].startswith("warm") else "cold")
        key = None
        for k in mreg:
            if k[0] == r["part"] and k[1] == pch and float(r["dT_K"]) == k[3] and r["floor"] == k[4]:
                if k[2] == lev or (lev.startswith("KIT 10 to 25 C") and k[2].startswith("KIT 10 to 25 C median")) or \
                   (r["lever_name"].startswith("central cold KIT") and k[2].startswith("KIT 10 to 25 C median")):
                    key = k
                    break
        if key is None:
            unmatched.append((r["part"], r["chamber"], lev, r["dT_K"], r["floor"]))
            continue
        (iln, xl), (ino, xn) = mreg[key]["ln"], mreg[key]["no"]
        w = f"{r['part']} {r['chamber']} {lev} dT={r['dT_K']} {r['floor']}"
        A.num("experiment s_per_K", xl["s_per_K"], r["s_per_K"], w)
        A.num("experiment sd_e registered (log-normal)", xl["sd_e"], r["sd_e_registered"], w)
        A.num("experiment sd_e normal shape", xn["sd_e"], r["sd_e_normal_shape"], w)
        A.same("experiment n calculated (log-normal)", int(xl["n"]), int(float(r["n_pairs"])), w)
        A.num("experiment power at n", xl["power_at_n"], r["power_at_n"], w)
        A.same("experiment n calculated (normal shape)", int(xn["n"]), int(float(r["n_pairs_normal_shape"])), w)
        A.same("experiment central flag", b(xl["central"]), b(r["central"]), w)
        A.same("experiment conditional flag", b(xl["conditional"]), b(r["conditional"]), w)
        sl, sn = sim_ln.get(iln), sim_no.get(ino)
        if sl:
            A.num("V4 simulated n (log-normal), Monte Carlo", sl["n_sim"], r["n_sim_lognormal"], w)
            if b(sl["V4_pass_5pct"]) != b(r["v4_pass"]):
                v4_disagree.append({"case": w, "n_calc": int(xl["n"]), "checker_n_sim": int(sl["n_sim"]),
                                    "pipeline_n_sim": int(float(r["n_sim_lognormal"])), "checker_pass": b(sl["V4_pass_5pct"]),
                                    "pipeline_pass": b(r["v4_pass"])})
        if sn:
            A.num("checker's normal-shape simulated n vs pipeline's calculated normal-shape n", sn["n_sim"], r["n_pairs_normal_shape"], w)
        ss = sim_same.get(iln)
        if ss:
            A.num("V4 normal check, same sd(e) as registered (post-freeze), Monte Carlo", ss["n_sim"], r["n_sim_normal_same_sd"], w)
            A.same("V4 normal check pass (post-freeze reading)", b(ss["V4_pass_5pct"]), b(r["normal_check_pass"]), w)
        A.num("sim power at n_sim (log-normal), Monte Carlo", sl["sim_power_at_n_sim"] if sl else float("nan"),
              r["sim_power_at_n_sim_lognormal"] if sl else float("nan"), w)
    # the figure curves on the checker's 0.5 K grid
    pcur = pcsv("results/stage4_experiment_curves.csv")
    mg = mcsv("experiment_grid.csv")
    for r in mg:
        ch = "warm" if r["chamber"].startswith("warm") else "cold"
        y = [p for p in pcur if p["chamber"] == ch and p["floor"] == r["floor"] and abs(float(p["dT_K"]) - float(r["dT_K"])) < 1e-9]
        if not y:
            continue
        col = "n_pairs" if r["shape"] == "lognormal" else "n_pairs_normal_shape"
        A.same(f"figure curve n ({r['shape']})", int(r["n"]), int(float(y[0][col])), f"{ch} {r['floor']} dT={r['dT_K']}")
    # ---------------- scorecard ----------------
    ps = pjson("results/stage4_scorecard.json")["_summary"]
    ms = mjson("scorecard.json")
    for k in ("P1", "P2", "P3", "P4", "P5", "P6", "P7"):
        A.same("verdict " + k, ms[k]["held"], ps[k], k)
    # ---------------- the pipeline's checks ----------------
    checks = {}
    for v in ("v1", "v2", "v3", "v4"):
        d = pjson(f"verification/stage4/{v}.json")
        checks[v] = {k: d[k] for k in d if k in ("pass", "pass_all", "failures", "max_abs_rel_diff", "max_abs_diff_K",
                                                 "all_bitwise_identical", "normal_check_pass_all", "normal_check_max_abs_rel_diff")}
    my_ln = mcsv("simulation_lognormal_all.csv")
    my_no = mcsv("simulation_normal_all.csv")
    checks["checker_V4_lognormal_pass_count"] = sum(b(r["V4_pass_5pct"]) for r in my_ln)
    checks["checker_V4_normal_pass_count"] = sum(b(r["V4_pass_5pct"]) for r in my_no)
    checks["pipeline_V4_lognormal_pass_count"] = sum(b(r["v4_pass"]) for r in pe)
    checks["pipeline_V4_normal_pass_count"] = sum(b(r["normal_check_pass"]) for r in pe)
    # ---------------- post-hoc D66 (cold bracket) recomputed with the checker's solver ----------------
    sys.path.insert(0, HERE)
    import module as MOD
    import levers as LV
    warm, cold, floors = LV.load()
    med_cold = float(np.median([abs(lv["s"]) for lv in cold]))
    med_warm = float(np.median([abs(lv["s"]) for lv in warm]))
    ph = pjson("results/stage4_posthoc.json", POST)
    for c in ph["cold_bracket"]["cases"]:
        Ta = c["ambient_C"] + 273.15
        q = ({"C/4": 0.25, "C/2": 0.5, "1C": 1.0, "2C": 2.0}[c["rate"]] * 5.0) ** 2 * c["R_ohm"]
        T = MOD.solve_symmetric(q, c["h"], Ta)
        sp_ = float(T.max() - T.min())
        A.num("post-hoc D66 spread_K (checker's solver)", sp_, c["spread_K"], f"{c['ambient_C']:g} C R({c['resistance_of_C']:g}) {c['rate']} h={c['h']:g}")
        mband = med_cold if c["ambient_C"] == 10.0 else med_warm
        A.num("post-hoc D66 median divergence", mband * sp_, c["median_divergence"], f"{c['ambient_C']:g} C R({c['resistance_of_C']:g}) {c['rate']} h={c['h']:g}")
    out = {"pipeline_commits": {"registered": REG, "posthoc": POST}, "checker_freeze_commit": "fd5d1e7",
           "quantities": A.q, "mismatches": A.mismatch, "experiment_rows_unmatched": unmatched,
           "V4_pass_disagreements": v4_disagree, "checks": checks}
    with open(os.path.join(CMP, "compare.json"), "w") as f:
        json.dump(out, f, indent=1)
    for k, v in A.q.items():
        if "max_abs" in v:
            print(f"{k:70s} n={v['n']:5d} max_abs={v['max_abs']:.3e} max_rel={v['max_rel']:.3e}  {v['where_max_abs']}")
        else:
            print(f"{k:70s} n={v['n']:5d} disagree={v['disagree']}")
    print("mismatches:", len(A.mismatch))
    for m in A.mismatch[:60]:
        print("  ", m)
    print("unmatched experiment rows:", unmatched)
    print("V4 pass disagreements:", len(v4_disagree))
    for m in v4_disagree:
        print("  ", m)
    print(json.dumps(checks, indent=1))


if __name__ == "__main__":
    main()
