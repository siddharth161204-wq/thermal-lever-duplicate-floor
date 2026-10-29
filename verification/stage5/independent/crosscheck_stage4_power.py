#!/usr/bin/env python3
"""Independent Stage 5 check: after the checker's own power code (partc.py) was written and run, the C1 calculated
detectable |s| and the C2 calculated pairs are recomputed with Stage 4's archived functions experiment.power_t and
experiment.n_required (src/stage4/experiment.py at the base commit; Stage 4 code, not Stage 5), as the task allows.
The detectable |s| uses Brent's method on |s| with power_t(n, |s| dT_rms, sd(e)), to 1e-12 relative.

Output: out/crosscheck_stage4_power.json.
Usage: PYTHONDONTWRITEBYTECODE=1 python3 verification/stage5/independent/crosscheck_stage4_power.py
"""
import importlib.util, json, math, os, sys

from scipy import optimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ind_common as C   # noqa: E402


def load_stage4_experiment():
    fn = C.path("src/stage4/experiment.py")
    spec = importlib.util.spec_from_file_location("stage4_experiment_archived", fn)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod, C.sha256_file(fn)


def main():
    E, h = load_stage4_experiment()
    Cj = json.load(open(os.path.join(C.OUT, "C.json")))
    out = {"stage4_experiment_py_sha256": h, "sha256_at_v1.3.0": C.sha256_at("v1.3.0", "src/stage4/experiment.py"),
           "C1": [], "C2": []}
    worst = 0.0
    for r in Cj["C1"]:
        n, rms, sd = r["n_pairs"], r["dT_rms_K"], r["sd_e"]
        f = lambda s: E.power_t(n, s * rms, sd) - 0.80
        s4 = optimize.brentq(f, 1e-9, 1.0, xtol=1e-16, rtol=1e-13, maxiter=500)
        mine = r["detectable_s_calc_pct_per_K"] / 100.0
        rel = abs(s4 - mine) / mine
        worst = max(worst, rel)
        out["C1"].append({"case": r["case"], "floor": r["floor"], "stage4_power_t_pct_per_K": 100 * s4,
                          "checker_pct_per_K": 100 * mine, "rel_diff": rel})
    out["C1_max_rel_diff"] = worst
    agree = True
    for r in Cj["C2"]:
        eff = r["band_median_abs_s_per_K"] * r["dT_rms_K"]
        n4 = E.n_required(eff, r["sd_e"])
        n4s = E.n_required(eff, r["sd_e"] * math.sqrt(r["F_shared_cell"]))
        p4 = E.power_t(r["n_pairs_available"], eff, r["sd_e"])
        ok = n4 == r["n_required_calc"] and n4s == r["n_required_shared_cell"]
        agree &= ok
        out["C2"].append({"band": r["band"], "dataset": r["dataset"], "floor": r["floor"], "stage4_n_required": n4,
                          "checker_n_required": r["n_required_calc"], "stage4_n_shared": n4s,
                          "checker_n_shared": r["n_required_shared_cell"], "stage4_power_t": p4,
                          "checker_power_calc": r["power_calc"], "power_abs_diff": abs(p4 - r["power_calc"]), "n_equal": ok})
    out["C2_all_n_equal"] = agree
    out["C2_max_power_abs_diff"] = max(x["power_abs_diff"] for x in out["C2"])
    C.write_json("crosscheck_stage4_power.json", out)
    print(json.dumps({k: v for k, v in out.items() if k not in ("C1", "C2")}, indent=1))


if __name__ == "__main__":
    main()
