#!/usr/bin/env python3
"""Independent Stage 5 check: Monte Carlo spread of Part C's simulated quantities (the checker's robustness study;
not a registered quantity). For every C1 case, 20 replicate 10,000-dataset simulations under the registered shape
(own seeds [BASE_SEED, 700 + case, floor, r]): the spread of the simulated detectable |s| and the share of replicates
in which V4 (within 5 % of the calculation) passes. For the two C2 cases D75 names (KIT warm and cold, KIT floor),
5 replicate bisections of the simulated pairs needed, and the spread of the simulated power at the band median.

Output: out/C_mc_spread.json.
Usage: PYTHONDONTWRITEBYTECODE=1 python3 verification/stage5/independent/mc_spread_c.py
"""
import json, math, os, sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ind_common as C   # noqa: E402
import partc as P        # noqa: E402

R1, R2 = 20, 5


def main():
    rows, _ = C.build_pairs()
    sh = C.shapes()
    med = C.band_medians()
    Cj = json.load(open(os.path.join(C.OUT, "C.json")))
    sel = {"KIT": C.select(rows, "KIT"), "Sandia": C.select(rows, "Sandia"), "pooled": C.select(rows),
           "KIT warm": C.select(rows, "KIT", "warm"), "KIT cold": C.select(rows, "KIT", "cold"),
           "Sandia warm": C.select(rows, "Sandia", "warm"), "Sandia cold": C.select(rows, "Sandia", "cold"),
           "pooled warm": C.select(rows, None, "warm"), "pooled cold": C.select(rows, None, "cold"),
           "KIT all 142": C.select(rows, "KIT", scored=False)}
    out = {"C1": [], "C2": []}
    for ci, r in enumerate(Cj["C1"]):
        dT = np.array([x["dT_K"] for x in sel[r["case"]]])
        shp = sh[r["floor"]]
        fi = 0 if r["floor"] == "Geslin" else 1
        calc = r["detectable_s_calc_pct_per_K"] / 100.0
        vals, pw_med = [], []
        for rep in range(R1):
            fd = P.FixedDesign(dT, shp, "lognormal", C.rng_for(700 + ci, fi, rep))
            vals.append(fd.detectable_bisection(calc))
            b = "warm" if r["case"].endswith("warm") else ("cold" if r["case"].endswith("cold") else None)
            if b:
                pw_med.append(fd.power(abs(med[b]["s"])))
        vals = np.array(vals)
        rel = (vals - calc) / calc
        rec = {"case": r["case"], "floor": r["floor"], "calc_pct_per_K": 100 * calc,
               "registered_run_sim_pct_per_K": r["detectable_s_sim_lognormal_pct_per_K"],
               "replicates": R1, "sim_mean_pct_per_K": 100 * float(vals.mean()), "sim_sd_pct_per_K": 100 * float(vals.std(ddof=1)),
               "rel_diff_mean": float(rel.mean()), "rel_diff_sd": float(rel.std(ddof=1)),
               "v4_pass_share": float(np.mean(np.abs(rel) <= 0.05))}
        if pw_med:
            rec["power_at_band_median_mean"] = float(np.mean(pw_med))
            rec["power_at_band_median_sd"] = float(np.std(pw_med, ddof=1))
        out["C1"].append(rec)
        print(json.dumps(rec), flush=True)
    for b in ("warm", "cold"):
        r = [x for x in Cj["C2"] if x["band"] == b and x["dataset"] == "KIT" and x["floor"] == "KIT"][0]
        dT = np.array([x["dT_K"] for x in sel[f"KIT {b}"]])
        ns = []
        for rep in range(R2):
            n, _ = P.n_simulated(abs(med[b]["s"]), dT, sh["KIT"], "lognormal", r["n_required_calc"], (800 + (b == "cold"), rep))
            ns.append(n)
            print(b, rep, n, flush=True)
        out["C2"].append({"band": b, "dataset": "KIT", "floor": "KIT", "n_calc": r["n_required_calc"],
                          "registered_run_n_sim": r["n_required_sim_lognormal"], "replicate_n_sim": ns,
                          "mean": float(np.mean(ns)), "sd": float(np.std(ns, ddof=1))})
    C.write_json("C_mc_spread.json", out)


if __name__ == "__main__":
    main()
