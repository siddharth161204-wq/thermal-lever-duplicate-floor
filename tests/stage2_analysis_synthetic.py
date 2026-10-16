#!/usr/bin/env python3
"""Synthetic test of src/stage2/analysis.py on a made-up per-cell table with the Sandia design (86 cells). The two
published laws are replaced by dummies here, so this test evaluates no law at the Sandia conditions.
Run: python3 tests/stage2_analysis_synthetic.py WORKDIR
"""
import json, os, sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src", "stage2"))
import analysis as AN  # noqa: E402
import laws as LW      # noqa: E402

LW.law_secants = lambda T1, T2, c, soh, c_charge=0.5: {"law_1_kucinskis2022": -999.0, "law_2_schimpe2018": -999.0}
LW.law_curves = lambda temps, c, soh, c_charge=0.5: []

DESIGN = {  # (T, window, c_dis): replicates, from the archive's cell list
    "LFP": {(15, "0-100", 1): 2, (15, "0-100", 2): 2, (25, "0-100", 0.5): 1, (25, "0-100", 1): 4, (25, "0-100", 2): 2,
            (25, "0-100", 3): 4, (25, "20-80", 0.5): 4, (25, "20-80", 3): 1, (25, "40-60", 0.5): 2, (25, "40-60", 3): 2,
            (35, "0-100", 1): 4, (35, "0-100", 2): 2},
    "NCA": {(15, "0-100", 1): 2, (15, "0-100", 2): 2, (25, "0-100", 0.5): 2, (25, "0-100", 1): 4, (25, "0-100", 2): 2,
            (25, "20-80", 0.5): 4, (25, "40-60", 0.5): 2, (35, "0-100", 1): 4, (35, "0-100", 2): 2},
    "NMC": {(15, "0-100", 1): 2, (15, "0-100", 2): 2, (25, "0-100", 0.5): 2, (25, "0-100", 1): 4, (25, "0-100", 2): 2,
            (25, "0-100", 3): 4, (25, "20-80", 0.5): 4, (25, "20-80", 3): 2, (25, "40-60", 0.5): 2, (25, "40-60", 3): 2,
            (35, "0-100", 1): 4, (35, "0-100", 2): 2}}
S_TRUE = {"LFP": 0.05, "NCA": 0.03, "NMC": 0.02}

if __name__ == "__main__":
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(1)
    rows = []
    for chem, des in DESIGN.items():
        for (T, win, c), n in des.items():
            lo, hi = map(int, win.split("-"))
            for k in range(n):
                letter = "abcd"[k]
                cid = f"SNL_18650_{chem}_{T}C_{win}_0.5/{c}C_{letter}"
                rise = 0.3 * c * c + rng.normal(0, 0.1 * c)
                L = 1000.0 * np.exp(-S_TRUE[chem] * (T - 25 + rise)) * np.exp(rng.normal(0, 0.01))
                rows.append({"cell_id": cid, "key": cid.replace("/", "-"), "chemistry": chem, "ambient_C": float(T),
                             "soc_min": float(lo), "soc_max": float(hi), "c_charge": 0.5, "c_discharge": float(c),
                             "replicate": letter, "rated_Ah": 1.0, "group": f"{chem}_{T}C_{win}_0.5/{c}C",
                             "L_efc": L, "resolution_efc": 2.0, "censored": False, "efc_last_check": L * 1.2,
                             "rise_K": rise, "rest_T_C": T + rng.normal(0, 0.5)})
    cells = pd.DataFrame(rows)
    assert len(cells) == 86, len(cells)
    p = os.path.join(out, "cells.csv")
    cells.to_csv(p, index=False)
    lev, flo, sa, sb, mod, sc = AN.main(p, os.path.join(ROOT, "results"), out, {"LFP": 0.95, "NCA": 0.90, "NMC": 0.90})
    ok = True
    for chem in S_TRUE:
        e = lev[chem]
        print(chem, "s =", round(e["s_25_35_per_K"], 4), "CI", np.round(e["s_25_35_ci95"], 4), "true", S_TRUE[chem],
              "| floor mean", round(flo[chem]["mean"], 4), "pairs", flo[chem]["pairs_with_difference"])
        ok &= e["s_25_35_ci95"][0] - 0.01 < S_TRUE[chem] < e["s_25_35_ci95"][1] + 0.01
    for k in ("all", "1C or below", "2C and 3C"):
        r = sa[k]
        print(k, {x: (round(r[x], 3) if isinstance(r[x], float) else r[x]) for x in
                  ("pairs", "warmer_shorter", "warmer_longer", "ties", "share", "sign_test_p_one_sided", "slope", "sd_ratio_pred_over_obs")},
              "slope CI", np.round(r.get("slope_ci95", [np.nan, np.nan]), 3))
    ok &= 0.6 < sa["all"]["slope"] < 1.4
    print("module:", round(mod["implied_divergence"], 5), "scorecard:", {k: v.get("held") for k, v in sc.items()})
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
