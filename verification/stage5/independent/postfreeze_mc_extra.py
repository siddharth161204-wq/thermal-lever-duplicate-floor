#!/usr/bin/env python3
"""Stage 5 independent check, written after the freeze (ea3a5b4) for the comparison: the Monte Carlo spread of the
quantities whose spread the frozen outputs do not hold, so that every Monte Carlo difference can be read in units of
that spread. Nothing frozen is changed.

(1) Sandia measure (b) beta: 200 replicate 10,000-draw bootstraps and one of 1,000,000 draws (seeds [BASE_SEED, 1100, r]
    and [BASE_SEED, 1101]).
(2) C1 with normal e of the registered sd(e): 20 replicate 10,000-dataset simulations per case ([BASE_SEED, 1200 + case,
    floor, r]).
(3) C2 simulated pairs needed for the ten cases without frozen replicates: 5 replicate bisections each ([BASE_SEED,
    1300 + case, r, n]); the two KIT cases under the KIT floor have theirs in out/C_mc_spread.json.

Output: out_compare/postfreeze_mc_extra.json.
Usage: PYTHONDONTWRITEBYTECODE=1 python3 verification/stage5/independent/postfreeze_mc_extra.py
"""
import json, os, sys, time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ind_common as C   # noqa: E402
import partb as B        # noqa: E402
import partc as P        # noqa: E402

OUTC = os.path.join(C.HERE, "out_compare")


def main():
    t0 = time.time()
    rows, _ = C.build_pairs()
    sh = C.shapes()
    med = C.band_medians()
    out = {"note": "post-freeze Monte Carlo spread for the comparison; no frozen output is changed"}
    # (1) measure (b)
    S = C.select(rows, "Sandia")
    St = [B.stratum(B.as_xy(S, "b"))]
    reps = [B.bootstrap(St, (1100, r)) for r in range(200)]
    big = B.bootstrap(St, (1101,), n_draws=1_000_000)
    arr = np.array([x["beta_ci95"] for x in reps])
    out["Sandia measure b beta"] = {"replicates": 200, "lo_mean": float(arr[:, 0].mean()), "lo_sd": float(arr[:, 0].std(ddof=1)),
                                    "hi_mean": float(arr[:, 1].mean()), "hi_sd": float(arr[:, 1].std(ddof=1)),
                                    "big_ci95": big["beta_ci95"]}
    print("measure b", out["Sandia measure b beta"], flush=True)
    # (2) normal-e simulation spread
    Cj = json.load(open(os.path.join(C.OUT, "C.json")))
    sel = {"KIT": C.select(rows, "KIT"), "Sandia": C.select(rows, "Sandia"), "pooled": C.select(rows),
           "KIT warm": C.select(rows, "KIT", "warm"), "KIT cold": C.select(rows, "KIT", "cold"),
           "Sandia warm": C.select(rows, "Sandia", "warm"), "Sandia cold": C.select(rows, "Sandia", "cold"),
           "pooled warm": C.select(rows, None, "warm"), "pooled cold": C.select(rows, None, "cold"),
           "KIT all 142": C.select(rows, "KIT", scored=False)}
    out["C1 normal same sd"] = []
    for ci, r in enumerate(Cj["C1"]):
        dT = np.array([x["dT_K"] for x in sel[r["case"]]])
        fi = 0 if r["floor"] == "Geslin" else 1
        calc = r["detectable_s_calc_pct_per_K"] / 100.0
        vals = np.array([P.FixedDesign(dT, sh[r["floor"]], "normal", C.rng_for(1200 + ci, fi, rep)).detectable_bisection(calc)
                         for rep in range(20)])
        out["C1 normal same sd"].append({"case": r["case"], "floor": r["floor"], "replicates": 20,
                                         "mean_pct_per_K": 100 * float(vals.mean()), "sd_pct_per_K": 100 * float(vals.std(ddof=1))})
    print("normal-e done", f"{time.time() - t0:.0f} s", flush=True)
    # (3) C2 simulated pairs, the ten cases without frozen replicates
    out["C2 n_sim"] = []
    for ci, r in enumerate(Cj["C2"]):
        if r["dataset"] == "KIT" and r["floor"] == "KIT":
            continue
        name = f"{r['dataset']} {r['band']}"
        dT = np.array([x["dT_K"] for x in sel[name]])
        ns = []
        for rep in range(5):
            n, _ = P.n_simulated(abs(med[r["band"]]["s"]), dT, sh[r["floor"]], "lognormal", r["n_required_calc"], (1300 + ci, rep))
            ns.append(n)
        out["C2 n_sim"].append({"dataset": r["dataset"], "band": r["band"], "floor": r["floor"], "n_calc": r["n_required_calc"],
                                "registered_run_n_sim": r["n_required_sim_lognormal"], "replicate_n_sim": ns,
                                "mean": float(np.mean(ns)), "sd": float(np.std(ns, ddof=1))})
        print(name, r["floor"], ns, f"{time.time() - t0:.0f} s", flush=True)
    os.makedirs(OUTC, exist_ok=True)
    with open(os.path.join(OUTC, "postfreeze_mc_extra.json"), "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
