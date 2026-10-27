#!/usr/bin/env python3
"""POST-FREEZE: D62's normal check read literally, "normal differences of the same sd(e)": normal noise with the
registered (log-normal) sd(e), simulated n against the calculated registered n. The frozen simulate.py used the normal
shape's own sd(e) = m sqrt(pi/2) instead; both check the noncentral-t calculation for normal noise. Same bisection,
10,000 experiments per n, seeds SeedSequence([20261027, case index, n, 3]).
Output: out/postfreeze_simulation_normal_same_sd.csv
"""
import csv, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
sys.path.insert(0, HERE)
import simulate as S                                # noqa: E402

S.SHAPE_CODE["normal_same_sd"] = 3


def main():
    par_all = S.shapes()
    with open(os.path.join(OUT, "experiment.csv")) as f:
        rows = list(csv.DictReader(f))
    res = []
    for idx, r in enumerate(rows):
        if r["shape"] != "lognormal":
            continue
        dT, s, n_calc = float(r["dT_K"]), float(r["s_per_K"]), int(r["n"])
        p = dict(par_all[r["floor"]])
        p["sd_norm"] = p["sd_ln"]                   # the registered sd(e), with normal noise
        t0 = time.time()
        # bisect_n draws normal noise for shape "normal"; the seed's shape code is changed to 3 below
        code = S.SHAPE_CODE["normal"]
        S.SHAPE_CODE["normal"] = 3
        n_sim, cache = S.bisect_n(n_calc, -s * dT, "normal", p, idx, 10000)
        S.SHAPE_CODE["normal"] = code
        res.append({"case_index": idx, "group": r["group"], "lever": r["lever"], "dT_K": dT, "floor": r["floor"],
                    "n_calc": n_calc, "n_sim": n_sim, "rel_diff": (n_sim - n_calc) / n_calc,
                    "V4_pass_5pct": abs(n_sim - n_calc) <= 0.05 * n_calc, "sim_power_at_n_sim": cache[n_sim][0],
                    "seconds": time.time() - t0})
        print(res[-1], flush=True)
    with open(os.path.join(OUT, "postfreeze_simulation_normal_same_sd.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(res[0].keys()))
        w.writeheader()
        for x in res:
            w.writerow(x)


if __name__ == "__main__":
    main()
