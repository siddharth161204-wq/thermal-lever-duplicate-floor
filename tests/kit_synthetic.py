#!/usr/bin/env python3
"""V4: synthetic EOC files for the 144 cyclic cells of inputs/kit/cyclic_cells.csv, with planted lifetimes, cell
temperatures and a planted lever per protocol; src/stage2/kit.py must recover them. No KIT value is used.
Run: python3 tests/kit_synthetic.py WORKDIR
"""
import json, math, os, sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src", "stage2"))
import kit as K  # noqa: E402

work = sys.argv[1]
os.makedirs(work, exist_ok=True)
cfg = pd.read_csv(os.path.join(ROOT, "inputs", "kit", "cyclic_cells.csv"))
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 3
rng = np.random.default_rng(SEED)
prots = sorted(cfg.protocol.unique())
S_TRUE = {p: 0.04 - 0.006 * i for i, p in enumerate(prots)}      # 0.040 down to -0.026 per K
truth = {}
for _, m in cfg.iterrows():
    T = m.age_temp + rng.normal(0.3, 0.3)
    L = 600.0 * math.exp(-S_TRUE[m.protocol] * (T - 25.0)) * math.exp(rng.normal(0, 0.03))
    q0 = 3.0 * (1 + rng.normal(0, 0.005))
    rows, t, tot = [], 1665598800.0, 0.0
    for cu in range(30):
        efc = tot / 3.0
        q = q0 * (1 - 0.1 * efc / L)
        tot += q
        rows.append((t, 1, 2, 0, 2.5, -1.0, q, -q, 25.0, 25.0, 12000.0, cu, tot))
        t += 3600 * 30
        for k in range(80):                                         # regular runs between check-ups
            dq = 3.0 * 0.9
            tot += dq if k % 2 else 0.0
            rows.append((t, 1, 1, int(k % 2 == 0), 4.2, -3.0, q, -dq if k % 2 else 0.0, T - 0.2, T + 0.2, 3600.0, cu + 1, tot))
            t += 3600
    cols = K.COLS
    df = pd.DataFrame(rows, columns=["timestamp_s", "timestamp_origin", "cyc_condition", "cyc_charged", "v_min_target_V",
                                     "i_dischg_min_A", "cap_aged_est_Ah", "delta_q_dischg_Ah", "t_start_degC", "t_end_degC",
                                     "cyc_duration_s", "num_cycles_checkup", "total_q_dischg_sum_Ah"])
    df.to_csv(os.path.join(work, f"cell_eocv2_P{m.parameter_id:03d}_{m.parameter_nr}_S00_C00.csv"), sep=";", index=False)
    truth[m.cell] = (L, T)
res = K.main(work, os.path.join(work, "out"))
cells = pd.read_csv(os.path.join(work, "out", "stage2_cells.csv")).set_index("cell")
errL = max(abs(cells.loc[c, "L90_efc"] - v[0]) / v[0] for c, v in truth.items() if np.isfinite(cells.loc[c, "L90_efc"]))
errT = max(abs(cells.loc[c, "T_C"] - v[1]) for c, v in truth.items())
inside = sum(r["ci95"][0] <= S_TRUE[r["protocol"]] <= r["ci95"][1] for r in res["lever_25_40"] if r["determinable"])
ndet = sum(r["determinable"] for r in res["lever_25_40"])
print(f"censored {res['censored90']}; max L error {errL:.2e}; max T error {errT:.2e} K; planted lever inside 95 % interval in {inside} of {ndet}")
ok = errL < 1e-9 and errT < 1e-9   # exact recovery; lever coverage is judged over several seeds
print("PASS" if ok else "FAIL")
json.dump({"max_rel_L_error": float(errL), "max_T_error_K": float(errT), "levers_inside": int(inside), "levers": int(ndet), "pass": bool(ok)},
          open(os.path.join(work, "v4.json"), "w"), indent=1)
sys.exit(0 if ok else 1)
