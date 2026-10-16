#!/usr/bin/env python3
"""Synthetic end-to-end test of the Stage 2 pipeline (no Sandia data): writes time series in the archive's standard
column layout with planted lifetimes and self-heating, then checks that steps.py, lifetimes.py and process_cells.py
recover them. Run: python3 tests/stage2_synthetic.py WORKDIR
"""
import json, math, os, sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src", "stage2"))
import steps as ST          # noqa: E402
import process_cells as PC  # noqa: E402

DT = 30.0
V_HI, V_CUT = 4.2, 2.5


def make_cell(path, rated, c_dis, soc_lo, soc_hi, alpha, heat_K_per_C2, offset_K, T_amb=25.0,
              n_blocks=6, per_block=60, tau=600.0, seed=0, rest_s=1800.0):
    rng = np.random.default_rng(seed)
    rows, t, cyc = [], 0.0, 0
    state = {"T": T_amb + offset_K, "efc": 0.0, "qc": 0.0, "qd": 0.0}

    def cap():
        return rated * 0.98 * (1.0 - alpha * state["efc"])

    def emit(seg_dur, I, v0, v1, kind):
        nonlocal t
        n = max(1, int(math.ceil(seg_dur / DT)))
        ts = np.linspace(0.0, seg_dur, n + 1)[1:]
        c = abs(I) / rated
        target = T_amb + offset_K + (heat_K_per_C2 * c * c if kind != "rest" else 0.0)
        prev_t = 0.0
        for k, tt in enumerate(ts):
            dt = tt - prev_t
            prev_t = tt
            state["T"] = target + (state["T"] - target) * math.exp(-dt / tau)
            if kind == "charge":
                state["qc"] += abs(I) * dt / 3600.0
            elif kind == "discharge":
                state["qd"] += abs(I) * dt / 3600.0
                state["efc"] += abs(I) * dt / 3600.0 / rated
            v = v0 + (v1 - v0) * (k + 1) / len(ts)
            rows.append((t + tt, cyc, I, v, state["qc"], state["qd"], T_amb + rng.normal(0, 0.02),
                         state["T"] + rng.normal(0, 0.03)))
        t += seg_dur

    def cycle(c_ch, c_d, lo, hi, rest_s):
        nonlocal cyc
        cyc += 1
        state["qc"] = state["qd"] = 0.0
        q = cap() * (hi - lo)
        emit(q / (c_ch * rated) * 3600.0, c_ch * rated, V_CUT + 0.3 if lo == 0 else 3.4, V_HI if hi == 1 else 4.0, "charge")
        emit(rest_s, 0.0, 4.1, 4.1, "rest")
        q = cap() * (hi - lo)
        emit(q / (c_d * rated) * 3600.0, -c_d * rated, V_HI - 0.1, V_CUT if lo == 0 else 3.5, "discharge")
        emit(rest_s, 0.0, 3.3, 3.3, "rest")

    emit(600.0, 0.0, 3.0, 3.0, "rest")
    for b in range(n_blocks + 1):
        for _ in range(3):
            cycle(0.5, 0.5, 0.0, 1.0, 600.0)
        if b < n_blocks:
            for _ in range(per_block):
                cycle(0.5, c_dis, soc_lo, soc_hi, rest_s)
    a = np.array(rows)
    df = pd.DataFrame({"Date_Time": "", "Test_Time (s)": a[:, 0], "Cycle_Index": a[:, 1].astype(int), "Current (A)": a[:, 2],
                       "Voltage (V)": a[:, 3], "Charge_Capacity (Ah)": a[:, 4], "Discharge_Capacity (Ah)": a[:, 5],
                       "Charge_Energy (Wh)": 0.0, "Discharge_Energy (Wh)": 0.0,
                       "Environment_Temperature (C)": a[:, 6], "Cell_Temperature (C)": a[:, 7]})
    df.to_csv(path, index=False)
    return df


def by_hand_rise(df, rated):
    """Independent recomputation of measure (a) on the synthetic file (V2-style): label, runs, time-weighted means."""
    t = df["Test_Time (s)"].to_numpy(); I = df["Current (A)"].to_numpy(); T = df["Cell_Temperature (C)"].to_numpy()
    lab = np.where(I > 0.005 * rated, 1, np.where(I < -0.005 * rated, -1, 0))
    b = np.flatnonzero(np.diff(lab)) + 1
    starts, ends = np.r_[0, b], np.r_[b, lab.size]
    runs = []
    for s, e in zip(starts, ends):
        tt, TT = t[s:e], T[s:e]
        dur = tt[-1] - tt[0]
        m = np.trapezoid(TT, tt) / dur if dur > 0 else TT[0]
        runs.append((lab[s], dur, m, abs(np.trapezoid(np.abs(I[s:e]), tt) / dur) / rated if dur > 0 else 0))
    return runs


if __name__ == "__main__":
    work = sys.argv[1]
    os.makedirs(work, exist_ok=True)
    rated = 3.2
    spec = [("SNL_18650_NCA_25C_0-100_0.5-1C_a", 1.0, 0, 1, 4.0e-4, 1.2, 0.4),
            ("SNL_18650_NCA_25C_0-100_0.5-1C_b", 1.0, 0, 1, 4.4e-4, 1.5, -0.3),
            ("SNL_18650_NCA_25C_20-80_0.5-0.5C_a", 0.5, 0.2, 0.8, 5.0e-4, 1.0, 0.0),
            ("SNL_18650_NCA_25C_0-100_0.5-0.5C_a", 0.5, 0, 1, 1.0e-3, 1.0, 0.2)]
    planted = {}
    for k, (key, c, lo, hi, alpha, heat, off) in enumerate(spec):
        p = os.path.join(work, key + "_timeseries.csv")
        per = 60 if c != 0.5 or hi - lo < 1 else 20
        df = make_cell(p, rated, c, lo, hi, alpha, heat, off, per_block=per, n_blocks=8 if hi - lo < 1 else 6, seed=k)
        st, stats = ST.segment(p, 0.005 * rated, -1, False, chunksize=7919)
        st.to_csv(os.path.join(work, key + "_steps.csv"), index=False)
        stT, _ = ST.segment(p, 0.005 * rated, -1, True, chunksize=10007)
        stT.to_csv(os.path.join(work, key + "_steps_T.csv"), index=False)
        planted[key] = {"L_true": 0.1 / alpha, "rows": stats["rows"], "steps": len(st)}
        # chunk independence: segment again in one chunk and compare
        one, _ = ST.segment(p, 0.005 * rated, -1, True, chunksize=10 ** 9)
        assert len(one) == len(stT), (len(one), len(stT))
        assert np.allclose(one["ah_int"], stT["ah_int"]) and np.allclose(one["t_cell_mean"], stT["t_cell_mean"], equal_nan=True)
        runs = by_hand_rise(df, rated)
        assert len(runs) == len(stT), (len(runs), len(stT))
        assert np.allclose([r[2] for r in runs], stT["t_cell_mean"].to_numpy())
    params = {"ah_column": "ah_ctr", "v_cut_low": {"NCA": V_CUT, "LFP": 2.0, "NMC": 2.5}, "tol_v": 0.02, "tol_c": 0.05,
              "thresholds": {"NCA": 0.90, "LFP": 0.95, "NMC": 0.90}, "rest_baseline": "mean"}
    cells, checks = PC.process(work, params, with_temperature=True)
    ok = True
    for _, r in cells.iterrows():
        tru = planted[r.key]["L_true"]
        err = abs(r.L_efc - tru)
        print(f"{r.cell_id}: L={r.L_efc:.2f} true={tru:.2f} err={err:.3f} res={r.resolution_efc:.2f} "
              f"checks={r.n_checks} rise={r.rise_K:.3f} K (end-baseline {r.rise_alt_K:.3f}) rest={r.rest_T_C:.3f} C units={r.n_units} "
              f"Lsens={r.get('L_sens_efc', float('nan')):.2f}")
        ok &= err < 1.0
    a, b = cells.set_index("key").loc[["SNL_18650_NCA_25C_0-100_0.5-1C_a", "SNL_18650_NCA_25C_0-100_0.5-1C_b"]].itertuples()
    ok &= b.rise_K > a.rise_K                     # planted: b heats more (1.5 vs 1.2 K per C^2)
    ok &= abs((a.rest_T_C - b.rest_T_C) - 0.7) < 0.1   # planted offsets +0.4 and -0.3 K
    print("checks per cell:", checks.groupby("cell_id").size().to_dict())
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
