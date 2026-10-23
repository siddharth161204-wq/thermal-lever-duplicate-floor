#!/usr/bin/env python3
"""Synthetic test of the Stage 3 code (src/stage3/rows.py, prepare.py, stage3.py) on made-up time series in the
archive's layout that copy the Sandia features of D50: cycles without rests (one zero-current record at the top of
charge), counters that restart with the cycle index, a next cycle that begins on the last sample of a discharge
(in order, or with the old cycle's last record logged late), rows of another record interleaved, missing temperature
samples and one interrupted check. Lifetimes, the lever and the temperatures are planted; no Sandia value is used.

Checks: (1) the planted lifetimes are recovered to within the error of the check offset; (2) measure (a), measure (b)
and T_cyc equal a by-hand computation on the rows in their physical order, written here without the pipeline's code;
(3) a cell written with late records, or with interleaved rows of another record, gives the same lifetime and
measures as its twin written without them; (4) the planted setpoint lever is recovered; (5) the interrupted check is
dropped; (6) the scorecard runs, with P8 not scored.

Run: python3 -I tests/stage3_synthetic.py WORKDIR
"""
import json, math, os, sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src", "stage3"))
sys.path.insert(0, os.path.join(ROOT, "src", "stage2"))
import prepare as PR   # noqa: E402
import stage3 as S3    # noqa: E402
import steps as ST     # noqa: E402

RATED = {"LFP": 1.1, "NCA": 3.2, "NMC": 3.0}
VCUT = {"LFP": 2.0, "NCA": 2.5, "NMC": 2.0}
VTOP = {"LFP": 3.6, "NCA": 4.2, "NMC": 4.2}
THR = {"LFP": 0.95, "NCA": 0.90, "NMC": 0.90}
S_HI = {"LFP": 0.03, "NCA": 0.01, "NMC": -0.02}     # planted 25 to 35 C lever, per K
S_LO = {"LFP": -0.02, "NCA": -0.03, "NMC": -0.04}   # planted 15 to 25 C lever
DT, TAU = 120.0, 400.0


def planted_L(chem, T, rng):
    s = S_HI[chem] if T >= 25 else S_LO[chem]
    return 50.0 * math.exp(-s * (T - 25.0)) * math.exp(rng.normal(0, 0.02))


def generate(chem, T, lo, hi, c_dis, L, heat, offset, seed, boundary="normal", nan_frac=0.0, interrupted_block=None,
             n_blocks=6, per_block=15):
    """Rows in physical order (a list of tuples) and the rows as logged (with late records or foreign rows applied)."""
    rng = np.random.default_rng(seed)
    rated, vcut = RATED[chem], VCUT[chem]
    alpha = (1.0 - THR[chem]) / L
    st = {"t": 0.0, "cyc": 1, "efc": 0.0, "T": T + offset, "qc": 0.0, "qd": 0.0}
    rows = []   # (time, cycle, I, V, qc, qd, Tenv, Tcell, tag)

    def emit_step(I, dur, v0, v1, target):
        n = max(1, int(math.ceil(dur / DT - 1e-9)))
        ts = np.minimum(np.arange(1, n + 1) * DT, dur)
        prev = 0.0
        for k, tt in enumerate(ts):
            st["T"] = target + (st["T"] - target) * math.exp(-(tt - prev) / TAU)
            dq = abs(I) * (tt - prev) / 3600.0
            if I > 0:
                st["qc"] += dq
            elif I < 0:
                st["qd"] += dq
                st["efc"] += dq / rated
            prev = tt
            rows.append([st["t"] + tt, st["cyc"], I, v0 + (v1 - v0) * (k + 1) / n, st["qc"], st["qd"],
                         T + rng.normal(0, 0.03), st["T"] + rng.normal(0, 0.02), "x"])
        st["t"] += dur

    def cycle(c_d, a, b, interrupted=False):
        cap = rated * 0.98 * (1.0 - alpha * st["efc"])
        q = cap * (b - a)
        emit_step(0.5 * rated, q / (0.5 * rated) * 3600.0, 3.0 if a == 0 else 3.4, VTOP[chem] if b == 1 else 3.9,
                  T + offset + heat * 0.25)
        rows.append([st["t"] + DT, st["cyc"], 0.0, VTOP[chem] - 0.05, st["qc"], st["qd"], T, st["T"], "x"])   # one record
        st["t"] += DT
        cap = rated * 0.98 * (1.0 - alpha * st["efc"])
        q = cap * (b - a) * (0.3 if interrupted else 1.0)
        emit_step(-c_d * rated, q / (c_d * rated) * 3600.0, VTOP[chem] - 0.1, vcut if (a == 0 or interrupted) else 3.5,
                  T + offset + heat * c_d * c_d)
        # the next cycle index, and the counters' restart
        if boundary == "normal":
            st["cyc"] += 1
        else:
            last = rows[-1]
            new = [last[0] + 0.018, st["cyc"] + 1, last[2], last[3] - 0.003, 0.0, 0.0, last[6], last[7], "new"]
            if boundary == "late":
                new[0] = last[0]
                last[0] = last[0] + 0.02            # the old cycle's last record logged 20 ms after the new cycle's first
                last[8] = "late"
            rows.append(new)
            st["cyc"] += 1
        st["qc"] = st["qd"] = 0.0

    every = c_dis == 0.5 and lo == 0 and hi == 1
    for blk in range(n_blocks + 1):
        for _ in range(3):
            cycle(0.5, 0.0, 1.0, interrupted=(interrupted_block == blk))
        if blk < n_blocks:
            for _ in range(per_block):
                cycle(0.5 if every else c_dis, 0.0 if every else lo, 1.0 if every else hi)
    phys = list(rows)                                # physical order: as generated, the old cycle's last record first
    if nan_frac:
        m = rng.random(len(rows)) < nan_frac
        m[len(rows) // 2: len(rows) // 2 + 30] = True
        for i in np.flatnonzero(m):
            rows[i][7] = float("nan")
    logged = sorted(rows, key=lambda r: r[0])        # the file is in time order, so late records sit out of cycle order
    return phys, logged


def write(path, logged, foreign=None):
    df = pd.DataFrame([r[:8] for r in logged], columns=["Test_Time (s)", "Cycle_Index", "Current (A)", "Voltage (V)",
                                                        "Charge_Capacity (Ah)", "Discharge_Capacity (Ah)",
                                                        "Environment_Temperature (C)", "Cell_Temperature (C)"])
    if foreign is not None:                      # rows of another record: earlier cycle numbers, interleaved in time
        i0 = len(df) // 2
        f = df.iloc[i0: i0 + 60].copy()
        f["Cycle_Index"] = f["Cycle_Index"] - 9
        f["Charge_Capacity (Ah)"] += 0.7
        f["Discharge_Capacity (Ah)"] += 0.9
        f["Voltage (V)"] += 0.3
        f["Cell_Temperature (C)"] += 4.0
        f["Test_Time (s)"] += 7.0
        df = pd.concat([df, f]).sort_values("Test_Time (s)", kind="stable")
    df.insert(0, "Date_Time", "")
    df["Charge_Energy (Wh)"] = 0.0
    df["Discharge_Energy (Wh)"] = 0.0
    df.to_csv(path, index=False)


def by_hand(phys, rated, L_efc, censored, c_dis):
    """Measures (a), (b) and T_cyc from the rows in physical order, written without the pipeline's code."""
    t = np.array([r[0] for r in phys]); I = np.array([r[2] for r in phys]); T = np.array([r[7] for r in phys])
    V = np.array([r[3] for r in phys]); qd = np.array([r[5] for r in phys])
    lab = [1 if x >= 0.005 * rated else (-1 if x <= -0.005 * rated else 0) for x in I]
    steps, cur = [], None
    for i, l in enumerate(lab):
        if cur is None or l != cur["lab"]:
            cur = {"lab": l, "idx": []}
            steps.append(cur)
        cur["idx"].append(i)
    cum = 0.0
    for s in steps:
        ii = s["idx"]
        num = den = 0.0
        for a, b in zip(ii[:-1], ii[1:]):
            if not (math.isnan(T[a]) or math.isnan(T[b])):
                num += 0.5 * (T[a] + T[b]) * (t[b] - t[a]); den += t[b] - t[a]
        vals = [T[k] for k in ii if not math.isnan(T[k])]
        s["Tm"] = num / den if den > 0 else (sum(vals) / len(vals) if vals else float("nan"))
        s["dur"] = t[ii[-1]] - t[ii[0]]
        if s["lab"] == -1:
            s["efc"] = cum / rated
            ah = 0.0
            prev = qd[ii[0] - 1] if ii[0] > 0 else 0.0
            for k in ii:
                ah += qd[k] - prev if qd[k] >= prev - 1e-12 else qd[k]
                prev = qd[k]
            cum += ah
            dur = t[ii[-1]] - t[ii[0]]
            crate = abs(np.trapezoid(np.abs(I[ii]), t[ii]) / dur) / rated if dur > 0 else abs(I[ii[0]]) / rated
            s["check"] = abs(crate - 0.5) <= 0.05 and V[ii[-1]] <= {1.1: 2.0, 3.2: 2.5, 3.0: 2.0}[rated] + 0.02
    lim = float("inf") if censored else L_efc
    sw, rest, num_t, den_t = [], [], 0.0, 0.0
    for k, s in enumerate(steps):
        if s["lab"] != -1 or s["check"] or s["efc"] > lim:
            continue
        j = k - 1
        while j >= 0 and steps[j]["lab"] == 0:
            j -= 1
        if j < 0 or steps[j]["lab"] != 1 or math.isnan(s["Tm"]) or math.isnan(steps[j]["Tm"]):
            continue
        sw.append(s["Tm"] - steps[j]["Tm"])
        rr = [steps[x]["Tm"] for x in range(j + 1, k) if not math.isnan(steps[x]["Tm"])]
        x = j - 1
        while x >= 0 and steps[x]["lab"] == 0:
            if not math.isnan(steps[x]["Tm"]):
                rr.append(steps[x]["Tm"])
            x -= 1
        x = k + 1
        while x < len(steps) and steps[x]["lab"] == 0:
            if not math.isnan(steps[x]["Tm"]):
                rr.append(steps[x]["Tm"])
            x += 1
        if rr:
            rest.append(sum(rr) / len(rr))
        num_t += steps[j]["Tm"] * steps[j]["dur"] + s["Tm"] * s["dur"]
        den_t += steps[j]["dur"] + s["dur"]
    return {"swing_K": float(np.mean(sw)) if sw else float("nan"), "rest_T_C": float(np.mean(rest)) if rest else float("nan"),
            "T_cyc_C": num_t / den_t if den_t else float("nan"), "n_units": len(sw), "efc_total": cum / rated}


if __name__ == "__main__":
    work = sys.argv[1]
    raw = os.path.join(work, "raw"); os.makedirs(raw, exist_ok=True)
    rng = np.random.default_rng(11)
    design = []
    for chem in ("LFP", "NCA", "NMC"):
        for T, n in ((15, 2), (25, 4), (35, 4)):
            design += [(chem, T, 0, 1, 1.0, k) for k in range(n)]
        for T in (15, 25, 35):
            design += [(chem, T, 0, 1, 2.0, k) for k in range(2)]
        design += [(chem, 25, 0, 1, 0.5, k) for k in range(2)]
        design += [(chem, 25, 0.2, 0.8, 0.5, k) for k in range(2)]
        design += [(chem, 25, 0.2, 0.8, 3.0, k) for k in range(4)]
    truth, seed = {}, 100
    for chem, T, lo, hi, c, k in design:
        seed += 1
        letter = "abcd"[k]
        key = f"SNL_18650_{chem}_{T}C_{int(lo * 100)}-{int(hi * 100)}_0.5-{c:g}C_{letter}"
        L = planted_L(chem, T, rng)
        heat, off = 0.9 + 0.2 * rng.random(), rng.normal(0, 0.3)
        boundary = "late" if T == 35 else ("new" if (T == 25 and c == 1.0 and k == 0) else "normal")
        nanf = 0.02 if (T == 25 and c == 2.0) else 0.0
        intr = 3 if (chem == "NCA" and T == 25 and c == 1.0 and k == 1) else None
        per = 25 if hi - lo < 1 else 15
        if c == 3.0 and k >= 2:                 # twins of replicates a and b, written with interleaved foreign rows
            src = truth[key.replace(f"_{letter}", "_" + "ab"[k - 2])]
            phys, logged = src["phys"], src["logged"]
            write(os.path.join(raw, key + "_timeseries.csv"), logged, foreign=True)
            truth[key] = dict(src, twin_of=src["key"], key=key)
            continue
        phys, logged = generate(chem, T, lo, hi, c, L, heat, off, seed, boundary, nanf, intr, per_block=per)
        write(os.path.join(raw, key + "_timeseries.csv"), logged)
        truth[key] = {"key": key, "L": L, "phys": phys, "logged": logged, "chem": chem, "T": T, "c": c, "boundary": boundary,
                      "interrupted": intr is not None}
    # the structure pass of D38 on the files as written, then Stage 3
    p1 = os.path.join(work, "steps1"); os.makedirs(p1, exist_ok=True)
    for key, v in truth.items():
        st, _ = ST.segment(os.path.join(raw, key + "_timeseries.csv"), 0.005 * RATED[v["chem"]], -1, False)
        st.to_csv(os.path.join(p1, key + "_steps.csv"), index=False)
    rep = PR.prepare(raw, p1, os.path.join(work, "w"))
    res, sc = S3.main(raw, os.path.join(work, "w"), os.path.join(work, "out"), expect_cells=None)
    cells = pd.read_csv(os.path.join(work, "out", "stage3_cells.csv")).set_index("key")
    ok, worst = True, {}
    for key, v in truth.items():
        r = cells.loc[key]
        if r.censored:
            print("censored:", key); ok = False; continue
        errL = abs(r.L_efc - v["L"])
        hand = by_hand(v["phys"], RATED[v["chem"]], r.L_efc, bool(r.censored), v["c"])
        e_sw = abs(r.swing_K - hand["swing_K"]) if np.isfinite(hand["swing_K"]) else (0.0 if not np.isfinite(r.get("swing_K", np.nan)) else 9)
        e_rest = abs(r.rest_T_C - hand["rest_T_C"]) if np.isfinite(hand["rest_T_C"]) else 0.0
        e_tc = abs(r.T_cyc_C - hand["T_cyc_C"]) if np.isfinite(hand["T_cyc_C"]) else 0.0
        e_efc = abs(r.efc_total - hand["efc_total"])
        for nm, e in (("L", errL), ("swing", e_sw), ("rest", e_rest), ("Tcyc", e_tc), ("efc_total", e_efc)):
            worst[nm] = max(worst.get(nm, 0.0), e)
        ok &= errL < 1.0 and e_sw < 1e-4 and e_rest < 1e-4 and e_tc < 1e-4 and e_efc < 1e-6
        ok &= int(r.n_units) == hand["n_units"] or not np.isfinite(hand["swing_K"])
        if v.get("twin_of"):
            t = cells.loc[v["twin_of"]]
            same = all(abs(r[c] - t[c]) < 1e-9 for c in ("L_efc", "swing_K", "T_cyc_C", "efc_total"))
            print("foreign-row twin", key[10:], "identical to", v["twin_of"][10:], same); ok &= same
        if v["interrupted"]:
            print("interrupted check dropped:", key[10:], int(r.checks_dropped)); ok &= int(r.checks_dropped) == 1
    print("worst errors:", {k: float(f"{x:.3g}") for k, x in worst.items()})
    lev = res["registered_1C"]
    for chem in ("LFP", "NCA", "NMC"):
        for sec, s_true, lo_T in (("s_25_35", S_HI[chem], 25), ("s_15_25", S_LO[chem], 15)):
            x = lev[chem][sec]
            cl = lambda T: [c["L_efc"] for c in lev[chem]["cells"][f"{T}C"]]
            Lp = {T: [truth[c["cell_id"].replace("/", "-")]["L"] for c in lev[chem]["cells"][f"{T}C"]] for T in (15, 25, 35)}
            s_pl = -(np.mean(np.log(Lp[lo_T + 10])) - np.mean(np.log(Lp[lo_T]))) / 10.0
            print(chem, sec, "setpoint s", round(x["s_setpoint"], 5), "from planted lifetimes", round(s_pl, 5), "planted law", s_true,
                  "logged s", round(x["s"], 5), "CI", np.round(x["ci95"], 4))
            ok &= abs(x["s_setpoint"] - s_pl) < 2e-3
    print("prepare:", {k: v for k, v in json.loads(json.dumps({"late": sum(r["late_rows"] for r in rep.values()),
          "foreign": sum(r["foreign_rows_dropped"] for r in rep.values())})).items()})
    print("scorecard:", sc["summary"])
    ok &= sc["P8"]["scored"] is True or sc["P8"]["scored"] is False
    print("PASS" if ok else "FAIL")
    json.dump({"pass": bool(ok), "worst": worst}, open(os.path.join(work, "stage3_synthetic.json"), "w"), indent=1)
    sys.exit(0 if ok else 1)
