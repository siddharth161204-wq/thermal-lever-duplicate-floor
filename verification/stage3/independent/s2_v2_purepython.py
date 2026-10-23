#!/usr/bin/env python3
"""Stage 3 independent check, step 2: V2 of DECISIONS.md D49 as registered, plus four extra cells (unregistered),
recomputed from the raw time-series CSV with the csv module and plain loops only (no numpy, no pandas), as a second
implementation beside s1_cells.py. Uses none of the pipeline's functions.

Usage: python3 -I s2_v2_purepython.py <out_json> <cell_id> [<cell_id> ...]
"""
import csv, json, math, os, sys

RAW = {"LFP": "/tmp/sandia/raw_LFP/SNL LFP", "NCA": "/tmp/sandia/raw_NCA/SNL NCA", "NMC": "/tmp/sandia/raw_NMC/SNL NMC"}
RATED = {"LFP": 1.1, "NCA": 3.2, "NMC": 3.0}
CUT = {"LFP": 2.0, "NCA": 2.5, "NMC": 2.0}
THR = {"LFP": 0.95, "NCA": 0.90, "NMC": 0.90}


def fnum(s):
    return float(s) if s != "" else None


def read_rows(cid):
    p = cid.split("_")
    chem = p[2]
    fname = "_".join(p[:5] + [p[5].replace("/", "-"), p[6]]) + "_timeseries.csv"
    path = os.path.join(RAW[chem], fname)
    rows = []
    with open(path, newline="") as fh:
        rd = csv.reader(fh)
        hdr = next(rd)
        ix = {h: i for i, h in enumerate(hdr)}
        for rec in rd:
            rows.append((float(rec[ix["Test_Time (s)"]]), float(rec[ix["Cycle_Index"]]), float(rec[ix["Current (A)"]]),
                         float(rec[ix["Voltage (V)"]]), float(rec[ix["Charge_Capacity (Ah)"]]),
                         float(rec[ix["Discharge_Capacity (Ah)"]]), fnum(rec[ix["Cell_Temperature (C)"]]),
                         fnum(rec[ix["Environment_Temperature (C)"]])))
    return chem, path, rows


def analyse(cid, unit_to_show=0):
    chem, path, rows = read_rows(cid)
    C = RATED[chem]
    p = cid.split("_")
    window = p[4]
    c_ch, c_dis = (float(x) for x in p[5].rstrip("C").split("/"))
    # D50 (1): running maximum, late records, rows of another record; drop the latter, stable sort by cycle index
    runmax = -1e300
    late = other = 0
    kept = []
    for i, r in enumerate(rows):
        runmax = max(runmax, r[1])
        if r[1] == runmax - 1:
            late += 1
        if r[1] < runmax - 1:
            other += 1
            continue
        kept.append((r[1], i, r))
    kept.sort(key=lambda x: x[0])            # Python's sort is stable: equal cycle indices keep logged order
    R = [x[2] for x in kept]
    n = len(R)
    # D40 labels and steps
    lab = [0 if abs(r[2]) < 0.005 * C else (1 if r[2] > 0 else -1) for r in R]
    steps = []
    s = 0
    for j in range(1, n + 1):
        if j == n or lab[j] != lab[s]:
            steps.append([s, j - 1, lab[s]])
            s = j
    # per-sample counter rule, D50 (2)
    def ah(col, s, e):
        tot = R[s][col] if s == 0 else 0.0
        for j in range(max(s, 1), e + 1):
            d = R[j][col] - R[j - 1][col]
            tot += d if d >= 0 else R[j][col]
        return tot

    def reg_ah(col, s, e):
        before = R[s - 1][col] if s > 0 else 0.0
        last = R[e][col]
        return last - before if last >= before else last

    # D50 (3) step temperature
    def tmean(s, e, col=6):
        num = den = 0.0
        npair = 0
        for j in range(s, e):
            a, b = R[j][col], R[j + 1][col]
            if a is not None and b is not None:
                dt = R[j + 1][0] - R[j][0]
                num += 0.5 * (a + b) * dt
                den += dt
                npair += 1
        if npair and den != 0:
            return num / den
        # no valid pair, or valid pairs whose durations sum to zero (two records at one instant): mean of valid samples
        v = [R[j][col] for j in range(s, e + 1) if R[j][col] is not None]
        return sum(v) / len(v) if v else None

    # discharges, checks (D41)
    dis = [k for k, st in enumerate(steps) if st[2] == -1]
    dah, cum, qual, imean = {}, {}, {}, {}
    run = 0.0
    tot_dis = tot_reg = 0.0
    for k in dis:
        s, e, _ = steps[k]
        dah[k] = ah(5, s, e)
        cum[k] = run
        run += dah[k]
        tot_dis += dah[k]
        tot_reg += reg_ah(5, s, e)
        im = sum(R[j][2] for j in range(s, e + 1)) / (e - s + 1)
        imean[k] = im
        qual[k] = abs(-im - 0.5 * C) <= 0.05 * C and abs(R[e][3] - CUT[chem]) <= 0.02
    checks = []
    if window == "0-100" and c_dis == 0.5:
        for k in dis:
            if qual[k]:
                checks.append([k])
    else:
        cur = []
        for k in dis:
            if qual[k]:
                cur.append(k)
            else:
                if cur:
                    checks.append(cur)
                cur = []
        if cur:
            checks.append(cur)
    cap = [sum(dah[k] for k in g) / len(g) for g in checks]
    efc = [cum[g[0]] / C for g in checks]
    q0 = cap[0]
    keep = [i for i in range(len(checks)) if cap[i] >= 0.5 * q0]
    # D42 lifetime
    thr = THR[chem] * q0
    L = res = None
    br = None
    for a, b in zip(keep[:-1], keep[1:]):
        if cap[b] < thr:
            L = efc[a] + (thr - cap[a]) * (efc[b] - efc[a]) / (cap[b] - cap[a])
            res = efc[b] - efc[a]
            br = (a, b)
            break
    censored = L is None
    # units (D43, D50 (5)): every non-check discharge with the last charge step before it
    units = []
    last_charge = None
    for k, st in enumerate(steps):
        if st[2] == 1:
            last_charge = k
        if st[2] == -1 and not qual[k]:
            units.append((last_charge, k))
    used = []
    for c, d in units:
        efc_d = cum[d] / C
        if not censored and efc_d > L:
            continue
        if c is None:
            continue
        cm, dm = tmean(steps[c][0], steps[c][1]), tmean(steps[d][0], steps[d][1])
        if cm is None or dm is None:
            continue
        rests = []
        for j in (c - 1, c + 1, d - 1, d + 1):
            if 0 <= j < len(steps) and steps[j][2] == 0 and j not in rests:
                rests.append(j)
        rm = [tmean(steps[j][0], steps[j][1]) for j in rests]
        rm = [x for x in rm if x is not None]
        used.append({"charge_step": c, "discharge_step": d, "efc_at_discharge_start": efc_d, "charge_mean_C": cm,
                     "discharge_mean_C": dm, "swing_K": dm - cm,
                     "charge_duration_s": R[steps[c][1]][0] - R[steps[c][0]][0],
                     "discharge_duration_s": R[steps[d][1]][0] - R[steps[d][0]][0],
                     "rest_steps": rests, "rest_means_C": rm, "rest_mean_C": (sum(rm) / len(rm)) if rm else None})
    a_meas = sum(u["swing_K"] for u in used) / len(used) if used else None
    bl = [u["rest_mean_C"] for u in used if u["rest_mean_C"] is not None]
    b_meas = sum(bl) / len(bl) if bl else None
    if used:
        num = sum(u["charge_mean_C"] * u["charge_duration_s"] + u["discharge_mean_C"] * u["discharge_duration_s"] for u in used)
        den = sum(u["charge_duration_s"] + u["discharge_duration_s"] for u in used)
        T_cyc = num / den
    else:
        T_cyc = None
    out = {"cell_id": cid, "file": path, "rows_logged": len(rows), "late_records": late, "rows_other_record_dropped": other,
           "steps": len(steps), "steps_charge": sum(1 for x in steps if x[2] == 1),
           "steps_discharge": len(dis), "steps_rest": sum(1 for x in steps if x[2] == 0),
           "discharge_ah_total": tot_dis, "registered_counter_ah_total": tot_reg,
           "registered_share": tot_reg / tot_dis if tot_dis else None,
           "checks": len(checks), "checks_dropped_half": len(checks) - len(keep), "first_check_Ah": q0,
           "L_efc": L, "resolution_efc": res, "censored": censored,
           "bracketing_checks": None if br is None else [{"index": i, "n_discharges": len(checks[i]), "capacity_Ah": cap[i],
                                                         "fraction_of_first": cap[i] / q0, "efc": efc[i],
                                                         "discharge_ah": [dah[k] for k in checks[i]]} for i in br],
           "units_used": len(used), "measure_a_K": a_meas, "measure_b_C": b_meas, "T_cyc_C": T_cyc}
    # where no unit is used: the first unit with both step means, whatever its EFC, to show where temperature begins
    if not used:
        for c, d in units:
            if c is None:
                continue
            cm, dm = tmean(steps[c][0], steps[c][1]), tmean(steps[d][0], steps[d][1])
            if cm is not None and dm is not None:
                out["first_unit_with_temperature"] = {"cycle_index": R[steps[d][0]][1], "efc_at_discharge_start": cum[d] / C,
                                                      "charge_mean_C": cm, "discharge_mean_C": dm}
                break
        out["units_in_lifetime_window"] = sum(1 for c, d in units if c is not None and (censored or cum[d] / C <= L))
        out["non_check_discharges"] = len(units)
    # units by hand: the first used unit and the one in the middle of the used units, with their raw samples
    out["example_units"] = []
    for which in sorted({0, len(used) // 2}) if used else []:
        u = dict(used[which])
        u["index_among_used_units"] = which
        cs, ce = steps[u["charge_step"]][0], steps[u["charge_step"]][1]
        ds, de = steps[u["discharge_step"]][0], steps[u["discharge_step"]][1]
        u["charge_samples"] = ce - cs + 1
        u["discharge_samples"] = de - ds + 1
        u["charge_first_last_time_s"] = [R[cs][0], R[ce][0]]
        u["discharge_first_last_time_s"] = [R[ds][0], R[de][0]]
        u["discharge_samples_t_T_I"] = [[R[j][0], R[j][6], R[j][2]] for j in range(ds, de + 1)]
        u["charge_first_last_T"] = [R[cs][6], R[ce][6]]
        u["rest_records"] = [[[R[j][0], R[j][6], R[j][2], R[j][3]] for j in range(steps[k][0], steps[k][1] + 1)]
                             for k in u["rest_steps"]]
        u["cycle_index_charge_discharge"] = [R[cs][1], R[ds][1]]
        u["discharge_ah"] = dah[u["discharge_step"]]
        u["charge_samples_t_T_I"] = [[R[j][0], R[j][6], R[j][2]] for j in range(cs, ce + 1)]
        out["example_units"].append(u)
    return out


def main(out_path, cids):
    res = [analyse(c) for c in cids]
    with open(out_path, "w") as fh:
        json.dump(res, fh, indent=1)
    for r in res:
        print(r["cell_id"], "L", r["L_efc"], "res", r["resolution_efc"], "a", r["measure_a_K"], "b", r["measure_b_C"],
              "Tcyc", r["T_cyc_C"], "units", r["units_used"], "dis_ah", r["discharge_ah_total"], "reg_share", r["registered_share"])


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
