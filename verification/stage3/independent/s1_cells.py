#!/usr/bin/env python3
"""Stage 3 independent check, step 1: the per-cell table for all 86 cells, rebuilt from the raw time series.

Rules from DECISIONS.md D39 to D43 and D50, read as stated in v2_independent_report.md ("Readings"). Nothing from
src/stage2, src/stage3 or tests is used or was opened. Run with `python3 -I s1_cells.py <out_dir>`.

Primary readings (alternatives are computed beside them where the wording leaves a choice):
  steps     maximal runs of one label in the ordered rows (no break at a cycle-index change)
  checks    a discharge qualifies when |mean sample current| is within 0.05C of 0.5C and its last voltage within
            0.02 V of the cut-off; consecutive qualifying discharges in the sequence of discharge steps form one check
            (any other discharge step between them splits the check); 0-100 % cells at a 0.5C discharge: every
            qualifying discharge its own check; capacity = mean ampere-hours; EFC = discharge before / rated;
            checks below half the first check's capacity dropped
  lifetime  first kept check strictly below threshold x first check capacity; linear interpolation in EFC between
            it and the previous kept check; resolution = their EFC spacing; censored at the last kept check
  units     every discharge step that is not a qualifying check discharge, with the last charge step before it;
            EFC at the discharge's start <= L (all units when censored); a unit is skipped when its charge or
            discharge step mean is missing
  (a)       mean over units of (discharge step mean - charge step mean)
  (b)       per unit, the mean of the step means of the rest steps directly before or after its charge or its
            discharge (each rest step once); cell (b) = mean over units that have one
  T_cyc     sum over the units' charge and discharge steps of step mean x step duration, over summed duration
            (duration = time of last sample - time of first sample)
"""
import json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib


def _runs(mask):
    """Start and end positions of runs of True in a boolean array."""
    d = np.diff(np.r_[0, mask.astype(np.int8), 0])
    return np.flatnonzero(d == 1), np.flatnonzero(d == -1) - 1


def process(cid, idx, keep_detail=False):
    m = lib.parse_id(cid)
    chem, C = m["chemistry"], m["rated_Ah"]
    a = lib.load(cid, idx)
    n_logged = len(a["t"])
    r, n_late, n_other = lib.ordered_rows(a)
    other_rows = None
    if n_other:
        _, _, oth = lib.row_structure(a["cyc"])
        oc = a["cyc"][oth]
        w = np.flatnonzero(oth)
        own = a["cyc"][w[0]:w[-1] + 1][~oth[w[0]:w[-1] + 1]]
        other_rows = {"n": int(oth.sum()), "cycles": [float(oc.min()), float(oc.max())],
                      "span_rows": int(w[-1] - w[0] + 1), "own_cycles_in_span": [float(own.min()), float(own.max())]}
    late_rows = None
    if n_late:
        _, late, _ = lib.row_structure(a["cyc"])
        wl = np.flatnonzero(late)
        dtl = a["t"][wl] - a["t"][wl - 1]
        late_rows = {"n": int(late.sum()), "max_time_after_previous_row_s": float(dtl.max()),
                     "min_time_after_previous_row_s": float(dtl.min())}
    t, cyc, I, V, Qc, Qd, T, Te = (r[k] for k in ("t", "cyc", "I", "V", "Qc", "Qd", "Tcell", "Tenv"))
    n = len(t)
    lab = lib.labels(I, C)
    st, en = lib.steps(lab)
    kind = lab[st]
    ns = en - st + 1
    dur = t[en] - t[st]
    csI = np.r_[0.0, np.cumsum(I)]
    imean = (csI[en + 1] - csI[st]) / ns
    vlast = V[en]
    ah_d = lib.step_ah_per_sample(Qd, st, en)
    ah_c = lib.step_ah_per_sample(Qc, st, en)
    reg_d = lib.counter_ah_registered(Qd, st, en)
    tm = lib.step_temperature(t, T, st, en)
    te = lib.step_temperature(t, Te, st, en)
    Tm = tm["mean"]
    Em = te["mean"]

    # ---- counter decreases in the ordered rows (D50 (2)) ----
    dq = np.diff(Qd)
    dec = np.flatnonzero(dq < 0) + 1                     # row of the new (lower) value
    dec_same_cycle = dec[cyc[dec] == cyc[dec - 1]]
    sid_row = np.repeat(np.arange(len(st)), ns)
    dec_same_cycle_kinds = [int(kind[sid_row[j]]) for j in dec_same_cycle]

    # ---- discharges ----
    dis = np.flatnonzero(kind == lib.DISCHARGE)
    ah_dis = ah_d[dis]
    cum_before = np.r_[0.0, np.cumsum(ah_dis)][:-1]
    efc_start = cum_before / C
    total_dis_ah = float(ah_dis.sum())
    reg_dis = reg_d[dis]
    n_differ = int(np.sum(np.abs(reg_dis - ah_dis) > 1e-9))
    # trapezoid of the step's own samples (D50 (2) bound), discharges of five samples or more
    big = dis[ns[dis] >= 5]
    if len(big):
        sid = np.repeat(np.arange(len(st)), ns)
        same = sid[:-1] == sid[1:]
        tr = np.where(same, 0.5 * (np.abs(I[:-1]) + np.abs(I[1:])) * np.diff(t), 0.0) / 3600.0
        trap = np.bincount(sid[:-1], weights=tr, minlength=len(st))
        ratio = ah_d[big] / trap[big]
        trap_ratio = [float(ratio.min()), float(ratio.max())]
    else:
        trap_ratio = [None, None]

    # ---- checks (D41) ----
    qual = (np.abs(-imean[dis] - 0.5 * C) <= 0.05 * C) & (np.abs(vlast[dis] - lib.CUTOFF[chem]) <= 0.02)
    # alternative mean current: ampere-hours over duration (only to count disagreements)
    with np.errstate(divide="ignore", invalid="ignore"):
        i_alt = ah_dis * 3600.0 / dur[dis]
    qual_alt = (np.abs(i_alt - 0.5 * C) <= 0.05 * C) & (np.abs(vlast[dis] - lib.CUTOFF[chem]) <= 0.02)
    every_own = (m["window"] == "0-100") and (m["c_discharge"] == 0.5)
    if every_own:
        groups = [[k] for k in np.flatnonzero(qual)]
    else:
        s_, e_ = _runs(qual)
        groups = [list(range(s, e + 1)) for s, e in zip(s_, e_)]
    checks = []
    for g in groups:
        checks.append({"first_discharge_step": int(dis[g[0]]), "n_discharges": len(g),
                       "capacity_Ah": float(np.mean(ah_dis[g])), "efc": float(efc_start[g[0]]),
                       "cycle_index": float(cyc[st[dis[g[0]]]]), "t_start_s": float(t[st[dis[g[0]]]])})
    q0 = checks[0]["capacity_Ah"] if checks else float("nan")
    for c in checks:
        c["dropped_half"] = bool(c["capacity_Ah"] < 0.5 * q0)
    kept = [c for c in checks if not c["dropped_half"]]

    # ---- lifetime (D42) ----
    thr_frac = lib.THRESHOLD[chem]
    thr = thr_frac * q0
    L = res = float("nan")
    censored = True
    bracket = None
    for k in range(1, len(kept)):
        if kept[k]["capacity_Ah"] < thr:
            c0, c1 = kept[k - 1], kept[k]
            L = c0["efc"] + (thr - c0["capacity_Ah"]) * (c1["efc"] - c0["efc"]) / (c1["capacity_Ah"] - c0["capacity_Ah"])
            res = c1["efc"] - c0["efc"]
            censored = False
            bracket = [k - 1, k]
            break
    censor_efc = kept[-1]["efc"] if (censored and kept) else float("nan")

    # alternative check grouping (unregistered, diagnostic): consecutive qualifying discharges where only discharge
    # steps of at most two samples lie between them
    alt = None
    if not every_own:
        tiny = ns[dis] <= 2
        qq = qual.copy()
        # bridge: a run of tiny non-qualifying discharges between two qualifying ones joins them
        groups2 = []
        cur = []
        k = 0
        while k < len(dis):
            if qual[k]:
                cur.append(k)
            elif tiny[k] and cur:
                # look ahead past tiny ones
                j = k
                while j < len(dis) and tiny[j] and not qual[j]:
                    j += 1
                if j < len(dis) and qual[j]:
                    k = j
                    continue
                groups2.append(cur)
                cur = []
                k = j
                continue
            else:
                if cur:
                    groups2.append(cur)
                cur = []
            k += 1
        if cur:
            groups2.append(cur)
        ch2 = [{"capacity_Ah": float(np.mean(ah_dis[g])), "efc": float(efc_start[g[0]])} for g in groups2]
        ch2 = [c for c in ch2 if c["capacity_Ah"] >= 0.5 * ch2[0]["capacity_Ah"]] if ch2 else []
        L2 = r2 = float("nan")
        for k in range(1, len(ch2)):
            if ch2[k]["capacity_Ah"] < thr_frac * ch2[0]["capacity_Ah"]:
                c0, c1 = ch2[k - 1], ch2[k]
                L2 = c0["efc"] + (thr_frac * ch2[0]["capacity_Ah"] - c0["capacity_Ah"]) * (c1["efc"] - c0["efc"]) / (c1["capacity_Ah"] - c0["capacity_Ah"])
                r2 = c1["efc"] - c0["efc"]
                break
        alt = {"n_checks": len(groups2), "L": L2, "resolution": r2}

    # ---- sensitivity lifetime (D42, unscored; reading S1 in the report) ----
    sens = {"applies": m["window"] == "0-100"}
    if sens["applies"]:
        cyc_full = (np.abs(-imean[dis] - m["c_discharge"] * C) <= 0.05 * C) & (np.abs(vlast[dis] - lib.CUTOFF[chem]) <= 0.02)
        if not every_own:
            cyc_full &= ~qual
        caps = ah_dis[cyc_full]
        efcs = efc_start[cyc_full]
        for ref_name, ref in (("ref_first_cycling_capacity", caps[0] if len(caps) else np.nan), ("ref_first_check", q0)):
            Ls = np.nan
            below = np.flatnonzero(caps < thr_frac * ref)
            below = below[below > 0]
            if len(below):
                k = below[0]
                Ls = efcs[k - 1] + (thr_frac * ref - caps[k - 1]) * (efcs[k] - efcs[k - 1]) / (caps[k] - caps[k - 1])
            sens[ref_name] = float(Ls)
        sens["n_cycling_discharges"] = int(cyc_full.sum())

    # ---- units (D43, D50 (5), (6)) ----
    is_chk_dis = np.zeros(len(st), bool)
    is_chk_dis[dis[qual]] = True
    pos_charge = np.where(kind == lib.CHARGE, np.arange(len(st)), -1)
    last_charge = np.maximum.accumulate(pos_charge)            # last charge step at or before each step
    unit_d = dis[~qual]
    unit_c = np.where(unit_d > 0, last_charge[np.maximum(unit_d - 1, 0)], -1)
    unit_efc = efc_start[~qual]
    in_life = np.ones(len(unit_d), bool) if censored else (unit_efc <= L)
    has_c = unit_c >= 0
    cm = np.where(has_c, Tm[np.maximum(unit_c, 0)], np.nan)
    dm = Tm[unit_d]
    ok = in_life & has_c & ~np.isnan(cm) & ~np.isnan(dm)
    # adjacent rest steps (reading B1)
    nst = len(st)

    def rest_at(j):
        return (j >= 0) & (j < nst) & (kind[np.clip(j, 0, nst - 1)] == lib.REST)

    cand = np.stack([unit_c - 1, unit_c + 1, unit_d - 1, unit_d + 1], axis=1)
    cand_ok = rest_at(cand) & has_c[:, None]
    # de-duplicate (unit_c+1 can equal unit_d-1)
    dup = np.zeros_like(cand_ok)
    for i in range(4):
        for j in range(i):
            dup[:, i] |= (cand[:, i] == cand[:, j]) & cand_ok[:, j]
    cand_ok &= ~dup
    rest_means = np.where(cand_ok, Tm[np.clip(cand, 0, nst - 1)], np.nan)
    with np.errstate(invalid="ignore"):
        cnt = np.sum(~np.isnan(rest_means), axis=1)
        unit_b = np.where(cnt > 0, np.nansum(rest_means, axis=1) / np.maximum(cnt, 1), np.nan)
    swing = dm - cm
    n_units_total = int(len(unit_d))
    n_units_life = int(in_life.sum())
    n_ok = int(ok.sum())
    meas_a = float(np.mean(swing[ok])) if n_ok else float("nan")
    okb = ok & ~np.isnan(unit_b)
    meas_b = float(np.mean(unit_b[okb])) if okb.any() else float("nan")
    if n_ok:
        cc, dd = unit_c[ok], unit_d[ok]
        num = np.sum(Tm[cc] * dur[cc] + Tm[dd] * dur[dd])
        den = np.sum(dur[cc] + dur[dd])
        T_cyc = float(num / den)
        numE = np.sum((Tm[cc] - Em[cc]) * dur[cc] + (Tm[dd] - Em[dd]) * dur[dd])
        cell_minus_env = float(numE / den)
    else:
        T_cyc = cell_minus_env = float("nan")

    # alternatives for the unit definition (diagnostic, unscored readings)
    def alt_units(mask):
        o2 = ok & mask
        if not o2.any():
            return {"n": 0, "a": float("nan"), "T_cyc": float("nan"), "b": float("nan")}
        cc, dd = unit_c[o2], unit_d[o2]
        ob = o2 & ~np.isnan(unit_b)
        return {"n": int(o2.sum()), "a": float(np.mean(swing[o2])),
                "T_cyc": float(np.sum(Tm[cc] * dur[cc] + Tm[dd] * dur[dd]) / np.sum(dur[cc] + dur[dd])),
                "b": float(np.mean(unit_b[ob])) if ob.any() else float("nan")}

    prev_nonrest_is_charge = np.zeros(len(unit_d), bool)
    nonrest = np.flatnonzero(kind != lib.REST)
    p = np.searchsorted(nonrest, unit_d) - 1
    prev_nonrest_is_charge = (p >= 0) & (kind[nonrest[np.maximum(p, 0)]] == lib.CHARGE)
    nominal = np.abs(-imean[unit_d] - m["c_discharge"] * C) <= 0.05 * C
    multi = ns[unit_d] >= 3
    # first unit per charge only (each charge paired with the first non-check discharge after it)
    first_per_charge = np.r_[True, unit_c[1:] != unit_c[:-1]] if len(unit_c) else np.zeros(0, bool)
    alts = {"U2_nominal_rate_only": alt_units(nominal), "U3_previous_nonrest_step_is_charge": alt_units(prev_nonrest_is_charge),
            "U4_three_samples_or_more": alt_units(multi), "U5_first_discharge_after_each_charge": alt_units(first_per_charge)}
    # alternative (b) readings
    # B2: only the rest steps after the charge and after the discharge (the two switches of the unit)
    candB2 = np.stack([unit_c + 1, unit_d + 1], axis=1)
    okB2 = rest_at(candB2) & has_c[:, None]
    okB2[:, 1] &= ~((candB2[:, 1] == candB2[:, 0]) & okB2[:, 0])
    rmB2 = np.where(okB2, Tm[np.clip(candB2, 0, nst - 1)], np.nan)
    with np.errstate(invalid="ignore"):
        cB2 = np.sum(~np.isnan(rmB2), axis=1)
        bB2 = np.where(cB2 > 0, np.nansum(rmB2, axis=1) / np.maximum(cB2, 1), np.nan)
    okB2u = ok & ~np.isnan(bB2)
    meas_b_B2 = float(np.mean(bB2[okB2u])) if okB2u.any() else float("nan")

    out = {
        **m,
        "rows_logged": n_logged, "late_records": n_late, "rows_other_record_dropped": n_other, "rows_ordered": n,
        "late_detail": late_rows, "other_detail": other_rows,
        "steps_total": int(len(st)), "steps_charge": int(np.sum(kind == lib.CHARGE)),
        "steps_discharge": int(len(dis)), "steps_rest": int(np.sum(kind == lib.REST)),
        "discharge_ah_total": total_dis_ah, "charge_ah_total": float(ah_c[kind == lib.CHARGE].sum()),
        "efc_total": total_dis_ah / C,
        "registered_counter_ah_total": float(reg_dis.sum()),
        "registered_counter_share": float(reg_dis.sum() / total_dis_ah) if total_dis_ah else float("nan"),
        "discharges_registered_rule_differs": n_differ,
        "per_sample_over_trapezoid_min_max": trap_ratio,
        "discharge_counter_decreases": int(len(dec)),
        "discharge_counter_decreases_same_cycle": int(len(dec_same_cycle)),
        "discharge_counter_decreases_same_cycle_step_kinds": dec_same_cycle_kinds,
        "qualifying_discharges": int(qual.sum()),
        "qualifying_discharges_alt_current": int(qual_alt.sum()),
        "checks": int(len(checks)), "checks_dropped_half": int(sum(c["dropped_half"] for c in checks)),
        "first_check_Ah": q0, "threshold": thr_frac,
        "L_efc": L, "resolution_efc": res, "censored": censored, "censor_efc": censor_efc,
        "bracket_checks": bracket,
        "alt_grouping_bridge_tiny": alt,
        "sensitivity": sens,
        "units_total": n_units_total, "units_in_life": n_units_life, "units_used": n_ok,
        "units_used_b": int(okb.sum()),
        "measure_a_K": meas_a, "measure_b_C": meas_b, "T_cyc_C": T_cyc, "cell_minus_env_K": cell_minus_env,
        "measure_b_B2_C": meas_b_B2,
        "unit_alternatives": alts,
        "zero_duration_pairs_steps": int(np.sum(tm["zero_den"])),
        "steps_partial_temperature": int(np.sum((tm["nvalid"] > 0) & (tm["nvalid"] < ns))),
    }
    out["checks_list"] = checks
    if keep_detail:
        out["_detail"] = {"r": r, "st": st, "en": en, "kind": kind, "Tm": Tm, "dur": dur, "imean": imean,
                          "ah_d": ah_d, "ah_c": ah_c, "dis": dis, "qual": qual, "unit_c": unit_c, "unit_d": unit_d,
                          "unit_efc": unit_efc, "ok": ok, "unit_b": unit_b, "cand": cand, "cand_ok": cand_ok,
                          "efc_start": efc_start}
    return out


def main(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    idx = lib.index()
    rows = []
    for cid in sorted(idx):
        o = process(cid, idx)
        rows.append(o)
        print(f"{cid:40s} L={o['L_efc']:.3f} res={o['resolution_efc']:.3f} cens={o['censored']} a={o['measure_a_K']:.4f} "
              f"b={o['measure_b_C']:.4f} Tcyc={o['T_cyc_C']:.4f} units={o['units_used']}", flush=True)
    with open(os.path.join(out_dir, "cells_full.json"), "w") as fh:
        json.dump(rows, fh, indent=0, default=float)


if __name__ == "__main__":
    main(sys.argv[1])
