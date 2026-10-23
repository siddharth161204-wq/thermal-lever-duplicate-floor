#!/usr/bin/env python3
"""Stage 3 independent check, step 1: the per-cell table for all 86 cells, rebuilt from the raw time series.

Rules from DECISIONS.md D39 to D43 and D50. Nothing from src/stage2, src/stage3 or tests is used or was opened.
Usage: python3 -I s1_cells.py <out_dir> [option=value ...]

Every reading the text leaves open is an option; the defaults are the primary readings stated in the report:
  check_v     two_sided   last voltage within 0.02 V of the cut-off (|V - cut| <= 0.02); one_sided: V <= cut + 0.02
  grouping    G1          consecutive qualifying discharges in the sequence of discharge steps form one check (any
                          other discharge step between them splits it); G2: discharge steps of one or two samples
                          between qualifying discharges do not split a check
  units       U1          every discharge step that is not a qualifying (check) discharge, with the last charge step
                          before it; U2: only those at the cell's cycling rate (within 0.05C); U3: only those whose
                          previous non-rest step is a charge; U5: only the first such discharge after each charge
  zero_dt     sample_mean a step whose valid pairs have zero summed duration takes the mean of its valid samples;
                          nan: its mean is missing
  rests       B1          rest steps directly before or after the unit's charge or discharge, each once; B2: only the
                          rest steps directly after the charge and directly after the discharge
"""
import json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

DEFAULTS = {"check_v": "two_sided", "grouping": "G1", "units": "U1", "zero_dt": "sample_mean", "rests": "B1"}


def _runs(mask):
    d = np.diff(np.r_[0, mask.astype(np.int8), 0])
    return np.flatnonzero(d == 1), np.flatnonzero(d == -1) - 1


def check_groups(qual, tiny, every_own, grouping):
    """Lists of discharge positions (in the sequence of discharge steps) forming each check."""
    if every_own:
        return [[k] for k in np.flatnonzero(qual)]
    if grouping == "G1":
        s_, e_ = _runs(qual)
        return [list(range(s, e + 1)) for s, e in zip(s_, e_)]
    groups, cur = [], []
    k, n = 0, len(qual)
    while k < n:
        if qual[k]:
            cur.append(k)
            k += 1
        elif tiny[k] and cur:
            j = k
            while j < n and tiny[j] and not qual[j]:
                j += 1
            if j < n and qual[j]:
                k = j
                continue
            groups.append(cur)
            cur = []
            k = j
        else:
            if cur:
                groups.append(cur)
            cur = []
            k += 1
    if cur:
        groups.append(cur)
    return groups


def lifetime(checks, thr_frac):
    q0 = checks[0]["capacity_Ah"] if checks else float("nan")
    for c in checks:
        c["dropped_half"] = bool(c["capacity_Ah"] < 0.5 * q0)
    kept = [c for c in checks if not c["dropped_half"]]
    thr = thr_frac * q0
    for k in range(1, len(kept)):
        if kept[k]["capacity_Ah"] < thr:
            c0, c1 = kept[k - 1], kept[k]
            L = c0["efc"] + (thr - c0["capacity_Ah"]) * (c1["efc"] - c0["efc"]) / (c1["capacity_Ah"] - c0["capacity_Ah"])
            return {"L_efc": L, "resolution_efc": c1["efc"] - c0["efc"], "censored": False, "censor_efc": float("nan"),
                    "bracket": [checks.index(c0), checks.index(c1)], "first_check_Ah": q0}
    return {"L_efc": float("nan"), "resolution_efc": float("nan"), "censored": True,
            "censor_efc": kept[-1]["efc"] if kept else float("nan"), "bracket": None, "first_check_Ah": q0}


def process(cid, idx, opt=None, keep_detail=False):
    opt = dict(DEFAULTS, **(opt or {}))
    m = lib.parse_id(cid)
    chem, C = m["chemistry"], m["rated_Ah"]
    a = lib.load(cid, idx)
    n_logged = len(a["t"])
    r, n_late, n_other = lib.ordered_rows(a)
    t, cyc, I, V, Qc, Qd, T, Te = (r[k] for k in ("t", "cyc", "I", "V", "Qc", "Qd", "Tcell", "Tenv"))
    n = len(t)
    lab = lib.labels(I, C)
    st, en = lib.steps(lab)
    kind = lab[st]
    nst = len(st)
    ns = en - st + 1
    dur = t[en] - t[st]
    csI = np.r_[0.0, np.cumsum(I)]
    imean = (csI[en + 1] - csI[st]) / ns          # arithmetic mean of the step's sample currents
    vlast = V[en]
    ah_d = lib.step_ah_per_sample(Qd, st, en)
    ah_c = lib.step_ah_per_sample(Qc, st, en)
    reg_d = lib.counter_ah_registered(Qd, st, en)
    Tm = lib.step_temperature(t, T, st, en, zero_dt=opt["zero_dt"])["mean"]
    Em = lib.step_temperature(t, Te, st, en, zero_dt=opt["zero_dt"])["mean"]

    out = {**m, "options": opt, "rows_logged": n_logged, "late_records": n_late, "rows_other_record_dropped": n_other,
           "rows_ordered": n, "steps_total": int(nst), "steps_charge": int(np.sum(kind == lib.CHARGE)),
           "steps_discharge": int(np.sum(kind == lib.DISCHARGE)), "steps_rest": int(np.sum(kind == lib.REST))}

    # structure details (D50 (1), (2))
    if n_late or n_other:
        _, late, oth = lib.row_structure(a["cyc"])
        if n_late:
            wl = np.flatnonzero(late)
            # lag behind the first logged record of the next cycle (the record before it in logged order is in the
            # next cycle; the first record of that cycle is found by searching backwards)
            lag = []
            for j in wl:
                k = j - 1
                while k > 0 and a["cyc"][k - 1] == a["cyc"][j - 1]:
                    k -= 1
                lag.append(a["t"][j] - a["t"][k])
            out["late_detail"] = {"max_lag_behind_first_record_of_next_cycle_s": float(max(lag)),
                                  "min_lag_s": float(min(lag))}
        if n_other:
            w = np.flatnonzero(oth)
            own = a["cyc"][w[0]:w[-1] + 1][~oth[w[0]:w[-1] + 1]]
            out["other_detail"] = {"n": int(oth.sum()), "cycles": [float(a["cyc"][oth].min()), float(a["cyc"][oth].max())],
                                   "span_rows": int(w[-1] - w[0] + 1), "own_cycles_in_span": [float(own.min()), float(own.max())]}
    dq = np.diff(Qd)
    dec = np.flatnonzero(dq < 0) + 1
    dec_same = dec[cyc[dec] == cyc[dec - 1]]
    sid_row = np.repeat(np.arange(nst), ns)
    out["discharge_counter_decreases"] = int(len(dec))
    out["discharge_counter_decreases_same_cycle"] = int(len(dec_same))
    out["discharge_counter_decreases_same_cycle_step_kinds"] = [int(kind[sid_row[j]]) for j in dec_same]

    # discharges and ampere-hours (D40, D50 (2))
    dis = np.flatnonzero(kind == lib.DISCHARGE)
    ah_dis = ah_d[dis]
    efc_start = np.r_[0.0, np.cumsum(ah_dis)][:-1] / C
    total = float(ah_dis.sum())
    out.update({"discharge_ah_total": total, "charge_ah_total": float(ah_c[kind == lib.CHARGE].sum()), "efc_total": total / C,
                "registered_counter_ah_total": float(reg_d[dis].sum()),
                "registered_counter_share": float(reg_d[dis].sum() / total) if total else float("nan"),
                "discharges_registered_rule_differs": int(np.sum(np.abs(reg_d[dis] - ah_dis) > 1e-9))})
    big = dis[ns[dis] >= 5]
    if len(big):
        same = sid_row[:-1] == sid_row[1:]
        tr = np.where(same, 0.5 * (np.abs(I[:-1]) + np.abs(I[1:])) * np.diff(t), 0.0) / 3600.0
        trap = np.bincount(sid_row[:-1], weights=tr, minlength=nst)
        ratio = ah_d[big] / trap[big]
        out["per_sample_over_trapezoid_min_max"] = [float(ratio.min()), float(ratio.max())]
    else:
        out["per_sample_over_trapezoid_min_max"] = [None, None]

    # checks (D41)
    vok = (np.abs(vlast[dis] - lib.CUTOFF[chem]) <= 0.02) if opt["check_v"] == "two_sided" else (vlast[dis] <= lib.CUTOFF[chem] + 0.02)
    qual = (np.abs(-imean[dis] - 0.5 * C) <= 0.05 * C) & vok
    every_own = (m["window"] == "0-100") and (m["c_discharge"] == 0.5)
    groups = check_groups(qual, ns[dis] <= 2, every_own, opt["grouping"])
    checks = [{"first_discharge_step": int(dis[g[0]]), "n_discharges": len(g), "capacity_Ah": float(np.mean(ah_dis[g])),
               "efc": float(efc_start[g[0]]), "cycle_index": float(cyc[st[dis[g[0]]]]), "t_start_s": float(t[st[dis[g[0]]]]),
               "discharge_ah": [float(x) for x in ah_dis[g]]} for g in groups]
    lt = lifetime(checks, lib.THRESHOLD[chem])
    out.update({"qualifying_discharges": int(qual.sum()), "checks": len(checks),
                "checks_dropped_half": int(sum(c["dropped_half"] for c in checks)), "threshold": lib.THRESHOLD[chem], **lt})
    out["checks_list"] = checks
    L, censored = lt["L_efc"], lt["censored"]

    # sensitivity lifetime (D42; unscored): per-cycle capacity at the cycling rate = full discharges (within 0.05C of the
    # cycling rate, last voltage as for checks) that are not check discharges (all qualifying discharges for 0.5C cells),
    # threshold x first check capacity (primary) or x first cycling capacity, linear interpolation in EFC
    sens = {"applies": m["window"] == "0-100"}
    if sens["applies"]:
        cf = (np.abs(-imean[dis] - m["c_discharge"] * C) <= 0.05 * C) & vok
        if not every_own:
            cf &= ~qual
        caps, efcs = ah_dis[cf], efc_start[cf]
        for name, ref in (("L_sens_ref_first_check", lt["first_check_Ah"]), ("L_sens_ref_first_cycling", caps[0] if len(caps) else np.nan)):
            Ls, rs = float("nan"), float("nan")
            below = np.flatnonzero(caps < lib.THRESHOLD[chem] * ref)
            below = below[below > 0]
            if len(below):
                k = below[0]
                thr = lib.THRESHOLD[chem] * ref
                Ls = float(efcs[k - 1] + (thr - caps[k - 1]) * (efcs[k] - efcs[k - 1]) / (caps[k] - caps[k - 1]))
                rs = float(efcs[k] - efcs[k - 1])
            sens[name] = Ls
            sens[name + "_spacing"] = rs
        sens["n_cycling_discharges"] = int(cf.sum())
    out["sensitivity"] = sens

    # units (D43, D50 (4)-(6))
    pos_charge = np.where(kind == lib.CHARGE, np.arange(nst), -1)
    last_charge = np.maximum.accumulate(pos_charge)
    ud = dis[~qual]
    uc = np.where(ud > 0, last_charge[np.maximum(ud - 1, 0)], -1)
    uefc = efc_start[~qual]
    sel = np.ones(len(ud), bool)
    if opt["units"] == "U2":
        sel = np.abs(-imean[ud] - m["c_discharge"] * C) <= 0.05 * C
    elif opt["units"] == "U3":
        nonrest = np.flatnonzero(kind != lib.REST)
        p = np.searchsorted(nonrest, ud) - 1
        sel = (p >= 0) & (kind[nonrest[np.maximum(p, 0)]] == lib.CHARGE)
    elif opt["units"] == "U5":
        sel = np.r_[True, uc[1:] != uc[:-1]] if len(uc) else np.zeros(0, bool)
    ud, uc, uefc = ud[sel], uc[sel], uefc[sel]
    in_life = np.ones(len(ud), bool) if censored else (uefc <= L)
    has_c = uc >= 0
    cm = np.where(has_c, Tm[np.maximum(uc, 0)], np.nan)
    dm = Tm[ud]
    ok = in_life & has_c & ~np.isnan(cm) & ~np.isnan(dm)

    def rest_at(j):
        return (j >= 0) & (j < nst) & (kind[np.clip(j, 0, nst - 1)] == lib.REST)

    cand = (np.stack([uc - 1, uc + 1, ud - 1, ud + 1], axis=1) if opt["rests"] == "B1" else np.stack([uc + 1, ud + 1], axis=1))
    if len(ud) == 0:
        cand = cand.reshape(0, cand.shape[1] if cand.ndim == 2 else 4)
    cok = rest_at(cand) & has_c[:, None]
    dup = np.zeros_like(cok)
    for i in range(cand.shape[1]):
        for j in range(i):
            dup[:, i] |= (cand[:, i] == cand[:, j]) & cok[:, j]
    cok &= ~dup
    rmeans = np.where(cok, Tm[np.clip(cand, 0, nst - 1)], np.nan)
    with np.errstate(invalid="ignore"):
        cnt = np.sum(~np.isnan(rmeans), axis=1)
        ub = np.where(cnt > 0, np.nansum(rmeans, axis=1) / np.maximum(cnt, 1), np.nan)
    swing = dm - cm
    okb = ok & ~np.isnan(ub)
    out.update({"units_total": int(len(ud)), "units_in_life": int(in_life.sum()), "units_used": int(ok.sum()),
                "units_used_b": int(okb.sum()),
                "measure_a_K": float(np.mean(swing[ok])) if ok.any() else float("nan"),
                "measure_b_C": float(np.mean(ub[okb])) if okb.any() else float("nan")})
    if ok.any():
        cc, dd = uc[ok], ud[ok]
        den = np.sum(dur[cc] + dur[dd])
        out["T_cyc_C"] = float(np.sum(Tm[cc] * dur[cc] + Tm[dd] * dur[dd]) / den)
        out["cell_minus_env_K"] = float(np.sum((Tm[cc] - Em[cc]) * dur[cc] + (Tm[dd] - Em[dd]) * dur[dd]) / den)
        out["first_unit_cycle"] = float(cyc[st[dd[0]]])
        out["last_unit_efc"] = float(uefc[ok][-1])
    else:
        out["T_cyc_C"] = out["cell_minus_env_K"] = float("nan")
        out["first_unit_cycle"] = out["last_unit_efc"] = None
    # units in the lifetime window that were skipped for a missing step mean
    out["units_in_life_skipped_missing_mean"] = int(np.sum(in_life & has_c & (np.isnan(cm) | np.isnan(dm))))
    if keep_detail:
        out["_detail"] = {"r": r, "st": st, "en": en, "kind": kind, "Tm": Tm, "dur": dur, "imean": imean, "ah_d": ah_d,
                          "dis": dis, "qual": qual, "uc": uc, "ud": ud, "uefc": uefc, "ok": ok, "ub": ub, "cand": cand,
                          "cok": cok, "efc_start": efc_start}
    return out


def main(out_dir, opt):
    os.makedirs(out_dir, exist_ok=True)
    idx = lib.index()
    rows = [process(cid, idx, opt) for cid in sorted(idx)]
    tag = "_".join(f"{k}-{v}" for k, v in sorted(opt.items())) or "primary"
    path = os.path.join(out_dir, f"cells_{tag}.json")
    with open(path, "w") as fh:
        json.dump(rows, fh, indent=0, default=float)
    print(path)


if __name__ == "__main__":
    o = dict(kv.split("=", 1) for kv in sys.argv[2:])
    main(sys.argv[1], o)
