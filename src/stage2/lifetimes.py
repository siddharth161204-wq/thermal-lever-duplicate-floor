#!/usr/bin/env python3
"""Stage 2: from a cell's step table to its capacity checks, its lifetime and its two temperature measures.

Every rule here is one written in DECISIONS.md (Stage 2) before the computation it governs; the parameters are passed
in from the decisions file's values in analysis.py, not chosen here.

Checks. A discharge step is check-like when its mean current is within tol_c of 0.5C (rated capacity) and it ends
at the lower voltage cut-off (last voltage within tol_v of it). Consecutive check-like discharges with no other
discharge step between them form one capacity check. For a cell whose own cycling is a full 0.5C discharge
(0 to 100 % window, 0.5C discharge), every check-like discharge is its own check. A check's capacity is the mean of
its discharges' ampere-hours; its equivalent full cycles (EFC) are the cumulative discharge ampere-hours of the cell,
all discharge steps included, before its first discharge, over the rated capacity.

Lifetime. EFC at which the check capacity first falls below a fraction f of the first check's capacity, by linear
interpolation between the two bracketing checks; the EFC spacing of those two checks is the resolution. A cell
whose checks never fall below is censored at its last check.

Sensitivity lifetime (0 to 100 % windows). The same fraction of the median of the first ten cycling discharges after
the first check, applied to a centred running median of eleven cycling discharges at the cycling rate, with linear
interpolation between consecutive cycling discharges.

Temperature (only after the predictions are public). For each cycling discharge that is not part of a check (every
discharge, for a cell whose cycling is itself a full 0.5C discharge), its unit is the charge step before it and the discharge itself, with the rest steps that touch either (the rests
between the previous non-rest step and the charge, between the charge and the discharge, and between the discharge
and the next non-rest step). The unit's rise is the duration-weighted mean of the step-mean skin temperatures over
its charge and discharge minus the duration-weighted mean over its rests. Measure (a), the self-heating rise, is the
mean rise over the units up to the cell's lifetime (all units of a censored cell); measure (b), the rest-period skin
temperature, is the mean over the same units of their rest temperature.
"""
import numpy as np
import pandas as pd


def discharges(steps, rated_Ah, v_cut_low, tol_v, tol_c, every_discharge_is_check, ah_col="ah"):
    d = steps[steps["kind"] == "discharge"].copy()
    d["crate"] = d["mean_abs_I_A"] / rated_Ah
    d["full"] = d["v_last"] <= v_cut_low + tol_v
    d["check_like"] = d["full"] & ((d["crate"] - 0.5).abs() <= tol_c)
    d["cum_ah_before"] = d[ah_col].cumsum().shift(fill_value=0.0)
    d["efc_before"] = d["cum_ah_before"] / rated_Ah
    cid, k, prev_check = [], -1, False
    for is_c in d["check_like"].to_numpy():
        if is_c and (every_discharge_is_check or not prev_check):
            k += 1
        cid.append(k if is_c else -1)
        prev_check = bool(is_c)
    d["check_id"] = cid
    return d


def checks_table(d, ah_col="ah"):
    c = d[d["check_id"] >= 0]
    if c.empty:
        return pd.DataFrame(columns=["check_id", "efc", "t_start_s", "q_Ah", "n_discharges", "q_each"])
    g = c.groupby("check_id", sort=True)
    out = pd.DataFrame({"efc": g["efc_before"].first(), "t_start_s": g["t_start_s"].first(), "q_Ah": g[ah_col].mean(),
                        "n_discharges": g[ah_col].size(),
                        "q_each": g[ah_col].apply(lambda s: ";".join(f"{x:.5f}" for x in s))}).reset_index()
    return out


def crossing(efc, q, q_thr):
    """First crossing of q below q_thr, interpolated linearly in efc; returns (L, resolution, k) or (nan, nan, None)."""
    efc, q = np.asarray(efc, float), np.asarray(q, float)
    for k in range(1, q.size):
        if q[k] < q_thr:
            if q[k - 1] >= q_thr:
                w = (q[k - 1] - q_thr) / (q[k - 1] - q[k])
                return efc[k - 1] + w * (efc[k] - efc[k - 1]), efc[k] - efc[k - 1], k
            return efc[k], np.nan, k
    return np.nan, np.nan, None


def check_lifetime(ch, frac):
    if len(ch) < 2:
        return {"L_efc": np.nan, "resolution_efc": np.nan, "censored": True,
                "efc_last_check": float(ch["efc"].iloc[-1]) if len(ch) else np.nan,
                "q_first_Ah": float(ch["q_Ah"].iloc[0]) if len(ch) else np.nan, "n_checks": len(ch)}
    q0 = float(ch["q_Ah"].iloc[0])
    L, res, k = crossing(ch["efc"], ch["q_Ah"], frac * q0)
    return {"L_efc": float(L), "resolution_efc": float(res), "censored": k is None,
            "efc_last_check": float(ch["efc"].iloc[-1]), "q_first_Ah": q0,
            "q_last_frac": float(ch["q_Ah"].iloc[-1] / q0), "q_min_frac": float(ch["q_Ah"].min() / q0),
            "n_checks": len(ch), "efc_first_check": float(ch["efc"].iloc[0]),
            "bracket": [float(ch["efc"].iloc[k - 1]), float(ch["efc"].iloc[k])] if k else None}


def sensitivity_lifetime(d, frac, c_dis, tol_c, ah_col="ah", n_ref=10, window=11):
    """Lifetime from the per-cycle discharge capacity at the cycling rate (0 to 100 % windows only)."""
    first = d.loc[d["check_id"] >= 0, "efc_before"]
    if first.empty:
        return {"L_sens_efc": np.nan, "censored_sens": True}
    cyc = d[(d["check_id"] < 0) & d["full"] & ((d["crate"] - c_dis).abs() <= tol_c * max(1.0, c_dis / 0.5))
            & (d["efc_before"] >= first.iloc[0])]
    if c_dis == 0.5:                                      # every full 0.5C discharge is a check for these cells
        cyc = d[(d["check_id"] >= 0)]
    if len(cyc) < n_ref + window:
        return {"L_sens_efc": np.nan, "censored_sens": True, "n_cycling": len(cyc)}
    q = cyc[ah_col].to_numpy()
    efc = cyc["efc_before"].to_numpy()
    ref = float(np.median(q[:n_ref]))
    rm = pd.Series(q).rolling(window, center=True, min_periods=1).median().to_numpy()
    L, res, k = crossing(efc, rm, frac * ref)
    return {"L_sens_efc": float(L), "resolution_sens_efc": float(res), "censored_sens": k is None,
            "q_ref_sens_Ah": ref, "n_cycling": len(cyc), "efc_last_cycling": float(efc[-1])}


def temperature_measures(steps, d, L_efc, efc_last, rest_baseline="mean", every_discharge_is_check=False):
    """Measures (a) and (b) over the cycling units up to the lifetime (DECISIONS.md, Stage 2). Needs the temperature pass.

    rest_baseline "mean": the duration-weighted mean of the rest steps' time-weighted mean temperatures (the brief's
    wording). "end": the mean of each adjacent rest's last sample, the point of each rest nearest the chamber."""
    kinds = steps["kind"].to_numpy()
    tmean = steps["t_cell_mean"].to_numpy()
    tlast = steps["t_cell_last"].to_numpy()
    dur = steps["duration_s"].to_numpy()
    idx_of = {s: i for i, s in enumerate(steps["step"].to_numpy())}
    limit = L_efc if np.isfinite(L_efc) else efc_last
    units = []
    for _, r in (d if every_discharge_is_check else d[d["check_id"] < 0]).iterrows():
        if r["efc_before"] > limit:
            break
        i_d = idx_of[r["step"]]
        j = i_d - 1
        while j >= 0 and kinds[j] == "rest":
            j -= 1
        if j < 0 or kinds[j] != "charge":
            continue
        i_c = j
        rest_idx = list(range(i_c + 1, i_d))                    # between charge and discharge
        j = i_c - 1
        while j >= 0 and kinds[j] == "rest":
            rest_idx.append(j); j -= 1                          # before the charge
        j = i_d + 1
        while j < len(kinds) and kinds[j] == "rest":
            rest_idx.append(j); j += 1                          # after the discharge
        cyc_idx = [i_c, i_d]
        w_c, w_r = dur[cyc_idx], dur[rest_idx]
        if not rest_idx or w_r.sum() <= 0 or w_c.sum() <= 0:
            continue
        t_cyc = float(np.sum(tmean[cyc_idx] * w_c) / w_c.sum())
        if rest_baseline == "end":
            t_rest = float(np.mean(tlast[rest_idx]))
        else:
            t_rest = float(np.sum(tmean[rest_idx] * w_r) / w_r.sum())
        if not (np.isfinite(t_cyc) and np.isfinite(t_rest)):
            continue
        units.append((r["efc_before"], t_cyc, t_rest, t_cyc - t_rest, w_r.sum()))
    if not units:
        return {"rise_K": np.nan, "rest_T_C": np.nan, "n_units": 0}
    u = np.array(units)
    return {"rise_K": float(u[:, 3].mean()), "rest_T_C": float(u[:, 2].mean()), "cycling_T_C": float(u[:, 1].mean()),
            "n_units": int(len(u)), "rise_sd_K": float(u[:, 3].std(ddof=1)) if len(u) > 1 else np.nan,
            "rest_duration_median_s": float(np.median(u[:, 4]))}
