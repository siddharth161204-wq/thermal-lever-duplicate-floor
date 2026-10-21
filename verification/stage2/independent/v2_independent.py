"""Independent V2 recomputation of Stage 2 (DECISIONS.md D27 to D33) and scoring (D34).

Written from the wording of DECISIONS.md D25 to D35 and predictions_stage2.md only,
without reading the pipeline's code (src/stage2, tests). Reads the raw KIT EOC files
and the cell configuration table, writes its own per-cell table and summary JSON.

Usage: python3 -I v2_independent.py EOC_DIR CELL_TABLE OUT_DIR
"""
import sys
import os
import re
import glob
import json
import itertools
import datetime as dt

import numpy as np
import pandas as pd
from scipy import stats

Q_NOM = 3.0                      # Ah, D26
R_GAS = 8.314462618              # J/(mol K), D30
SEED = 20261021                  # D31
N_MODULES = 200_000              # D31
N_CELLS_MODULE = 8               # D31
PELTIER_CUTOFF = dt.datetime(2023, 4, 23, 0, 0, tzinfo=dt.timezone.utc).timestamp()  # D30
OFFER_REG = 0.391                # K, D33 (registered)
OFFER_DYN = 1.17                 # K, D33 (post hoc, unscored)
GESLIN_MEAN = 0.030              # D33
GESLIN_P95 = 0.087               # D33


# ----------------------------------------------------------------------------- reading
def read_eoc(path):
    d = pd.read_csv(path, sep=';')
    d = d[d['timestamp_origin'] != 5].reset_index(drop=True)      # D27: never used
    return d


def checkups(d):
    """D27: one capacity per check-up (cyc_condition 2 grouped by num_cycles_checkup)."""
    cu = d[d['cyc_condition'] == 2]
    out = []
    for k, g in cu.groupby('num_cycles_checkup', sort=True):
        c = g[(g['cyc_charged'] == 0)
              & np.isclose(g['v_min_target_V'], 2.5, atol=1e-9)
              & np.isclose(g['i_dischg_min_A'], -1.0, atol=1e-9)]
        if len(c) == 0:
            continue                                   # no capacity measured in this check-up
        aq = c['delta_q_dischg_Ah'].abs()
        r = c.loc[aq.idxmax()]
        cap = abs(r['delta_q_dischg_Ah'])
        efc = (r['total_q_dischg_sum_Ah'] - cap) / Q_NOM
        out.append(dict(cu=int(k), cap_Ah=cap, efc=efc, timestamp_s=r['timestamp_s']))
    return pd.DataFrame(out).sort_values('cu').reset_index(drop=True)


def crossing(efc, val, thr):
    """First fall below thr, linear interpolation between the bracketing points.
    Returns (L, censored). Censored -> L is the last point's EFC."""
    efc = np.asarray(efc, float)
    val = np.asarray(val, float)
    below = np.nonzero(val < thr)[0]
    if len(below) == 0:
        return float(efc[-1]), True
    j = below[0]
    if j == 0:
        return float(efc[0]), False                   # cannot happen for check-ups (reference is point 0)
    x0, x1, y0, y1 = efc[j - 1], efc[j], val[j - 1], val[j]
    return float(x0 + (thr - y0) * (x1 - x0) / (y1 - y0)), False


def running_estimate_L(d, cap0, frac=0.9, nan_mode='drop', edge='strict'):
    """D28 sensitivity lifetime: centred running median (11 rows) of cap_aged_est_Ah on
    cyc_condition 1 rows, first falling below frac*cap0, interpolated linearly in EFC
    (total_q_dischg_sum_Ah / 3.0) between rows.
    nan_mode 'drop': rows with NaN estimate removed before the median (primary reading).
    nan_mode 'keep_skip': NaN rows kept in the window, ignored inside it.
    edge 'strict': median only where all eleven rows exist; 'shrink': partial windows at the ends."""
    reg = d[d['cyc_condition'] == 1]
    x = reg['total_q_dischg_sum_Ah'].to_numpy(float) / Q_NOM
    y = reg['cap_aged_est_Ah']
    if nan_mode == 'drop':
        keep = y.notna().to_numpy()
        x, y = x[keep], y[keep].reset_index(drop=True)
        minp = 11 if edge == 'strict' else 1
        med = y.rolling(11, center=True, min_periods=minp).median().to_numpy(float)
    elif nan_mode == 'keep_skip':
        y = y.reset_index(drop=True)
        minp = 11 if edge == 'strict' else 1
        if edge == 'strict':
            # all eleven rows must exist in the series; NaNs inside are skipped
            arr = y.to_numpy(float)
            med = np.full(len(arr), np.nan)
            for i in range(5, len(arr) - 5):
                w = arr[i - 5:i + 6]
                if np.isfinite(w).any():
                    med[i] = np.nanmedian(w)
        else:
            med = y.rolling(11, center=True, min_periods=1).median().to_numpy(float)
    elif nan_mode == 'keep_nan':                      # a NaN anywhere in the window -> no median
        y = y.reset_index(drop=True)
        med = y.rolling(11, center=True).median().to_numpy(float)
    else:
        raise ValueError(nan_mode)
    ok = np.isfinite(med)
    xs, ms = x[ok], med[ok]
    if len(ms) == 0:
        return np.nan, True
    return crossing(xs, ms, frac * cap0)


def cell_temperature(d, L, censored, t_first_cu, cutoff=None):
    """D29: weighted mean of (t_start+t_end)/2 by cyc_duration_s over cyc_condition 1 rows
    with duration > 0, from the first check-up up to the 90 % crossing (EFC <= L), all rows
    when censored. cutoff: drop rows with timestamp after it (D30 sensitivity)."""
    reg = d[(d['cyc_condition'] == 1) & (d['cyc_duration_s'] > 0)]
    reg = reg[reg['timestamp_s'] >= t_first_cu]
    if not censored:
        reg = reg[reg['total_q_dischg_sum_Ah'] / Q_NOM <= L]
    if cutoff is not None:
        reg = reg[reg['timestamp_s'] <= cutoff]
    w = reg['cyc_duration_s'].to_numpy(float)
    t = 0.5 * (reg['t_start_degC'].to_numpy(float) + reg['t_end_degC'].to_numpy(float))
    if w.sum() <= 0:
        return np.nan, 0
    return float(np.sum(w * t) / np.sum(w)), int(len(w))


# ----------------------------------------------------------------------------- per cell
def per_cell(eoc_dir, table):
    files = glob.glob(os.path.join(eoc_dir, 'cell_eocv2_P*.csv'))
    index = {}
    for f in files:
        m = re.match(r'cell_eocv2_P(\d{3})_(\d+)_S(\d+)_C(\d+)\.csv$', os.path.basename(f))
        if m:
            index[(int(m.group(1)), int(m.group(2)))] = (f, int(m.group(3)), int(m.group(4)))
    rows, cu_all = [], []
    for _, c in table.iterrows():
        pid, pnr = int(c['parameter_id']), int(c['parameter_nr'])
        f, sl, ch = index[(pid, pnr)]
        s_tab, c_tab = c['slave_channel'].split(':')
        assert int(s_tab[1:]) == sl and int(c_tab[1:]) == ch, (pid, pnr)
        d = read_eoc(f)
        cu = checkups(d)
        cap0 = cu['cap_Ah'].iloc[0]
        L90, cens90 = crossing(cu['efc'], cu['cap_Ah'], 0.9 * cap0)
        L80, cens80 = crossing(cu['efc'], cu['cap_Ah'], 0.8 * cap0)
        t_first = cu['timestamp_s'].iloc[0]
        T, nrows = cell_temperature(d, L90, cens90, t_first)
        Tpre, nrows_pre = cell_temperature(d, L90, cens90, t_first, cutoff=PELTIER_CUTOFF)
        Lest, cest = running_estimate_L(d, cap0, 0.9, 'drop', 'strict')
        variants = {}
        for nm, ed in [('drop', 'shrink'), ('keep_skip', 'strict'), ('keep_skip', 'shrink'), ('keep_nan', 'strict')]:
            lv, cv = running_estimate_L(d, cap0, 0.9, nm, ed)
            variants[f'L90_est_{nm}_{ed}'] = lv
            variants[f'cens_est_{nm}_{ed}'] = cv
        # timing of the lifetime window (for reference)
        reg = d[(d['cyc_condition'] == 1) & (d['cyc_duration_s'] > 0)]
        win = reg if cens90 else reg[reg['total_q_dischg_sum_Ah'] / Q_NOM <= L90]
        rows.append(dict(
            cell=c['cell'], parameter_id=pid, parameter_nr=pnr, age_temp=int(c['age_temp']),
            window=c['window'], chg=float(c['age_chg_rate']), dchg=float(c['age_dischg_rate']),
            protocol=c['protocol'], n_checkups=len(cu), cap0_Ah=cap0,
            cap_last_Ah=cu['cap_Ah'].iloc[-1], efc_last=cu['efc'].iloc[-1],
            L90_efc=L90, censored90=cens90, L80_efc=L80, censored80=cens80,
            T_C=T, n_rows_T=nrows, T_pre_peltier_C=Tpre, n_rows_T_pre=nrows_pre,
            L90_est_efc=Lest, censored_est=cest,
            window_end_utc=(dt.datetime.fromtimestamp(win['timestamp_s'].max(), dt.timezone.utc).isoformat() if len(win) else ''),
            **variants))
        for _, r in cu.iterrows():
            cu_all.append(dict(cell=c['cell'], **r.to_dict()))
    return pd.DataFrame(rows), pd.DataFrame(cu_all)


# ----------------------------------------------------------------------------- set and lever
def set_stats(cells, Lcol='L90_efc', ccol='censored90'):
    out = []
    for pid, g in cells.groupby('parameter_id', sort=True):
        u = g[~g[ccol]]
        n = len(u)
        lnL = np.log(u[Lcol].to_numpy(float))
        out.append(dict(
            parameter_id=pid, age_temp=int(g['age_temp'].iloc[0]), protocol=g['protocol'].iloc[0],
            window=g['window'].iloc[0], chg=g['chg'].iloc[0], dchg=g['dchg'].iloc[0],
            n_uncensored=n, enters=n >= 2,
            mean_lnL=float(lnL.mean()) if n else np.nan,
            var_lnL=float(lnL.var(ddof=1)) if n >= 2 else np.nan,
            T_group=float(u['T_C'].mean()) if n else np.nan,
            T_group_pre=float(u['T_pre_peltier_C'].mean()) if n else np.nan,
            lower_bound_efc=(float(g.loc[g[ccol], Lcol].min()) if (n < 2 and g[ccol].any()) else np.nan)))
    return pd.DataFrame(out)


def pooled_sigma(sets):
    e = sets[sets['n_uncensored'] >= 2]
    num = ((e['n_uncensored'] - 1) * e['var_lnL']).sum()
    den = (e['n_uncensored'] - 1).sum()
    return float(np.sqrt(num / den)), int(den), int(len(e))


def levers(sets, sigma, t_lo, t_hi):
    out = []
    for prot, g in sets.groupby('protocol', sort=False):
        a = g[g['age_temp'] == t_lo].iloc[0]
        b = g[g['age_temp'] == t_hi].iloc[0]
        det = bool(a['enters'] and b['enters'])
        r = dict(protocol=prot, window=a['window'], chg=a['chg'], dchg=a['dchg'],
                 t_lo=t_lo, t_hi=t_hi, determinable=det,
                 n_lo=int(a['n_uncensored']), n_hi=int(b['n_uncensored']))
        if det:
            dlnL = b['mean_lnL'] - a['mean_lnL']
            dT = b['T_group'] - a['T_group']
            s = -dlnL / dT
            se = sigma * np.sqrt(1 / a['n_uncensored'] + 1 / b['n_uncensored']) / dT
            # D30 sensitivities
            s_set = -dlnL / float(t_hi - t_lo)
            se_set = sigma * np.sqrt(1 / a['n_uncensored'] + 1 / b['n_uncensored']) / float(t_hi - t_lo)
            T_lo_pre = a['T_group_pre'] if t_lo == 25 else a['T_group']
            T_hi_pre = b['T_group_pre'] if t_hi == 25 else b['T_group']
            dT_pre = T_hi_pre - T_lo_pre
            s_pre = -dlnL / dT_pre
            se_pre = sigma * np.sqrt(1 / a['n_uncensored'] + 1 / b['n_uncensored']) / dT_pre
            Tk = 0.5 * (a['T_group'] + b['T_group']) + 273.15
            Ea = s * R_GAS * Tk ** 2
            r.update(T_lo=a['T_group'], T_hi=b['T_group'], dT=dT, mean_lnL_lo=a['mean_lnL'],
                     mean_lnL_hi=b['mean_lnL'], s=s, se=se, ci_lo=s - 1.96 * se, ci_hi=s + 1.96 * se,
                     s_setpoint=s_set, se_setpoint=se_set,
                     T_lo_pre=T_lo_pre, T_hi_pre=T_hi_pre, dT_pre=dT_pre, s_pre_peltier=s_pre, se_pre_peltier=se_pre,
                     Ea_J_per_mol=Ea, Ea_kJ_per_mol=Ea / 1000.0, Tmean_K=Tk,
                     module_reg=OFFER_REG * s, module_dyn=OFFER_DYN * s,
                     module_reg_below_mean=bool(OFFER_REG * s < GESLIN_MEAN),
                     module_reg_below_p95=bool(OFFER_REG * s < GESLIN_P95),
                     module_dyn_below_mean=bool(OFFER_DYN * s < GESLIN_MEAN),
                     module_dyn_below_p95=bool(OFFER_DYN * s < GESLIN_P95))
        out.append(r)
    return pd.DataFrame(out)


# ----------------------------------------------------------------------------- floor, warmer
def pair_diff(a, b):
    return abs(a - b) / (0.5 * (a + b))


def floor(cells, sigma):
    pre, allp = [], []
    for pid, g in cells.groupby('parameter_id', sort=True):
        g = g.sort_values('parameter_nr')
        u = g[~g['censored90']]
        c1 = g[g['parameter_nr'] == 1].iloc[0]
        c2 = g[g['parameter_nr'] == 2].iloc[0]
        if (not c1['censored90']) and (not c2['censored90']):
            pre.append(dict(parameter_id=pid, age_temp=int(g['age_temp'].iloc[0]), protocol=g['protocol'].iloc[0],
                            cell_a=c1['cell'], cell_b=c2['cell'], L_a=c1['L90_efc'], L_b=c2['L90_efc'],
                            diff=pair_diff(c1['L90_efc'], c2['L90_efc'])))
        for (_, x), (_, y) in itertools.combinations(list(u.iterrows()), 2):
            allp.append(dict(parameter_id=pid, age_temp=int(g['age_temp'].iloc[0]), cell_a=x['cell'], cell_b=y['cell'],
                             diff=pair_diff(x['L90_efc'], y['L90_efc'])))
    pre, allp = pd.DataFrame(pre), pd.DataFrame(allp)
    rng = np.random.default_rng(SEED)
    lnL = rng.normal(0.0, sigma, size=(N_MODULES, N_CELLS_MODULE))
    L = np.exp(lnL)
    m = L.mean(axis=1)
    shortfall = (m - L.min(axis=1)) / m
    res = dict(n_pairs=int(len(pre)), mean=float(pre['diff'].mean()), p95=float(np.percentile(pre['diff'], 95)),
               median=float(np.median(pre['diff'])), max=float(pre['diff'].max()),
               worst_of_8_median_shortfall=float(np.median(shortfall)),
               all_pairs_n=int(len(allp)), all_pairs_mean=float(allp['diff'].mean()),
               all_pairs_p95=float(np.percentile(allp['diff'], 95)),
               per_temp={int(t): dict(n=int(len(g)), mean=float(g['diff'].mean()), p95=float(np.percentile(g['diff'], 95)))
                         for t, g in pre.groupby('age_temp')})
    # alternative layouts of the same draws (for the report only)
    rng2 = np.random.default_rng(SEED)
    L2 = np.exp(rng2.normal(0.0, sigma, size=(N_CELLS_MODULE, N_MODULES))).T
    m2 = L2.mean(axis=1)
    res['worst_of_8_alt_layout_8xN'] = float(np.median((m2 - L2.min(axis=1)) / m2))
    return res, pre, allp


def warmer(cells):
    out = []
    for pid, g in cells.groupby('parameter_id', sort=True):
        u = g[~g['censored90']].sort_values('parameter_nr')
        for (_, x), (_, y) in itertools.combinations(list(u.iterrows()), 2):
            differ = bool(x['T_C'] != y['T_C'])
            if differ:
                w, c = (x, y) if x['T_C'] > y['T_C'] else (y, x)
                hit = bool(w['L90_efc'] < c['L90_efc'])
            else:
                w, c, hit = x, y, None
            out.append(dict(parameter_id=pid, age_temp=int(g['age_temp'].iloc[0]), protocol=g['protocol'].iloc[0],
                            cell_a=x['cell'], cell_b=y['cell'], T_a=x['T_C'], T_b=y['T_C'], L_a=x['L90_efc'], L_b=y['L90_efc'],
                            dT=abs(x['T_C'] - y['T_C']), T_differ=differ, warmer=w['cell'] if differ else '',
                            warmer_shorter=hit))
    wp = pd.DataFrame(out)

    def share(df):
        dd = df[df['T_differ']]
        n = int(len(dd)); k = int(dd['warmer_shorter'].sum())
        p = float(stats.binomtest(k, n, 0.5, alternative='greater').pvalue) if n else np.nan
        return dict(n_pairs=n, k_warmer_shorter=k, share=(k / n if n else np.nan), p_one_sided=p)

    res = share(wp)
    res['n_pairs_total'] = int(len(wp))
    res['per_temp'] = {int(t): share(g) for t, g in wp.groupby('age_temp')}
    return res, wp


# ----------------------------------------------------------------------------- scoring (D34)
def score(lev_hi, lev_lo, fl, wr):
    hi = lev_hi.set_index('protocol')
    lo = lev_lo.set_index('protocol')
    sc = {}
    p1set = hi[np.isclose(hi['chg'], 0.33)]
    assert len(p1set) == 6
    if p1set['determinable'].all():
        mean6 = float(p1set['s'].mean())
        ok = bool((p1set['s'] > 0).all() and 0.01 <= mean6 <= 0.06)
        sc['P1'] = dict(verdict='holds' if ok else 'fails', mean_s=mean6, all_positive=bool((p1set['s'] > 0).all()),
                        s=p1set['s'].to_dict())
    else:
        sc['P1'] = dict(verdict='not scorable', not_determinable=p1set.index[~p1set['determinable']].tolist())
    # P2
    win = {}
    for w in ['0-100', '10-100', '10-90']:
        a = hi[(hi['window'] == w) & np.isclose(hi['chg'], 1.67) & np.isclose(hi['dchg'], 1.0)].iloc[0]
        b = hi[(hi['window'] == w) & np.isclose(hi['chg'], 0.33) & np.isclose(hi['dchg'], 1.0)].iloc[0]
        c = hi[(hi['window'] == w) & np.isclose(hi['chg'], 0.33) & np.isclose(hi['dchg'], 0.33)].iloc[0]
        if a['determinable'] and b['determinable']:
            win[w] = dict(verdict='holds' if a['s'] < b['s'] else 'fails', s_167_100=float(a['s']), s_033_100=float(b['s']),
                          s_033_033=(float(c['s']) if c['determinable'] else None))
        else:
            win[w] = dict(verdict='not scorable', s_167_100=(float(a['s']) if a['determinable'] else None),
                          s_033_100=(float(b['s']) if b['determinable'] else None))
    vs = [v['verdict'] for v in win.values()]
    if 'fails' in vs:
        p2 = 'fails'
    elif all(v == 'holds' for v in vs):
        p2 = 'holds'
    else:
        p2 = 'not scorable (in part)'
    sc['P2'] = dict(verdict=p2, windows=win)
    det = hi[hi['determinable']]
    sc['P3'] = dict(verdict='holds' if ((det['s'] > 0).any() and (det['s'] < 0).any()) else 'fails',
                    n_determinable=int(len(det)), n_pos=int((det['s'] > 0).sum()), n_neg=int((det['s'] < 0).sum()))
    sc['P4'] = dict(verdict='holds' if (0.015 <= fl['mean'] <= 0.060 and 0.044 <= fl['p95'] <= 0.174) else 'fails',
                    mean=fl['mean'], p95=fl['p95'],
                    mean_in=bool(0.015 <= fl['mean'] <= 0.060), p95_in=bool(0.044 <= fl['p95'] <= 0.174))
    sc['P5'] = dict(verdict='holds' if 0.45 <= wr['share'] <= 0.65 else 'fails', share=wr['share'])
    sc['P6'] = dict(verdict='holds' if not (det['s'] > 0.077).any() else 'fails', max_s=float(det['s'].max()),
                    max_protocol=det['s'].idxmax())
    if p1set['determinable'].all():
        m = float(p1set['s'].mean())
        sc['P7'] = dict(verdict='holds' if (-0.0099 < m < 0.052) else 'fails', mean_s=m)
    else:
        sc['P7'] = dict(verdict='not scorable')
    p8set = lo[lo['chg'] >= 1.0 - 1e-9]
    assert len(p8set) == 6
    neg = [(p, bool(r['determinable'] and r['s'] < 0)) for p, r in p8set.iterrows()]
    nneg = sum(v for _, v in neg)
    sc['P8'] = dict(verdict='holds' if nneg >= 3 else 'fails', n_negative=int(nneg),
                    per_protocol={p: dict(determinable=bool(p8set.loc[p, 'determinable']),
                                          s=(float(p8set.loc[p, 's']) if p8set.loc[p, 'determinable'] else None),
                                          negative=v) for p, v in neg})
    return sc


# ----------------------------------------------------------------------------- main
def main(eoc_dir, table_path, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    table = pd.read_csv(table_path)
    table = table[(table['parameter_id'] >= 17) & (table['parameter_id'] <= 64)]
    assert len(table) == 144
    cells, cus = per_cell(eoc_dir, table)
    cells.to_csv(os.path.join(out_dir, 'my_cells.csv'), index=False)
    cus.to_csv(os.path.join(out_dir, 'my_checkups.csv'), index=False)

    sets = set_stats(cells)
    sigma, dof, nsets = pooled_sigma(sets)
    sets.to_csv(os.path.join(out_dir, 'my_sets.csv'), index=False)
    lev_hi = levers(sets, sigma, 25, 40)
    lev_lo = levers(sets, sigma, 10, 25)
    lev_hi.to_csv(os.path.join(out_dir, 'my_lever_25_40.csv'), index=False)
    lev_lo.to_csv(os.path.join(out_dir, 'my_lever_10_25.csv'), index=False)

    # secondary threshold 80 % (same machinery, T of D29 unchanged)
    sets80 = set_stats(cells, 'L80_efc', 'censored80')
    sigma80, dof80, nsets80 = pooled_sigma(sets80)
    lev80 = levers(sets80, sigma80, 25, 40)

    fl, pre, allp = floor(cells, sigma)
    pre.to_csv(os.path.join(out_dir, 'my_pairs_prespecified.csv'), index=False)
    allp.to_csv(os.path.join(out_dir, 'my_pairs_all.csv'), index=False)
    wr, wp = warmer(cells)
    wp.to_csv(os.path.join(out_dir, 'my_warmer_pairs.csv'), index=False)
    sc = score(lev_hi, lev_lo, fl, wr)

    def lev_dict(df):
        return {r['protocol']: {k: (v.item() if hasattr(v, 'item') else v) for k, v in r.items()} for _, r in df.iterrows()}

    summary = dict(
        n_cells=int(len(cells)), n_censored90=int(cells['censored90'].sum()), n_censored80=int(cells['censored80'].sum()),
        n_sets_entering=int(sets['enters'].sum()),
        sets_not_entering={int(r['parameter_id']): dict(protocol=r['protocol'], age_temp=int(r['age_temp']),
                                                        n_uncensored=int(r['n_uncensored']), lower_bound_efc=r['lower_bound_efc'])
                           for _, r in sets[~sets['enters']].iterrows()},
        pooled_sigma_lnL=sigma, pooled_dof=dof, pooled_n_sets=nsets,
        lever_25_40=lev_dict(lev_hi), lever_10_25=lev_dict(lev_lo),
        lever_25_40_L80=dict(sigma=sigma80, n_sets=nsets80, per_protocol=lev_dict(lev80)),
        floor=fl, warmer=wr, scorecard=sc,
        peltier_cutoff_unix=PELTIER_CUTOFF)
    with open(os.path.join(out_dir, 'my_summary.json'), 'w') as fh:
        json.dump(summary, fh, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    return summary


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
