"""Compare the independent V2 recomputation (frozen outputs of v2_independent.py) with the
pipeline's Stage 2 outputs, and score P1 to P8 from the independent numbers (D34).

Usage: python3 -I compare.py MY_DIR RESULTS_DIR OUT_JSON
MY_DIR holds my_*.csv / my_summary.json written by v2_independent.py before any pipeline
output was opened. The 'extras' block (80 % floor, check-up resolution at the crossing,
median within-pair dT) is computed here from MY_DIR after the pipeline's output structure
was seen; it uses the same definitions and is labelled as such.
"""
import sys
import os
import json
import itertools

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2_independent as v2  # noqa: E402


def diff_stats(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    both_nan = np.isnan(a) & np.isnan(b)
    one_nan = np.isnan(a) ^ np.isnan(b)
    ok = ~(np.isnan(a) | np.isnan(b))
    ad = np.abs(a[ok] - b[ok])
    rd = ad / np.maximum(np.abs(b[ok]), 1e-300)
    i = int(np.argmax(ad)) if ad.size else None
    j = int(np.argmax(rd)) if rd.size else None
    return dict(n_compared=int(ok.sum()), n_nan_both=int(both_nan.sum()), n_nan_mismatch=int(one_nan.sum()),
                max_abs=(float(ad.max()) if ad.size else None), max_rel=(float(rd.max()) if rd.size else None),
                n_exact=int((ad == 0).sum()), argmax_abs=i, argmax_rel=j)


def num(x, y):
    x = float(x); y = float(y)
    return dict(mine=x, pipeline=y, abs_diff=abs(x - y), rel_diff=(abs(x - y) / abs(y) if y != 0 else (0.0 if x == y else float('inf'))))


def main(my_dir, res_dir, out_json):
    mine = pd.read_csv(os.path.join(my_dir, 'my_cells.csv'))
    cus = pd.read_csv(os.path.join(my_dir, 'my_checkups.csv'))
    summ = json.load(open(os.path.join(my_dir, 'my_summary.json')))
    lev_hi = pd.read_csv(os.path.join(my_dir, 'my_lever_25_40.csv')).set_index('protocol')
    lev_lo = pd.read_csv(os.path.join(my_dir, 'my_lever_10_25.csv')).set_index('protocol')
    my_pre = pd.read_csv(os.path.join(my_dir, 'my_pairs_prespecified.csv'))
    my_wp = pd.read_csv(os.path.join(my_dir, 'my_warmer_pairs.csv'))

    pc = pd.read_csv(os.path.join(res_dir, 'stage2_cells.csv'))
    pl = json.load(open(os.path.join(res_dir, 'stage2_lever.json')))
    ps = json.load(open(os.path.join(res_dir, 'stage2_scorecard.json')))
    pp = pd.read_csv(os.path.join(res_dir, 'stage2_pairs.csv'))
    pw = pd.read_csv(os.path.join(res_dir, 'stage2_warmer_pairs.csv'))

    out = {}
    # ------------------------------------------------------------------ per cell
    m = mine.merge(pc, on='cell', suffixes=('_mine', '_pipe'), how='outer', indicator=True)
    assert (m['_merge'] == 'both').all() and len(m) == 144
    cell_cmp = {}
    for mc, pcn in [('L90_efc', 'L90_efc'), ('L80_efc', 'L80_efc'), ('T_C', 'T_C'),
                    ('T_pre_peltier_C', 'T_pre_peltier_C'), ('L90_est_efc', 'L90_est_efc'),
                    ('n_checkups', 'n_checkups'), ('cap0_Ah', 'q_first_Ah'), ('efc_last', 'efc_last_checkup'),
                    ('n_rows_T', 'runs_in_T')]:
        a = m[mc + '_mine'] if mc + '_mine' in m else m[mc]
        b = m[pcn + '_pipe'] if pcn + '_pipe' in m else m[pcn]
        st = diff_stats(a, b)
        for key in ('argmax_abs', 'argmax_rel'):
            if st[key] is not None:
                okmask = ~(np.isnan(np.asarray(a, float)) | np.isnan(np.asarray(b, float)))
                st[key] = m.loc[okmask, 'cell'].iloc[st[key]]
        cell_cmp[f'{mc} vs {pcn}'] = st
    for col in ['censored90', 'censored80']:
        a = m[col + '_mine'].astype(bool)
        b = m[col + '_pipe'].astype(str).str.lower().map({'true': True, 'false': False})
        cell_cmp[col] = dict(n_mismatch=int((a != b).sum()), n_censored_mine=int(a.sum()), n_censored_pipe=int(b.sum()),
                             censored_cells_mine=m.loc[a, 'cell'].tolist(), censored_cells_pipe=m.loc[b, 'cell'].tolist())
    # L90_est variants
    est_var = {}
    for col in ['L90_est_efc_mine'] + [c for c in mine.columns if c.startswith('L90_est_') and c != 'L90_est_efc']:
        st = diff_stats(m[col], m['L90_est_efc_pipe'])
        est_var[col] = dict(max_abs=st['max_abs'], max_rel=st['max_rel'], n_exact=st['n_exact'])
    cell_cmp['L90_est_efc_variants_vs_pipeline'] = est_var
    out['per_cell'] = cell_cmp

    # ------------------------------------------------------------------ levers
    def lever_cmp(mydf, plist, has_pre):
        rows = {}
        worst = {}
        for e in plist:
            p = e['protocol']
            r = mydf.loc[p]
            d = {}
            if bool(r['determinable']) != bool(e['determinable']):
                d['determinable_mismatch'] = True
            for mk, pk in [('T_lo', 'T_lo_C'), ('T_hi', 'T_hi_C'), ('dT', 'dT_K'), ('s', 's_per_K'), ('se', 'se'),
                           ('ci_lo', ('ci95', 0)), ('ci_hi', ('ci95', 1)), ('s_setpoint', 's_setpoint_per_K'),
                           ('Ea_kJ_per_mol', 'Ea_kJ_per_mol'), ('n_lo', 'n_lo'), ('n_hi', 'n_hi')] + \
                    ([('s_pre_peltier', 's_25_without_after_peltier')] if has_pre else []):
                if isinstance(pk, tuple):
                    pv = e[pk[0]][pk[1]] if e.get(pk[0]) is not None else None
                else:
                    pv = e.get(pk)
                if pv is None or (not bool(r['determinable']) and mk not in ('n_lo', 'n_hi')):
                    continue
                d[pk if isinstance(pk, str) else f'{pk[0]}[{pk[1]}]'] = num(r[mk], pv)
            # module translation: D33 says 0.391 K x s; pipeline stores offers in module_offers_K
            off_reg = pl['module_offers_K']['registered']
            off_dyn = pl['module_offers_K']['dynamic_post_hoc']
            rows[p] = d
            if not (bool(r['determinable']) and e.get('module_reg') is not None):
                d['present_in_pipeline'] = sorted(e.keys())
                continue
            d['module_reg'] = num(abs(r['s']) * off_reg, e['module_reg'])
            d['module_reg_signed_check'] = dict(pipeline_equals_abs_s_times_offer=bool(np.isclose(e['module_reg'], abs(e['s_per_K']) * off_reg, rtol=1e-12)),
                                                mine_with_0p391=float(0.391 * r['s']))
            d['module_dyn'] = num(abs(r['s']) * off_dyn, e['module_dyn'])
            rows[p] = d
            for k, v in d.items():
                if isinstance(v, dict) and 'abs_diff' in v:
                    w = worst.setdefault(k, dict(max_abs=0.0, max_rel=0.0))
                    w['max_abs'] = max(w['max_abs'], v['abs_diff'])
                    w['max_rel'] = max(w['max_rel'], v['rel_diff'])
        return dict(per_protocol=rows, worst_by_field=worst)

    out['lever_25_40'] = lever_cmp(lev_hi, pl['lever_25_40'], True)
    out['lever_10_25'] = lever_cmp(lev_lo, pl['lever_10_25'], False)
    l80 = pd.DataFrame(summ['lever_25_40_L80']['per_protocol']).T
    l80.index.name = 'protocol'
    for c in l80.columns:
        try:
            l80[c] = pd.to_numeric(l80[c])
        except (ValueError, TypeError):
            pass
    l80['determinable'] = l80['determinable'].astype(bool)
    out['lever_25_40_at_80pct'] = lever_cmp(l80, pl['lever_25_40_at_80pct'], True)
    out['pooled'] = dict(sigma2_lnL=num(summ['pooled_sigma_lnL'] ** 2, pl['sigma2_lnL']),
                         sigma_lnL=num(summ['pooled_sigma_lnL'], pl['floor']['sigma_lnL']),
                         sets_entering=dict(mine=summ['n_sets_entering'], pipeline=pl['sets_entering']),
                         censored90=dict(mine=summ['n_censored90'], pipeline=pl['censored90']),
                         censored80=dict(mine=summ['n_censored80'], pipeline=pl['censored80']),
                         sigma_lnL_80=num(summ['lever_25_40_L80']['sigma'], pl['floor_at_80pct']['sigma_lnL']))

    # ------------------------------------------------------------------ floor
    fl = summ['floor']
    out['floor'] = dict(pairs=dict(mine=fl['n_pairs'], pipeline=pl['floor']['pairs']),
                        mean=num(fl['mean'], pl['floor']['mean']), median=num(fl['median'], pl['floor']['median']),
                        p95=num(fl['p95'], pl['floor']['p95']),
                        worst_of_8=num(fl['worst_of_8_median_shortfall'], pl['floor']['worst_of_8']),
                        worst_of_8_alt_layout_8xN=num(fl['worst_of_8_alt_layout_8xN'], pl['floor']['worst_of_8']),
                        all_pairs_n=dict(mine=fl['all_pairs_n'], pipeline=pl['floor']['all_pairs']['pairs']),
                        all_pairs_mean=num(fl['all_pairs_mean'], pl['floor']['all_pairs']['mean']),
                        all_pairs_p95=num(fl['all_pairs_p95'], pl['floor']['all_pairs']['p95']))
    # per pair (pre-specified)
    mp = my_pre.merge(pp, on='parameter_id', how='outer', indicator=True)
    out['floor']['pair_rows'] = dict(n_mine=int(len(my_pre)), n_pipeline=int(len(pp)),
                                     same_sets=bool((mp['_merge'] == 'both').all()),
                                     L1=diff_stats(mp['L_a'], mp['L1']), L2=diff_stats(mp['L_b'], mp['L2']),
                                     rel_diff=diff_stats(mp['diff'], mp['rel_diff']))

    # ------------------------------------------------------------------ warmer
    wr = summ['warmer']
    w = dict(all=dict(pairs=dict(mine=wr['n_pairs'], pipeline=pl['warmer']['all']['pairs']),
                      warmer_shorter=dict(mine=wr['k_warmer_shorter'], pipeline=pl['warmer']['all']['warmer_shorter']),
                      share=num(wr['share'], pl['warmer']['all']['share']),
                      p_one_sided=num(wr['p_one_sided'], pl['warmer']['all']['p_one_sided'])))
    for t in ['0', '10', '25', '40']:
        a = wr['per_temp'][t]; b = pl['warmer'][f'{t} C']
        w[f'{t} C'] = dict(pairs=dict(mine=a['n_pairs'], pipeline=b['pairs']),
                           warmer_shorter=dict(mine=a['k_warmer_shorter'], pipeline=b['warmer_shorter']),
                           share=num(a['share'], b['share']), p_one_sided=num(a['p_one_sided'], b['p_one_sided']))
    # per pair rows: match on unordered cell pair
    my_wp['key'] = my_wp.apply(lambda r: tuple(sorted([r['cell_a'], r['cell_b']])), axis=1)
    pw['key'] = pw.apply(lambda r: tuple(sorted([r['warm'], r['cool']])), axis=1)
    mw = my_wp.merge(pw, on='key', how='outer', indicator=True, suffixes=('_mine', '_pipe'))
    both = mw[mw['_merge'] == 'both']
    pipe_ws = both['warmer_shorter_pipe'].astype(str).str.lower().map({'true': True, 'false': False})
    mine_ws = both['warmer_shorter_mine'].astype(str).str.lower().map({'true': True, 'false': False})
    my_Lw = np.where(both['warmer'] == both['cell_a'], both['L_a'], both['L_b'])
    my_Lc = np.where(both['warmer'] == both['cell_a'], both['L_b'], both['L_a'])
    w['pair_rows'] = dict(n_mine=int(len(my_wp)), n_pipeline=int(len(pw)), n_matched=int(len(both)),
                          warmer_identity_mismatch=int((both['warmer'] != both['warm']).sum()),
                          warmer_shorter_mismatch=int((mine_ws != pipe_ws).sum()),
                          dT=diff_stats(both['dT'], both['dT_K']),
                          L_warm=diff_stats(my_Lw, both['L_warm']), L_cool=diff_stats(my_Lc, both['L_cool']))
    out['warmer'] = w

    # ------------------------------------------------------------------ extras (after seeing structure)
    ex = {}
    # 80 % floor with the same machinery
    c80 = mine.copy()
    c80['L90_efc'] = c80['L80_efc']; c80['censored90'] = c80['censored80']
    s80 = summ['lever_25_40_L80']['sigma']
    fl80, _, _ = v2.floor(c80, s80)
    ex['floor_at_80pct'] = dict(pairs=dict(mine=fl80['n_pairs'], pipeline=pl['floor_at_80pct']['pairs']),
                                mean=num(fl80['mean'], pl['floor_at_80pct']['mean']),
                                median=num(fl80['median'], pl['floor_at_80pct']['median']),
                                p95=num(fl80['p95'], pl['floor_at_80pct']['p95']),
                                worst_of_8=num(fl80['worst_of_8_median_shortfall'], pl['floor_at_80pct']['worst_of_8']),
                                all_pairs_n=dict(mine=fl80['all_pairs_n'], pipeline=pl['floor_at_80pct']['all_pairs']['pairs']),
                                all_pairs_mean=num(fl80['all_pairs_mean'], pl['floor_at_80pct']['all_pairs']['mean']),
                                all_pairs_p95=num(fl80['all_pairs_p95'], pl['floor_at_80pct']['all_pairs']['p95']))
    # median within-pair dT
    dd = my_wp[my_wp['T_differ']]
    ex['warmer_dT_median_K'] = {'all': num(np.median(dd['dT']), pl['warmer']['all']['dT_median_K'])}
    for t in [0, 10, 25, 40]:
        ex['warmer_dT_median_K'][f'{t} C'] = num(np.median(dd.loc[dd['age_temp'] == t, 'dT']), pl['warmer'][f'{t} C']['dT_median_K'])
    # check-up spacing at the crossing (res = EFC gap of the bracketing check-ups)
    res = []
    for cell, g in cus.groupby('cell'):
        g = g.sort_values('cu')
        row = dict(cell=cell)
        for frac, tag, ccol in [(0.9, '90', 'censored90'), (0.8, '80', 'censored80')]:
            thr = frac * g['cap_Ah'].iloc[0]
            below = np.nonzero(g['cap_Ah'].to_numpy() < thr)[0]
            row['res' + tag] = (float(g['efc'].iloc[below[0]] - g['efc'].iloc[below[0] - 1]) if len(below) else np.nan)
        res.append(row)
    res = pd.DataFrame(res).merge(pc[['cell', 'res90_efc', 'res80_efc', 'L90_efc', 'L80_efc', 'censored90', 'censored80']], on='cell')
    ex['res90_efc'] = {k: v for k, v in diff_stats(res['res90'], res['res90_efc']).items() if not k.startswith('argmax')}
    ex['res80_efc'] = {k: v for k, v in diff_stats(res['res80'], res['res80_efc']).items() if not k.startswith('argmax')}
    r90 = (res['res90'] / res['L90_efc']).dropna()
    r80 = (res['res80'] / res['L80_efc']).dropna()
    ex['resolution_over_L90'] = dict(median=num(np.median(r90), pl['resolution_over_L90']['median']),
                                     p95=num(np.percentile(r90, 95), pl['resolution_over_L90']['p95']))
    ex['resolution_over_L80'] = dict(median=num(np.median(r80), pl['resolution_over_L80']['median']),
                                     p95=num(np.percentile(r80, 95), pl['resolution_over_L80']['p95']))
    out['extras_after_structure_seen'] = ex

    # ------------------------------------------------------------------ scorecard
    sc = summ['scorecard']
    key_map = {'holds': True, 'fails': False}
    s_cmp = {}
    for p in ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8']:
        mv = sc[p]['verdict']
        pv = ps[p]['held']
        s_cmp[p] = dict(mine=mv, pipeline=('holds' if pv is True else 'fails' if pv is False else str(pv)),
                        agree=(key_map.get(mv) == pv))
    # numbers inside the scorecard
    nums = {}
    nums['P1.mean'] = num(sc['P1']['mean_s'], ps['P1']['mean'])
    for k, v in ps['P1']['values'].items():
        nums[f'P1.values[{k}]'] = num(sc['P1']['s'][k], v)
    nums['P7.mean'] = num(sc['P7']['mean_s'], ps['P7']['mean'])
    for wname, wv in ps['P2']['windows'].items():
        mw_ = sc['P2']['windows'][wname]
        nums[f'P2.{wname}.s_1.67/1.00'] = num(mw_['s_167_100'], wv['s_1.67/1.00'])
        nums[f'P2.{wname}.s_0.33/1.00'] = num(mw_['s_033_100'], wv['s_0.33/1.00'])
        nums[f'P2.{wname}.s_0.33/0.33'] = num(mw_['s_033_033'], wv['s_0.33/0.33 (unscored)'])
        s_cmp['P2'].setdefault('windows', {})[wname] = dict(mine=mw_['verdict'], pipeline=('holds' if wv['held'] else 'fails'))
    nums['P3.n_determinable'] = num(sc['P3']['n_determinable'], ps['P3']['n_determinable'])
    nums['P3.min'] = num(lev_hi.loc[lev_hi['determinable'], 's'].min(), ps['P3']['min'])
    nums['P3.max'] = num(lev_hi.loc[lev_hi['determinable'], 's'].max(), ps['P3']['max'])
    nums['P4.mean'] = num(sc['P4']['mean'], ps['P4']['mean'])
    nums['P4.p95'] = num(sc['P4']['p95'], ps['P4']['p95'])
    nums['P4.pairs'] = num(summ['floor']['n_pairs'], ps['P4']['pairs'])
    nums['P5.pairs'] = num(summ['warmer']['n_pairs'], ps['P5']['pairs'])
    nums['P5.warmer_shorter'] = num(summ['warmer']['k_warmer_shorter'], ps['P5']['warmer_shorter'])
    nums['P5.share'] = num(sc['P5']['share'], ps['P5']['share'])
    nums['P5.p_one_sided'] = num(summ['warmer']['p_one_sided'], ps['P5']['p_one_sided'])
    nums['P5.dT_median_K'] = num(np.median(dd['dT']), ps['P5']['dT_median_K'])
    nums['P6.max'] = num(sc['P6']['max_s'], ps['P6']['max'])
    for k, v in ps['P8']['values'].items():
        nums[f'P8.values[{k}]'] = num(sc['P8']['per_protocol'][k]['s'], v)
    mine_neg = sorted(k for k, v in sc['P8']['per_protocol'].items() if v['negative'])
    s_cmp['P8']['negative_lists_agree'] = (mine_neg == sorted(ps['P8']['negative']))
    s_cmp['P8']['of_lists_agree'] = (sorted(sc['P8']['per_protocol']) == sorted(ps['P8']['of']))
    worst = max(nums.items(), key=lambda kv: kv[1]['rel_diff'])
    out['scorecard'] = dict(verdicts=s_cmp, numbers=nums,
                            all_verdicts_agree=all(v['agree'] for v in s_cmp.values()),
                            largest_rel_diff_in_scorecard=dict(field=worst[0], **worst[1]))
    # scorecard inside stage2_lever.json identical to stage2_scorecard.json?
    out['scorecard']['lever_json_scorecard_identical_to_scorecard_json'] = (json.dumps(pl['scorecard'], sort_keys=True) == json.dumps(ps, sort_keys=True))
    out['simple_counts'] = dict(cells=dict(mine=int(len(mine)), pipeline=pl['cells']),
                                module_offers_K=dict(pipeline=pl['module_offers_K'], D33_wording=dict(registered=0.391, dynamic_post_hoc=1.17)))

    # ------------------------------------------------------------------ largest differences over everything
    scan = []

    def walk(o, path):
        if isinstance(o, dict):
            if 'abs_diff' in o and 'rel_diff' in o:
                scan.append((path, o['abs_diff'], o['rel_diff']))
                return
            if 'max_abs' in o and 'max_rel' in o and o['max_abs'] is not None:
                scan.append((path, o['max_abs'], o['max_rel']))
                return
            for k, v in o.items():
                walk(v, f'{path}/{k}')
    walk(out, '')
    sens = [s for s in scan if ('L90_est' in s[0] or 'alt_layout' in s[0])]
    reg = [s for s in scan if s not in sens]
    a = max(reg, key=lambda s: s[1]); r = max(reg, key=lambda s: s[2])
    out['largest_differences'] = dict(
        registered_and_reported_numbers=dict(n_fields_compared=len(reg), max_abs=dict(field=a[0], value=a[1]),
                                             max_rel=dict(field=r[0], value=r[2])),
        sensitivity_lifetime_L90_est_efc=dict(primary_reading=out['per_cell']['L90_est_efc vs L90_est_efc'],
                                              matching_reading='L90_est_keep_skip_strict (NaN estimate rows kept as rows, skipped inside the 11-row window): 144 of 144 exact'),
        worst_of_8_alternative_layout=out['floor']['worst_of_8_alt_layout_8xN'])

    # ------------------------------------------------------------------ my per-cell values
    keep = ['cell', 'parameter_id', 'parameter_nr', 'age_temp', 'protocol', 'n_checkups', 'cap0_Ah', 'efc_last',
            'L90_efc', 'censored90', 'L80_efc', 'censored80', 'T_C', 'T_pre_peltier_C', 'n_rows_T',
            'L90_est_efc', 'L90_est_keep_skip_strict', 'window_end_utc']
    out['my_cells'] = mine[keep].rename(columns={'L90_est_efc': 'L90_est_efc_drop_nan_first',
                                                 'L90_est_keep_skip_strict': 'L90_est_efc_nan_skipped_in_window'}).to_dict(orient='records')
    out['my_levers'] = dict(
        lever_25_40={p: dict(s=r['s'], se=r['se'], ci=[r['ci_lo'], r['ci_hi']], T_25=r['T_lo'], T_40=r['T_hi'], n_25=int(r['n_lo']), n_40=int(r['n_hi']),
                             s_setpoint=r['s_setpoint'], s_pre_peltier=r['s_pre_peltier'], Ea_kJ_per_mol=r['Ea_kJ_per_mol'])
                     for p, r in lev_hi.iterrows()},
        lever_10_25={p: dict(s=r['s'], se=r['se'], ci=[r['ci_lo'], r['ci_hi']], T_10=r['T_lo'], T_25=r['T_hi'], n_10=int(r['n_lo']), n_25=int(r['n_hi']),
                             s_setpoint=r['s_setpoint'], s_pre_peltier=r['s_pre_peltier'], Ea_kJ_per_mol=r['Ea_kJ_per_mol'])
                     for p, r in lev_lo.iterrows()},
        pooled_sigma_lnL=summ['pooled_sigma_lnL'], pooled_dof=summ['pooled_dof'], pooled_n_sets=summ['pooled_n_sets'])
    out['my_floor'] = summ['floor']
    out['my_warmer'] = summ['warmer']
    out['my_scorecard'] = summ['scorecard']

    with open(out_json, 'w') as fh:
        json.dump(out, fh, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    return out


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
