"""Assemble the machine-readable V2 record: comparison.json plus provenance and flags.

Usage: python3 -I assemble.py COMPARISON_JSON FROZEN_DIR RESULTS_DIR OUT_JSON
"""
import sys
import os
import json
import hashlib
import platform

import numpy as np
import pandas as pd
import scipy


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as fh:
        h.update(fh.read())
    return h.hexdigest()


def main(cmp_json, frozen, res_dir, out_json):
    c = json.load(open(cmp_json))
    frozen_at = open(os.path.join(frozen, 'FROZEN_AT_UTC.txt')).read().strip()
    meta = dict(
        check='V2 (DECISIONS.md D35): lifetimes, temperatures, levers and the floor recomputed from the raw EOC files by separate code; plus an independent scoring of P1 to P8 by D34',
        date='2026-10-21',
        independence=('Code written from DECISIONS.md D25 to D35, predictions_stage2.md and inputs/kit/cyclic_cells.csv only. '
                      'Nothing under src/stage2 or tests was opened. Pipeline outputs were opened only after the independent outputs '
                      'were frozen and hashed; of the results, only stage2_cells.csv, stage2_lever.json, stage2_scorecard.json, '
                      'stage2_pairs.csv and stage2_warmer_pairs.csv were read.'),
        independent_outputs_frozen_at_utc=frozen_at,
        frozen_output_sha256={f: sha(os.path.join(frozen, f)) for f in sorted(os.listdir(frozen))
                              if f.endswith(('.csv', '.json', '.py'))},
        pipeline_outputs_compared_sha256={f: sha(os.path.join(res_dir, f)) for f in
                                          ['stage2_cells.csv', 'stage2_lever.json', 'stage2_scorecard.json',
                                           'stage2_pairs.csv', 'stage2_warmer_pairs.csv']},
        data=dict(eoc_dir='/tmp/kit/eoc', files_used='cell_eocv2_P017 to P064, 144 files (parameter sets 17 to 64)',
                  rows_with_timestamp_origin_5_in_these_files=0),
        software=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__),
        scripts=['independent/v2_independent.py (the recomputation, frozen before any comparison)',
                 'independent/compare.py (comparison and extras)', 'independent/assemble.py (this record)'],
    )
    flags = [
        dict(id='F1', severity='none for any registered number or verdict (unscored sensitivity)',
             item='D28 sensitivity lifetime L90_est_efc: handling of NaN cap_aged_est_Ah rows',
             finding=('D28 does not say how rows with a missing running estimate enter the eleven-row centred median. Under my primary '
                      'reading (rows with NaN estimate removed before the median) 10 of 144 cells differ from the pipeline, by at most '
                      '1.58 EFC (3.6 %, P031-1). Keeping those rows as rows and skipping the NaN inside each window reproduces the pipeline '
                      'exactly for all 144 cells, with or without partial windows at the ends. Both readings fit the wording.'),
             suggestion='State the NaN rule in one clause of D28 or the note (the pipeline counts NaN rows in the window and ignores their value).'),
        dict(id='F2', severity='none (unscored numbers)',
             item='D33 module translation',
             finding=('D33 words the offer as 0.391 K x s; the pipeline multiplies |s| by 0.391279 K, the unrounded value in Stage 1 results/module_spread.csv '
                      'that D33 names as the source. The offer differs by a relative 7.1e-4, and the absolute value only changes the sign for negative s. '
                      'My recomputation with the pipeline offer and |s| agrees to 3e-14 relative.'),
             suggestion='Either write 0.391279 K and |s| in D33 or the note, or say that the stated 0.391 K is rounded.'),
        dict(id='F3', severity='none (no numerical effect)',
             item='D30 post-23-April-2023 sensitivity for the 10 to 25 C lever',
             finding=('The pipeline reports s_25_without_after_peltier for the 25 to 40 C lever only; D30 says the low-temperature lever is "the same", which can be read to '
                      'include its sensitivities. Every cell lifetime window (all temperatures) closes by 27 March 2023 (latest P049-2), so no row is removed and the sensitivity '
                      'equals the main estimate exactly for both levers; T_pre_peltier_C equals T_C for all 144 cells.'),
             suggestion='Say in the note that the sensitivity is vacuous at 90 % because no 25 C window reaches 23 April 2023.'),
        dict(id='F4', severity='none (representation)',
             item='Censored lifetimes in stage2_cells.csv',
             finding='For censored cells the pipeline leaves L90_efc / L80_efc empty; the censoring EFC is in efc_last_checkup and equals my value exactly (1 cell at 90 %, 14 at 80 %).',
             suggestion='None needed; optionally say so in the column description.'),
        dict(id='F5', severity='none',
             item='D31 worst-of-eight layout',
             finding=('With default_rng(20261021).normal(0, sigma, size=(200000, 8)) my median shortfall equals the pipeline bit for bit (0.10890735181666544). '
                      'Drawing the same stream as (8, 200000) gives 0.10885962, a Monte Carlo difference of 4.8e-5. D31 does not fix the layout; the pipeline uses the natural one.'),
             suggestion='Optionally record the array layout next to the seed.'),
    ]
    observations = [
        'P044-2 is the only cell censored at 90 %: two check-ups (0.26 and 53.0 EFC, 90.9 % at the second) and no check-up afterwards although it kept cycling to 199 EFC. The 0-100 1.67/1.00C 25 C set enters with two cells, so all 48 sets enter and the lower-bound rule of D28 is not used.',
        'P1 fails on two counts: 0-100 0.33/1.00C has s = -0.0062 per K (interval -0.0151 to 0.0027), and the mean of the six, 0.00952 per K, is below the 0.01 edge. The verdict is therefore not decided by the 0.0005 margin to the band edge.',
        'The check-up spacing at the 90 % crossing is large relative to L (median 0.61 L, 95th percentile 1.94 L, as the pipeline also reports), so L90 rests on linear interpolation across wide gaps; this is relevant when reading the P4 failure (floor mean 7.3 %, 95th percentile 20.5 %).',
        'Per temperature (my numbers, not in the pipeline outputs read): pre-specified pair floor mean 5.4 % at 0 C, 9.0 % at 10 C, 10.9 % at 25 C, 4.3 % at 40 C.',
        'On the 10 to 25 C lever, 0.391279 K x |s| for 10-100 1.67/1.00C is 4.5 %, above the 3.0 % mean floor; P6 and D33 concern the 25 to 40 C lever, where the largest |s| x 0.391279 K is 1.35 %.',
    ]
    rec = dict(meta=meta,
               summary=dict(
                   all_eight_verdicts_agree=c['scorecard']['all_verdicts_agree'],
                   verdicts={p: v['mine'] for p, v in c['scorecard']['verdicts'].items()},
                   largest_differences=c['largest_differences'],
                   per_cell_headline={k: {kk: vv for kk, vv in v.items() if kk in ('n_compared', 'max_abs', 'max_rel', 'n_exact', 'n_mismatch', 'n_nan_mismatch')}
                                      for k, v in c['per_cell'].items() if k != 'L90_est_efc_variants_vs_pipeline'}),
               flags=flags, observations=observations)
    rec.update({k: v for k, v in c.items()})
    with open(out_json, 'w') as fh:
        json.dump(rec, fh, indent=1)


if __name__ == '__main__':
    main(*sys.argv[1:5])
