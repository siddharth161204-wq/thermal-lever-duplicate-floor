# Stage 2, V2: independent recomputation and independent read

Date: 21 October 2026. Machine-readable record: `verification/stage2/v2_independent.json`. Scripts: `verification/stage2/independent/` (`v2_independent.py` is the recomputation, `compare.py` the comparison, `assemble.py` the record).

## Verdict

The independent code reproduces every per-cell lifetime, censoring flag and temperature, both levers with their standard errors, intervals and sensitivities, the pooled variance, the floor, the worst-of-eight simulation, the all-pairs sensitivity and the warmer-replicate share and p value, all to floating-point rounding. There are 583 comparison entries, where a per-cell or per-pair column counts as one entry. Across them, the largest absolute difference is 4.5e-13 (in EFC) and the largest relative difference is 2.3e-13. All eight verdicts agree: P1 fails, P2 holds, P3 holds, P4 fails, P5 fails, P6 holds, P7 holds, P8 holds. Every number in `stage2_scorecard.json` matches my recomputation to 1.6e-15 relative or better.

The only value that differs beyond rounding is the unscored sensitivity lifetime `L90_est_efc`. D28 does not say how rows with a missing running estimate enter the median, and the two readings differ in 10 of 144 cells, by at most 1.58 EFC (3.6 %). Four further points are matters of wording or representation, and none changes a number that is scored (flags F1 to F5 below).

## How the check was done

I wrote the code from DECISIONS.md D25 to D35, `predictions_stage2.md` and `inputs/kit/cyclic_cells.csv` only. I opened nothing under `src/stage2` or `tests`. I read the 144 raw EOC files in `/tmp/kit/eoc/` (parameter sets 17 to 64). None of these files has a row with timestamp_origin 5, so that rule removes nothing. My outputs were frozen and hashed at 06:47:59 UTC, and only then did I open the five permitted pipeline outputs. Those files carry a modification time of 06:43:17 UTC, and their SHA-256 hashes are in the JSON. Software: Python 3.13.16, numpy 2.5.3, pandas 3.0.5, scipy 1.18.1. A fresh run of the copied script reproduces the frozen outputs byte for byte.

Where the wording left a choice, I read it as follows.

- D27. Within each check-up (cyc_condition 2, grouped by num_cycles_checkup), I took the discharge rows with v_min_target_V 2.5 and i_dischg_min_A -1.0, and of those the row with the largest |delta_q_dischg_Ah|. Every check-up has at most one such row. One check-up (P041-3, CU 20) holds only a charge row, so it has no capacity and is not counted. EFC = (total_q_dischg_sum_Ah - |delta_q_dischg_Ah|) / 3.0.
- D28. The reference is the first check-up's capacity. "Below" means strictly below. L is interpolated linearly between the last check-up at or above the threshold and the first one below it. A censored cell is censored at the EFC of its last check-up that has a capacity. The sensitivity lifetime uses the cyc_condition 1 rows, a centred eleven-row median of cap_aged_est_Ah and row EFC = total_q_dischg_sum_Ah / 3.0. My primary reading drops rows with a missing estimate before taking the median. I computed three alternatives beside it (F1).
- D29. T is the cyc_duration_s-weighted mean of (t_start_degC + t_end_degC) / 2 over cyc_condition 1 rows with a duration above zero, at or after the first check-up, with row EFC at or below L. A censored cell uses all its rows. No regular-operation row precedes the first check-up in any cell.
- D30. The cutoff is 1682208000 s (23 April 2023, 00:00 UTC), and rows with a later timestamp_s are dropped. The pooled variance uses ddof 1 over the 48 sets, which gives 95 degrees of freedom. The activation energy uses the mean of the two group temperatures plus 273.15 K.
- D31. Cells 1 and 2 are parameter_nr 1 and 2. I used np.percentile with its default linear rule. The simulation draws default_rng(20261021).normal(0, sigma, size=(200000, 8)), exponentiates, takes (mean - min) / mean per module and then the median.
- D32. A pair counts when the two T values differ at all. The p value is scipy binomtest with alternative "greater".
- D34. The bands of P1, P4 and P5 are closed, the P7 interval is open, and "exceeds" in P6 means strictly greater.

## Per-cell comparison (144 cells)

| Quantity | Compared | Largest abs. difference | Largest rel. difference |
|---|---|---|---|
| L90_efc | 143 uncensored | 1.1e-13 EFC | 3.8e-16 |
| L80_efc | 130 uncensored | 1.1e-13 EFC | 4.4e-16 |
| censored90, censored80 | 144 | 0 mismatches (1 and 14 censored) | |
| T_C | 144 | 0 (identical) | 0 |
| T_pre_peltier_C | 144 | 0 (identical) | 0 |
| L90_est_efc, my primary reading | 144 | 1.58 EFC (P031-1) | 3.6 % (P031-1) |
| L90_est_efc, NaN rows kept and skipped | 144 | 0 (identical) | 0 |

For a censored cell, the pipeline leaves L90_efc or L80_efc empty and puts the censoring EFC in efc_last_checkup. That value equals mine exactly for the one cell censored at 90 % (P044-2, 53.0 EFC) and the fourteen censored at 80 %. The check-up count, the first check-up capacity, the last check-up EFC and the number of runs in T also agree exactly for every cell. T_pre_peltier_C equals T_C in all 144 cells because no cell's lifetime window reaches 23 April 2023. The latest window to close is P049-2's, on 27 March 2023.

## Summary numbers

| Number | Mine | Pipeline | Largest rel. difference |
|---|---|---|---|
| Pooled sigma^2 of ln L (48 sets, 95 dof) | 0.0068018081 | 0.0068018081 | 0 |
| Lever 25 to 40 C: s, SE, interval, 15 K s, post-April s, Ea, module values (12 protocols) | see below | same | 4.0e-14 |
| Lever 10 to 25 C: the same fields (12 protocols) | see below | same | 2.0e-14 |
| Lever 25 to 40 C at 80 % (9 determinable, the same 3 not) | | same | 3.5e-16 |
| Floor, 47 pre-specified pairs: mean / median / 95th percentile | 0.073323 / 0.040902 / 0.204751 | same | 2.7e-16 |
| Worst of eight, median shortfall | 0.10890735 | 0.10890735 | 0 |
| All-pairs sensitivity, 142 pairs: mean / 95th percentile | 0.078459 / 0.228151 | same | 4.9e-16 |
| Warmer share, all temperatures | 57 / 142 = 0.40141 | same | 0 |
| One-sided sign test p | 0.99267 | same | 0 |
| Shares at 0 / 10 / 25 / 40 C | 12/36, 15/36, 14/34, 16/36 | same, and the p values | 0 |
| Per-pair rows (pairs and warmer pairs files) | 47 and 142 | same pairs, same warmer cell, same outcome | 2.3e-13 |

I added three items after seeing which secondary numbers the pipeline reports: the 80 % floor (40 pairs, mean 0.081801, 95th percentile 0.220968, worst of eight 0.142736), the check-up spacing at the crossing relative to L (median 0.6116, 95th percentile 1.9373 at 90 %) and the median within-pair temperature difference (0.3278 K overall). They use the same definitions. The summary values agree to 7.5e-16 relative or better, and the per-cell spacings behind them to 4.1e-15.

The 25 to 40 C lever, my values (identical to the pipeline's to rounding):

| Protocol | n 25/40 | s per K | 95 % interval |
|---|---|---|---|
| 0-100 0.33/0.33C | 3/3 | +0.01303 | +0.00412 to +0.02194 |
| 0-100 0.33/1.00C | 3/3 | -0.00624 | -0.01514 to +0.00267 |
| 0-100 1.00/1.00C | 3/3 | -0.02685 | -0.03588 to -0.01782 |
| 0-100 1.67/1.00C | 2/3 | -0.03455 | -0.04454 to -0.02455 |
| 10-100 0.33/0.33C | 3/3 | +0.01360 | +0.00459 to +0.02261 |
| 10-100 0.33/1.00C | 3/3 | +0.00672 | -0.00221 to +0.01566 |
| 10-100 1.00/1.00C | 3/3 | -0.00422 | -0.01310 to +0.00465 |
| 10-100 1.67/1.00C | 3/3 | -0.01915 | -0.02803 to -0.01026 |
| 10-90 0.33/0.33C | 3/3 | +0.01910 | +0.01004 to +0.02815 |
| 10-90 0.33/1.00C | 3/3 | +0.01090 | +0.00187 to +0.01992 |
| 10-90 1.00/1.00C | 3/3 | -0.00194 | -0.01070 to +0.00683 |
| 10-90 1.67/1.00C | 3/3 | -0.00995 | -0.01880 to -0.00111 |

The 10 to 25 C lever is negative at all twelve protocols. For the six P8 protocols, s is -0.0334 (0-100 1.00/1.00C), -0.0358 (0-100 1.67/1.00C), -0.0417 (10-100 1.00/1.00C), -0.1159 (10-100 1.67/1.00C), -0.0388 (10-90 1.00/1.00C) and -0.0524 (10-90 1.67/1.00C) per K.

## Independent scoring by D34

| | Rule | My number | Mine | Pipeline |
|---|---|---|---|---|
| P1 | the six 0.33 C charge protocols all have s > 0, and their mean is in [0.01, 0.06] | 0-100 0.33/1.00C is -0.0062; mean 0.00952 | fails | fails |
| P2 | s at 1.67/1.00 is below s at 0.33/1.00 in each window | -0.0345 < -0.0062; -0.0191 < +0.0067; -0.0100 < +0.0109 | holds | holds |
| P3 | the determinable estimates take both signs | 12 determinable: 5 positive, 7 negative | holds | holds |
| P4 | mean in [1.5, 6.0] %, 95th percentile in [4.4, 17.4] % | 7.33 % and 20.48 % over 47 pairs | fails | fails |
| P5 | share in [45, 65] % | 40.1 % (57 of 142) | fails | fails |
| P6 | no determinable s above 0.077 per K | largest 0.0191 (10-90 0.33/0.33C) | holds | holds |
| P7 | mean of P1's six strictly between -0.0099 and 0.052 | 0.00952 | holds | holds |
| P8 | at least three of the six 1C-and-above protocols negative from 10 to 25 C | 6 of 6 | holds | holds |

Every protocol is determinable at both levers (all 48 sets enter), so no part of P1, P2 or P8 depends on the rules for protocols that are not determinable. The P2 windows agree one by one, and so do the P8 list of negative protocols and its denominator. The scorecard block inside `stage2_lever.json` is identical to `stage2_scorecard.json`. I found no rule that the pipeline applies differently from D34. No verdict sits on a band edge: P1 fails on the sign of 0-100 0.33/1.00C regardless of its mean, which misses the 0.01 edge by 0.0005.

## Flags

F1 (unscored sensitivity). D28 does not fix how a row with a missing cap_aged_est_Ah (NaN; a few per cell, about one per check-up interval) enters the eleven-row centred median. Dropping those rows first, my primary reading, gives 10 cells that differ from the pipeline: P031-1 by +1.58 EFC, P043-2 by +0.50 and P041-2 by -0.38, with the rest at 0.33 EFC or less. Keeping those rows as rows and skipping their value inside each window reproduces the pipeline exactly in all 144 cells, whether or not partial windows are allowed at the ends. A third reading, where any NaN in the window voids the median, differs by up to 0.95 EFC. The wording supports both the primary and the pipeline's reading. I suggest one clause in D28 or the note that states the pipeline's rule. The effect is negligible next to the gap between the check-up and running-estimate crossings themselves (for example 510.6 against 41.1 EFC in P049-1), which is V3's subject.

F2 (unscored). D33 writes the offer as 0.391 K x s. The pipeline multiplies |s| by 0.391279 K, the unrounded value in Stage 1's `results/module_spread.csv`, which D33 names as the source. The offer differs by 7.1e-4 relative. The absolute value matters only for negative s, where a signed product would compare trivially below the floor. With the pipeline's offer and |s|, my values agree to 3e-14 relative. I suggest the note state 0.391279 K and |s|, or say that 0.391 K is a rounding.

F3 (no numerical effect). The pipeline reports the post-23-April sensitivity (s_25_without_after_peltier) for the 25 to 40 C lever only. D30's "the low-temperature lever is the same" can be read to include that sensitivity. It is vacuous for both levers: no lifetime window at any temperature reaches 23 April 2023, so the sensitivity equals the main estimate exactly and T_post_peltier_C is empty for every cell. The note should say that this sensitivity removes nothing at 90 %.

F4 (representation). Censored cells have an empty L in `stage2_cells.csv`, with the censoring EFC in efc_last_checkup. This is correct, but it is easy to misread without a sentence that says so.

F5 (none). D31 does not fix the array layout of the simulation. The (200000, 8) layout reproduces the pipeline exactly. Drawing the same stream as (8, 200000) gives 0.10886 instead of 0.10891, a Monte Carlo difference of 4.8e-5. Recording the layout next to the seed would make the number reproducible from the text alone.

## Observations for the note

P044-2 is the only cell censored at 90 %. It has two check-ups (0.26 and 53.0 EFC, at 90.9 % at the second) and no later one, although it kept cycling to 199 EFC. Its set enters with two cells, so the lower-bound rule of D28 is never used at 90 %.

At 90 %, the check-up spacing at the crossing has a median of 0.61 L and a 95th percentile of 1.94 L. The lifetimes therefore rest on linear interpolation across wide gaps, which belongs in any reading of the P4 failure. By my computation, the pre-specified floor's mean is 5.4 % at 0 C, 9.0 % at 10 C, 10.9 % at 25 C and 4.3 % at 40 C.

On the 10 to 25 C lever, 0.391279 K x |s| for 10-100 1.67/1.00C is 4.5 %, above the 3.0 % mean floor. P6 and D33 concern the 25 to 40 C lever, where the largest value is 1.35 %.
