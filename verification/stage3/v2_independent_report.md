# Stage 3, V2: independent recomputation and independent read

Date: 23 October 2026 (UTC). Machine-readable record: `verification/stage3/v2_independent.json`. Scripts, the reading and flag lists and my frozen outputs: `verification/stage3/independent/`. Results checked: commit 9c6bb84 on branch `stage3`.

## Verdict

Every registered number and all ten verdicts reproduce.

When I ran my separate code with the two readings the pipeline uses where D41 and D43 leave a choice, it reproduced the pipeline's per-cell table to floating-point rounding. Lifetimes agree to 1.1e-14 relative. Swing, rest temperature, unit counts and check lists are identical, and T_cyc agrees to 1.9e-16. All 1,306 group-level fields I compared also reproduce. They cover levers with their bootstrap intervals, the laws, the floor, all 86 pairs, both oriented-pair files, the thermal share with its slope intervals, the module numbers and the scorecard. The largest absolute difference is 1.3e-13 and the largest relative difference 7.0e-13.

My own primary readings were frozen at 20:58:11 UTC, before any pipeline output was opened. Under them the ten verdicts are the same:

- P1, P2, P6 and P9 hold.
- P3, P5, P7 and P10 do not hold.
- P4 and P8 are not scored.

Every number a verdict rests on agrees to 6.5e-5 relative or better (the NCA lever and P9). The reported intervals agree to 2.3e-4, and P6's reported median difference to 7.3e-3. Larger differences occur only in unscored 2C levers and D47 pair averages. Two reading choices explain them (F1, F2).

V2 as registered (D49) passes. In the three V2 cells (LFP 15 C 0-100 0.5/1C a, LFP 15 C 0-100 0.5/2C a and LFP 25 C 0-100 0.5/3C a), lifetime, measure (a), measure (b) and T_cyc agree with the pipeline to 7.1e-15. A second implementation in plain Python (csv module and loops) agrees with my numpy code to 5.6e-15. The four extra cells agree to the same precision. Every D50 count reproduces exactly, cell by cell.

On the three points I was asked to judge:

- (i) Reporting P4 as "not scored" is a defensible literal reading, but it is not forced. The reason given names the wrong rule, and the two natural repairs give opposite verdicts. Taking the mean T_cyc over the three cells that have one, P4 holds (ratio to Law 2 2.92, against a limit of 3). With the 10 K setpoint, P4 does not hold (ratio 3.09).
- (ii) P7 "did not hold" on one untied pair is what D48 says, but the verdict carries no information. Without the tie rule, which D46 does not allow, the share would be 24 of 60 = 0.40, inside the band.
- (iii) The 0.5C exclusion is applied consistently.

Flags that change or could change a scored number or verdict:

- F5 (high): the P4 verdict depends on a reading the rules did not foresee.
- F6 (medium): P7 rests on one pair, untied by 13.4 EFC.
- F3 (medium): spurious checks double the NMC floor mean that P5 scores (0.375 against 0.178 without them), but P5 fails either way.
- F4 (low to medium): one-record discharge steps split capacity checks. Under the other reading, NMC's registered lever would move by 6.4 % and NCA's by 1.8 %, with no verdict change.
- F2 and F1 (low): the two reading differences, worth at most 6.5e-5 relative on a scored number.
- F7 (low to medium): a skin-temperature channel that does not follow its cell; no verdict change.

F8 to F10 change nothing that is scored.

## How the check was done

I wrote the code from these sources only:

- DECISIONS.md (read in full, D1 to D50) and `predictions_stage3.md`;
- `inputs/sandia/inventory.csv`, `inputs/floor.json` and `results/module_spread.csv`;
- `results/stage2_lever.csv` and `.json`, read unchanged per D44 and D50 (9);
- `src/thermal_lever.py`, imported for Law 1 (`ln_L_kucinskis`) and Law 2 (`schimpe_cycle_coeffs`), as the brief allows.

Before the freeze I did not open:

- anything under `src/stage2`, `src/stage3`, `tests` or `/tmp/sandia/steps`;
- the pipeline's own check outputs in `verification/stage3`;
- any `results/stage3_*` file.

I read the 86 raw time series in `/tmp/sandia/raw_*` read-only, with `python3 -I`, from scripts kept outside the data folders. Software: Python 3.13.16, numpy 2.5.3, pandas 3.0.5, scipy 1.18.1.

The scripts are:

- `s0_cache.py` parses each CSV once and hashes it.
- `lib.py` holds the row, step, ampere-hour and temperature rules.
- `s1_cells.py` builds the per-cell table, with every reading I identified as an option.
- `s2_v2_purepython.py` is the second implementation, without numpy, for the V2 and extra cells.
- `s3_groups.py` computes the levers, laws, floor, thermal share, module numbers and verdicts.
- `s4_variants.py` runs every reading option end to end.
- `s5_freeze.py` writes `independent/frozen/`, records the UTC time and the SHA-256 of every output, script and input in `frozen/freeze_manifest.json`, and makes no other change.

Only after the freeze did I hash the pipeline files (at 20:58:59 UTC; they were last modified between 20:02 and 20:13 UTC and are unchanged since 9c6bb84) and open them. The comparison scripts `s6_compare_cells.py`, `s7_compare_groups.py`, `s8_assemble.py` and `s9_tables.py` were written after that.

The per-cell comparison showed the pipeline uses two readings that differ from my primary ones. Both were already options in the frozen code (one-sided check voltage; units whose previous non-rest step is a charge). Re-running the frozen code with those two options is the "pipeline readings" run below.

Two things in that run were written after seeing the pipeline, and are labelled as such:

- the thermal-share layout: which subsets carry a slope, and "one pair per group" taken as the warmest and coolest cells;
- the search for the unscored sensitivity-lifetime rule.

During the check, at 20:42:06 UTC, commit 767bb9c (not mine) added early versions of `lib.py`, `s0_cache.py`, `s1_cells.py` and `s3_groups.py` from this folder. The files were changed after that, and the frozen versions are the ones hashed in the manifest. The same commit added `results/stage3_posthoc.json` and `src/stage3/posthoc.py`, and it changed none of the inputs or results compared here. I never opened `src/stage3/posthoc.py`. I opened `results/stage3_posthoc.json` only after the freeze, with the coordinator's leave, to compare my point (i) numbers; nothing of it entered my recomputation. I made no commit.

### Readings where the wording left a choice

The full list, R1 to R21, is in `independent/readings.json`. These are the ones that matter:

- **Rows (D50 (1)).** A late record has a cycle index exactly one below the running maximum on the rows as logged. A row more than one below is another record's and is dropped. Then the rows get a stable sort by cycle index.
- **Steps (D40).** A step is a maximal run of one label (rest when |I| < 0.005C). There is no break at a cycle-index change. A step's duration is the time of its last sample minus that of its first.
- **Ampere-hours (D50 (2)).** Each step sums its own counter transitions, a decrease counting the new value in full. The rule is applied per step, not through a running sum over the file. My first version used a running sum, and it put a 1e-9 relative drift on lifetimes until I changed it.
- **Checks (D41).** A discharge qualifies on the arithmetic mean of its sample currents. My primary reading of "within 0.02 V of the cut-off" is two-sided; the pipeline is one-sided (F1). Consecutive qualifying discharges in the sequence of discharge steps form a check, so any discharge step between them splits it (F4). A discharge in a check later dropped by the half rule still counts as part of a check.
- **Lifetime (D42).** A check "falls below" the threshold strictly. The lifetime is interpolated linearly between the last kept check at or above the threshold and the first below, and the resolution is their spacing. No cell is censored.
- **Units (D43).** My primary reading pairs every non-check discharge with the last charge before it. The pipeline takes a non-check discharge only when the step before it, rests aside, is a charge (F2).
- **Step temperature (D50 (3)).** 1,086 steps have valid pairs whose durations sum to zero: two records at one instant, 1,074 of them rests. Such a step takes the mean of its valid samples, the limit of the trapezoidal mean. The pipeline does the same.
- **Measure (b).** For each unit, I average the step means of the rest steps directly before or after its charge or its discharge, each counted once. The cell's value is the mean over its units. The pipeline does the same.
- **Group temperature (D43).** My primary reading is strict: undefined if an uncensored cell has no T_cyc. The repairs are computed beside it (point (i)).
- **Bootstrap (D44).** Each iteration calls integers(0, n_lo, n_lo) and then integers(0, n_hi, n_hi). The order is the brief's. A secant with an undefined group temperature draws nothing, and intervals use np.percentile with its default linear rule.
- **Ties (D45, D46).** A tie is |L_a - L_b| < max(res_a, res_b) in EFC. The sign test uses scipy binomtest, "greater" and two-sided. Bands are closed, and "exceeds" and "less than" are strict.

## V2 as registered (D49)

The selection rule takes, in cell-ID order, the first cell with a 1C discharge, the first with 2C and the first with 3C whose lifetime is uncensored. No cell is censored, so these are the three cells expected.

| Cell | L (EFC) | resolution | bracketing checks (fraction of first, EFC) | measure (a), K | measure (b), C | T_cyc, C | units | largest abs. difference to pipeline |
|---|---|---|---|---|---|---|---|---|
| LFP 15 C 0-100 0.5/1C a | 2518.6538798701613 | 486.66364 | 6 (0.954443, 2279.667) and 7 (0.945396, 2766.331) | 0.5557414393299501 | 15.493074145975092 | 16.113602134693167 | 2718 | 0 |
| LFP 15 C 0-100 0.5/2C a | 2137.927689393939 | 222.01818 | 5 (0.956746, 2041.257) and 6 (0.941252, 2263.276) | 1.8565158314759256 | 14.954118485946838 | 15.859771392949181 | 2321 | 3.6e-15 (b) |
| LFP 25 C 0-100 0.5/3C a | 1257.2879679144346 | 699.68182 | 6 (0.955174, 1124.760) and 7 (0.927860, 1824.442) | 2.0421262029297766 | 24.676875759572667 | 25.323956411287174 | 1348 | 7.1e-15 (b) |

The first lifetime can be checked by hand. The threshold is 0.95 x 1.0316667 = 0.9800833 Ah. Interpolating between 0.9846667 Ah at 2279.6673 EFC and 0.9753333 Ah at 2766.3309 EFC gives 2279.6673 + 0.491071 x 486.6636 = 2518.654 EFC.

The table below shows one unit per cell, from the middle of each cell's units. The full sample lists are in the JSON under `v2_registered`. In each of these cells the only rest next to a unit is a single zero-current record between charge and discharge. Recomputing the trapezoidal means from the listed samples gives the same values to the last digit shown.

| Cell, unit | charge: samples, duration, mean | rest record (t s, T C, V) | discharge: samples, duration, mean, Ah | swing |
|---|---|---|---|---|
| LFP 15 C 1C a, unit 1359 (cycle 1375, 1275.26 EFC) | 59, 6867.706 s, 15.751278 C | 15261064.969, 15.217, 3.519 V | 28, 3176.619 s, 16.408970 C, 1.0070 Ah | 0.657691 K |
| LFP 15 C 2C a, unit 1160 (cycle 1176, 1086.11 EFC) | 58, 6774.031 s, 14.842766 C | 11033153.314, 14.088, 3.515 V | 14, 1519.069 s, 16.701180 C, 1.0010 Ah | 1.858414 K |
| LFP 25 C 3C a, unit 674 (cycle 693, 636.74 EFC) | 58, 6679.959 s, 24.998234 C | 6095260.269, 24.648, 3.529 V | 10, 969.771 s, 27.046326 C, 0.9990 Ah | 2.048092 K |

The 3C discharge's ten samples (t in s, T in C) are:

- (6095380.315, 25.182), (6095500.339, 25.915), (6095620.363, 26.296), (6095740.387, 26.540), (6095860.411, 26.844);
- (6095980.435, 27.362), (6096100.460, 27.758), (6096220.484, 28.291), (6096340.508, 29.203), (6096350.086, 29.203).

Their trapezoidal mean over the 969.771 s is 27.046326 C.

The four extra cells (unregistered) exercise D50:

| Cell | rows logged | D50 feature | total discharge Ah (registered rule's share) | checks (dropped) | L (EFC) | measure (a), K | measure (b), C | T_cyc, C | units |
|---|---|---|---|---|---|---|---|---|---|
| NMC 35 C 0-100 0.5/1C a | 103,615 | 213 late records | 1955.243 (4.53 %) | 7 (0) | 162.18508096405245 | 1.2478123915974648 | 35.12291205268629 | 36.74001411246633 | 179 |
| LFP 35 C 0-100 0.5/1C c | 382,258 | 1,217 late records; no cell temperature for cycles 1 to 257 (the index then jumps to 259) | 3186.120 (2.80 %) | 16 (0) | 222.08962934987738 | none | none | none | 0 of 237 in the window |
| NMC 25 C 20-80 0.5/3C a | 332,662 | 4,313 rows of another record dropped | 5150.135 (100 %) | 13 (0) | 493.29412365591236 | 7.738506361271583 | 24.886931140049207 | 27.843582611792176 | 881 |
| NCA 25 C 0-100 0.5/0.5C b | 228,011 | capacity collapses after about 385 EFC; one interrupted check at 200 EFC | 1429.422 (100 %) | 914 (436) | 87.20171875000004 | none (no units) | none | none | 0 |

All four agree with the pipeline to 2.8e-14 EFC and 7.1e-15 K. LFP 35 C 1C c has 237 units inside its lifetime, and none has a temperature. Its first unit with one starts at 241.41 EFC (cycle 262), after L = 222.09 EFC.

## D50 counts

| Count | mine | pipeline (stage3_prepare.json) | D50 text |
|---|---|---|---|
| late records | 9,776 | 9,776 | 9,776 |
| rows of another record dropped | 4,313 | 4,313 | 4,313 |
| discharges | 279,973 | 279,973 | 279,973 |
| discharges where counter_ah differs from the per-sample sum | 25,412 | 25,412 | 25,412 |
| discharge-counter decreases in the ordered rows | 275,861 | 275,861 | 275,861 |
| of which without a cycle-index change | 2 | 2 | 2 |
| cells where counter_ah differs | 18 | 18 | 18 |

Per-cell counts agree in all 86 cells: no mismatch in any of the seven counts above. Discharge ampere-hours agree to 9.1e-13 Ah per cell, and the counter_ah totals agree exactly.

The late records sit in all 18 cells at 35 C and in the two NCA cells at 15 C and 2C. There are 1,099 to 1,370 per LFP cell at 35 C and 143 to 271 in each NCA and NMC cell at 35 C. They lag the first record of the next cycle by at most 0.417 s at 35 C, and by 60.09 and 60.08 s in the two 15 C cells.

The 4,313 rows of cycles 1211 to 1229 lie within 8,750 rows that carry the cell's own cycles 1253 to 1267. No other file has a row more than one cycle below its running maximum.

The two decreases without a cycle-index change are both in rests of the two NCA 15 C 2C cells. They come from the 60 s late records, whose counters read zero.

The undercount claim holds. On the ordered rows, counter_ah counts the following shares of the discharge ampere-hours:

- 2.78 to 2.89 % in the LFP cells at 35 C;
- 4.48 to 4.56 % in the NMC cells at 35 C and 1C;
- 41.94 and 45.07 % in the NMC cells at 35 C and 2C;
- 65.58 and 58.03 % in NCA 35 C 1C c and d.

These match D50's 2.8 to 2.9, 4.5 to 4.6, 42 to 45 and 58 to 66 % to rounding. In LFP 35 C 1C a, for example, the per-sample sum is 3083.961 Ah and counter_ah gives 89.023 Ah. On the rows as logged, counter_ah gives 1150.516 Ah, which is the inventory's efc_total of 1045.92 EFC; the true value is 2803.6 EFC. NMC 35 C 1C a has 1955.243 Ah against 88.581 Ah. counter_ah agrees with the per-sample sum in every discharge of the other 68 cells. The per-sample sum lies between 0.9949 and 1.4783 times the trapezoidal integral of the step's own samples in every discharge of five samples or more, against D50's 0.995 to 1.48.

## Comparison tables

"Frozen primary" is my frozen output under my own readings. "Pipeline readings" is my frozen code run after the freeze with the pipeline's two readings.

### Per cell, 86 cells

| Column (mine vs pipeline) | Frozen primary: max abs (cell) | max rel | Pipeline readings: max abs | max rel |
|---|---|---|---|---|
| L_efc | 1.176e-02 (NMC 25 C 20-80 0.5/0.5C a, F1) | 2.843e-04 | 2.024e-11 | 1.072e-14 |
| resolution_efc | 5.684e-14 | 3.073e-16 | 5.684e-14 | 3.073e-16 |
| censored | 0 | 0 | 0 | 0 |
| first check capacity | 8.882e-16 | 3.008e-16 | 8.882e-16 | 3.008e-16 |
| checks found | 6 (NMC 25 C 0-100 0.5/0.5C a, F1) | 1.087e-02 | 0 | 0 |
| checks dropped by the half rule | 0 | 0 | 0 | 0 |
| efc_total | 6.421e-10 | 2.376e-13 | 6.421e-10 | 2.376e-13 |
| measure (a), swing_K | 3.136 (NCA 25 C 2C a, F2) | 2.383e-01 | 0 | 0 |
| measure (b), rest_T_C | 6.374e-01 (NCA 25 C 2C a, F2) | 2.043e-02 | 0 | 0 |
| T_cyc_C | 7.760e-01 (NCA 25 C 2C a, F2) | 2.662e-02 | 7.105e-15 | 1.938e-16 |
| cell minus environment | 7.953e-01 (NCA 25 C 2C a, F2) | 1.918e-01 | 1.388e-17 | 6.594e-16 |
| units used | 300 (NCA 35 C 2C a, F2) | 1 | 0 | 0 |
| units skipped for a missing step mean | 1 (LFP 35 C 1C c) | 4.237e-03 | 0 | 0 |

In the pipeline-readings run, every check list agrees in all 86 cells: same count, EFC to 1.1e-13, capacity to 8.9e-16, and the same kept flags. The pipeline leaves units_with_rest_record blank where I write 0, in the six cells without units.

### Levers (D44)

| Secant | s mine (frozen) | s pipeline | rel. difference | 95 % interval, mine | 95 % interval, pipeline | pipeline-readings run, max rel over all lever fields |
|---|---|---|---|---|---|---|
| 1C LFP 25-35 C | not formed (T_35 undefined) | not determinable | | | | |
| 1C LFP 15-25 C | 0.0337184 | 0.0337184 | 9.5e-15 | 0.02949 to 0.03827 | 0.02949 to 0.03827 | 1.7e-14 |
| 1C NCA 25-35 C | -0.0203719 | -0.0203732 | 6.5e-05 | -0.02947 to -0.01201 | -0.02947 to -0.01201 | 1.5e-14 |
| 1C NCA 15-25 C | -0.0131514 | -0.0131514 | 7.8e-15 | -0.05220 to 0.00742 | -0.05220 to 0.00742 | 1.0e-14 |
| 1C NMC 25-35 C | -0.0244324 | -0.0244322 | 7.3e-06 | -0.04331 to -0.00445 | -0.04331 to -0.00445 | 3.8e-14 |
| 1C NMC 15-25 C | -0.0562271 | -0.0562271 | 0 | -0.08435 to -0.03100 | -0.08435 to -0.03100 | 3.7e-15 |
| 2C LFP 25-35 C | 0.1100562 | 0.1100582 | 1.8e-05 | 0.08302 to 0.13729 | 0.08302 to 0.13730 | 6.7e-15 |
| 2C LFP 15-25 C | 0.0204987 | 0.0204987 | 5.1e-14 | 0.01184 to 0.02908 | 0.01184 to 0.02908 | 1.1e-13 |
| 2C NCA 25-35 C | 0.0133076 | 0.0129948 | 2.4e-02 | 0.00723 to 0.01888 | 0.00737 to 0.01775 | 1.3e-14 |
| 2C NCA 15-25 C | -0.0247593 | -0.0229705 | 7.8e-02 | -0.03347 to -0.01645 | -0.03347 to -0.01428 | 7.3e-15 |
| 2C NMC 25-35 C | -0.0740104 | -0.0714182 | 3.6e-02 | -0.09763 to -0.05074 | -0.09494 to -0.04860 | 7.8e-16 |
| 2C NMC 15-25 C | -0.0440014 | -0.0440014 | 0 | -0.05069 to -0.03689 | -0.05069 to -0.03689 | 2.7e-15 |

The pipeline-readings run reproduces every lever field to 1.1e-13 relative or better. That covers s, the group temperatures and mean ln L, the setpoint s, n, both laws at the logged temperatures and at the setpoints, and both interval ends. The intervals match draw for draw, which confirms the bootstrap order, including that the LFP 1C 25-35 C secant draws nothing.

Where the inputs are identical, as in the four 15-25 C secants and LFP 2C, my frozen intervals match too. The pipeline's values, which both of my runs reproduce for the 1C series, are:

| Secant | T_lo C | T_hi C | setpoint s | Law 1 | Law 2 | Law 1, setpoints | Law 2, setpoints |
|---|---|---|---|---|---|---|---|
| 1C LFP 25-35 C | 24.3460 | undefined (35.0863 over three cells) | 0.13819 | | | | |
| 1C LFP 15-25 C | 15.9611 | 24.3460 | 0.02827 | -0.05108 | 0.01259 | -0.05108 | 0.01091 |
| 1C NCA 25-35 C | 27.2610 | 38.1952 | -0.02228 | -0.01452 | 0.04863 | -0.02102 | 0.04453 |
| 1C NCA 15-25 C | 18.4527 | 27.2610 | -0.01158 | -0.04794 | 0.02381 | -0.05724 | 0.00717 |
| 1C NMC 25-35 C | 26.1927 | 36.7726 | -0.02585 | -0.01729 | 0.04692 | -0.02102 | 0.04453 |
| 1C NMC 15-25 C | 17.9901 | 26.1927 | -0.04612 | -0.05197 | 0.02047 | -0.05724 | 0.00717 |

The KIT levers beside NMC are read unchanged: -0.006236 and -0.026850 per K from 25 to 40 C, and -0.031041 and -0.033410 from 10 to 25 C. So is the KIT floor (47 pairs, mean 0.073323, median 0.040902, 95th percentile 0.204751).

### Floor (D45)

| Chemistry | pairs (groups) | mean, mine | mean, pipeline | median | 95th percentile, mine | 95th percentile, pipeline | below resolution | pipeline-readings run, max rel |
|---|---|---|---|---|---|---|---|---|
| LFP | 30 (10) | 0.175359 | 0.175359 | 0.080220 | 0.89471 | 0.89471 | 25 of 30 | 4.1e-15 |
| NCA | 24 (9) | 0.107545 | 0.107545 | 0.090254 | 0.25548 | 0.25548 | 23 of 24 | 4.3e-15 |
| NMC | 32 (12) | 0.374767 | 0.374782 | 0.137908 | 1.53971 | 1.53980 | 28 of 32 | 2.8e-15 |

No pair is censored. All 86 pairs agree in difference, resolution and below-resolution flag (pipeline-readings run: 7.0e-13 relative).

### Thermal share (D46)

Sign test. Every count, share and p value below is the same in both of my runs and in the pipeline.

| Subset | measure | pairs | ties | untied | warmer shorter | share | p one-sided | p two-sided |
|---|---|---|---|---|---|---|---|---|
| all | (a) | 60 | 59 | 1 | 0 | 0 | 1 | 1 |
| 1C or below | (a) | 36 | 36 | 0 | 0 | n/a | n/a | n/a |
| 2C and 3C | (a) | 24 | 23 | 1 | 0 | 0 | 1 | 1 |
| LFP | (a) | 20 | 19 | 1 | 0 | 0 | 1 | 1 |
| NCA | (a) | 16 | 16 | 0 | 0 | n/a | n/a | n/a |
| NMC | (a) | 24 | 24 | 0 | 0 | n/a | n/a | n/a |
| one pair per group | (a) | 23 | 22 | 1 | 0 | 0 | 1 | 1 |
| all | (b) | 60 | 59 | 1 | 1 | 1 | 0.5 | 1 |
| 1C or below | (b) | 36 | 36 | 0 | 0 | n/a | n/a | n/a |
| 2C and 3C | (b) | 24 | 23 | 1 | 1 | 1 | 0.5 | 1 |
| LFP | (b) | 20 | 19 | 1 | 1 | 1 | 0.5 | 1 |
| NCA | (b) | 16 | 16 | 0 | 0 | n/a | n/a | n/a |
| NMC | (b) | 24 | 24 | 0 | 0 | n/a | n/a | n/a |
| one pair per group | (b) | 23 | 22 | 1 | 1 | 1 | 0.5 | 1 |

Slopes and ratios, which the pipeline reports for NCA and NMC only (F9):

| Subset | measure | slope, pipeline | interval, pipeline | slope, pipeline-readings run | interval, pipeline-readings run | ratio (c), pipeline | my frozen slope |
|---|---|---|---|---|---|---|---|
| NCA | (a) | -0.808597 | -2.99039 to 1.85358 | -0.808597 | -2.99039 to 1.85358 | 0.29240 | -0.82099 |
| NMC | (a) | 2.300041 | -1.09996 to 8.19030 | 2.300041 | -1.09996 to 8.19030 | 0.20167 | 2.30573 |
| NCA | (b) | 0.665671 | -3.09537 to 0.78601 | 0.665671 | -3.09537 to 0.78601 | 0.37528 | 0.65279 |
| NMC | (b) | -1.085390 | -12.20012 to 4.86860 | -1.085390 | -12.20012 to 4.86860 | 0.05162 | -0.36462 |

The intervals match draw for draw only if the bootstrap draws for these four subsets alone, measure (a) before (b). If every subset drew, the NMC (b) interval would move by up to 0.69. All 60 oriented pairs per measure agree in warmer cell, dT, predicted and observed ln ratio and outcome (pipeline-readings run: 6.9e-13 relative). In my frozen run, one measure-(b) orientation differs: NCA 35 C 2C a and b, a tie (F2).

### Module (D47) and P9

| Quantity | mine (frozen) | pipeline | rel. difference | pipeline-readings run |
|---|---|---|---|---|
| abs(s_NCA) x 0.391279 K | 0.007971 | 0.007972 | 6.5e-05 | 1.5e-14 or better |
| interval | 0.004697 to 0.011531 | 0.004699 to 0.011531 | 2.3e-04 | |
| NCA pairs (16), mean / max | 0.023838 / 0.081064 | 0.026819 / 0.081633 | 0.11 / 7.0e-03 (F2) | |
| all pairs (60), mean / max | 0.016568 / 0.108161 | 0.017373 / 0.108168 | 4.6e-02 / 6.5e-05 (F2) | |

The offer, 0.3912788825035136 K, is the largest dT_module_K with crate at or below 0.5 in `results/module_spread.csv` (h = 10).

### Scorecard (D48)

| Prediction | number it rests on, mine (frozen) | pipeline | my verdict | pipeline verdict |
|---|---|---|---|---|
| P1 | NMC 25-35 C s = -0.024432, in [-0.035, 0.005] | -0.024432 | holds | holds |
| P2 | NMC 15-25 C s = -0.056227, negative, larger than 0.024432 | -0.056227 | holds | holds |
| P3 | NCA 25-35 C s = -0.020372, outside [0, 0.05] | -0.020373 | did not hold | did not hold |
| P4 | LFP 25-35 C s undefined | "a 25 to 35 C lever not determinable" | not scored | not scored |
| P5 | means LFP 0.17536, NCA 0.10754, NMC 0.37477, all above 0.10 | 0.17536, 0.10754, 0.37478 | did not hold | did not hold |
| P6 | 30 of 36 below 1 K = 0.8333; median 0.2913 K | 0.8333 of 36; median 0.2934 K | holds | holds |
| P7 | 0 of 1 untied, 59 ties | 0 of 1, 59 ties | did not hold | did not hold |
| P8 | 0 untied pairs at 1C or below, 1 at 2C and 3C | 0 and 1 | not scored | not scored |
| P9 | 0.007971, below 0.015 | 0.007972 | holds | holds |
| P10 | LFP +0.03372, NCA -0.01315, NMC -0.05623 | the same | did not hold | did not hold |

Every reading option of mine, run end to end, gives the same ten verdicts (`frozen/variants.json`). The options are: one-sided check voltage, no splitting of checks by one-record steps, four unit definitions, missing means for zero-duration steps, and two rest definitions. Across them, the registered NCA lever moves by at most 1.8 % and the NMC lever by at most 6.4 % (both F4).

## The three points

**(i) P4.** LFP 35 C 0-100 0.5/1C c is uncensored, with L = 222.09 EFC. Its file has no cell temperature for cycles 1 to 257, after which the cycle index jumps to 259, so none of the 237 units inside its lifetime has a step mean. D43 makes a group's temperature "the mean T_cyc of its uncensored cells", so T_35 is undefined, and the logged lever of D44 cannot be formed. D48 then says a prediction resting on a lever that is not determinable is not scored. Read in the ordinary sense of "cannot be determined", the pipeline's verdict follows.

The reading is not forced, though. D44 defines "not determinable" by censoring alone, and no cell is censored, so the scorecard's reason names the wrong rule. D43 also handles missing temperature elsewhere by skipping, not by voiding ("a unit with a missing step mean is skipped"). The case was not foreseen, and the two natural repairs give opposite verdicts:

- **Mean T_cyc over the cells that have one.** The pipeline's own lever file lists this T_35, 35.0863 C over cells a, b and d, beside the undefined lever. With it, s = 0.128661 per K and Law 2 over the same secant at the 95 % threshold is 0.044007 per K, a ratio of 2.924. LFP exceeds NCA (-0.0204) and NMC (-0.0244), and the signs agree, so P4 holds. Its margin to the factor-3 edge is 0.08, and an unregistered bootstrap of this s (0.099 to 0.176 per K) spans the edge.
- **The 10 K setpoint.** s = 0.138187 per K, and the ratio to Law 2 over 25 to 35 C is 3.094 (3.140 against Law 2 at the logged temperatures), so P4 does not hold. LFP still exceeds NCA (-0.0223) and NMC (-0.0258) at their setpoints.
- **Dropping the cell from both means** gives s = 0.105382 per K and a ratio of 2.395, so P4 holds.

These values come from the pipeline-readings run. My frozen values differ by under 1e-5, with the same verdicts. `results/stage3_posthoc.json`, opened only after the freeze, gives the same three readings: s = 0.1286613775935446, ratio 2.923658665640619 and "would hold"; setpoint s = 0.13818713692762047, ratio 3.0941492452282153 and "would not hold"; and 0.10538204857447925 without the cell. These agree with my pipeline-readings run to 1.3e-15 relative and with my frozen primary to 9.6e-6.

I agree that "not scored" is a faithful literal reading. It should be reported with the correct cause and with both repairs, since the verdict depends on this reading (F5).

**(ii) P7.** D48 says "P7 holds if the share of D46 (a) over all pairs lies in [0.35, 0.65]". It sets a minimum pair count for P8 only, and its "not scored" clause covers only levers that are not determinable, which P7 does not use. With 59 ties and one untied pair the share is 0 of 1, outside the band, so "did not hold" is what D48 says, and I agree with that reading.

The verdict is degenerate, however. With one pair the share can only be 0 or 1, so P7 could not have held whatever the physics; the sign test gives p = 1.0 both ways. The pair is LFP 35 C 0-100 0.5/2C a and b, at 667.81 and 419.15 EFC, so 248.66 EFC apart against a resolution of 235.25 EFC. Measure (b) gives the opposite share, 1 of 1, on the same pair. Had its resolution been 13.4 EFC larger there would be no untied pair at all. I recomputed the ties independently: they arise because the check spacing (about 100 to 700 EFC) exceeds the replicate differences. Without the tie rule, which D46 does not allow, 24 of the 60 pairs have the warmer cell shorter-lived (share 0.40, two-sided p = 0.155), inside the band. Both of my runs and the post hoc file give this number, so the registered tie rule alone turns the verdict (F6).

**(iii) The 0.5C exclusion.** It is applied consistently. No cell discharged at 0.5C appears in either oriented-pair file (60 pairs each, all at 1, 2 or 3C), in P6 (36 pairs, all 1C), in the D47 pair numbers (60 and 16 pairs) or in any slope. The floor keeps all 86 pairs, as D50 (5) says, and the 0.5C cells' measures are reported with measures_used = False.

Including the 21 partial-window 0.5C pairs that do have measure (a), as D48's P6 wording alone would, raises P6's fraction to 0.895 and changes no verdict. Two small points:

- D43's "carry no swing" holds only for the 0-100 % cells. The partial-window 0.5C cells have a measured swing of -0.15 to -0.75 K, which D50 (5) rightly sets aside.
- `stage3_pairs.csv` labels the 26 pairs at 0.5C "1C or below", although no thermal number uses them.

## Flags

**F1 (low; scored effect at most 1.4e-5 absolute on the NMC floor mean, no verdict).**

- What: D41's "last voltage within 0.02 V of the cut-off" is implemented one-sided (at or below cut-off + 0.02 V); "within" reads two-sided. D41's parenthetical "lowest discharge end voltages 1.99 V" is not true of the files: end voltages reach 1.912 V in NMC 25 C 0-100 0.5/0.5C a and 1.94 to 1.95 V in several NMC and LFP cells, as the inventory's own column shows.
- Where: four cells. NMC 25 C 0-100 0.5/0.5C a has 552 against 546 checks, all six extra after its lifetime. In the three NMC 25 C 20-80 cells of F3, two or three more 1.8 Ah discharges join the spurious checks. Lifetimes move by at most 0.012 EFC.
- Effect: the NMC floor mean is 0.374767 against 0.374782, and its 95th percentile 1.53971 against 1.53980.
- Post hoc file: `results/stage3_posthoc.json`, opened after my freeze, lists the same two-sided lifetimes for these cells as my frozen primary, digit for digit.
- Suggestion: state the one-sided rule in D41 and correct the parenthetical.

**F2 (low for verdicts, medium for unscored numbers).**

- What: D43's unit is "a cycling discharge that is not part of a check and the charge step before it". The pipeline takes a non-check discharge only when the step before it, rests aside, is a charge. My literal primary reading took every non-check discharge with the last charge before it.
- Where they differ: on the one-record discharge tails that end many 2C discharges after a short rest. Inside the lifetime window these number 215 in NCA 25 C 2C a, 300 and 279 in NCA 35 C 2C a and b, and 156 and 199 in NMC 35 C 2C a and b. They also differ on the cycling discharge that follows a one-record "discharge" at the end of a CV charge in the 35 C cells, one to four per cell.
- Registered levers: NCA 25-35 C -0.0203719 against -0.0203732 per K, NMC -0.0244324 against -0.0244322. P9 moves by 6.5e-5 relative. P6's median moves by 0.0021 K and its fraction is unchanged.
- Unscored numbers: the 2C levers move by up to 7.8 % (NCA 15-25 C -0.02476 against -0.02297). NCA 25 C 2C a's swing is 10.02 against 13.16 K. The D47 mean over the NCA pairs is 2.38 % against 2.68 %. One measure-(b) orientation flips, on a tie.
- Judgement: the pipeline's reading is physically the better one, since a tail is the end of the same discharge, but D43 should say so.

**F3 (medium; changes the NMC floor mean that P5 scores by a factor of two, not the P5 verdict).**

- What: in NMC 25 C 20-80 0.5/0.5C a, b and c, coulomb-limited cycling discharges count as capacity checks. These cells cycle a fixed 1.8 Ah, 60 % of rated, at 0.5C. As the window drifts down, that discharge ends near 2.0 V for one to five cycles. Those discharges meet D41 (0.5C mean current, end voltage at the cut-off) and the half rule (0.61 of the first check). Each becomes a check at 0.61, below the 90 % threshold.
- Effect on lifetimes: they become 68.3, 36.7 and 310.6 EFC. Without those checks (my unregistered sensitivity) they are 406.7, 489.1 and 392.1 EFC. Both codes follow D41 as written.
- Effect on the floor: the NMC mean is 0.375 with these checks and 0.178 without. P5 fails either way, since LFP (0.175) and NCA (0.108) are above 0.10 as well. These cells do not enter the thermal share.
- Post hoc file: `results/stage3_posthoc.json`, opened after my freeze, reports the same sensitivity. Its lifetimes without these checks and its NMC floor mean of 0.17823 agree with mine to 1.6e-15.
- Suggestion: report the NMC floor both ways, and exclude discharges that end at the programmed cycling capacity from the check rule.

**F4 (low to medium; no verdict change under either reading).**

- What: one- and two-record discharge steps, about -0.08 A at the end of a CV charge (2 to 10 per cell at 35 C), split capacity checks, because D41 groups consecutive qualifying discharges. Both codes split them, as the inventory did. This leaves checks of one or two discharges 1 to 5 EFC before the rest of their check.
- Where it bites: in NCA 35 C 1C b the crossing falls inside such a split, and the resolution is 1.78 EFC instead of 64.9 EFC.
- Without the splitting, the lifetimes would be:
  - NCA 35 C 1C a: 227.14 instead of 229.69 EFC;
  - NCA 35 C 1C b: 233.63 instead of 234.82 EFC;
  - NMC 35 C 1C b: 154.36 instead of 164.87 EFC;
  - NMC 35 C 2C b: 167.99 instead of 186.15 EFC;
  - NCA 35 C 2C a: 267.23 instead of 264.94 EFC.
- Effect: NCA's registered lever would be -0.02000 instead of -0.02037 per K, and NMC's -0.02287 instead of -0.02443 (6.4 %). P1, P3 and P9 still hold or fail as scored, and the ties are still 59 of 60.
- Suggestion: mention it in the note, or ignore discharge steps of two records or fewer when grouping checks.

**F5 (high; the P4 verdict depends on an unforeseen reading).**

- What: see point (i). Strict reading: not scored. Mean T_cyc over the cells that have one: holds, with a ratio of 2.924. 10 K setpoint: does not hold, with a ratio of 3.094.
- Labelling: the scorecard's reason ("a 25 to 35 C lever not determinable") and `stage3_lever.csv` (determinable = False) cite D44's term, but D44's criterion, censoring, is not what failed. `stage3_lever.json` lists T_group_C at 35 C as 35.0863 C, a three-cell mean, beside the undefined lever without saying so.
- Suggestion: keep "not scored" with the cause stated correctly, and give both repairs in the note.

**F6 (medium; the verdict is literal but carries no information).**

- What: see point (ii). P7 rests on one untied pair, 13.4 EFC beyond its resolution. Measure (b) gives the opposite share on it.
- Without the tie rule: 24 of 60, share 0.40, inside the band. That reading is not allowed by D46.
- Suggestion: say so in the note, and add that the data show no sign signal either way.

**F7 (low to medium; data quality, no verdict change).**

- What: NCA 15 C 0-100 0.5/1C b's skin temperature shows no within-cycle swing over its whole lifetime window: 0.01 to 0.02 K per unit for cycles 1 to 200, against 3.9 to 4.2 K in its replicate a and 1.6 to 4.1 K in every other NCA 1C cell. After cycle 250 it shows a normal 4 to 5 K swing, while a's falls to 0.3 to 1.7 K after cycle 300. A detached or swapped sensor is likely.
- Effect: its swing of 0.019 K and T_cyc of 15.98 C set the NCA 15 C group temperature at 18.45 C, and so the registered NCA 15-25 C lever, -0.0132 per K. That lever enters P10, which fails on LFP anyway. The NCA 15 C pair's 3.98 K is one of the six P6 pairs at 1 K or more, and P6 holds with or without it.
- Related observations: in LFP 25 C 1C a and c, 70 to 76 % of units swing by less than 0.1 K, against 3 to 5 % in b and d. In NMC 25 C 3C c, 18 % do, against none in its replicates.
- Suggestion: one sentence in the note, and an NCA 15-25 C sensitivity without that cell.

**F8 (none for scored numbers; unscored, V1 only).**

- What: D42 says only "the per-cycle capacity at the cycling rate". The pipeline's sensitivity lifetime follows this rule:
  - capacities of non-check discharges within 0.05C of the cycling rate that end at or below cut-off + 0.02 V;
  - a centred 11-point running median of those capacities;
  - a reference equal to the median of the first ten;
  - linear interpolation in EFC at the discharge start.
- Reproduction: with that rule my code reproduces 60 of 61 cells exactly. NMC 25 C 1C a differs by one discharge spacing (0.867 EFC), because its running median equals the threshold to within 4.4e-16. My primary reading, with no smoothing and the first check as reference, differs by up to 327 EFC.
- Suggestion: state the rule.

**F9 (none for scored numbers; unscored).**

- What: the pipeline gives slopes (b) and ratios (c) only for NCA and NMC, since LFP pairs have no prediction. Its "one pair per group" takes the warmest and coolest cells; D46 fixes neither choice.
- Over the 40 pairs that do have a prediction, the all-pairs slope for measure (a) would be 1.36 under the pipeline's unit reading, with the signs agreeing in 24 pairs (the post hoc file gives the same), and 1.57 under mine (interval -2.05 to 2.71). The sign-test counts are the same either way.
- Suggestion: say in the note which subsets carry a slope.

**F10 (none; wording and representation).**

- D50 (1)'s "60 s for the single such record in each 15 C cell" is the lag behind the first record of the next cycle. That record follows the row before it by 0.064 s, lies inside a rest and reads 34.4 C against 47.6 C around it.
- D42's reason for the 90 % thresholds, that the NCA and NMC temperature-series cells "all reached 0.80 or lower", is not quite right. In 9 of the 32 such cells the lowest check is between 0.8001 and 0.8427 of the first. All fall below 0.90, so the threshold and every lifetime are unaffected.
- The pairs file labels 0.5C pairs "1C or below".
- The cells file leaves units_with_rest_record blank for cells without units.

## Reproduction

From `verification/stage3/independent/`:

1. `python3 -I s0_cache.py <cache folder>` (set S3V2_CACHE to that folder for the later steps).
2. `python3 -I s5_freeze.py` runs steps 1 to 4 and writes `frozen/`; a fresh run from a copy of the scripts reproduced all ten frozen outputs byte for byte.
3. After the freeze: `python3 -I s8_assemble.py` writes `../v2_independent.json`, running `s6_compare_cells.py` and `s7_compare_groups.py`.
4. `python3 -I s9_tables.py` prints the tables above from the record.

`s2_v2_purepython.py <out.json> <cell ids>` is the plain-Python recomputation of any cell.
