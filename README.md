# The thermal lever against the duplicate-cell floor

How large a temperature spread a lifetime divergence needs, and how large a module offers.

Author: Siddharth Satte, independent researcher, Pune, India (ORCID 0009-0001-3584-1008).

Nominally identical cells cycled under one protocol at one temperature still drift apart. In the Stanford dynamic-cycling dataset (Geslin et al., Nature Energy 2025; 92 SiOx-graphite/NCA 21700 cells, 47 protocols in duplicate, all at 35 C) the duplicate cells' lifetimes to 90 % state of health differ by 3.0 % on average and by 8.7 % at the 95th percentile, and the worst of eight identical cells falls 4.1 % below the module mean by selection alone (the floor, doi.org/10.5281/zenodo.23205413). A thermal explanation of cell-to-cell divergence in a module has to exceed that floor. This repository asks two questions before any module experiment or model is built:

1. The lever: how large a steady temperature difference between two cells does a published temperature-dependent ageing law need to produce a lifetime difference as large as the floor?
2. The availability: how large a temperature difference does a module at 35 C actually offer, cell to cell in a lumped network of eight 21700 cells, and inside one resolved cell?

## Result

At the dataset's own conditions (35 C, average rates C/16 to C/2, convective cooling) the thermal lever of a module stays inside the floor under both ageing laws used.

- The lever. Law 2, the full LiFePO4/graphite model of Schimpe et al. 2018, needs 1.68 K of steady temperature difference at 35 C for the 95th percentile (8.7 %), 0.79 K for the worst of eight (4.1 %) and 0.58 K for the mean (3.0 %). Law 1, the printed cycle counts of Kucinskis et al. 2022 for an NCA+LNO / Si-graphite 21700, the closest published match to the dataset's cells, needs 8.8 K with the opposite sign: at a C/2 charge its cycle life to 90 % SoH is 102, 126 and 124 cycles at 25, 35 and 45 C, so near 35 C a hotter cell lives slightly longer.
- The availability. An eight-cell LG M50 module in 35 C air offers at most 0.39 K between cells at the dataset's rates, and a resolved cell 0.09 K inside itself. The implied lifetime divergence is at most 2.0 % under Law 2 and 0.4 % under Law 1, below the 3.0 % mean floor. The module first offers Law 2's 1.68 K at 2C, at every convection coefficient from 5 to 25 W/(m2 K) (4.8 to 5.4 K, with the cells 20 to 52 K above the air); 1C gives 1.4 to 1.5 K; no rate up to 2C reaches Law 1's 8.8 K.
- What it means. This does not show that the dataset's own divergence is not thermal: its cells sat in separate test fixtures whose temperatures are not reported, and 0.6 K between two fixtures would suffice under Law 2. It shows that in a module at the dataset's conditions a thermal share would stay inside the floor, where a module test could not separate it from cell-to-cell chance, and that the measurement the thesis needs first is the dataset cell's own cycle life against temperature near 35 C at its C/2 charge, because the two laws disagree there in sign. The result is the lever and the availability, not a measurement of the thermal share.

The one-page note is `note/thermal_lever_note.pdf`; the figure is `results/thermal_lever_figure.png`.

## Pre-registration

The predictions in [predictions.md](predictions.md) were committed in `24ab8d31cd18aa4fe2a010d89f82271652c2d024` (10 October 2026, 07:25 UTC), before any script computed a lever at 35 C, a module spread at a registered rate and cooling, or any field of the resolved cell. This repository was created and that commit pushed at 08:01:59 UTC the same day; every run of V3 and of the registered steps came after the push. One commit lies between the predictions and the push, `7791192`, which changed the linear solver of the resolved cell's steady state after a direct solve ran out of memory on the refined mesh (DECISIONS.md D21), and added the comparison and figure scripts and the release metadata; it was tested with zero heat only. The first commit holds the brief's script as received.

After the push, one change was made to a check and none to a registered calculation: the first run of V3 returned a diverged solution in its analytic radial check (e), a breakdown of that check's linear solve at round-off, so the check was moved to the solver of the main steady solve and V3 was rerun in full (D22; the first run's output is kept in `verification/`). The figure layout and the note were written after the runs.

## Prediction scorecard

| | Prediction | Result | Computed (central estimate) |
|--|--|--|--|
| P1 | Delta T_floor(8.7 %) at 35 C lies between 1 and 3 K for every law | Failed | Law 2: 1.68 K (1.7 K); Law 1: -8.81 K (-8.8 K) |
| P2 | At C/2 with h = 10 W/(m2 K) the eight-cell spread is below 1 K | Held | 0.391 K (0.39 K); 0.368 K at h = 5, 0.348 K at h = 25 |
| P3 | At C/2 the resolved cell's spread is below 0.5 K | Held | 0.092 K at h = 5, 10 and 25 (0.12 K) |
| P4 | The module first exceeds Delta T_floor(8.7 %) at 1C or higher with h = 10, and one rate lower with h = 5 | Failed | First clause held: 2C for Law 2, no rate up to 2C for Law 1 (as estimated); second clause failed: h = 5 gives the same rate as h = 10 for both laws (as estimated) |
| P5 | The two laws' Delta T_floor(8.7 %) differ by less than a factor of 2 | Failed | A factor of 5.24 to 5.26, with opposite signs (5.2) |

P1 and P5 failed because Law 1, the chemistry nearest the dataset's, sits near its optimum temperature at 35 C. The predictions file expected both failures, and the failure of P4's second clause, before any lever was computed, because its central estimates were hand calculations from the copied law values and a two-node reduction of the network; as it says, their agreement with the scripts tests the implementation more than the physics. Scores are in `results/prediction_scorecard.json`.

## Results

Table 1. The lever at 35 C. s = -d ln L/dT; Delta T_floor = f/s (D17; the exact pair inversion, in `results/lever.json`, differs by less than 0.1 %). Law 1 has no discharge-rate dependence (the paper tested one discharge rate), so its lever is the same at every rate.

| Law | Average rate | s (1/K) | L at 35 C (cycles to 90 %) | Delta T_floor, 3.0 % (K) | 4.1 % (K) | 8.7 % (K) |
|---|---|---|---|---|---|---|
| Law 1 | C/16 to C/2 | -0.0099 | 126 | -3.02 | -4.14 | -8.81 |
| Law 2 | C/16 | +0.0520 | 1196 | +0.58 | +0.79 | +1.68 |
| Law 2 | C/10 | +0.0520 | 1629 | +0.58 | +0.79 | +1.68 |
| Law 2 | C/5 | +0.0519 | 2395 | +0.58 | +0.79 | +1.68 |
| Law 2 | C/2 | +0.0519 | 3458 | +0.58 | +0.79 | +1.68 |

Table 2. Heat per cell, (c x 5 A)^2 x 18.71 mOhm, and the spreads (K): hottest minus coldest cell of the lumped module, hottest minus coldest point of the resolved cell. In the module the hottest cell is always an inner cell and the coldest a corner cell.

| Rate | Heat per cell (mW) | Module, h = 5 | Module, h = 10 | Module, h = 25 | Cell, h = 5 | Cell, h = 10 | Cell, h = 25 | Cell hottest minus mean, h = 10 |
|---|---|---|---|---|---|---|---|---|
| C/16 | 1.83 | 0.0058 | 0.0062 | 0.0055 | 0.0014 | 0.0014 | 0.0014 | 0.0007 |
| C/10 | 4.68 | 0.0150 | 0.0158 | 0.0140 | 0.0037 | 0.0037 | 0.0037 | 0.0019 |
| C/5 | 18.7 | 0.060 | 0.063 | 0.056 | 0.015 | 0.015 | 0.015 | 0.0075 |
| C/2 | 117 | 0.368 | 0.391 | 0.348 | 0.092 | 0.092 | 0.092 | 0.047 |
| 1C | 468 | 1.41 | 1.52 | 1.38 | 0.368 | 0.368 | 0.367 | 0.188 |
| 2C | 1871 | 4.77 | 5.41 | 5.25 | 1.47 | 1.47 | 1.47 | 0.754 |

At 1C and 2C the module's mean temperature is 41 to 49 C and 58 to 84 C, so the lever taken at 35 C no longer describes those cells, and at 2C every cell is above 55 C, beyond both laws' fitted ranges. Those rows show what it would take, not a prediction for a working module.

Table 3. Implied lifetime divergence s x Delta T at h = 10 (%), cell to cell in the module and within the resolved cell.

| Rate | Law 1, module | Law 1, cell | Law 2, module | Law 2, cell |
|---|---|---|---|---|
| C/16 | 0.006 | 0.001 | 0.032 | 0.007 |
| C/10 | 0.016 | 0.004 | 0.082 | 0.019 |
| C/5 | 0.062 | 0.015 | 0.33 | 0.076 |
| C/2 | 0.39 | 0.091 | 2.03 | 0.48 |
| 1C | 1.5 | 0.36 | 7.9 | 1.9 |
| 2C | 5.3 | 1.5 | 28 | 7.6 |

Table 4. Step 3: the lowest average rate at which the spread exceeds the magnitude of Delta T_floor. The answer is the same at h = 5, 10 and 25 W/(m2 K).

| Law | Floor | Module | Resolved cell |
|---|---|---|---|
| Law 1 | 3.0 % mean | 2C | none up to 2C |
| Law 1 | 4.1 % worst of 8 | 2C | none up to 2C |
| Law 1 | 8.7 % 95th percentile | none up to 2C | none up to 2C |
| Law 2 | 3.0 % mean | 1C | 2C |
| Law 2 | 4.1 % worst of 8 | 1C | 2C |
| Law 2 | 8.7 % 95th percentile | 2C | none up to 2C |

Secondary results, none of them scored (`results/comparison.json`):

- Law 1 by charge rate. At the printed charge rates the lever changes sign: s = +0.036 per K at 0.2C (Delta T_floor(8.7 %) = 2.4 K), +0.007 at 0.4C (12.8 K) and -0.027 at 0.6C (-3.3 K). The registered C/2 value lies between the last two, so the lever of this chemistry at 35 C is not known even in sign.
- Law 1 with assumed noise (D16, an assumption, not a stated error). If each printed count carried the floor's per-cell noise averaged over two cells, s at a C/2 charge lies between -0.0117 and -0.0081 per K (95 %), and |Delta T_floor(8.7 %)| between 7.4 and 10.8 K. The hand figure in D16 (0.0013 per K) left out the averaging over the two printed charge rates that make up the C/2 value; the Monte Carlo standard error is 0.0009 per K.
- Law 2 by term. The calendar term, part of the paper's own model of a cycling cell, carries 57 to 78 % of the loss at 35 C from 2C down to C/16. The calendar term alone gives s = 0.052 per K and the cycling terms alone 0.051, so the lever does not hinge on the calendar share. The more precise values printed in the paper's text give the same s within 0.1 %.
- The dataset's C/2 charge (D19). The cycle-averaged heat of a dataset cell, which charges at C/2, is 8.0, 5.0 and 2.5 times the heat at the average rate for C/16, C/10 and C/5. At C/16 the module spread is then 0.044 to 0.049 K and Law 2's implied divergence at most 0.26 %, still below every floor.
- Cold plate (D12). One row of eight on a plate at 35 C under the bases: the C/2 module spread falls to 0.016 to 0.036 K, while a base-cooled resolved cell at h = 10 holds 0.29 K inside itself (1.16 K at 1C), because most of its heat now leaves along the jelly roll's axis.
- Sensitivities at C/2 and h = 10: no strip 0.448 K; two strips per pair 0.347 K; the received radiation form (the four-fold understatement corrected in D10) 0.655 K; the 10 s and 30 s pulse resistances (15.47 and 20.37 mOhm) 0.324 and 0.425 K.
- Exploratory: at 1C and h = 10 the module's cells sit at 44.5 to 46.1 C, on Law 1's steep 45 to 55 C segment, where Law 1 evaluated at each cell's own temperature gives a 14.8 % difference between the hottest and coldest cell. Law 1's lever is local, and small only near 35 C; the segment itself is coarsely resolved, since the printed counts at 55 C reach 95 % and 90 % SoH at the same 31 cycles.

## Verification

- V1, before the predictions. Law 1: the transcription of the printed table reproduces every value and statement quoted in the paper's text. Law 2: 8.6 % capacity loss after 2800 full cycles at 1C and 25 C (the paper: less than 9 %) and 10.4 % at 45 C (the paper: 12 %), inside the paper's 2 point model error. Both laws kept.
- V2, before the predictions. The module network's energy balance closes to round-off (1e-12 relative or better), and with the strips at 1e3 W/K the eight cells agree within 0.0004 K.
- V3, after the push (`verification/v3_resolved_cell.json`). (a) Energy closes within 8e-13 of the source. (b) The resolved cell's surface-mean temperature, which sets its losses, matches a lumped cell with the same losses within 2e-7 of the rise. (c) The brief's literal check, the volume-mean temperature against the lumped cell within 1 % of the rise, fails, as D15 stated before the run: the volume mean lies 1.8, 2.6 and 5.2 % of the rise above the lumped value at h = 5, 10 and 25 (hand estimate 2.4, 3.6 and 7.1 %), because a lumped cell has no internal gradient. (d) One mesh refinement (39,990 to 297,057 nodes) changes the cell's spread by 0.83 %. (e) The roll operator reproduces the analytic steady radial solution within 0.13 %. Numerical acceptance, (a), (b), (d) and (e), passes; every conjugate-gradient solve converged.

## What changed from the brief

- The ageing laws. The brief's examples, Wang et al. 2011 and Waldmann et al. 2014, could not be read (subscription full texts behind a CAPTCHA, which was not bypassed), so two open-access laws whose full texts were read replace them, chosen by a rule fixed before any of their numbers were used (D2 to D4).
- The resistance. Neither Chen et al. 2020 nor a readable copy of O'Regan et al. 2022 reports a full-cell DC resistance of the LG M50, so the resistance is derived with PyBaMM's Doyle-Fuller-Newman model and the ORegan2022 parameter set (D8). It is a model value, not a measurement.
- The radiation term follows the two-cell study's corrected form (D10); the received form is reported as a sensitivity.
- The note does not say that the dataset's divergence is not thermal in origin, which the brief allowed "if it is the case". It is not shown: the dataset's cells were cycled in separate fixtures, not in a module, their temperatures are not reported, and under Law 2 a 0.6 K difference between fixtures would reach the mean floor.

## Inputs

- The floor: copied from the published re-analysis (`inputs/floor.json`, with the file hash).
- Two ageing laws, both open access and read in full (`inputs/SOURCES.md`): Kucinskis et al. 2022 (J. Power Sources 549, 232129), a 5 Ah NCA+LNO / Si-graphite 21700, the closest published chemistry to the dataset's cells, whose printed cycle counts to 90 % SoH against temperature and charge rate are the law (`inputs/kucinskis2022_fig7_cycles.csv`); and Schimpe et al. 2018 (J. Electrochem. Soc. 165, A181), the full temperature-dependent LiFePO4/graphite model with its Table IV parameters (`inputs/schimpe2018_parameters.json`).
- The proxy cell, an LG M50, with the geometry, jelly-roll conductivities and can of the two-cell study (doi.org/10.5281/zenodo.23253738), whose resolved scikit-fem cell is reused unchanged in `src/twocell/`.
- The resistance for the I^2 R heat is derived, not measured: 18.71 mOhm, the value whose I^2 R reproduces the steady C/2 cycling heat of PyBaMM's Doyle-Fuller-Newman model with the ORegan2022 parameter set at 35 C (D8; `inputs/lgm50_resistance.json`).

The dataset's cells are a different 21700; the proxy sets magnitudes, not the dataset's exact values.

## Files

- `results/lever.json`: the floors, the laws, s, Delta T_floor (and the exact pair inversion) per law and rate, the lifetime curves from 15 to 55 C, and every module case.
- `results/module_spread.csv`, `results/module_spread_coldplate.csv`: heat per cell, every cell's temperature, the spreads, the hottest and coldest positions and the energy balance, per rate and h.
- `results/resolved_cell.json`, `results/resolved_cell.csv`: the resolved cell per rate and h, with the hottest and coldest points.
- `results/implied_divergence.csv`, `results/comparison.json`: Step 2(c), Step 3 and the secondary results.
- `results/prediction_scorecard.json`, `results/thermal_lever_figure.png` (and `.pdf`), `note/thermal_lever_note.pdf`.
- `verification/`: V1 and V2 (`v1_v2.json`), V3 (`v3_resolved_cell.json`) and V3's first run (`v3_first_run_diverged_check_e.json`).
- `DECISIONS.md`: D1 to D20 fixed before the predictions, D21 before the push, D22 after it.

## How to run

```
python3 -m pip install -r requirements.txt
python3 src/lgm50_resistance.py --out inputs/lgm50_resistance.json   # the derived resistance (minutes)
python3 src/verify_pre.py --out verification                         # V1 and V2
python3 src/resolved_cell_steady.py --verify --out verification      # V3
python3 src/thermal_lever.py --out results                           # Step 1 and Step 2(a)
python3 src/resolved_cell_steady.py --out results                    # Step 2(b)
python3 src/compare.py --results results                             # Step 2(c), Step 3, scoring
python3 src/make_figure.py --results results                         # the figure
cd note && ./render_note.sh                                          # the note (pandoc and Chromium)
```

## Licence

Code: MIT (LICENSE). Note, figure and results: CC BY 4.0. The two-cell study's solver files in `src/twocell/` are MIT, by the same author.

## References

Geslin, A., Xu, L., Ganapathi, D., Moy, K., Chueh, W. C. & Onori, S. Dynamic cycling enhances battery lifetime. Nature Energy 10, 172-180 (2025). https://doi.org/10.1038/s41560-024-01675-8

Kucinskis, G., Bozorgchenani, M., Feinauer, M., Kasper, M., Wohlfahrt-Mehrens, M. & Waldmann, T. Arrhenius plots for Li-ion battery ageing as a function of temperature, C-rate, and ageing state: An experimental study. Journal of Power Sources 549, 232129 (2022). https://doi.org/10.1016/j.jpowsour.2022.232129

Schimpe, M., von Kuepach, M. E., Naumann, M., Hesse, H. C., Smith, K. & Jossen, A. Comprehensive modeling of temperature-dependent degradation mechanisms in lithium iron phosphate batteries. Journal of The Electrochemical Society 165, A181-A193 (2018). https://doi.org/10.1149/2.1181714jes

O'Regan, K., Brosa Planella, F., Widanage, W. D. & Kendrick, E. Thermal-electrochemical parameters of a high energy lithium-ion cylindrical battery. Electrochimica Acta 425, 140700 (2022). https://doi.org/10.1016/j.electacta.2022.140700

Chen, C.-H., Brosa Planella, F., O'Regan, K., Gastol, D., Widanage, W. D. & Kendrick, E. Development of experimental techniques for parameterization of multi-scale lithium-ion battery models. Journal of The Electrochemical Society 167, 080534 (2020). https://doi.org/10.1149/1945-7111/ab9050

Sulzer, V., Marquis, S. G., Timms, R., Robinson, M. & Chapman, S. J. Python Battery Mathematical Modelling (PyBaMM). Journal of Open Research Software 9, 14 (2021). https://doi.org/10.5334/jors.309

Satte, S. (2026). The cell-to-cell lifetime floor in the Geslin et al. dynamic-cycling dataset (v1.0.1). Zenodo. https://doi.org/10.5281/zenodo.23205413

Satte, S. (2026). Where and when a neighbour cell first reaches thermal-runaway onset temperatures: resolved against lumped, a pre-registered two-cell LG M50 pilot (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.23253738
