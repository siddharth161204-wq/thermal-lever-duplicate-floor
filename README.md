# The thermal lever against the duplicate-cell floor

How large a temperature spread a lifetime divergence needs, and how large a module offers.

Author: Siddharth Satte, independent researcher, Pune, India (ORCID 0009-0001-3584-1008).

Nominally identical cells cycled under one protocol at one temperature still drift apart. In the Stanford dynamic-cycling dataset (Geslin et al., Nature Energy 2025; 92 commercial cylindrical SiOx-graphite/NCA cells, 47 protocols in duplicate, all at 35 C) the duplicate cells' lifetimes to 90 % state of health differ by 3.0 % on average and by 8.7 % at the 95th percentile, and the worst of eight identical cells falls 4.1 % below the module mean by selection alone (the floor, doi.org/10.5281/zenodo.23205413). A thermal explanation of cell-to-cell divergence in a module has to exceed that floor. This repository asks two questions before any module experiment or model is built:

1. The lever: how large a steady temperature difference between two cells does a published temperature-dependent ageing law need to produce a lifetime difference as large as the floor?
2. The availability: how large a temperature difference does a module at 35 C actually offer, cell to cell in a lumped network of eight 21700 cells, and inside one resolved cell?

## Result

With the brief's constant-current heat, a module at the dataset's rates (C/16 to C/2, 35 C, convective cooling) stays inside the floor under both ageing laws used. A post-hoc estimate with the dataset's own discharge profiles and CC-CV charge lifts six of its eight dynamic C/2 protocols above the 3.0 % mean floor under Law 2, though not under Law 1's registered slope. Whether a module test at the dataset's conditions could show a thermal share therefore turns on the dataset cell's own temperature sensitivity near 35 C, which neither law fixes.

- The lever. Law 2, the full LiFePO4/graphite model of Schimpe et al. 2018, needs a steady 1.68 K at 35 C for the 95th percentile (8.7 %), 0.79 K for the worst of eight (4.1 %) and 0.58 K for the mean (3.0 %). Law 1, the printed cycle counts of Kucinskis et al. 2022 for an NCA+LNO / Si-graphite 21700, the closest published match found to the dataset's cells, gives 102, 126 and 124 cycles to 90 % SoH at 25, 35 and 45 C at a C/2 charge, so its cycle life peaks near 35 C. Its registered centred slope needs 8.8 K of the opposite sign (a hotter cell lives longer), but between 35 and 45 C, where its 0.4C and 0.6C rows slope oppositely (s = +0.032 and -0.029 per K), a hotter cell lives slightly shorter, and that is where the module's cells sit at the dataset's rates with constant-current heat (35.0 to 38.7 C). Near 35 C its lever is weak and its sign is not fixed.
- The availability (registered). An eight-cell LG M50 module in 35 C air offers at most 0.39 K between cells at the dataset's rates, and a resolved cell 0.09 K inside itself, an implied lifetime divergence of at most 2.0 % under Law 2 and 0.4 % under Law 1, below the 3.0 % mean floor. The module first offers Law 2's 1.68 K at 2C, at every convection coefficient from 5 to 25 W/(m2 K) (4.8 to 5.4 K), with every cell above 55 C, outside both laws' fitted ranges; 1C gives 1.4 to 1.5 K; no rate up to 2C gives Law 1's 8.8 K.
- The dataset's own profiles (post hoc, D24). The brief's heat is that of a constant current at the average rate. The dataset's dynamic discharges carry more: from the mean current and current variance its authors publish for every cell, with the dataset's CC-CV charge at C/2, the cycle-averaged heat of its C/2 protocols is up to 3.05 times the brief's C/2 heat. The module then offers up to 1.17 K (h = 10), and six of the eight dynamic C/2 protocols (Periodic A to E and Synthetic 2a) give 3.5 to 6.1 % under Law 2, above the 3.0 % mean floor; the other two give 2.8 %. The largest values come from Periodic C, whose cells averaged 0.60C and 0.61C rather than C/2; without it the maximum is 5.3 %. Against the hottest-minus-mean spread, D17's match for a shortfall below the module mean, no protocol reaches the 4.1 % worst of eight, and none reaches the 8.7 % 95th percentile. Law 1's registered slope keeps every protocol below 1.2 %; at h = 5, where Periodic C's cells sit just above 45 C on Law 1's steep segment, Law 1 at each cell's own temperature gives 15 %, above every floor (exploratory). At C/5 and below every protocol stays below 2.3 % under both laws.
- What it means. The measurement the thesis needs first is the dataset cell's own cycle life against temperature near 35 C at its C/2 charge: under an Arrhenius-type sensitivity a module at the dataset's dynamic C/2 conditions could show a thermal share above the mean floor; with the registered slope of the closest published law it could not, while that law at the cells' own temperatures could where they pass 45 C (exploratory); near 35 C it does not resolve even its own sign. The dataset's own divergence is neither shown nor ruled out as thermal: its cells sat in separate fixtures in a 35 C chamber, and their own temperatures are not reported. The result is the lever and the availability, not a measurement of the thermal share.

The one-page note is `note/thermal_lever_note.pdf`; the figure is `results/thermal_lever_figure.png`.

## Pre-registration

The predictions in [predictions.md](predictions.md) were committed in `24ab8d31cd18aa4fe2a010d89f82271652c2d024` (10 October 2026, 07:25 UTC), before any script computed a lever at 35 C, a module spread at a registered rate and cooling, or any field of the resolved cell. This repository was created at 08:01:27 UTC and that commit pushed at 08:01:59 UTC the same day; every run of V3 and of the registered steps came after the push. The tag and GitHub pre-release `pre-registration` mark that commit. One commit lies between the predictions and the push, `7791192`, which changed the linear solver of the resolved cell's steady state after a direct solve ran out of memory on the refined mesh (DECISIONS.md D21), and added the comparison and figure scripts and the release metadata; it was tested with zero heat only. The first commit holds the brief's script as received.

After the push:

- D22. The first run of V3 returned a diverged solution in its analytic radial check (e), a breakdown of that check's linear solve at round-off, so the check was moved to the solver of the main steady solve and V3 was rerun in full (the first run's output is kept in `verification/`).
- D23. Two independent checks were run on the results (below). Their corrections change no registered number: Law 2's exact forms are left blank where the cells lie outside its fitted range, as D18 says; the exact pair difference is also given with its sign; V1 for Law 2 is also given with all four terms of the model; the scorecard's keys say "computed"; the worst-of-8 floor is also compared with the hottest-minus-mean spread; the lifetime curves are written at every rate.
- D24. Post-hoc analyses, not pre-registered: the heat of the dataset's own discharge profiles with its CC-CV charge, the heat at which the module reaches Law 2's mean floor, the module's time constants, V1 with all four terms and the air-conductivity sensitivity (`src/posthoc.py`).
- A third check, of these post-hoc analyses and the rewritten text, led to the CC-CV charge being included and to the wording of the note and this README as they now stand.
- The figure layout and the note were written after the runs.

## Prediction scorecard

| | Prediction | Result | Computed (central estimate) |
|--|--|--|--|
| P1 | Delta T_floor(8.7 %) at 35 C lies between 1 and 3 K for every law | Failed, as predictions.md expected | Law 2: 1.68 K (1.7 K); Law 1: -8.81 K (-8.8 K) |
| P2 | At C/2 with h = 10 W/(m2 K) the eight-cell spread is below 1 K | Held | 0.391 K (0.39 K); 0.368 K at h = 5, 0.348 K at h = 25 |
| P3 | At C/2 the resolved cell's spread is below 0.5 K | Held | 0.092 K at h = 5, 10 and 25 (0.12 K) |
| P4 | The module first exceeds Delta T_floor(8.7 %) at 1C or higher with h = 10, and one rate lower with h = 5 | Failed, as expected | First clause held: 2C for Law 2, no rate up to 2C for Law 1 (as estimated); second clause failed: h = 5 gives the same rate as h = 10 for both laws (as estimated) |
| P5 | The two laws' Delta T_floor(8.7 %) differ by less than a factor of 2 | Failed, as expected | A factor of 5.24 to 5.26, of opposite signs (8.8 against 1.7 K, "a factor of about 5") |

P1 and P5 failed because Law 1, the chemistry nearest the dataset's, peaks near 35 C. The predictions file expected both failures, and the failure of P4's second clause, before any lever was computed, because its central estimates were hand calculations from the copied law values and a two-node reduction of the network; as it says, their agreement with the scripts tests the implementation more than the physics. Scores are in `results/prediction_scorecard.json`.

## Results

Table 1. The lever at 35 C. s = -d ln L/dT; Delta T_floor = f/s (D17; the exact pair inversion, in `results/lever.json`, differs by less than 0.1 %). Law 1 has no discharge-rate dependence (the paper tested one discharge rate), so its lever is the same at every rate. The lifetime curves from 15 to 55 C at every rate are in `results/lifetime_curves.csv`.

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

At 1C and 2C the module's mean temperature is 41 to 49 C and 58 to 84 C, so the lever taken at 35 C no longer describes those cells, and at 2C every cell is above 55 C, beyond both laws' fitted ranges. Those rows show what it would take, not a prediction for a working module. An independent finite-volume solve of the resolved cell gives spreads 1.2 % larger (hottest minus mean 2.3 % larger), from the registered base mesh's discretisation (D23).

Table 3. Implied lifetime divergence |s| x Delta T at h = 10 (%), cell to cell in the module and within the resolved cell.

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

The worst-of-8 floor is a shortfall below the module mean, so D17 also compares it with the hottest-minus-mean spread: the module then first offers it at 2C under Law 2 and at no rate up to 2C under Law 1, and the resolved cell at none.

Secondary results, none of them scored (`results/comparison.json`):

- Law 1 by charge rate. At the printed charge rates the lever changes sign: s = +0.036 per K at 0.2C (Delta T_floor(8.7 %) = 2.4 K), +0.007 at 0.4C (12.8 K) and -0.027 at 0.6C (-3.3 K). The registered C/2 value lies between the last two.
- Law 1 at each cell's own temperature (D18, D23). Every module cell is above 35 C, where Law 1's lifetime falls with temperature, slightly up to 45 C and steeply beyond, so the hotter cell lives shorter in every case, against the registered slope: by 0.05 % at C/2 (h = 10). At 1C the cells sit at 40 to 50 C, and the pair difference is 0.2 % (h = 25), 14.8 % (h = 10) and 19.5 % (h = 5), on the steep 45 to 55 C segment, whose 55 C counts reach 95 % and 90 % SoH at the same 31 cycles and rest on one cell. This is exploratory and fragile.
- Law 1 with assumed noise (D16, an assumption, not a stated error). If each printed count carried the floor's per-cell noise averaged over two cells, s at a C/2 charge lies between -0.0117 and -0.0081 per K (95 %), and |Delta T_floor(8.7 %)| between 7.4 and 10.8 K. The standard error is 0.0009 per K; the hand figure in D16 (0.0013) left out the averaging over the two printed charge rates that make up the C/2 value (D23).
- Law 2 by term. The calendar term, part of the paper's own model of a cycling cell, carries 57 to 78 % of the loss at 35 C from 2C down to C/16. The calendar term alone gives s = 0.052 per K and the cycling terms alone 0.051, so the lever does not hinge on the calendar share. The more precise values printed in the paper's text give the same s within 0.1 %. Reading D4's time average of the calendar factor as a root mean square instead of a mean would change s by under 0.05 % and L at 35 C by about 12 %.
- The dataset's C/2 charge (D19). The cycle-averaged heat of a dataset cell, which charges at C/2, is 8.0, 5.0 and 2.5 times the heat at the average rate for C/16, C/10 and C/5 with a constant-current discharge. At C/16 the module spread is then 0.044 to 0.049 K and Law 2's implied divergence at most 0.26 %. With the dataset's CC-CV charge instead of this 2 h constant-current charge the heat is 3 to 8 % lower (D24, post hoc). The steady spread at the cycle-averaged heat is used because the lever is linear to first order; the cells follow each phase (time constants 0.19 to 0.46 h, D24, post hoc), so during each C/2 charge, of about 2 h, the spread approaches its steady C/2 value, 0.35 to 0.39 K, still below every floor with a constant-current discharge.
- The dataset's own discharge profiles (post hoc, D24; `results/dataset_protocol_heat.csv`, `results/posthoc.json`). The authors' released features give each cell's mean discharge current c and normalized current variance NCV, which is the variance over the mean in C-rate units: for the constant-current protocols with a storage rest, two-level profiles, c + NCV equals the measured current of the same rate without storage within 3e-6, and for storage at full charge the rest fraction equals NCV / (c + NCV) within 3e-6; for constant current without rest NCV is zero. The discharge heat is then c^2 + c NCV, 1.7 to 4.3 times that of a constant current for the dynamic C/2 protocols. The dataset's CC-CV charge, from the proxy cell's Doyle-Fuller-Newman model, takes 2.23 h and carries 0.481 h of (I/I_1C)^2 per unit of discharged capacity (D19's constant-current form, 2 h and 0.5, would give 3 to 8 % more heat). The module's largest spreads are then 0.082 K (C/16), 0.17 K (C/10), 0.43 K (C/5) and 1.17 K (C/2) at h = 10, and Law 2's largest implied divergences 0.43, 0.88, 2.2 and 6.1 %. Six of the eight dynamic C/2 protocols exceed the 3.0 % mean floor under Law 2 at every h (Synthetic 2a by 0.1 point at h = 25); against the hottest-minus-coldest spread two or three exceed the 4.1 % worst of eight, against the hottest-minus-mean spread none; none exceeds the 8.7 % 95th percentile. Law 2 at each cell's own temperature, where the dynamic C/2 cells run at 39 to 43 C (h = 10), gives 4 to 7 % more. While each hot dynamic C/2 discharge lasts, the spread approaches its steady value at the discharge heat alone, up to 2.25 K (Periodic C, h = 10), above Law 2's 1.68 K. The open circles in panels (b) and (c) of the figure are these cells at h = 10.
- Cold plate (D12). One row of eight on a plate at 35 C under the bases: the C/2 module spread falls to 0.016 to 0.036 K, while a base-cooled resolved cell at h = 10 holds 0.29 K inside itself (1.16 K at 1C), because most of its heat now leaves along the jelly roll's axis.
- Sensitivities at C/2 and h = 10: no strip 0.448 K; two strips per pair 0.347 K; the received radiation form (D10) 0.655 K; the 10 s and 30 s pulse resistances (15.47 and 20.37 mOhm) 0.324 and 0.425 K; the air conductivity 10 % lower or higher, 0.406 or 0.378 K (post hoc, D24).

## Verification

- V1, before the predictions. Law 1: the transcription of the printed table reproduces every value and statement quoted in the paper's text. Law 2: 8.6 % capacity loss after 2800 full cycles at 1C and 25 C (the paper: less than 9 %) and 10.4 % at 45 C (the paper: 12 %), inside the paper's 2 point model error. The registered V1 left out the model's high-SOC term, as its script says; with all four terms it gives 8.9 % and 10.4 %, still inside (D23, computed post hoc). Both laws kept.
- V2, before the predictions. The module network's energy balance closes to round-off (1e-12 relative or better), and with the strips at 1e3 W/K the eight cells agree within 0.0004 K.
- V3, after the push (`verification/v3_resolved_cell.json`). (a) Energy closes within 8e-13 of the source. (b) The resolved cell's surface-mean temperature, which sets its losses, matches a lumped cell with the same losses within 2e-7 of the rise. (c) The brief's literal check, the volume-mean temperature against the lumped cell within 1 % of the rise, fails, as D15 stated before the run: the volume mean lies 1.8, 2.6 and 5.2 % of the rise above the lumped value at h = 5, 10 and 25 (hand estimate 2.4, 3.6 and 7.1 %), because a lumped cell has no internal gradient. (d) One mesh refinement (39,990 to 297,057 nodes) changes the cell's spread by 0.83 %. (e) The roll operator reproduces the analytic steady radial solution within 0.13 % of the maximum rise (1.6 % of the centre-to-surface difference). Numerical acceptance, (a), (b), (d) and (e), passes; every conjugate-gradient solve converged.
- Independent recomputation (`verification/independent_recomputation/`, with its report). A separate checker, who had not seen the study being built and did not use its code, recomputed 1,336 quantities from DECISIONS.md and the inputs at commit 0e4ee5e: the levers of both laws to 3e-9 or better, all 36 module cases (both arrangements) to 4e-10, every module-based implied divergence, every exceedance flag, Step 3's lowest rates and the scores; and the resolved cell, with an axisymmetric finite-volume solver checked against an analytic series, within 1.2 % (hottest minus coldest) and 2.3 % (hottest minus mean), the registered mesh being the coarser, a difference the cell-based implied divergences inherit. Its two discrepancies, V1's missing term and D16's standard error, are addressed in D23.
- Independent audit, of commit 0e4ee5e (report in `verification/reviews/`). A second checker traced every number then in this README, the note and the release metadata to the results files and found none wrong or mis-rounded; reran the pipeline in a scratch clone with byte-identical or round-off-identical outputs; matched the floor to its source (values and SHA-256, at geslin-floor-pilot v1.0.1); confirmed on GitHub that 24ab8d3 was pushed at 08:01:59 UTC; and checked the references against Crossref and DataCite. Its points on substance and wording (the sign of Law 1, the heat of dynamic profiles, the limits of the note) are taken up in the note and here; the text added for them (D23, D24) came after that audit. The Schimpe et al. Table IV values and the Kucinskis et al. Fig. 7 counts could not be fetched again for it; they were transcribed from the full texts when the inputs were built (D3, D4).
- Check of the post-hoc analysis and the rewritten text, at commit da17095 (report in `verification/reviews/`). A third checker recomputed the post-hoc analysis with its own code, matching it to round-off; confirmed the normalisation of the current variance by the two exact identities above; found that leaving out the CV part of the charge overstated the heat by 3 to 8 %, which is now included; and raised wording points on the sign of Law 1, the labelling of post-hoc numbers and the worst-of-8 comparison, which are taken up here and in the note.

## What changed from the brief

- The ageing laws. The brief's examples, Wang et al. 2011 and Waldmann et al. 2014, could not be read (subscription full texts behind a CAPTCHA, which was not bypassed), so two open-access laws whose full texts were read replace them, chosen by a rule fixed before any of their numbers were used (D2 to D4).
- The resistance. Chen et al. 2020 reports no full-cell DC resistance of the LG M50 and O'Regan et al. 2022 could not be read, so the resistance is derived with PyBaMM's Doyle-Fuller-Newman model and the ORegan2022 parameter set (D8). It is a model value, not a measurement.
- Some property values have no primary table: the air conductivity (0.03 W/(m K)) and the can's conductivity (15 W/(m K)), inherited from the two-cell study, the nickel conductivity behind the brief's strip value, and, for the post-hoc time constants only, the can's mass and specific heat (`inputs/SOURCES.md`).
- The radiation term follows the two-cell study's corrected form (D10); the received form is reported as a sensitivity.
- V3's acceptance. The brief's literal V3 test, the volume mean within 1 % of a lumped cell, cannot pass for a cell with internal heat, as D15 said before the run; acceptance rests on the surface mean, the energy balance, the refinement and the analytic radial check.
- The note does not say that the dataset's divergence is not thermal in origin, which the brief allowed "if it is the case". It is not shown: the dataset's cells were cycled in separate fixtures, not in a module, their own temperatures are not reported, and under Law 2 a 0.6 K difference between fixtures would reach the mean floor. The post-hoc profile estimate points the same way.

## Inputs

- The floor: copied from the published re-analysis (`inputs/floor.json`, with the file hash).
- Two ageing laws, both open access and read in full (`inputs/SOURCES.md`): Kucinskis et al. 2022 (J. Power Sources 549, 232129), a 5 Ah NCA+LNO / Si-graphite 21700, the closest published chemistry found to the dataset's cells, whose printed cycle counts to 90 % SoH against temperature and charge rate are the law (`inputs/kucinskis2022_fig7_cycles.csv`); and Schimpe et al. 2018 (J. Electrochem. Soc. 165, A181), the full temperature-dependent LiFePO4/graphite model with its Table IV parameters (`inputs/schimpe2018_parameters.json`).
- The proxy cell, an LG M50, with the geometry, jelly-roll conductivities and can of the two-cell study (doi.org/10.5281/zenodo.23253738), whose resolved scikit-fem cell is reused unchanged in `src/twocell/` (its comments cite that study's own decision numbers).
- The resistance for the I^2 R heat is derived, not measured: 18.71 mOhm, the value whose I^2 R reproduces the steady C/2 cycling heat of PyBaMM's Doyle-Fuller-Newman model with the ORegan2022 parameter set at 35 C (D8; `inputs/lgm50_resistance.json`).
- Post hoc only: the dataset authors' released per-cell discharge features (github.com/geslina/dynamic_cycling_Nature_Energy_2024 at commit 5b2f7f0), read from a local clone with hash checks. The upstream files are not included; the two features used, mean discharge current and normalized current variance, appear per cell, attributed, in `results/dataset_protocol_heat.csv`.

The LG M50 is a 21700 but a different cell from the dataset's; the proxy sets magnitudes, not the dataset's exact values.

## Files

- `results/lever.json`: the floors, the laws, s, Delta T_floor (and the exact pair inversion) per law and rate, the lifetime curves at C/2, and every module case. `results/lifetime_curves.csv`: the lifetime curves from 15 to 55 C at every rate.
- `results/module_spread.csv`, `results/module_spread_coldplate.csv`: heat per cell, every cell's temperature, the spreads, the hottest and coldest positions and the energy balance, per rate and h.
- `results/resolved_cell.json`, `results/resolved_cell.csv`: the resolved cell per rate and h, with the hottest and coldest points.
- `results/implied_divergence.csv`, `results/comparison.json`: Step 2(c), Step 3, the secondary results and the module network's constants.
- `results/prediction_scorecard.json`, `results/thermal_lever_figure.png` (and `.pdf`), `note/thermal_lever_note.pdf`.
- `results/posthoc.json`, `results/dataset_protocol_heat.csv`: the post-hoc analyses of D24.
- `verification/`: V1 and V2 (`v1_v2.json`), V3 (`v3_resolved_cell.json`), V3's first run (`v3_first_run_diverged_check_e.json`), the independent recomputation with its scripts and report, and the reports of the audit and of the check of the post-hoc analysis (`reviews/`).
- `DECISIONS.md`: D1 to D20 fixed before the predictions, D21 before the push, D22 to D24 after it.

## How to run

```
python3 -m pip install -r requirements.txt
python3 src/lgm50_resistance.py --out inputs/lgm50_resistance.json   # the derived resistance
python3 src/verify_pre.py --out verification                         # V1 and V2
python3 src/resolved_cell_steady.py --verify --out verification      # V3
python3 src/thermal_lever.py --out results                           # Step 1 and Step 2(a)
python3 src/resolved_cell_steady.py --out results                    # Step 2(b)
python3 src/compare.py --results results                             # Step 2(c), Step 3, scoring
git clone https://github.com/geslina/dynamic_cycling_Nature_Energy_2024 upstream
git -C upstream checkout 5b2f7f04d05072fe1f9bd8af664f23eadbde6317
python3 src/posthoc.py --upstream upstream --results results         # post hoc (D24)
python3 src/make_figure.py --results results                         # the figure
cd note && ./render_note.sh                                          # the note (pandoc and Chromium)
```

## Licence

Code: MIT (LICENSE). Note, figure and results: CC BY 4.0. The two-cell study's solver files in `src/twocell/` are MIT, by the same author. The two dataset features used post hoc (Geslin et al.) appear per cell in `results/dataset_protocol_heat.csv` with attribution; the dataset is CC BY 4.0 at the Stanford Digital Repository (doi.org/10.25740/td676xr4322), and the authors' code repository carries no licence file of its own.

## References

Geslin, A., Xu, L., Ganapathi, D., Moy, K., Chueh, W. C. & Onori, S. Dynamic cycling enhances battery lifetime. Nature Energy 10, 172-180 (2025). https://doi.org/10.1038/s41560-024-01675-8

Kucinskis, G., Bozorgchenani, M., Feinauer, M., Kasper, M., Wohlfahrt-Mehrens, M. & Waldmann, T. Arrhenius plots for Li-ion battery ageing as a function of temperature, C-rate, and ageing state: An experimental study. Journal of Power Sources 549, 232129 (2022). https://doi.org/10.1016/j.jpowsour.2022.232129

Schimpe, M., von Kuepach, M. E., Naumann, M., Hesse, H. C., Smith, K. & Jossen, A. Comprehensive modeling of temperature-dependent degradation mechanisms in lithium iron phosphate batteries. Journal of The Electrochemical Society 165, A181-A193 (2018). https://doi.org/10.1149/2.1181714jes

O'Regan, K., Brosa Planella, F., Widanage, W. D. & Kendrick, E. Thermal-electrochemical parameters of a high energy lithium-ion cylindrical battery. Electrochimica Acta 425, 140700 (2022). https://doi.org/10.1016/j.electacta.2022.140700

Chen, C.-H., Brosa Planella, F., O'Regan, K., Gastol, D., Widanage, W. D. & Kendrick, E. Development of experimental techniques for parameterization of multi-scale lithium-ion battery models. Journal of The Electrochemical Society 167, 080534 (2020). https://doi.org/10.1149/1945-7111/ab9050

Sulzer, V., Marquis, S. G., Timms, R., Robinson, M. & Chapman, S. J. Python Battery Mathematical Modelling (PyBaMM). Journal of Open Research Software 9, 14 (2021). https://doi.org/10.5334/jors.309

Satte, S. (2026). The cell-to-cell lifetime floor in the Geslin et al. dynamic-cycling dataset (v1.0.1). Zenodo. https://doi.org/10.5281/zenodo.23205413

Satte, S. (2026). Where and when a neighbour cell first reaches thermal-runaway onset temperatures: resolved against lumped, a pre-registered two-cell LG M50 pilot (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.23253738
