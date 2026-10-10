Filed by the author with this header. The text below is the checker's report as it returned it; its scripts stayed in its session and are not filed. It checked commit da17095. Its items were taken up in the commits that follow: the dataset's CC-CV charge is now part of the post-hoc heat (so the post-hoc ranges it quotes, which assume a 2 h constant-current charge, are superseded by `results/posthoc.json`; with the charge included the dynamic C/2 protocols give 3.5 to 6.1 % under Law 2 at h = 10, six of eight still above the mean floor at every h); post-hoc numbers are labelled as such; the statements on redistribution and on what the recomputation found are corrected; the worst-of-8 comparison is also made against the hottest-minus-mean spread; the qualifiers on Law 1 are added; the filed copy of the recomputation's README is marked as adapted. Its reading of the paper's Supplementary Table 3 was through a summarising tool and is not cited in the study; the normalisation of the current variance rests on the two exact identities in the data.

Report: independent check of what is new at da17095 (post-hoc D24, note, README)

**Summary.** I recomputed the post-hoc analysis with my own code and it matches the study to round-off. The normalisation the study inferred for NCV is confirmed by an exact identity in the data and by the authors' own definition. The note and README need these fixes:
* Five must-fix items, all cheap wording changes:
  * an unqualified Law 1 claim;
  * post-hoc numbers that are not labelled post hoc;
  * a false "not redistributed" statement;
  * "exactly" overstating what the recomputation found;
  * "its dynamic C/2 protocols" where only six of eight qualify.
* One substantive caveat that is not disclosed: ignoring the CV part of the charge biases the heat upward by about 6 %.

The read-only repo is untouched (git status clean, still da17095). My scripts are in /home/claude/check2/ (posthoc_recompute.py, ncv_alternatives.py, periodic_check.py, cv_phase.py, cv_effect.py, own_temperature.py, misc_checks.py, scan_text.py, load_upstream.py) and their outputs are in /home/claude/check2/out/.

## (1) Post-hoc recomputation (my code; module network from the earlier independent module.py)

| Quantity | Study | Mine | Verdict |
|---|---|---|---|
| Upstream file hashes and cell order | as in posthoc.json | same SHA-256; average current equals metadata avg_crate_exp for all 92 cells | ok |
| Cycle-averaged heat q, 92 cells | CSV | identical (relative difference 0) | ok |
| Heat factor 1 + NCV/c; ratio to D19 heat; ratio to brief heat | CSV | identical | ok |
| Module spread, 276 values | CSV | max relative difference 8.5e-14; hot minus mean 1.3e-12; T_max 1.4e-15 | ok |
| Resolved-cell spread scaled from C/2 | CSV | identical; linear to 0.02 % (2C/C2 ratio 15.997 against 16) | ok |
| Implied divergence | CSV | 8.6e-14 | ok |
| Floor-exceedance flags (1,656) | CSV | 0 differ | ok |
| By-rate summary in posthoc.json | as filed | every key identical | ok |
| Study's own src/posthoc.py rerun in a scratch clone | | byte-identical posthoc.json and CSV | ok |
| Up to 3.3 times the brief's C/2 heat | 3.3 | 3.265 | ok |
| Largest C/2 spread at h = 10 | 1.25 K | 1.2487 K | ok |
| Law 2 range for the six protocols at h = 10 | 3.7 to 6.5 % | 3.68 (Synthetic 2a) to 6.48 % (Periodic C) | ok |
| Synthetic 2b and 2c at h = 10 | (not stated) | 2.96 and 2.94 % | below floor |
| Law 1 maximum | below 1.3 % | 1.235 % | ok |
| C/5 and below, both laws | below 2.4 % | at most 2.32 % | ok |
| Largest spreads by rate | 0.085 / 0.18 / 0.45 K | 0.0848 / 0.1763 / 0.4471 K | ok |
| Law 2 maxima by rate | 0.44 / 0.92 / 2.3 % | 0.441 / 0.917 / 2.32 % | ok |
| Heat factor, dynamic C/2 | 1.7 to 4.3 | 1.73 to 4.30 | ok |
| Protocols above the worst-of-8 floor | "three or four" | 3 / 4 / 3 at h = 5 / 10 / 25 | ok |
| Mean-floor exceedances | six at every h | six at every h | ok |
| C_cell | 54.655 J/K | 54.655 J/K | ok |
| Secant loss conductance G | 0.03283 / 0.04462 / 0.07977 W/K | same | ok |
| Time constant tau | 0.4625 / 0.3402 / 0.1903 h | same; linearised-network slowest mode 0.459 / 0.342 / 0.196 h; difference mode 0.07 to 0.09 h | ok |
| V1 with four terms (my own implementation from the parameters JSON) | 8.91 % and 10.35 % | 8.909029 % and 10.353533 % | ok; inside the acceptance ranges [0, 11] and [10, 14] |
| Air-conductivity sensitivity | 0.4057 / 0.3913 / 0.3778 K | same; +3.69 % and -3.43 % | ok |
| D24's "1.48 times the heat" | 1.48 | 1.4807 at h = 10 only (1.58 at h = 5, 1.66 at h = 25) | add "at h = 10" |

**Algebra of q.** Correct. Per unit capacity:
* discharge window: 1/c h with mean(I^2) = c^2 + c NCV;
* charge: 2 h at 0.25;
* so q = R I1C^2 (c + NCV + 0.5) / (1/c + 2).

**Regenerative charge.** Correctly handled, provided "Average Current" is the signed net mean. My Periodic reconstruction supports that. The C/2 charge restores only the net charge, so this introduces no bias.

**CV omission biases q up.** The dataset charges CC-CV at C/2 to 4.2 V with a 0.05C cut-off over 4.2 to 3.1 V. For the proxy cell (PyBaMM DFN, ORegan2022, 35 C):
* the CV hold adds 0.36 h and 8.3 % of the charge;
* per unit capacity the charge then lasts 2.23 h and carries 0.481 instead of 0.5;
* q falls by 3.3 to 7.3 % across the 92 cells, and 6.1 to 6.9 % for the dynamic C/2 protocols;
* the h = 10 range becomes 3.5 to 6.1 %;
* six protocols still exceed 3.0 % at every h, but Synthetic 2a at h = 25 only just (3.09 %);
* worst-of-8 counts become 3 / 3 / 2.

Rests between steps (unknown) would also bias q up. The missing contact resistance biases it down.

## (2) The NCV normalisation

The normalisation (NCV = Var/mean in C-rate units) is confirmed:
* **Two-level storage profiles.** NCV + c equals the measured CC plateau current (0.099959, 0.199983, 0.499959 C, taken from the no-storage cells) within 5.4e-6 for all 18 storage cells. For the six 100 %-storage cells, "Rest Fraction at High SOC" equals NCV/(c + NCV) within 2e-6. The study's 4 % test against I_max is loose only because I_max includes a measured overshoot.
* **Other normalisations fail by factors** and are inconsistent across rates:
  * Var/mean^2: 8.8 to 16.6 times the published NCV;
  * Var: 0.06 to 0.11 times;
  * Var/I_max: 0.23 to 0.60 times;
  * Var/I_max^2: 0.45 to 6.0 times;
  * std/mean: 4.8 to 20 times;
  * std: 0.54 to 1.23 times;
  * Var in A^2 would need a capacity of 4.07, 3.41 or 2.97 Ah depending on the rate.
* **Authors' definition.** Supplementary Table 3 of Geslin et al. defines Normalized Current Variance as "Variance divided by average current", Average Current as the "overall protocol average current", and Relative Charge Fraction as the fraction of time charging. I read this through a web fetch tool that returns summaries, because direct download was blocked.
* **Periodic A to D.** I rebuilt them from the authors' Supplementary Tables 5 and 6: pulses at 3 or 5 times the average rate, 0.08C charge pulses, discharge-to-charge time ratio 5 or 10, rests otherwise. Taking the net mean:
  * A, B and C at C/10 and C/5 agree within 0.5 to 1 %; at C/2 within 1.5 to 2.6 %;
  * D is 3 to 5 % below the reconstruction at C/10 and C/5, and 13 % below at C/2, probably because logging under-resolves 10 s pulses;
  * a gross mean fits worse at C/10 (2 to 6 % off).

**Against.** Synthetic and Drive shapes cannot be reconstructed. The published feature may slightly understate mean(I^2) for very short pulses, which biases q down.

**Confidence.** Near-certain for the definition and units. High that mean(I^2) = c^2 + c NCV holds to a few per cent.

## (3) Problems

### Must fix

**M1. A Law 1 claim that only holds for the registered slope.**
* README, "What it means": "under the closest published law it could not".
* This holds only for Law 1's registered slope. At h = 5 the hottest dynamic C/2 protocols put the module cells at 44.3 to 46.9 C. Law 1 at each cell's own temperature then gives a pair difference of -16.0 to -16.1 % (Periodic C) and -4.7 to -4.9 % (Periodic D). This is the same exploratory D18 form the README uses at 1C.
* At h = 10 and 25 it gives at most 0.15 %, so the note's h = 10 sentence is fine.
* Fix: "with its registered slope it could not; at h = 5 (exploratory) it could".
* Add the same qualifier, as should-fix, to:
  * README Result intro: "though not under Law 1";
  * README Result: "Under Law 1 every protocol stays below 1.3 %";
  * Secondary results: "under Law 1 none exceeds any floor";
  * .zenodo.json.

**M2. Post-hoc items not labelled post hoc.**
* Note, Method 6: "(8.9 % with its high-SOC term)".
* Note, Limits 2: "time constants of 0.2 to 0.5 h".
* README sensitivities bullet: "the air conductivity ... 0.406 or 0.378 K".
* README V1 bullet: "with all four terms ... (D23)".
* README D19 bullet: "(D24)" only.
* inputs/SOURCES.md: "+3.7 % and -3.4 %".

**M3. "Not redistributed" is false.**
* Where it appears: README Inputs ("not redistributed"), README Licence ("are not redistributed here"), and SOURCES.md item 9.
* results/dataset_protocol_heat.csv reproduces c_mean (= Average Current) and normalized_current_variance at full precision for all 92 cells.
* Also, the upstream repository has no licence file. The CC BY 4.0 statement rests on the SDR record, which I could not open.

**M4. README, Verification: "the implied divergences, Step 3 and the scores exactly".**
* The recomputation report says the flags, lowest rates and scores were reproduced.
* But 72 cell-derived implied values, 2 Step 3 values and the 3 P3 values differ by up to 1.2 %.

**M5. README, Result intro: "lifts its dynamic C/2 protocols above the 3.0 % mean floor".**
* Only six of eight do. Write "six of its eight".

### Should fix

* **S1. Disclose the CV omission.** Say it in D24, the README and the note, with its direction and the size from Section 1.
* **S2. Worst-of-8 comparison.** "The hottest also pass the 4.1 % worst of eight" and "three or four" use hottest minus coldest. D17's matching quantity, hottest minus mean, is exactly half the spread here, giving at most 3.24 % (h = 10). None of the protocols reaches 4.1 %. The README gives this comparison for the registered results (after Table 4) but not for the post-hoc ones.
* **S3. "Six of eight" is margin-dependent.**
  * Synthetic 2b and 2c are 0.03 to 0.05 points below the floor at h = 10.
  * At h = 10 the module mean is 38 to 45 C, where Law 2's local s is 3 to 8 % higher than at 35 C. Evaluating Law 2 at each cell's own temperature therefore gives all eight above 3.0 % (2b 3.10 %, 2c 3.07 %; the six span 3.87 to 6.94 %).
  * The CV correction moves the numbers the other way.
* **S4. Note, Limits 2: the time-constant argument covers only the charge phase.**
  * For the dynamic C/2 protocols the discharge window (1.6 to 2.2 h) is hotter.
  * Its steady spread is 2.06 to 2.25 K for Periodic C and 1.62 to 1.78 K for Periodic D. That is above Law 2's 1.68 K.
  * "Still below every floor" is true only for the constant-current D19 case.
  * Say that the cycle average is used because the lever is linear.
* **S5. Wording on Law 1's sign above 35 C.**
  * Note headline 1: "above 35 C a hotter cell lives slightly shorter ... of unknown sign". "Slightly" holds only from 35 to 45 C (126.0 to 124.5 cycles); above 45 C the drop is steep.
  * The sign is unknown because the 0.4C and 0.6C rows have opposite 35 to 45 C slopes (+0.032 and -0.029 per K); say that.
  * README "which is where every module cell sits" is true only at the dataset's rates with constant-current heat (35.0 to 38.8 C). It fails at 1C (40.1 to 49.5 C) and for the post-hoc protocols at h = 5.
* **S6. NCV evidence.** The README states the normalisation as fact, while SOURCES.md says it was inferred. Cite Supplementary Table 3 once the author has read it, and replace the 4 % test with the 5e-6 and 2e-6 identities.
* **S7. The 6.5 % and 3.3x come from a 0.6C protocol.** Both come from Periodic C "C/2", whose cells average 0.60C and 0.61C. The authors set them apart and drop cell_045 from their Fig. 3 analysis. A dataset co-author will notice.
* **S8. Filed checker README edited without marking.**
  * "Who had not seen the study being built" is not in the checker's REPORT.md.
  * The author added it to the filed copy of the checker's README, and removed its run-time estimates, without marking either edit. Mark the copy as adapted.
* **S9. Note, P5: "a factor of 5.2".** The computed values are 5.245 to 5.258, so three of the four rates round to 5.3. The README's "5.24 to 5.26" is right.
* **S10. Two more inherited values.** The time constants use the can's "about 9 g" and 500 J/(kg K), from a code comment in src/twocell/lumped_two_cell.py. List them with the three values that have no primary table.
* **S11. O'Regan wording.** "Nor a readable copy of O'Regan et al. 2022 reports" implies a copy was read. D8 says it could not be read.
* **S12. .zenodo.json.** It says "Two independent checks are included", but only the recomputation is filed. Its post-hoc range also lacks "(h = 10)".

### Fine (verified against the results and verification files)
* Every other number in the note and in the listed README sections, including:
  * 1.68 / 0.79 / 0.58 K;
  * 102, 126 and 124 cycles;
  * 4.8 to 5.4 K and 1.4 to 1.5 K;
  * every cell above 55 C at 2C (minimum 55.18 C);
  * mean temperatures of 41 to 49 C and 58 to 84 C;
  * cells at 40 to 50 C at 1C;
  * Law 1 at 1C: 0.2 % / 14.8 % / 19.5 % (note: "0.2 % to 20 %");
  * 0.05 % at C/2;
  * Law 1 by charge rate: +0.036 / +0.007 / -0.027 per K, 2.4 / 12.8 / -3.3 K;
  * the Monte Carlo interval and the 0.0009 per K standard error;
  * 57 to 78 %, and 0.051 against 0.052;
  * cold plate 0.016 to 0.036 K, 0.29 K and 1.16 K;
  * all sensitivities;
  * V1, V2 and V3 figures;
  * the scorecard apart from S9;
  * 1,336 quantities, 3e-9, 4e-10, 1.2 % and 2.3 %.
* The sign statements match comparison.json: Law 1's signed pair difference is negative in every case and blank at 2C.
* Pre-registration facts checked through the GitHub API:
  * repository created at 08:01:27;
  * main created at 08:01:59 with 7791192 on top of 24ab8d3;
  * tag and pre-release point to 24ab8d3.
* The floor's values and SHA-256 match the audit's evidence file.
* Note PDF: one page (pdfinfo), and its content matches the .md.
* No em or en dashes in any tracked file, either PDF or any commit message (only U+2212 minus signs).
* No schedule language; no letter or pilot label for this study.
* The result is not called a measurement of the thermal share.

## (4) Could not check
* The audit's own findings ("found none wrong or mis-rounded", the Crossref and DataCite checks): its report is not filed.
* Run times beyond commit times. They are consistent with runs after the push.
* The brief's phrase "if it is the case": the brief is not in the repository.
* Verbatim Supplementary text: I saw it only through a summarising tool.
* The SDR dataset (blocked):
  * the raw shapes of the Synthetic and Drive profiles;
  * its licence page;
  * whether per-cell temperatures were logged. The paper and SI report none, but the raw files might; a dataset co-author can answer.
* The dataset cell's own CV duration and any rests: I used the proxy cell instead.
* The Methods page number cited (p. 177).
