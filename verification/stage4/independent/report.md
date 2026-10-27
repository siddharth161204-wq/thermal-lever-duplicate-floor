# Stage 4 independent check

(The checker's report, saved by the main session from the checker's hand-back text, since the checker could not write report files itself; the checker's own files, record and comparison are in this folder.)

Own code from D56 to D63, the brief, the predictions and the Stage 1 to 3 inputs, built on `c1e58a1`. Nothing from `src/stage4`, `results/stage4_*`, `verification/stage4` or `tests/stage4_synthetic.py` was opened, and no later commit was logged, before the freeze: `fd5d1e7`, 2026-10-27T10:33:27Z, pushed 10:33:36Z to `stage4-independent`. Post-freeze files are named `compare*` or `postfreeze*`.

## Method

- (a) `lgm50_resistance.py` imported unchanged, `T_K` set per ambient, every run recorded. A separate constant-current DFN integration on a 1 s grid, run raw and with a spike screen.
- (b) My own network from D5 to D12, solved by nested `brentq` on the symmetric reduction. Reference: an 8-node Newton at 50 digits. Cross-check: `module_steady`.
- (c) `SteadyCell(1.0)` as D59 states, with the analytic q/(4π k_r L) bound.
- (d) Levers checked against D56; |s| × spread with intervals, ln offsets, numpy median and 95th percentile, map. A plain-Python second implementation of the two V3 cases.
- (e) Noncentral-t power by quadrature over the chi distribution (matches mpmath to 1e-16 and the textbook sizes 34, 15 and 199). Bisection simulation with 10,000 experiments per n over all 144 cases and both shapes; 200,000 for P7's two cases.
- (f) D63 scoring.

## Agreement

Against the registered results in `6f658f1`, no deterministic number, count, flag, label, map entry or verdict differs, across 81 quantity classes (`out_compare/compare.json`).

| Quantity | Largest difference |
|--|--|
| Resistance: every run, R at 35, 25 and 10 C, ratios, the same two failed runs | identical to the last bit |
| Module cell temperatures, 36 cases (own solver) | 1.7e-13 K |
| Module spreads | 1.1e-13 K (relative 8.8e-13) |
| Module via `module_steady` (V1 included) | identical to the last bit |
| Resolved cell, within-cell spread | 2.8e-13 K |
| Implied divergences, medians, 95th percentiles | 1.4e-14 absolute, 8.8e-13 relative |
| Median levers, offers and detectable, map | identical |
| Experiment: floor shapes, calculated n (all 144 rows, both shapes), figure curves | identical |
| Power at n | 1.7e-12 |

Only the Monte Carlo results differ, because the random seeds differ: the simulated n is off by up to 43 pairs at n of about 5,700 (0.8 %), and by one pair at n of 3 to 9. V4 passes in 84 of 144 cases here and 82 in the pipeline; 10 cases fall on opposite sides of the 5 % line.

## Verdicts (the same as the pipeline's)

- P1 held: at most 1.98 % against 2.99 %.
- P2 did not hold: 5 of 13 warm levers exceed 8.72 % (up to 16.5 %).
- P3 held: ratios 2.65 to 2.72.
- P4 held: 13.4 to 14.7 %.
- P5 did not hold: every 10 C, 1C case exceeds 7.33 %.
- P6 held: ratio 0.229.
- P7 did not hold, and is scored: cold chamber n = 3; warm chamber n = 25, and 23 when simulated with 200,000 experiments per n.

## Choices the rules left

| Point | Reading (both sides, unless stated) | Effect |
|--|--|--|
| "Fails to solve" | a raised exception | the same two failures |
| Resistance taken at the ambient | as registered | decides P2, P4, P5 (flag 2) |
| Two-sided power | quadrature (checker); lower tail set to zero where scipy gives nan (D64) | same n |
| V4 seeds | different | 10 boundary differences |
| Normal check sd(e) | registered sd(e) (pipeline; checker after the freeze) | none on verdicts |
| Band membership | set points (checker after the freeze) | 17 of 36 outside |
| P2 lever test | point estimates | none |
| P7's central cold lever | cold band median, 3.5800 | the C1 value 3.6139 also gives n = 3 |
| Figure grid | 0.5 K (checker), 0.1 K (pipeline) | agree where they overlap |

## Flags

1. The registered 10 C resistance is right only by the luck of the time grid. At default solver tolerances the C/2 charge at 10 C carries a sub-second heat spike (13.8 kW at t = 3954.70 s, about 3.4 kJ against about 1.1 kJ for the rest of the run) while the voltage is smooth, so it is a solver artifact. The procedure's 10 s output grid steps over it; on a 1 s grid the checker's separate code gives 55.21 against 47.95 mOhm (+15 %). At tight tolerance (rtol 1e-8) the spike vanishes and the registered charge value is confirmed to 0.05 %. D65's "no isolated spike" therefore describes the samples, not the solution, and V2's 0.02 % agreement rests on its tight tolerances.
2. Taking the resistance at the ambient decides three verdicts, not only those D66 names. With the heat at each cell's own temperature (the checker's sensitivity, with the same procedure's R(T)): P2 would hold (3.10 %, and the largest lever gives 8.23 % against 8.72 %), P4 would fail (5.54 % at h = 5), and P5 would hold (5.54 to 7.27 % against 7.33 %). D66 does not cover P2.
3. Levers are applied far outside the temperatures they were measured at: 17 of 36 cases have cells outside their band's set points, and P4 rests on 10 C, 1C cells at 24 to 50 C under levers measured from 10 to 25 C. At 10 C and 2C the 95th percentiles reach 75 to 98 %, where |s| × spread no longer approximates a pair difference.
4. V4's 5 % rule is coarse at small n and depends on the seed. P7's warm case simulates to n = 23 (power 0.8024 at 23, 0.7897 at 22), 8 % below 25, so V4 fails there; with 10,000 experiments the checker's run passed it at 24. P7's verdict is unaffected either way, since all values are below 30.
5. P7 turns on the registered median-of-|s| reading: under the signed median the warm chamber needs 221 pairs and P7 would hold. The note's 25 pairs is an upper estimate (23 simulated).
6. P5's 35 C clause holds by only 0.29 points (7.046 against 7.332 %).
7. The 10 C C/2 charge reaches 4.2 V after 2.75 Ah (55 % of 5 Ah), so R(10 C) averages a full-window discharge with a half-window charge; the discharge alone would be 2.9 % higher. No verdict moves.
8. Wording: "a 10 C module" should be "a module in 10 C air", since its 1C cells run at 24 to 50 C; "offers a share from C/2" holds only under the upper heat bracket; the note's "KIT and Sandia, 25 to 40 C" for the median lever is KIT's alone (the Sandia lever is 25 to 35 C).
9. The unscored 2C within-cell ratio is 0.306, close to one third; an analytic infinite-cylinder bound gives 0.339.
10. Where the checker agrees with the pipeline: the C/16 charge values at 35 and 25 C are solver spikes (D65), and scipy returns nan at n = 3 in 12 cases (D64), where the checker's quadrature gives the same n.

## Files

Frozen scripts and `out/` with `manifest.json`; `compare.py`, `postfreeze_*.py`, `out_compare/compare.json`, `record.json`.
