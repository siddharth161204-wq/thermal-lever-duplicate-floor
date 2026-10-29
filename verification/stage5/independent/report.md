# Stage 5 independent check

(The checker's two reports, saved by the main session from the checker's hand-back text; the checker's own files, record and comparison are in this folder.)

Own code from D71 to D78, the brief, the predictions and the archived inputs, built on `da613a1`. Nothing under `src/stage5/`, `results/stage5_*`, `verification/stage5/` outside this folder or `tests/stage5_synthetic.py` was opened or listed, and no later commit was logged, shown or fetched, before the freeze: `ea3a5b4`, committed 2026-10-29T04:40:07Z and pushed 04:40:14Z to `stage5-independent`. Every input has the SHA-256 it has at v1.3.0, and a clean rerun regenerates every frozen output bit for bit. The comparison (`4e41831`, pushed 04:58:38Z) read the pipeline with `git show` at `101ca54`, `157827f`, `67685ae` and `ee4a403`. Post-freeze files are named `compare*` or `postfreeze*`; the frozen `record.json` is untouched.

## Frozen results

- Pairs: KIT 142 in the file, none censored, the 36 at 0 C not scored, 106 scored in 36 sets (warm 70 in 24, cold 36 in 12; the warm interval is the narrower in all twelve protocols); 104 of the 106 are ties, the two untied (P049, P052) in the warm band. Sandia 86 pairs, 60 with measure (a), the four LFP pairs at 35 C without a determinable lever, 56 scored in 21 groups (warm 36 in 11, cold 20 in 10), all ties; measure (b) orients 21 of them the other way. Pooled 162 in 57.
- Part B (its own seeds): KIT beta -0.0862 (-1.9135 to +1.4274), Sandia +1.5496 (-0.7070 to +2.4816), pooled +1.0956 (-0.5574 to +1.9466); s_w KIT warm -5.191 (-16.01 to +8.20), cold -0.793 (-9.32 to +4.05) % per K; signs 45 of 106 (p = 0.1448), 32 of 56, 77 of 162; residual sd KIT 0.1309, Sandia 0.1581, pooled 0.1417, with the same floor verdicts; sd(e) rebuilt from the floors exactly.
- Part C: detectable |s| KIT 6.5599 (KIT floor) and 2.8061 (Geslin) % per K, Sandia 2.5818 and 1.1044, pooled 2.3771 and 1.0168, by its own noncentral-t quadrature, which agrees with Stage 4's `power_t` to 3.1e-15; every C2 pair count matches.
- Verdicts: P1, P2, P4, P5, P6 and P8 hold, P3 and P7 fail; the same in 200 replicate bootstraps and at 1,000,000 draws.

## Agreement with the registered results (`out_compare/compare.json`)

| Class | Largest difference |
|--|--|
| Pair table, 228 rows: scored, band, tie, orientation by measures (a) and (b) | no mismatch |
| Pair table numbers (dT, d, p, lever and interval, lifetimes, resolution, measure (b)) | 6.9e-18 absolute, 2.2e-16 relative |
| beta (six values) | 2.2e-16 absolute |
| s_w (eight values) | 6.9e-18 absolute |
| Residual standard deviations | identical |
| Counts, sign counts, sign p-values and shares | identical |
| C1 calculated, 20 cases (with normal shape, F and shared cells) | 1.5e-14 absolute, 2.1e-13 relative (root-finding tolerances) |
| C2 calculated power; pairs needed | 2.8e-17; identical |
| P1 to P8, P6's scoring, B4's floor verdicts | identical |

Monte Carlo quantities differ by their seeds only: across 92 bootstrap ends, simulated detectable |s|, simulated powers and simulated pair counts, every registered value lies within 2.23 standard deviations of the mean of the checker's independent replicate runs (four beyond 2, against 4.2 expected). P7's width ratio is 0.947 registered, 0.954 in the checker's run and 0.953 with 1,000,000 draws. The one Monte Carlo verdict that differs is V4 in the unscored Sandia cold band under the KIT floor: the registered run fails it (5.15 % below), the checker's passes (4.70 %), and 15 of its 20 replicate runs pass.

## Flags (taken up in DECISIONS.md D79 to D81)

1. The KIT untied bracket holds two warm pairs and no cold pair; it carries no information.
2. A 10,000-draw interval end moves by its Monte Carlo spread (about 0.04 on KIT's beta); compare ends within several of those spreads, not to rounding.
3. The pooled warm within-set lever's registered interval (-7.14 to -0.10 % per K) excludes zero by chance: 1,000,000 draws give -7.11 to +0.03, and 65.5 % of 200 independent runs include zero.
4. B4's pooled interval clears the KIT sd(e) by 0.0014; for KIT the comparison is partly circular, since the KIT floor comes from cells 1 and 2 of these same sets.
5. The KIT warm band runs against the between-set lever: its per-band beta (not registered) is -7.75 (-17.81 to -2.00 with 1,000,000 draws; -11.02 to -5.83 leaving out any one set), the cold band's +0.32 (-0.78 to +2.71); with one pair per set the warm band's sign test gives p = 0.15.
6. P7 fails on the sum of p^2 (Sandia 0.077 against KIT 0.030), set by Sandia's larger swing differences, not on the ties, which both datasets share.
7. One set carries each fit: KIT P036 holds 52 % of KIT's sum of p^2 (beta without it -0.37; leaving out any one set, -0.43 to +0.17); the Sandia group NMC 25 C 0-100 % 0.5/3C holds 57 % of Sandia's and 41 % of the pooled weight (without it Sandia 0.68, pooled 0.32; P7's ratio 1.30); in that group cell c swings 5.99 K against 10.12 to 11.30 K and lives 107.4 EFC against 135.6 to 145.6.
8. Without the NCA 15 C 1C pair whose sensor did not follow its cell, Sandia's interval excludes zero in every replicate.
9. Sandia's 2C and 3C groups carry the 1C lever (D72), though Stage 3's 2C levers differ (NCA 25 to 35 C +1.30 against -2.04, NMC -7.14 against -2.44 % per K, and others); these 23 pairs in 13 groups hold 84 % of Sandia's weight; beta over the 1C groups alone is 0.04 (-4.40 to +3.92), over the 2C and 3C groups alone 1.83 (-0.47 to +3.81). The pooled analysis mixes KIT's mean-temperature differences with Sandia's swing differences.
10. Some V4 verdicts in the unscored Sandia bands depend on the seed; P6's case passes in all 20 replicates.
11. P6 holds against the KIT floor; against the Geslin floor the KIT detectable |s| is 2.81 % per K, below the cold median.
12. C1's calculation treats pairs as independent; the shared-cell factor F (1.5 to 1.72) is given beside it.
13. C2's simulation, which resamples dT, needs more pairs than the calculation at small n (Sandia cold under the Geslin floor 19 against 13).
14. The checker's simulated KIT cold count (198) is low by chance; its replicates give 203 to 207 (registered 207).
15. D77 names V3's Sandia pair twice; the parenthesis's pair is the one both computed.
16. The brief's P1 expectation cites 142 pairs at a median 0.6 K; 106 are scored, at a median pair dT of 0.26 K.
17. Resampling whole protocols instead of sets changes no verdict: with 1,000,000 draws KIT -1.66 to +1.25, Sandia -0.15 to +2.50, pooled -0.24 to +1.93 (P7's ratio 0.91).
18. scipy's noncentral t returned nan at 36 of Brent's trial points far from the root (nc 11.8 and 23.7); its quadrature filled them, and no reported value is affected.

## Text corrections (`out_compare/text_review.json`)

Of 96 quoted numbers, 95 round correctly from the pipeline. The exception: the KIT warm band's spread is 0.3447 K root mean square, so 0.34 K, not 0.35 K (D80 (5), the README, the note and the Zenodo text). Further points on wording: "the within-set noise is the floor" and "the residual spread is the KIT floor's" overstate a comparison that is partly circular and whose pooled interval excludes the floor ("consistent with the KIT floor"); the note's "show no thermal lever" should read "cannot show a thermal lever"; flags 5, 7 and 9 belong in D80's account of P3 and P7; V4's Sandia cold failure holds "in this run"; a Sandia group's F is 1 (two cells) or 2 (four cells), the values 1.13 to 1.85 being those of the bands and of all 56 pairs; the 6.6 % per K, 4,221 and 128 are against the KIT floor (Geslin: 2.81, 774 and 25); Figure 4's KIT curve is drawn at all 106 pairs' 0.44 K spread, where the warm median needs 2,640 pairs; s_w interval ends are best quoted to 0.1 % per K; the README's Stage 5 run lines lacked `-I`.
