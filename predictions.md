# Predictions

The thermal lever against the duplicate-cell floor: how large a temperature spread a lifetime divergence needs, and how large a module offers.

Written and committed before any script has computed a lever at 35 C, a module spread at a registered rate and cooling, or any field of the resolved cell. The order of work after this commit is fixed: V3, then Step 1 (the lever), Step 2(a) (the lumped module), Step 2(b) (the resolved cell), the comparison and the scoring.

## What was known when these were written

- The inputs and every decision in `DECISIONS.md` (D1 to D20).
- V1 and V2 (`verification/v1_v2.json`): both pass. V1 for Law 1 confirms the transcription of the printed table; V1 for Law 2 gives 8.6 % capacity loss after 2800 full cycles at 1C and 25 C (the paper: less than 9 %) and 10.4 % at 45 C (the paper: 12 %), inside the paper's 2 point model error. V2 closes the network's energy balance to 1e-12 and holds eight nodes within 0.0004 K with the strips at 1e3 W/K.
- The law parameters, because copying them is Step 1 of the brief. For Law 1 the parameters are printed cycle counts, so its lever follows from them by a secant; for Law 2 it follows in closed form. The central estimates below are hand calculations from the copied values, which is the brief's reasoning without data, not outputs of the scripts.
- For the module, the central estimates come from a hand reduction of the eight-node network to its two node types (the four corner cells and the four inner cells are alike by symmetry), with radiation linearised. The script solves the full nonlinear network, so agreement there tests the implementation more than the physics. This is said here so that the scorecard is not read as more than it is.

## The five predictions (from the brief), with central estimates

P1. Delta T_floor for f = 8.7 % at 35 C lies between 1 and 3 K for every law used.
Central estimates. Law 1 (NCA+LNO / Si-graphite 21700, charge at C/2): ln L at C/2 is the mean of the 0.4C and 0.6C rows, so L(25 C) = sqrt(141 x 74) = 102 cycles and L(45 C) = sqrt(123 x 126) = 124 cycles to 90 % SoH; s = -(ln 124.5 - ln 102.1) / 20 K = -0.0099 per K; Delta T_floor(8.7 %) = -8.8 K: hotter cells live longer, and 8.8 K of difference is needed. Law 2 (LiFePO4 / graphite, full model, charge at C/2): at 35 C and a C/2 discharge the three square-root terms are about 1.06e-3 (calendar), 5.5e-4 (high-temperature cycling) and 9.0e-5 (low-temperature cycling) per square root of a cycle, and the high-SOC term is negligible; with activation energies of +20.6, +32.7 and -55.5 kJ/mol the effective energy is 20.5 kJ/mol, and because the terms are square roots s = 2 x 20.5 kJ/mol / (R T^2) = 0.052 per K, so Delta T_floor(8.7 %) = 1.7 K. At the slower discharge rates the calendar share grows, but the effective energy stays near 20.5 kJ/mol and s near 0.052 per K. Expected: P1 fails, because Law 1 needs 8.8 K; it holds for Law 2.
For the record: Delta T_floor(3.0 %) is 0.58 K and Delta T_floor(4.1 %) 0.79 K for Law 2, and -3.0 K and -4.1 K for Law 1.

P2. At C/2 average rate with h = 10 W/(m2 K), the eight-cell lumped spread Delta T_module is below 1 K.
Central estimate: 0.39 K, the inner cells hotter than the corner cells, all about 2.4 to 2.8 K above the air. Hand reduction: with the corner and inner loss conductances G_c and G_i and the pair coupling G_p, the corner rise theta_c satisfies 2q = theta_c [G_i (1 + d) + G_c] and the spread is d theta_c, with d = (G_c - G_i) / (G_i + 2 G_p). At h = 10: G_c = 0.0525 W/K, G_i = 0.0377 W/K (convection plus linearised radiation to the room) and G_p = 0.0269 W/K (strip 6.2, gap 17.4 and radiation 3.2 mW/K), with q = 0.117 W per cell. The same reduction gives 0.37 K at h = 5 and 0.35 K at h = 25. Expected: P2 holds.

P3. At C/2 the within-cell spread Delta T_cell of the resolved cell is below 0.5 K.
Central estimate: 0.12 K at h = 10 (radial conduction q''' R^2 / (4 k_r) = 0.10 K with q''' = 4.77 kW/m3 and k_r = 1.3 W/mK, plus about 0.02 K between the mid-height surface and the rim), between 0.12 and 0.15 K across h; hottest minus volume mean about 0.05 K. P3 is scored at h = 10 and must also hold at h = 5 and 25. Expected: P3 holds.

P4. The lowest rate at which Delta T_module exceeds Delta T_floor(8.7 %) with h = 10 W/(m2 K) is 1C or higher, and with h = 5 W/(m2 K) it is one step lower.
Central estimates. Law 2: at h = 10 the spread at 1C is about 1.54 K against 1.68 K, a near miss, so the lowest rate is 2C; at h = 5 the spread at 1C is about 1.44 K, so the lowest rate is again 2C. Law 1: no rate up to 2C reaches 8.8 K (2C gives about 5 to 6 K). Expected: the first clause holds for both laws; the second clause fails for both.
Scoring rule, fixed now: scored per law on the rate grid C/16, C/10, C/5, C/2, 1C, 2C, comparing the spread with the magnitude of Delta T_floor; "no rate up to 2C" counts as "1C or higher"; the second clause holds only if the rate at h = 5 is exactly one grid step below the rate at h = 10, and fails if either is "no rate up to 2C". P4 holds only if both clauses hold for both laws.

P5. The two laws' Delta T_floor(8.7 %) differ by less than a factor of 2.
Central estimate: 8.8 K against 1.7 K in magnitude, a factor of about 5. Expected: P5 fails.

## Scoring

A prediction holds when the computed value satisfies the brief's statement. The central estimates are recorded with their errors but are not scored. Values are taken at 35 C, at the registered resistance (18.71 mOhm) and with the registered law settings (Law 1 at a C/2 charge; Law 2 with the Table IV values).

## Secondary results, stated now and not scored

- Implied lifetime divergence s x Delta T at C/2 and h = 10: about 2.0 % cell to cell for Law 2 (0.052 x 0.39 K), below the 3.0 % mean floor; about 0.4 % for Law 1. Within the cell, about 0.6 % for Law 2.
- At the dataset's own rates (C/16 to C/2) and every h, the implied cell-to-cell divergence stays below the 3.0 % mean floor for both laws; the largest is about 2 %, Law 2 at C/2.
- At the cycle-averaged heat of the dataset's protocol, which includes the C/2 charge (D19), the heat at C/16 is eight times the brief's, but the implied divergence still stays below the mean floor at every rate.
- The cold-plate arrangement: the cell-to-cell spread falls to about 0.015 K at C/2, while the within-cell spread of a base-cooled cell rises to about 0.3 to 0.45 K, because most of the heat now leaves through the base along the jelly roll's axis.
