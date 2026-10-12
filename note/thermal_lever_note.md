# The thermal lever against the duplicate-cell floor: how large a temperature spread a lifetime divergence needs, and how large a module offers

Siddharth Satte, Independent researcher, Pune, India, 10 October 2026. doi.org/10.5281/zenodo.23309131. Predictions commit 24ab8d3, public at 08:02 UTC before any lever or spread was computed. [github.com/siddharth161204-wq/thermal-lever-duplicate-floor](https://github.com/siddharth161204-wq/thermal-lever-duplicate-floor){.nobr}.

<div class="headline">

1. For a lifetime difference as large as the duplicate-cell floor's 95th percentile (8.7 %, doi.org/10.5281/zenodo.23205413), two cells at 35 C need a steady 1.68 K under Law 2, an Arrhenius-type LiFePO4 model (0.58 K for the 3.0 % mean). Law 1, the closest published match found to the dataset's cells (NCA+LNO / Si-graphite), peaks near 35 C (102, 126 and 124 cycles to 90 % SoH at 25, 35 and 45 C, C/2 charge): its registered centred slope needs 8.8 K of the opposite sign, yet from 35 to 45 C, where its 0.4C and 0.6C rows slope oppositely, a hotter cell lives slightly shorter: its lever there is weak and of unknown sign.
2. An eight-cell LG M50 module in 35 C air first offers Law 2's 1.68 K at 2C for convection coefficients of 5 to 25 W/(m^2^ K) (4.8 to 5.4 K), with every cell above 55 C, outside both laws' ranges; 1C gives 1.4 to 1.5 K, and no rate up to 2C gives 8.8 K.

</div>

![](../results/thermal_lever_figure.png)

Figure 1. (a) |ΔT~floor~| of each law against the floor f; floors dotted, the 95th percentile's bootstrap interval grey. (b) Hottest minus coldest cell of the lumped module (solid) and point of one resolved cell (dashed), with each law's ΔT~floor~ at 3.0 and 8.7 %. (c) The implied divergence |s| × ΔT. Open circles, post hoc: the module at h = 10 with each dataset cell's own cycle-averaged heat.

## Method

1. Floor: the published re-analysis of the 45 duplicate pairs of Geslin et al. (Nature Energy 2025), all at 35 C, at 90 % SOH: pair lifetime difference 2.99 % on average and 8.72 % at the 95th percentile; worst of eight cells 4.10 % below their mean.
2. Laws: Law 1 is the printed cycles to 90 % SoH of Kucinskis et al. (J. Power Sources 2022) at 15 to 55 C and charge rates 0.2C to 0.6C (0.5C discharge), interpolated to a C/2 charge; Law 2 is the full model of Schimpe et al. (J. Electrochem. Soc. 2018, Table IV) at a C/2 charge and each discharge rate; its calendar term, which the paper adds for a cycling cell, carries 57 to 78 % of the loss at 35 C (cycling terms alone: s = 0.051 against 0.052 per K). s = −d ln L/dT at 35 C; ΔT~floor~ = f/s.
3. Proxy cell: LG M50 (21.1 by 70.1 mm), a 21700 but not the dataset's cell, so it sets magnitudes, not values; heat (c&nbsp;×&nbsp;5&nbsp;A)^2^ × 18.71 mOhm, the resistance reproducing the C/2 cycling heat of PyBaMM's Doyle-Fuller-Newman model (ORegan2022) at 35 C.
4. Module: two rows of four at 23.1 mm pitch; between neighbours a 6.2 mW/K nickel strip, 17.4 mW/K across the 2 mm air gap and grey radiation (view factor 0.16); convection and radiation (emissivity 0.8) to 35 C air.
5. Resolved cell: the two-cell study's scikit-fem cell unchanged (roll 1.3 W/mK across and 25 W/mK along the layers, 0.25 mm can), uniform heat, no neighbour, same losses, steady state.
6. Checks: Law 2 gives 8.6 % (8.9 % with its high-SOC term, post hoc) and 10.4 % loss after 2800 cycles at 1C, 25 and 45 C (paper: under 9 % and 12 %); energy closes to 1e-12; the resolved cell's surface mean matches a lumped cell, and its volume mean lies 1.8 to 5.2 % of the rise above, failing the brief's 1 % test as foreseen; a mesh refinement moves the spread 0.8 %. An independent recomputation reproduces the levers and the module to round-off and puts the cell's spread 1.2 % higher.

## Limits

1. Neither law was fitted on the dataset's cells. Law 1 rests on a 10 K table: at each cell's own temperature, all above 35 C, it gives the hotter cell the shorter life, against its registered slope, and at 1C, with the cells at 40 to 50 C, pair differences from 0.2 % (h = 25) to 20 % (h = 5) (exploratory).
2. The proxy cell and a derived resistance without contact resistance set the heat (a 10 s or 30 s pulse value gives 0.32 or 0.43 K at C/2), held steady at its cycle average, which the linear lever allows to first order; the cells follow each phase (time constants 0.2 to 0.5 h, post hoc), so the spread nears 0.35 to 0.39 K in each C/2 charge and about 2 K in the hottest dynamic C/2 discharges.
3. One lumped geometry, one strip per pair, no diagonal coupling; the radiation form matters (the received form gives 0.65 K at C/2). A cold plate under the bases cuts the C/2 module spread to 0.02 to 0.04 K but raises the cell's own to 0.29 K.

## Prediction scorecard

| | Prediction | Result | Computed |
|--|-------------------------------------------------|------|---------------------------------------------|
| P1 | ΔT~floor~(8.7 %) lies between 1 and 3 K for every law | Failed* | Law 2: 1.68 K; Law 1: −8.81 K (central 1.7 and −8.8 K) |
| P2 | At C/2 and h = 10 the module spread is below 1 K | Held | 0.39 K (central 0.39 K); 0.37 and 0.35 K at h = 5 and 25 |
| P3 | At C/2 the resolved cell's spread is below 0.5 K | Held | 0.092 K at every h (central 0.12 K) |
| P4 | The module first exceeds ΔT~floor~(8.7 %) at 1C or higher with h = 10, one rate lower with h = 5 | Failed* | First clause held (Law 2 at 2C, Law 1 not up to 2C); h = 5 gives the same rates |
| P5 | The two laws' ΔT~floor~(8.7 %) differ by less than a factor of 2 | Failed* | A factor of 5.2 to 5.3, of opposite signs (central about 5) |

<p class="credit">&#42; A failure that predictions.md expected before any lever was computed.</p>

## What it means for the thesis

At the dataset's rates with constant-current heat a module stays inside the floor under both laws (at most 0.39 K, 2.0 % under Law 2), but with the current variance the dataset publishes for its dynamic discharges and its CC-CV charge, six of its eight dynamic C/2 protocols reach 3.5 to 6.1 % under Law 2 (post hoc, h = 10), above the 3.0 % mean, while Law 1's registered slope keeps every protocol below 1.2 %. Whether a module can show a thermal share above the floor therefore turns on the dataset cell's own cycle life against temperature near 35 C, which the closest published law cannot fix even in sign: that is the first measurement the thesis needs, and the dataset's own divergence, from cells in separate fixtures in a 35 C chamber whose own temperatures are not reported, is neither shown nor ruled out as thermal.
