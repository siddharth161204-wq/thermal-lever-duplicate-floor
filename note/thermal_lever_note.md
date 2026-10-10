# The thermal lever against the duplicate-cell floor: how large a temperature spread a lifetime divergence needs, and how large a module offers

Siddharth Satte, Independent researcher, Pune, India, 10 October 2026. Predictions commit 24ab8d3, public at 08:02 UTC before any lever or spread was computed. [github.com/siddharth161204-wq/thermal-lever-duplicate-floor](https://github.com/siddharth161204-wq/thermal-lever-duplicate-floor){.nobr}.

<div class="headline">

1. To give two cells a lifetime difference as large as the 95th percentile of the duplicate-cell floor (8.7 %, doi.org/10.5281/zenodo.23205413), Law 2, an Arrhenius-type LiFePO4 model, needs 1.68 K of steady temperature difference at 35 C (0.58 K for the 3.0 % mean). Law 1, the closest published match to the dataset's cells (an NCA+LNO / Si-graphite 21700 charged at C/2), needs 8.8 K with the opposite sign: at 35 C its cycle life is near its optimum (102, 126 and 124 cycles to 90 % SoH at 25, 35 and 45 C), so a hotter cell lives slightly longer.
2. An eight-cell LG M50 module in 35 C air first offers 1.68 K at 2C, for every convection coefficient from 5 to 25 W/(m^2^ K) (4.8 to 5.4 K, with the cells 20 to 52 K above the air); 1C gives 1.4 to 1.5 K, and Law 1's 8.8 K is not reached up to 2C. At the dataset's own rates, C/16 to C/2, the module offers at most 0.39 K and a resolved cell 0.09 K inside itself: a lifetime divergence of at most 2.0 % under Law 2 and 0.4 % under Law 1, below the 3.0 % mean floor.

</div>

![](../results/thermal_lever_figure.png)

Figure 1. (a) The temperature difference each law needs against the floor f; floors dotted, the 95th percentile's bootstrap interval grey. (b) Hottest minus coldest cell of the lumped module (solid) and hottest minus coldest point of one resolved cell (dashed), with each law's ΔT~floor~ at 3.0 and 8.7 %. (c) The implied divergence s × ΔT.

## Method

1. Floor: the published re-analysis of the 45 duplicate pairs of Geslin et al. (Nature Energy 2025), all at 35 C, at 90 % SOH: pair lifetime difference 2.99 % on average and 8.72 % at the 95th percentile; worst of eight cells 4.10 % below their mean.
2. Laws: Law 1 is the printed cycles to 90 % SoH of Kucinskis et al. (J. Power Sources 2022) at 15 to 55 C and charge rates 0.2C to 0.6C (0.5C discharge), taken at a C/2 charge; Law 2 is the full model of Schimpe et al. (J. Electrochem. Soc. 2018, Table IV) at a C/2 charge and each discharge rate; its calendar term, which the paper adds for a cycling cell, carries 57 to 78 % of the loss at 35 C (the cycling terms alone give s = 0.051 against 0.052 per K). s = −d ln L/dT at 35 C; ΔT~floor~ = f/s.
3. Proxy cell: an LG M50 (21.1 by 70.1 mm), a different 21700 from the dataset's that sets magnitudes, not its values; heat (c&nbsp;×&nbsp;5&nbsp;A)^2^ × 18.71 mOhm, the resistance that reproduces the C/2 cycling heat of PyBaMM's Doyle-Fuller-Newman model (ORegan2022 set) at 35 C.
4. Module: two rows of four at 23.1 mm pitch; between neighbours a 6.2 mW/K nickel strip, 17.4 mW/K across the 2 mm air gap and grey radiation (view factor 0.16); convection and radiation (emissivity 0.8) to 35 C air.
5. Resolved cell: the two-cell study's scikit-fem cell unchanged (roll 1.3 W/mK across and 25 W/mK along the layers, 0.25 mm can), uniform heat, no neighbour, same losses, steady state.
6. Checks: Law 2 gives 8.6 and 10.4 % loss after 2800 cycles at 1C, 25 and 45 C (paper: under 9 % and 12 %); energy closes to 1e-12; against a lumped cell the surface mean agrees and the volume mean is 1.8 to 5.2 % of the rise higher, failing the brief's 1 % test as foreseen (a lumped cell has no internal gradient); a mesh refinement moves the spread 0.8 %; the radial profile is within 0.13 % of the analytic one.

## Limits

1. Neither law was fitted on the dataset's cells; Law 1 is a 10 K table at one discharge rate, so its lever is local: at 1C (h = 10) the cells sit at 44.5 to 46.1 C, on its steep 45 to 55 C segment, where their 1.5 K spread gives a 15 % pair difference (exploratory).
2. Steady heat at the average current, from a derived resistance without contact resistance (a 10 s or 30 s pulse value gives 0.32 or 0.43 K at C/2). The dataset also charges at C/2, which raises the cycle-averaged heat up to eight-fold at C/16, where the spread is then 0.05 K, still below every floor.
3. One lumped geometry, one strip per pair, no diagonal coupling; the radiation form matters (the received form gives 0.65 K at C/2). A cold plate under the bases cuts the C/2 module spread to 0.02 to 0.04 K but raises the cell's own to 0.29 K.

## Prediction scorecard

| | Prediction | Result | Computed |
|--|-------------------------------------------------|------|---------------------------------------------|
| P1 | ΔT~floor~(8.7 %) lies between 1 and 3 K for every law | Failed | Law 2: 1.68 K; Law 1: −8.81 K (central 1.7 and −8.8 K) |
| P2 | At C/2 and h = 10 the module spread is below 1 K | Held | 0.39 K (central 0.39 K); 0.37 and 0.35 K at h = 5 and 25 |
| P3 | At C/2 the resolved cell's spread is below 0.5 K | Held | 0.092 K at every h (central 0.12 K, 30 % high) |
| P4 | The module first exceeds ΔT~floor~(8.7 %) at 1C or higher with h = 10, one rate lower with h = 5 | Failed | First clause held (Law 2 at 2C, Law 1 not up to 2C); h = 5 gives the same rates |
| P5 | The two laws' ΔT~floor~(8.7 %) differ by less than a factor of 2 | Failed | A factor of 5.2, with opposite signs (central 5.2) |

## What it means for the thesis

At the dataset's own conditions (35 C, C/16 to C/2, convective cooling) a module's thermal lever stays inside the duplicate-cell floor under both laws; only near 2C does it clear the 95th percentile, and under Law 2 alone. That does not show that the dataset's own divergence is not thermal (its cells sat in separate fixtures of unreported temperature, and 0.6 K would suffice under Law 2), but it does show that the thesis first needs the dataset cell's own cycle life against temperature near 35 C at its charge rate, where the two laws disagree in sign.
