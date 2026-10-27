# The thermal lever against the duplicate-cell floor

How large a temperature spread a lifetime divergence needs, and how large a module offers.

Author: Siddharth Satte, independent researcher, Pune, India (ORCID 0009-0001-3584-1008).

Archived at https://doi.org/10.5281/zenodo.23604443 (version v1.2.0: Stage 1, and Stages 2 and 3 below; all versions https://doi.org/10.5281/zenodo.23309130).

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

## Stage 2: the lever measured on public data with replicates (v1.1.0)

Stage 1 ended on one unmeasured number: the slope of cycle life against cell temperature near 35 C, which the closest published law cannot fix even in sign. Stage 2 measures it on the KIT ageing dataset of Luh and Blank (doi.org/10.35097/1969, version 2, CC BY 4.0; version 1 doi.org/10.35097/1947; data descriptor Scientific Data 11, 1004, 2024): LG INR18650HG2 cells, NMC with a SiOx-graphite anode, the anode family of the Geslin cells, cycled in twelve protocols (windows 0-100, 10-100 and 10-90 % by charge/discharge rates 0.33/0.33, 0.33/1.00, 1.00/1.00 and 1.67/1.00 C), each on three cells at 0, 10, 25 and 40 C in temperature-controlled oil pools, with capacity check-ups at 25 C and the cell temperature logged at every run. The lifetime L is the equivalent full cycles to 90 % of each cell's first check-up capacity; s = -d ln L / dT between group means, against the logged group temperatures (DECISIONS.md D25 to D35). The Sandia study planned first waits for its files, which batteryarchive.org sends only on request, and will run as Stage 3.

The eight predictions in [predictions_stage2.md](predictions_stage2.md) were committed in `f15ee119cc512a2a06bf28b39fd7705b6cbfd785` and pushed at 06:37 UTC on 21 October 2026, before any capacity, temperature or date of the cell data was read.

- The lever near 35 C is a property of the protocol. From 25 to 40 C (logged group temperatures 14.6 to 15.1 K apart) s runs from -3.45 % per K (0-100 %, 1.67/1.00 C) to +1.91 % per K (10-90 %, 0.33/0.33 C), each with a 95 % interval of about plus or minus 0.9 % per K (pooled within-set standard deviation of ln L 0.082). It is positive in five of the six protocols with a 0.33 C charge (their mean 0.95 % per K), negative as a point estimate in all six with a 1C or 1.67 C charge, and lower at 1.67 C than at 0.33 C charge in every window. From 10 to 25 C every protocol's lever is negative (-1.8 to -11.6 % per K). The setpoint temperature difference in place of the logged one moves s by at most 0.07 % per K. At the secondary 80 % threshold s moves by up to 5 % per K and changes sign in two of the nine protocols determinable there, so the lever also depends on how far the cells have aged.
- The replicate floor of this dataset is larger than the Geslin floor. The pre-specified pairs (cells 1 and 2 of each of 47 sets) differ by 7.3 % on average (median 4.1 %) and 20.5 % at the 95th percentile, and the worst of eight cells would fall 10.9 % below the mean; all 142 pairs give 7.9 % and 22.8 %, and the 80 % threshold 8.2 % and 22.1 %.
- The warmer replicate does not age faster. In 142 pairs within triplets, a median 0.33 K apart, the warmer cell reached 90 % first in 57 (40 %; one-sided sign test p = 0.99).
- Back to the module. At the largest measured magnitude, 3.45 % per K, Stage 1's registered offer of 0.391 K implies a 1.35 % lifetime divergence, below the 3.0 % mean floor at every protocol; the post-hoc dynamic offer of 1.17 K (D24) would give up to 4.0 %, above 3.0 % in two protocols and below 8.7 % in all.

| | Prediction (abridged) | Result | Computed |
|--|--|--|--|
| P1 | 0.33 C charge: s > 0 in every window, mean 1 to 6 % per K | Failed | 0-100 %, 0.33/1.00 C: -0.62 % per K; mean 0.95 % per K |
| P2 | s at 1.67 C charge below s at 0.33 C, every window | Held | -3.45 < -0.62; -1.91 < 0.67; -1.00 < 1.09 % per K |
| P3 | the twelve point estimates take both signs | Held | -3.45 to +1.91 % per K |
| P4 | floor mean 1.5 to 6.0 %, 95th percentile 4.4 to 17.4 % | Failed | 7.3 % and 20.5 % (47 pairs) |
| P5 | the warmer cell reaches 90 % first in 45 to 65 % of pairs | Failed | 40 % (57 of 142) |
| P6 | no protocol's s exceeds 7.7 % per K | Held | largest +1.91 % per K |
| P7 | the mean s of P1's protocols lies between -0.99 and 5.2 % per K | Held | 0.95 % per K |
| P8 | 10 to 25 C: s < 0 for at least half of the six 1C-and-faster protocols | Held | all six negative (-3.3 to -11.6 % per K) |

Checks (`verification/stage2/`). V1: the archive's MD5 equals the repository's; 228 main-experiment cells as the descriptor states, 48 cyclic sets of three, 144 EOC files, the protocol table from the configuration files. V2: a separate checker wrote its own code from D25 to D35 without seeing the pipeline and reproduced every cell's lifetimes, censoring and temperatures, both levers with their intervals and sensitivities, the floor, the worst of eight and the warmer share, across 583 compared values to 5e-13 at most, and reached the same eight verdicts (`verification/stage2/v2_independent_report.md`). V3: the EOC running capacity estimate crosses 90 % at a median 66 % fewer equivalent full cycles than the check-ups, because it is estimated at the operating temperature, so it does not serve as a sensitivity; every 25 C cell crossed 90 % by 27 March 2023, before the Peltier failure of 23 April 2023; the November 2022 pool excursion falls inside 63 of the 144 lifetime windows. V4: synthetic files with planted lifetimes, temperatures and levers: lifetimes and temperatures recovered to 1e-14, the planted levers inside their 95 % intervals in 65 of 72 cases over six seeds.

Limits. The 90 % crossing lies between check-ups three weeks apart that span a median 61 % of the lifetime (20 % at 80 %), so single lifetimes are coarse, though the floor at 80 % agrees. A logged cell temperature carries its sensor's offset, and cells in one pool differ by tenths of a kelvin. The HG2 is an 18650 and not the Geslin 21700; it shares the anode family. One cell (P044-2) stopped after its second check-up and is censored.

What it means for the thesis. Near 35 C a cell's lifetime changes by at most about 3.5 % per kelvin in this dataset, and with a sign set by the protocol, so the 0.39 K a module offers at the Geslin rates moves lifetimes by at most 1.4 %, well inside either floor; and among replicates in the same pool the warmer cell is not the one that ages first. The replicate spread itself, 2.5 times the Geslin floor here, is not ordered by the cells' logged temperatures.

The Stage 2 figure is `results/stage2_figure.png`; every number is in `results/stage2_*.csv` and `.json`.

## Stage 3: the lever in a second laboratory's data (v1.2.0)

Stage 2 measured the lever on one cell model in one laboratory. Stage 3 measures it on the Sandia National Laboratories cycling study of Preger et al. 2020 (Journal of The Electrochemical Society 167, 120532; doi.org/10.1149/1945-7111/abae37), whose files Battery Archive (www.batteryarchive.org) shared on request on 23 October 2026: 86 commercial 18650 cells, LFP (A123 APR18650M1A), NCA (Panasonic NCR18650B) and NMC (LG 18650HG2, the KIT cell), cycled at 15, 25 and 35 C with a 0.5C charge and discharge rates from 0.5C to 3C in windows of 0-100, 20-80 and 40-60 %, with repeated 0.5C capacity checks and the skin temperature logged. The lifetime L is the equivalent full cycles to 90 % of each cell's first check (95 % for LFP, whose 1C cells at 15 and 25 C never reached 90 %), and s = -d ln L / dT between the logged cycling temperatures of the 1C series at 15, 25 and 35 C (DECISIONS.md D38 to D50). The cells cycle without rests, so the self-heating measure is the within-cycle swing, the mean skin temperature over each discharge minus that over the charge before it, chosen by the author before any temperature was read (D43).

The ten predictions in [predictions_stage3.md](predictions_stage3.md) were committed in `bf4b99380c217aed7f875d06061a32b43271d622` and pushed at 19:34 UTC on 23 October 2026, before any lifetime or temperature value was computed. Structure checks after that push found that in all 18 cells at 35 C the cycler logs the last record of a cycle up to 0.42 s after the first record of the next, and that the next cycle often begins on the last sample of a discharge with the counters at zero: the registered ampere-hour rule would have counted under 3 % of the discharge of the LFP cells at 35 C. D50 puts the rows in cycle order, drops 4,313 rows of another record found in one file, and takes the counter difference sample by sample; it was pushed (`01e562b`, 20:11 UTC) before any lifetime or temperature was computed.

- The lever. On the 1C series the lifetimes of both layered-oxide chemistries rise from 25 to 35 C: s = -2.44 % per K for NMC (95 % interval -4.33 to -0.44), inside the band the same cell spans at KIT with a 1C discharge (-0.62 to -2.69 % per K from 25 to 40 C), and -2.04 % per K for NCA (-2.95 to -1.20). LFP's 25 to 35 C lever is not determinable as registered, because one of its 35 C cells logged no cell temperature before cycle 259, after its lifetime, which leaves the group temperature of D43 undefined; its lifetimes fall from a geometric mean of 1,873 to 470 cycles, and with the group temperature over the three cells that have one s = +12.9 % per K, about three times Law 2 (post hoc, D52). From 15 to 25 C s is -5.62 % per K for NMC, -1.32 for NCA (interval -5.22 to +0.74) and +3.37 for LFP: cold shortens NMC and NCA life and lengthens LFP's. The secondary 2C series gives +11.0 (LFP), +1.3 (NCA) and -7.1 % per K (NMC) from 25 to 35 C. The logged group temperatures lie 10.6 to 10.9 K apart from 25 to 35 C and 8.2 to 8.8 K apart from 15 to 25 C; the 10 K setpoint difference in their place moves s by at most 0.2 and 1.0 % per K respectively.
- The replicate floor here is mostly resolution. Except in the 0-100 % cells at 0.5C, which are checked at every cycle, the checks lie a median 60 to 470 equivalent full cycles apart, and the NCA and NMC cells of the 1C series reach 90 % within their first or second check interval, so a lifetime is an interpolation within one interval (the interval at the crossing is a median 26, 41 and 122 % of L for LFP, NCA and NMC). The pair differences average 17.5 % (LFP), 10.8 % (NCA) and 37.5 % (NMC), with 95th percentiles of 89, 26 and 154 %, but 76 of the 86 lie below their check spacing. The NMC value also carries an artefact: in its 20-80 % cells at 0.5C, runs of fixed 1.80 Ah cycling discharges reach the cut-off as the window drifts and count as checks at 0.61 of the first, setting three of the four lifetimes; without them the NMC floor is 17.8 % (post hoc, D53).
- The warmer replicate. 59 of the 60 pairs with a swing for both cells differ in lifetime by less than their check spacing, so the registered sign test rests on one pair, which clears its resolution by 13 cycles: P7's verdict carries no information. By the raw sign the warmer cell lived shorter in 24 of 60 (40 %, as in the KIT data; post hoc). In the 1C pairs the swing differs by less than 1 K in 83 % (median 0.29 K).
- Back to the module. The NCA lever puts Stage 1's 0.391 K offer at a 0.80 % lifetime divergence (interval 0.47 to 1.15 %), below the 3.0 % mean floor; with its negative sign, the hotter cell would live longer.

| | Prediction (abridged) | Result | Computed |
|--|--|--|--|
| P1 | NMC, 25 to 35 C: s between -3.5 and +0.5 % per K | Held | -2.44 % per K |
| P2 | NMC, 15 to 25 C: s < 0 and larger in magnitude than at 25 to 35 C | Held | -5.62 % per K |
| P3 | NCA, 25 to 35 C: s between 0 and 5 % per K | Failed | -2.04 % per K |
| P4 | LFP's s above NCA's and NMC's, with Law 2's sign, within a factor of 3 of it | Not scored | not determinable (a cell without temperature); post hoc 2.92 or 3.09 times Law 2 by reading |
| P5 | every chemistry's floor mean between 1.5 and 10 % | Failed | 17.5, 10.8 and 37.5 % |
| P6 | at 1C or below, at least 80 % of pairs differ in swing by less than 1 K | Held | 83 % (median 0.29 K, 36 pairs) |
| P7 | the warmer cell lives shorter in 35 to 65 % of pairs | Failed | 0 of 1 untied pair (59 ties) |
| P8 | that share larger at 2C and 3C than at 1C or below | Not scored | 1 and 0 untied pairs |
| P9 | 0.391 K under the NCA lever implies less than 1.5 % | Held | 0.80 % |
| P10 | 15 to 25 C: s < 0 in all three chemistries | Failed | LFP +3.37 % per K |

Checks (`verification/stage3/`). V1: two of the five cells agree within the resolution; in the two LFP cells at 15 C and 2C the lifetime from the 2C capacity is 401 and 435 cycles shorter than from the 0.5C checks, and in the LFP cell checked at every cycle the two crossings of 95 % lie 76 cycles apart (D54). V2 and the independent read: a separate checker wrote its own code from DECISIONS.md and the inputs and opened no pipeline code or Stage 3 result until its outputs were frozen. The three V2 cells agree to 7e-15 (a second implementation in plain Python to 6e-15), every count of D50 reproduces cell by cell, and with the pipeline's two readings of D41 and D43 every per-cell value and all 1,306 group-level fields reproduce to 7e-13; under its own readings the ten verdicts are the same (`verification/stage3/v2_independent_report.md`, D55). V3: the laws reproduce Stage 1's registered levers exactly, and the KIT values are read unchanged. V4: the Stage 2 synthetic tests pass on the Sandia reader, and a Stage 3 test with planted lifetimes and temperatures and the faults of D50 recovers them to 2e-12 (`tests/stage3_synthetic.py`).

Limits. The check spacing limits every lifetime, and with it the floor and the sign test. LFP's 95 % threshold lies on a flat fade, where one check can move the crossing by hundreds of cycles (one 35 C cell crosses at 222 cycles, its replicates at 558 to 655). Discharge steps of one or two records at the end of a CV charge split some checks; without that splitting NCA's lever would be -2.00 and NMC's -2.29 % per K (D55). A skin thermocouple's offset enters T_cyc but not the swing; one NCA cell at 15 C and 1C logs a swing of 0.02 K against 4.0 K for its replicate up to cycle 200, and without it the NCA 15 to 25 C lever is -4.1 % per K. Two cells at 35 C log no cell temperature before cycle 259. The data come from one test campaign of one laboratory, with one or two cells at 15 C.

What it means for the thesis. In a second laboratory's data, near 30 C, the lifetimes of both layered-oxide cells rise with temperature at a 1C discharge, by about 2 % per K, as for the same NMC cell at KIT; Stage 1's 0.39 K module offer then moves lifetimes by under 1 %, inside every floor measured. LFP alone ages markedly faster when warm, but its registered lever is not determinable. Neither dataset shows a thermal share among replicates: at KIT the warmer replicate aged first in 40 % of pairs, and here the checks are too far apart to order 59 of 60 pairs.

The Stage 3 figure is `results/stage3_figure.png`; every number is in `results/stage3_*.csv` and `.json`, the post-hoc numbers in `results/stage3_posthoc.json`.

## Stage 4: the measured levers in the module, and the experiment that would settle it (v1.3.0)

Stage 1 built the module that offers temperature differences between cells; Stages 2 and 3 measured the lever between chamber temperatures. Stage 4 joins the two and then sizes the experiment the thesis statements promise (DECISIONS.md D56 to D63). Stage 1's eight-cell LG M50 module is solved at 35, 25 and 10 C ambient, at average rates C/4, C/2, 1C and 2C and h = 5, 10 and 25 W/(m2 K), each cell's constant-current heat taken from the resistance that Stage 1's PyBaMM procedure (DFN, ORegan2022 parameters) gives at the ambient. The measured levers are applied to its cell-to-cell spread: at 35 and 25 C the warm band, the twelve KIT 25 to 40 C levers and the Sandia NMC 25 to 35 C lever of the same cell; at 10 C the cold band, the twelve KIT 10 to 25 C levers and the Sandia NMC, NCA and LFP 15 to 25 C levers, LFP on its own line. The implied divergence of a case under a lever is |s| times the spread, the lifetime difference between the hottest and the coldest cell, which is the quantity the duplicate floors measure; a case offers a thermal share when its median over the band's levers exceeds a floor's mean, and a detectable one when that median exceeds the floor's 95th percentile. The author fixed both readings on 27 October 2026, before the predictions were pushed: the median lever is the one at the median of |s|, and the median is taken over the levers, the eight cells entering through the spread (in this symmetric module the four corner cells share one temperature and the four inner cells another).

The seven predictions in [predictions_stage4.md](predictions_stage4.md) were committed in `c1e58a19d1653bb477de485c4fc9f68e956ee070` and pushed at 09:50 UTC on 27 October 2026, before any Stage 4 computation. The code followed (`e5ecf58`) and was corrected once before any result was committed (`aa74702`: scipy's noncentral t returns nan for a vanishing tail, D64); the results were pushed in `6f658f1` at 10:07 UTC.

- The module. PyBaMM's resistance rises from 18.71 mOhm at 35 C to 26.65 at 25 C and 47.95 at 10 C. At h = 10 the spread is 0.39, 1.52 and 5.41 K at C/2, 1C and 2C at 35 C (Stage 1's values, reproduced bitwise), 0.58, 2.23 and 7.57 K at 25 C, and 1.11, 4.12 and 12.34 K at 10 C; at 1C the 10 C spread is 2.65 to 2.72 times the 35 C spread over the three coolings. The resolved cell's own spread at 10 C is 0.94 K at 1C (0.23 of the module's) and 3.77 K at 2C (0.31).
- The implied divergence. The median levers are +1.30 % per K (warm; KIT 0-100 0.33/0.33C) and 3.58 % per K (cold; KIT 0-100 1.67/1.00C, negative). At 35 C no case up to 1C exceeds the Geslin mean floor (at most 1.98 %); 2C gives 6.2 to 7.0 %, above it, but below the 95th percentile, though 5 of the 13 warm levers exceed 8.7 % there (up to 18.7 %). At 25 C a share appears at 2C (8.5 to 9.9 %), detectable at h = 10 and 25. At 10 C a share appears from C/2 (3.5 to 4.0 %) and a detectable one from 1C (13.4 to 14.7 %). Against the KIT floor (7.3 and 20.5 %) a 35 C module never offers a share (largest 7.05 %), a 25 C module from 2C, and a 10 C module from 1C, detectable from 2C.
- The heat is an upper bracket. The resistance is that of the ambient, and the cells run warmer, where the model's resistance is lower (47.9 mOhm at 10 C, 18.7 at 35 C, 12.3 at 50 C). With each cell's heat at its own temperature (post hoc, D70, repeating the independent checker's sensitivity), the spreads fall to 0.96 to 1.09 K at 1C and 2.4 to 2.9 K at 2C in 35 C air, and to 1.5 to 2.0 K at 1C in 10 C air; the 1C divergence in 10 C air is then 5.5 to 7.3 %, a share against the Geslin mean but not a detectable one, and under this heat P2 and P5 would hold and P4 would fail. In 35 C air the map keeps its shape (a share at 2C, none detectable). With a single resistance of 25 or 35 C in 10 C air the 1C divergence is 7.6 to 8.6 % or 5.4 to 6.1 % (D66). At 2C under the registered heat the hottest cell reaches 60 to 133 C at every ambient, outside a cell's usual range and far outside the levers' secants.
- The experiment. With n split pairs, two duplicate cells held dT apart in each, a paired t-test of their log lifetimes at 5 % (two-sided) and 80 % power, and the floor's pair difference as the noise, log-normal through its mean and 95th percentile: a cold chamber at 10 C needs 3 pairs held 5 K apart under its median lever (6 at 2 K), and a warm chamber at 35 C 6 pairs at 5 K and 25 at 2 K, against the Geslin floor; against the KIT floor 5, 19, 23 and 128. The central KIT value at 1C (-0.42 % per K, 10-100 1.00/1.00C) needs 221 pairs at 2 K, the Sandia 1C value (-2.44 % per K) 9; across the warm band's levers, 2 K needs 6 to 1,042 pairs. A normal floor shape needs 0.71 to 1.00 times as many, and a simulation of the t-test under the log-normal shape up to 25 % fewer at small n (V4, D68).

| | Prediction (abridged) | Result | Computed |
|--|--|--|--|
| P1 | 35 C: no case up to 1C above the 3.0 % mean floor (median warm lever) | Held | at most 1.98 % |
| P2 | 35 C, 2C: above 3.0 % at h = 5, and no warm lever above 8.7 % | Failed | 6.22 % at h = 5, but 5 of 13 levers above 8.7 % (up to 16.5 %) |
| P3 | 10 C spread at 1C at least 1.3 times the 35 C spread | Held | 2.65 to 2.72 times |
| P4 | 10 C, 1C: above 3.0 % at every cooling, above 8.7 % at h = 5 | Held | 13.4 to 14.7 %; 13.4 % at h = 5 |
| P5 | KIT floor: no share at 35 C, and at 10 C only at 2C | Failed | none at 35 C (largest 7.05 % against 7.33 %), but 1C at 10 C (13.4 to 14.7 %) |
| P6 | within-cell spread at 1C and 10 C below a third of the module's | Held | 0.94 against 4.12 K (0.23) |
| P7 | cold chamber, 5 K: under 10 pairs; warm chamber, 2 K: over 30 pairs | Failed | 3 pairs; 25 pairs (23 simulated) |

Checks (`verification/stage4/`). V1: the nine 35 C cases at C/2, 1C and 2C reproduce Stage 1's module values bitwise, and the procedure rerun at 35 C gives Stage 1's resistance bitwise. V2: separate code that applies constant currents to the DFN without PyBaMM's experiment interface agrees with the procedure's resistance within 0.02 % at all three ambients. V3: two cases recomputed by hand agree exactly. V4: with normal differences the simulated t-test agrees with the calculated n within one pair up to n = 50 and within 2.3 % above 100; with the registered log-normal differences it agrees within 5 % in 82 of the 144 cases and within 3 % wherever n exceeds 100, and at small n needs up to 25 % fewer pairs, so the calculated numbers are an upper estimate there and the simulated ones stand beside them in `results/stage4_experiment.csv`. The independent read: a separate checker wrote its own code from D56 to D63 and the inputs, and opened no Stage 4 code or result until its outputs were frozen (`fd5d1e7`, 10:33 UTC). Its own network solver, resistance runs, resolved cell, divergences, power calculation (by quadrature) and scoring agree with the registered results to 2e-12 at most, with the same seven verdicts; only its seeded simulations differ, by up to 0.8 % in n. It found that the registered 10 C resistance is right only because the procedure's 10 s samples step over a sub-second solver spike in the 10 C charge, which tight tolerances remove (D69), and that the heat's bracket decides P2 and P5 as well as P4 (D70) (`verification/stage4/independent/report.md`).

Limits. The levers were measured between chamber set points, not between the cells of one module; applying them to a within-module spread assumes the lever is local and linear over a few kelvin, and 17 of the 36 cases have cells outside their band's set points. The heat at 25 and 10 C extrapolates a published parameter set and takes the resistance at the ambient, an upper bracket wherever the cells run warmer (D66, D70). The procedure's solver at its default tolerances produces isolated heat spikes: they set its C/16 charge values at 35 and 25 C (D65), and one in the 10 C C/2 charge falls between the 10 s samples, so the registered R(10 C) agrees with the tightly solved value (D69); two of its non-registered runs fail to solve. The module is Stage 1's lumped network with steady constant-current heat. The floor shape of the experiment is a choice, the pair difference is assumed to be held steady, and the duplicates are assumed to differ as the floors' cells did. Nothing here measures a thermal share.

What it means for the thesis. In 35 C air the measured levers produce a thermal share above the Geslin mean floor only at 2C, and a detectable one at no rate up to 2C, whichever way the heat is bracketed; in 10 C air, with more heat and the steeper cold levers, the module offers one from C/2 and a detectable one from 1C under the registered heat, and from 1C and 2C with each cell's heat at its own temperature. The experiment that would settle whether temperature differences drive duplicate divergence is small in the cold: three duplicate pairs held 5 K apart at 10 C detect the median cold lever against the Geslin floor. At 35 C, where the lever's sign depends on the protocol, 2 K needs about 25 pairs at the median lever and hundreds where the lever is small.

The Stage 4 figure is `results/stage4_figure.png`; every number is in `results/stage4_*.csv` and `.json`, the quoted ones in `results/stage4_paper_numbers.json` and the post-hoc ones in `results/stage4_posthoc.json`.

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

- Stage 2: the KIT dataset (`inputs/kit/MANIFEST.md`, `inputs/kit/cyclic_cells.csv`); the raw files are not redistributed and are fetched from doi.org/10.35097/1969.
- Stage 3: the Sandia files of Preger et al. 2020 (`inputs/sandia/MANIFEST.md`, every file's SHA-256 in `inputs/sandia/files_sha256.csv`, `inputs/sandia/inventory.csv`). Battery Archive shares them on request; the raw files are not redistributed, and their use asks for citations of the article and of www.batteryarchive.org. Per-cell derived values are in `results/stage3_cells.csv` and `stage3_checks.csv`.

- Stage 4: no new data. It reads Stage 1's module, resistance procedure and resolved cell, the KIT levers and floor (`results/stage2_lever.json`), the Sandia levers (`results/stage3_lever.csv`) and the Geslin floor (`inputs/floor.json`), every value listed with its file in D56.

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
- `DECISIONS.md`: D1 to D20 fixed before the predictions, D21 before the push, D22 to D24 after it; Stage 2's D25 to D35 before its predictions, D36 and D37 after; Stage 3's D38 to D49 before its predictions, D50 before its results, D51 to D55 after; Stage 4's D56 to D63 before its predictions, D64 to D70 after its results.
- Stage 2: `results/stage2_cells.csv` (every cell's lifetimes, censoring and temperatures), `stage2_checkups.csv`, `stage2_sets.csv`, `stage2_lever.csv` and `.json`, `stage2_pairs.csv`, `stage2_warmer_pairs.csv`, `stage2_scorecard.json`, `stage2_figure.png`; `verification/stage2/`. The Sandia code in `src/stage2/` (sandia_io, steps, lifetimes, laws, analysis, verify), written before the Sandia files existed, is used unchanged by Stage 3.
- Stage 3: `results/stage3_cells.csv` (every cell's checks, lifetime, resolution, sensitivity lifetime and temperature measures), `stage3_checks.csv`, `stage3_lever.csv` and `.json`, `stage3_pairs.csv`, `stage3_pairs_oriented_measure_a.csv` and `_b.csv`, `stage3_floor.csv`, `stage3_thermal_share.json`, `stage3_scorecard.json`, `stage3_prepare.json` (the counts of D50), `stage3_posthoc.json` (D51 to D54), `stage3_figure.png`; code in `src/stage3/`; `verification/stage3/` (V1, V3, V4 and the independent check).
- Stage 4: `results/stage4_resistance.json` (the procedure's every run at 35, 25 and 10 C), `stage4_module.csv` (every case: resistance, heat, the eight cell temperatures, spread, offsets), `stage4_resolved_cell.csv` and `.json`, `stage4_divergence.csv` (case by lever: offsets, implied divergence with its interval, the band's median and 95th percentile), `stage4_map.json` (verdicts per case and the map), `stage4_experiment.csv` and `.json` (pairs needed, calculated and simulated, with the floor shapes), `stage4_experiment_curves.csv`, `stage4_scorecard.json`, `stage4_posthoc.json` (D65 to D68), `stage4_posthoc_own_temperature.json` (D70), `stage4_paper_numbers.json`, `stage4_figure.png`; code in `src/stage4/`; `verification/stage4/` (V1 to V4, the heat traces of the procedure and the independent check); `tests/stage4_synthetic.py`.

## How to run

```
# Stage 2: unpack cfg.zip and cell_eocv2.zip of doi.org/10.35097/1969 into CFG_DIR and EOC_DIR
python3 src/stage2/kit.py --eoc EOC_DIR --out results                 # Stage 2 computations and scorecard
python3 src/stage2/kit_checks.py --cfg CFG_DIR --eoc EOC_DIR           # V1 and V3
python3 tests/kit_synthetic.py WORKDIR SEED                            # V4
python3 src/stage2/make_figure_kit.py                                  # the Stage 2 figure
# Stage 3: unpack the three Battery Archive zips (SNL LFP, NCA, NMC) into RAW_DIR, then for every time series
#   python3 -I src/stage2/steps.py FILE --i-rest-A <0.005 x rated Ah> --out STEPS_DIR/<cell>_steps.csv
python3 -I src/stage3/prepare.py --raw RAW_DIR --steps STEPS_DIR --work WORK --out results/stage3_prepare.json   # D50
python3 -I src/stage3/stage3.py --raw RAW_DIR --work WORK --out results            # Stage 3 computations and scorecard
python3 src/stage2/verify.py v1 --cells results/stage3_cells.csv --ids ... --out verification/stage3/v1.json   # V1
python3 -I tests/stage3_synthetic.py WORKDIR                                       # the Stage 3 synthetic test
python3 src/stage3/make_figure.py --results results                                # the Stage 3 figure
python3 -I src/stage3/posthoc.py --work WORK --results results                     # post hoc (D51 to D54)
# Stage 4 (reads only archived results; PyBaMM, scikit-fem and pyamg as in requirements.txt)
python3 -I src/stage4/resistance.py                     # Stage 1's procedure at 35, 25 and 10 C
python3 -I src/stage4/module.py                         # the module at three ambients, the resolved cell at 10 C, V1
python3 -I verification/stage4/v2_resistance.py         # V2
python3 -I src/stage4/divergence.py                     # the implied divergence and the map
python3 -I verification/stage4/v3_hand.py               # V3
python3 -I src/stage4/experiment.py                     # the pairs needed, with V4 (several minutes)
python3 -I src/stage4/score.py                          # the scorecard
python3 -I src/stage4/make_figure.py --results results  # the Stage 4 figure
python3 -I src/stage4/posthoc.py                        # post hoc (D65 to D68)
python3 -I src/stage4/posthoc_own_temperature.py        # post hoc: each cell's heat at its own temperature (D70)
python3 -I verification/stage4/procedure_heat_check.py  # the procedure's heat traces (D65)
python3 -I src/stage4/paper_numbers.py                  # the quoted numbers
python3 -I tests/stage4_synthetic.py                    # checks of the Stage 4 code on made-up inputs
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

## Tooling

The code and text in this repository were developed with the assistance of Claude (Anthropic) as a programming, writing and checking assistant, working from the author's brief; the study design and the five predictions are the author's, and the author is responsible for the decisions and interpretations.

## Licence

Code: MIT (LICENSE). Note, figure and results: CC BY 4.0. The two-cell study's solver files in `src/twocell/` are MIT, by the same author. The two dataset features used post hoc (Geslin et al.) appear per cell in `results/dataset_protocol_heat.csv` with attribution; the dataset is CC BY 4.0 at the Stanford Digital Repository (doi.org/10.25740/td676xr4322), and the authors' code repository carries no licence file of its own.

## References

Geslin, A., Xu, L., Ganapathi, D., Moy, K., Chueh, W. C. & Onori, S. Dynamic cycling enhances battery lifetime. Nature Energy 10, 172-180 (2025). https://doi.org/10.1038/s41560-024-01675-8

Kucinskis, G., Bozorgchenani, M., Feinauer, M., Kasper, M., Wohlfahrt-Mehrens, M. & Waldmann, T. Arrhenius plots for Li-ion battery ageing as a function of temperature, C-rate, and ageing state: An experimental study. Journal of Power Sources 549, 232129 (2022). https://doi.org/10.1016/j.jpowsour.2022.232129

Schimpe, M., von Kuepach, M. E., Naumann, M., Hesse, H. C., Smith, K. & Jossen, A. Comprehensive modeling of temperature-dependent degradation mechanisms in lithium iron phosphate batteries. Journal of The Electrochemical Society 165, A181-A193 (2018). https://doi.org/10.1149/2.1181714jes

O'Regan, K., Brosa Planella, F., Widanage, W. D. & Kendrick, E. Thermal-electrochemical parameters of a high energy lithium-ion cylindrical battery. Electrochimica Acta 425, 140700 (2022). https://doi.org/10.1016/j.electacta.2022.140700

Chen, C.-H., Brosa Planella, F., O'Regan, K., Gastol, D., Widanage, W. D. & Kendrick, E. Development of experimental techniques for parameterization of multi-scale lithium-ion battery models. Journal of The Electrochemical Society 167, 080534 (2020). https://doi.org/10.1149/1945-7111/ab9050

Sulzer, V., Marquis, S. G., Timms, R., Robinson, M. & Chapman, S. J. Python Battery Mathematical Modelling (PyBaMM). Journal of Open Research Software 9, 14 (2021). https://doi.org/10.5334/jors.309

Luh, M. & Blank, T. Comprehensive battery aging dataset: capacity and impedance fade measurements of a lithium-ion NMC/C-SiO cell [dataset, version 2: result data]. RADAR4KIT (2024). https://doi.org/10.35097/1969

Luh, M. & Blank, T. Comprehensive battery aging dataset: capacity and impedance fade measurements of a lithium-ion NMC/C-SiO cell [dataset]. RADAR4KIT (2024). https://doi.org/10.35097/1947

Luh, M. & Blank, T. Comprehensive battery aging dataset: capacity and impedance fade measurements of a lithium-ion NMC/C-SiO cell. Scientific Data 11, 1004 (2024). https://doi.org/10.1038/s41597-024-03831-x

Preger, Y., Barkholtz, H. M., Fresquez, A., Campbell, D. L., Juba, B. W., Romàn-Kustas, J., Ferreira, S. R. & Chalamala, B. Degradation of commercial lithium-ion cells as a function of chemistry and cycling conditions. Journal of The Electrochemical Society 167, 120532 (2020). https://doi.org/10.1149/1945-7111/abae37

Battery Archive. www.batteryarchive.org (source of the Sandia files, shared on request on 23 October 2026).

Satte, S. (2026). The cell-to-cell lifetime floor in the Geslin et al. dynamic-cycling dataset (v1.0.1). Zenodo. https://doi.org/10.5281/zenodo.23205413

Satte, S. (2026). Where and when a neighbour cell first reaches thermal-runaway onset temperatures: resolved against lumped, a pre-registered two-cell LG M50 pilot (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.23253738
