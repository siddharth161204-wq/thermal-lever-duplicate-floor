# Predictions, Stage 6

The lever within a fast-charging policy: in the fast-charging datasets of Severson et al. (Nature Energy, 2019) and Attia et al. (Nature, 2020), cells charged on the same policy in the same batch ran at different temperatures, and Stage 6 asks whether the warmer cell dies sooner, by how much per kelvin, whether that is the cell's own heat or its place in the chamber, and whether the ranking between policies follows the same lever.

Written before any file of these datasets was opened or listed. What is known of the data when this is committed is what the brief recalls and has not yet been verified at its source: that the lifetime is counted cycle by cycle, that charging at several C heats the cells by kelvins, that each cell's can temperature and internal resistance are logged every cycle with its charging policy, cycle life, barcode and channel, and that the cell is the A123 APR18650M1A, the LFP cell of the Sandia study. Known from this repository's archive (DECISIONS.md D82): Sandia's levers for this cell from Stage 3, s = +3.37 % per K from 15 to 25 C at 1C (2.95 to 3.83), +2.05 at 2C (1.18 to 2.91) and +11.0 from 25 to 35 C at 2C (8.3 to 13.7), with s = -d ln L / dT, positive when the warmer cell lives shorter; Stage 4's floors, the standard deviation of a pair's log-lifetime difference, 0.0445 for the Geslin duplicates and 0.1041 for the KIT replicates; and Stage 5's smallest detectable |s| for the KIT replicates, 6.56 % per K against the KIT floor. The choices the brief leaves open and the scoring rules are fixed in DECISIONS.md D83 onward before any lifetime is read, and those of P1 and P6 before the quantities they score are computed. Each prediction is followed by its expectation, the brief's, written with it and unscored. A prediction whose quantity the files turn out not to hold is recorded as not scorable.

P1. The pooled smallest detectable |s| under the KIT floor (A0.8) is below 3.4 % per K, Sandia's lever for this cell from 15 to 25 C.
Expectation: uncertain. Scored at the end of Part A0.

P2. s_w is positive and its 95 % interval excludes zero.
Expectation: uncertain.

P3. With the resistance held (B2), the temperature's slope keeps its sign and at least half its size.
Expectation: uncertain.

P4. s_b is positive and its interval excludes zero.
Expectation: likely (a policy that charges harder runs hotter and is expected to live shorter).

P5. The interval of s_b minus s_w excludes zero, with s_b the larger.
Expectation: likely.

P6. The reliability of the registered temperature (A0.5, odd against even cycles) exceeds 0.8.
Expectation: uncertain. Scored at the end of Part A0.

P7. s_w lies inside Sandia's interval for this cell from 25 to 35 C at 2C (8.3 to 13.7 % per K).
Expectation: unlikely.

P8. The within-set standard deviation of ln L after B1 lies between the Geslin floor (0.0445) and the KIT floor (0.1041).
Expectation: uncertain.
