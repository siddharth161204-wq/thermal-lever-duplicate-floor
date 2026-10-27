# Predictions, Stage 4

The measured levers run through Stage 1's module, and the experiment that would settle it: the KIT levers of Stage 2 and the Sandia levers of Stage 3, applied to the cell-to-cell spread of Stage 1's eight-cell LG M50 module at 35, 25 and 10 C ambient, against the Geslin and KIT floors; and the number of duplicate pairs, held a few kelvin apart, that would detect a thermal share against the floor.

Written from Stage 1's module and its archived 35 C values, the archived KIT and Sandia levers and the two floors, before any Stage 4 computation: no resistance had been derived at 25 or 10 C, no module case run at 25 or 10 C ambient or at C/4, no resolved-cell case run at 10 C, and no implied divergence or number of pairs computed. Known when committed: DECISIONS.md D56 to D63, including the author's two choices of 27 October 2026 (the median, or central, lever of a band is the lever at the median of |s| over the band; the distribution of the implied divergence is that of |s| times the spread over the band's levers). s is in % per K, positive when the warmer cells reach the threshold sooner; the scoring rules are D63. Each prediction is followed by its expectation, written with it and unscored.

P1. At 35 C ambient, no case up to 1C offers a thermal share above the 3.0 % mean floor under the median warm lever, at any cooling.
Expectation: likely (Stage 1's 1.4 to 1.5 K at 1C times 1.9 % per K at most is 2.9 %).

P2. At 35 C ambient and 2C, the module offers a thermal share above the mean floor at h = 5 under the median warm lever, and a detectable one (above 8.7 %) under none of the warm levers.
Expectation: uncertain.

P3. At 10 C ambient the module's spread at 1C exceeds its 35 C spread at 1C by at least 30 %, because the resistance rises in the cold.
Expectation: likely, conditional on PyBaMM's resistance at 10 C.

P4. At 10 C ambient and 1C, the module offers a thermal share above the mean floor under the median cold lever at every cooling, and a detectable one (above 8.7 %) at h = 5.
Expectation: likely for the first part, uncertain for the second.

P5. Under the KIT floors (7.3 and 20.5 %), no case at 35 C ambient offers a thermal share, and at 10 C ambient only 2C does.
Expectation: uncertain.

P6. The within-cell spread of the resolved cell at 1C and 10 C stays below one third of the module's cell-to-cell spread.
Expectation: likely.

P7. The experiment: at dT = 5 K in the cold chamber under the central cold lever, fewer than 10 duplicate pairs detect the thermal share against the Geslin floor; in the warm chamber at dT = 2 K under the central warm lever, more than 30 pairs are needed.
Expectation: uncertain; it depends on the floor's assumed shape.
