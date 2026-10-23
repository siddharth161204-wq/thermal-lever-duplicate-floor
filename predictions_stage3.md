# Predictions, Stage 3

The lever in a second laboratory's data: the Sandia National Laboratories 18650 cells of Preger et al. 2020 (LFP, NCA and NMC; the NMC is the LG 18650HG2, the KIT cell of Stage 2), shared by Battery Archive on 23 October 2026.

Written from the dataset's design, the Stage 1 laws and the Stage 2 KIT results, which the predictions use, before any Sandia lifetime or temperature value was read. Known when committed: the files' structure and the inventory (`inputs/sandia/inventory.csv`), and DECISIONS.md D38 to D49, including the author's choice of measure (a), the within-cycle swing, made because the cells cycle without rests. s is in % per K over the stated secant, positive when the warmer cells reach the threshold sooner; the scoring rules are D48. Each prediction is followed by the study's expectation in one sentence, written now and unscored.

P1 NMC, 25 to 35 C, 0 to 100 %, 0.5C charge, 1C discharge: s lies between minus 3.5 and plus 0.5 (central minus 1.5), inside the band the KIT protocols at 0.33/1.00 C and 1.00/1.00 C span for the same window (minus 0.62 to minus 2.69).
Expectation: more likely than not, since the same cell model at a 0.5C charge sits between the two KIT charge rates, though a 25 to 35 C secant need not match a 25 to 40 C one.

P2 NMC, 15 to 25 C: s is negative, and larger in magnitude than the 25 to 35 C value.
Expectation: likely, as every KIT protocol's 10 to 25 C lever was negative and larger than its 25 to 40 C lever.

P3 NCA, 25 to 35 C: s lies between 0 and 5 (central 2).
Expectation: uncertain; the nearest published law gives a weak lever of undetermined sign near 35 C, so the sign is close to even.

P4 LFP's 25 to 35 C s is larger than both NCA's and NMC's, and has the sign of Law 2 at the same conditions, within a factor of three of it.
Expectation: likely on the ordering and the sign, less so on the factor, since LFP ages slowly at these conditions and its lifetimes rest on a 95 % threshold.

P5 The replicate floor: for each chemistry with at least five pairs, the mean relative difference at the registered threshold lies between 1.5 and 10 % (central 4 %).
Expectation: likely, given the Geslin floor of 3.0 % and the KIT floor of 7.3 %.

P6 At discharge rates of 1C or below, the self-heating difference within a pair is below 1 K in at least 80 % of pairs (central 0.3 K).
Expectation: likely, since the swing at 1C is itself of the order of a kelvin and replicates differ by a fraction of it.

P7 The warmer cell, by self-heating, has the shorter life in between 35 and 65 % of pairs (central 50 %): at these offsets the thermal share among replicates is not detectable in the sign.
Expectation: likely, as in the KIT data, where the warmer replicate aged first in 40 % of pairs.

P8 That share is larger in the 2C and 3C pairs than in the pairs at 1C or below (scored only if each group holds at least five pairs).
Expectation: close to even; larger heating differences at high rates favour it, few pairs work against it.

P9 Under the measured NCA lever, Stage 1's offer of 0.391 K implies a lifetime divergence below 1.5 %, half the Geslin mean floor.
Expectation: likely, as 1.5 % needs a lever above 3.8 % per K.

P10 The 15 to 25 C lever is negative in all three chemistries.
Expectation: likely for NMC and NCA, less certain for LFP, which tolerates cold cycling at 0.5C charge better.
