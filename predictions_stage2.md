# Predictions, Stage 2

The lever measured on public data with replicates: the KIT ageing dataset of Luh and Blank (LG INR18650HG2, NMC with a SiOx-graphite anode; twelve protocols at 0, 10, 25 and 40 C, three cells each, cell temperatures logged).

Written from the dataset's design, its readme and the Stage 1 results, and committed and pushed before any capacity, cap_aged_est_Ah, temperature, date or other value of the cell data was read. Known when committed: the archive's file list and hashes, the readmes, the configuration files, the column definitions, the check-up schedule and the cell counts (DECISIONS.md D25 and D26). The definitions and scoring rules are D27 to D34. s = -d ln L / dT, positive when the hotter cell reaches the threshold sooner; L is the EFC to 90 % of the first check-up capacity. Each prediction is followed by the study's expectation, a judgement written now and not scored.

P1. For the protocols with the lowest charge rate in the design, the 25 to 40 C lever s is positive in every SOC window, and their mean lies between 1 and 6 % per K (central 3 % per K).
Expectation: likely holds on the sign, since at a 0.33 C charge the time-driven losses that grow with temperature should dominate; the magnitude band is the less certain part.

P2. Within each SOC window, s at the highest charge rate is lower than s at the lowest charge rate (all three windows).
Expectation: likely holds; a 1.67 C charge adds lithium plating, which a warmer cell suffers less, pulling s down; requiring all three windows makes a failure in one window plausible.

P3. Across the twelve protocols, the 25 to 40 C point estimates of s take both signs: the lever's sign near 35 C is a property of the protocol, not of the cell. (Uncertain; may fail.)
Expectation: close to even; it fails if plating at 25 C is too weak to outweigh the time-driven losses at any protocol.

P4. The replicate floor of this dataset is within a factor of two of the Geslin floor: mean pair difference at 90 % between 1.5 and 6.0 %, 95th percentile between 4.4 and 17.4 % (central 3.0 and 8.7 %).
Expectation: more likely than not; one commercial cell model in controlled pools should spread like the Geslin cells, but the pairs span four temperatures and twelve protocols, which may widen the tail.

P5. Within triplets, the warmer cell reaches 90 % first in 45 to 65 % of pairs (central 55 %): with the cells in oil pools, the thermal share among replicates is weak.
Expectation: likely holds; logged differences within a pool are a fraction of a kelvin to about 2 K and include sensor offsets.

P6. No protocol's s exceeds 7.7 % per K, so at every protocol Stage 1's largest module offer of 0.391 K implies a divergence below the 3.0 % mean floor.
Expectation: likely holds; 7.7 % per K is half again Law 2's Arrhenius value at 35 C.

P7. The mean s of P1's protocols lies between Law 1's registered secant (-0.99 % per K) and Law 2's value (5.2 % per K).
Expectation: likely holds if P1's central value is near 3 % per K; it fails if the low-rate cells are more temperature-sensitive than the LiFePO4 law.

P8. From 10 to 25 C, s is negative (the colder cells reach 90 % sooner) for at least half of the protocols whose charge rate is 1C or higher.
Expectation: likely holds; charging at 1C and above at 10 C favours plating, which the readme reports for fast charging at low temperature.
