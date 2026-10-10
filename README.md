# The thermal lever against the duplicate-cell floor

How large a temperature spread a lifetime divergence needs, and how large a module offers.

Author: Siddharth Satte, independent researcher, Pune, India (ORCID 0009-0001-3584-1008).

Nominally identical cells cycled under one protocol at one temperature still drift apart. In the Stanford dynamic-cycling dataset (Geslin et al., Nature Energy 2025; 92 SiOx-graphite/NCA 21700 cells, 47 protocols in duplicate, all at 35 C) the duplicate cells' lifetimes to 90 % state of health differ by 3.0 % on average and by 8.7 % at the 95th percentile, and the worst of eight identical cells falls 4.1 % below the module mean by selection alone (the floor, doi.org/10.5281/zenodo.23205413). A thermal explanation of cell-to-cell divergence in a module has to exceed that floor. This repository asks two questions before any module experiment or model is built:

1. The lever: how large a steady temperature difference between two cells does a published temperature-dependent ageing law need to produce a lifetime difference as large as the floor?
2. The availability: how large a temperature difference does a module at 35 C actually offer, cell to cell in a lumped network of eight 21700 cells, and inside one resolved cell?

## Pre-registration

The predictions in [predictions.md](predictions.md) were committed before any script computed a lever at 35 C, a module spread at a registered rate and cooling, or any field of the resolved cell, and the commit was pushed to this public repository before those runs. The pre-registration commit is recorded here once it exists.

## Inputs

- The floor: copied from the published re-analysis (`inputs/floor.json`, with the file hash).
- Two ageing laws, both open access and read in full (`inputs/SOURCES.md`): Kucinskis et al. 2022 (J. Power Sources 549, 232129), a 5 Ah NCA+LNO / Si-graphite 21700, the closest published chemistry to the dataset's cells, whose printed cycle counts to 90 % SoH against temperature and charge rate are the law; and Schimpe et al. 2018 (J. Electrochem. Soc. 165, A181), the full temperature-dependent LiFePO4/graphite model with its Table IV parameters. The brief's examples, Wang et al. 2011 and Waldmann et al. 2014, could not be read and are not used (DECISIONS.md D2).
- The proxy cell, an LG M50, with the geometry, jelly-roll conductivities and can of the two-cell study (doi.org/10.5281/zenodo.23253738), whose resolved scikit-fem cell is reused unchanged in `src/twocell/`.
- The resistance for the I^2 R heat is derived, not measured: 18.71 mOhm, the value whose I^2 R reproduces the steady C/2 cycling heat of PyBaMM's Doyle-Fuller-Newman model with the ORegan2022 parameter set at 35 C (D8).

The dataset's cells are a different 21700; the proxy sets magnitudes, not the dataset's exact values.

## How to run

```
python3 -m pip install numpy scipy matplotlib "scikit-fem==12.0.2" "pybamm==26.9.0.0"
python3 src/lgm50_resistance.py --out inputs/lgm50_resistance.json   # the derived resistance (minutes)
python3 src/verify_pre.py --out verification                         # V1 and V2
python3 src/resolved_cell_steady.py --verify --out verification      # V3
python3 src/thermal_lever.py --out results                           # Step 1 and Step 2(a)
python3 src/resolved_cell_steady.py --out results                    # Step 2(b)
```

## Licence

Code: MIT (LICENSE). Note, figure and results: CC BY 4.0. The two-cell study's solver files in `src/twocell/` are MIT, by the same author.

## References

Geslin, A., Xu, L., Ganapathi, D., Moy, K., Chueh, W. C. & Onori, S. Dynamic cycling enhances battery lifetime. Nature Energy 10, 172-180 (2025). https://doi.org/10.1038/s41560-024-01675-8

Kucinskis, G., Bozorgchenani, M., Feinauer, M., Kasper, M., Wohlfahrt-Mehrens, M. & Waldmann, T. Arrhenius plots for Li-ion battery ageing as a function of temperature, C-rate, and ageing state: An experimental study. Journal of Power Sources 549, 232129 (2022). https://doi.org/10.1016/j.jpowsour.2022.232129

Schimpe, M., von Kuepach, M. E., Naumann, M., Hesse, H. C., Smith, K. & Jossen, A. Comprehensive modeling of temperature-dependent degradation mechanisms in lithium iron phosphate batteries. Journal of The Electrochemical Society 165, A181-A193 (2018). https://doi.org/10.1149/2.1181714jes

O'Regan, K., Brosa Planella, F., Widanage, W. D. & Kendrick, E. Thermal-electrochemical parameters of a high energy lithium-ion cylindrical battery. Electrochimica Acta 425, 140700 (2022). https://doi.org/10.1016/j.electacta.2022.140700

Chen, C.-H., Brosa Planella, F., O'Regan, K., Gastol, D., Widanage, W. D. & Kendrick, E. Development of experimental techniques for parameterization of multi-scale lithium-ion battery models. Journal of The Electrochemical Society 167, 080534 (2020). https://doi.org/10.1149/1945-7111/ab9050

Sulzer, V., Marquis, S. G., Timms, R., Robinson, M. & Chapman, S. J. Python Battery Mathematical Modelling (PyBaMM). Journal of Open Research Software 9, 14 (2021). https://doi.org/10.5334/jors.309

Satte, S. (2026). The cell-to-cell lifetime floor in the Geslin et al. dynamic-cycling dataset (v1.0.1). Zenodo. https://doi.org/10.5281/zenodo.23205413

Satte, S. (2026). Where and when a neighbour cell first reaches thermal-runaway onset temperatures: resolved against lumped, a pre-registered two-cell LG M50 pilot (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.23253738
