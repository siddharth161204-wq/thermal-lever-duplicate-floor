# The fast-charging datasets: what was taken

Source: the Experimental Data Platform of the Toyota Research Institute, https://data.matr.io/1/, which holds the data of Severson et al., Nature Energy 4, 383 to 391 (2019), https://doi.org/10.1038/s41560-019-0356-8, and of Attia et al., Nature 578, 397 to 402 (2020), https://doi.org/10.1038/s41586-020-1994-5. Licence, as the site states it: "All data is released under CC BY 4" (https://creativecommons.org/licenses/by/4.0/). What the site, the papers and the authors' code say is in `SOURCES.md`.

The files were downloaded to the author's laptop on 16 and 17 November 2026 and read there (DECISIONS.md D83, D84); they are not redistributed here. The first downloads of the 2017-05-12 and 2018-04-12 structs stopped short when the laptop slept; those two were downloaded again, and the local copies read carry the browser's suffix " (1)" in their names. Each struct read has the size its HDF5 superblock records as its end-of-file address, so none is truncated.

| File | Address | Bytes | SHA-256 |
|---|---|---|---|
| 2017-05-12_batchdata_updated_struct_errorcorrect.mat | https://data.matr.io/1/api/v1/file/5c86c0b5fa2ede00015ddf66/download | 3,025,320,241 | 9d928ab978f0e3c70b31cb833a749fedd35094d01af76475d69b40aa3497f5ba |
| 2017-06-30_batchdata_updated_struct_errorcorrect.mat | https://data.matr.io/1/api/v1/file/5c86bf13fa2ede00015ddd82/download | 2,007,331,155 | 63ab200d09ecb237fee5ef3a5c5db76e3212e3206a0bd92f769e1427fed338b8 |
| 2018-04-12_batchdata_updated_struct_errorcorrect.mat | https://data.matr.io/1/api/v1/file/5c86bd64fa2ede00015ddbb2/download | 3,236,690,412 | 62c30e413b63e6144720e016deed3661fac8468641794a5807b123fe84717998 |
| 2019-01-24_batchdata_updated_struct_errorcorrect.mat | https://data.matr.io/1/api/v1/file/5dcef152110002c7215b2c90/download | 2,601,295,745 | 9604dfd58926a99d8aece76e1c407defb7359190edcc21e4d4dd217a726882af |

The batch pages: 2017-05-12, https://data.matr.io/1/projects/5c48dd2bc625d700019f3204/batches/5c86c0b5fa2ede00015ddf67; 2017-06-30, .../batches/5c86bf14fa2ede00015ddd83; 2018-04-12, .../batches/5c86bd64fa2ede00015ddbb3 (all under the project https://data.matr.io/1/projects/5c48dd2bc625d700019f3204); 2019-01-24, https://data.matr.io/1/projects/5d80e633f405260001c0b60a/batches/5dcef1fe110002c7215b2c94.

Also read, not redistributed:

| File | Address | Bytes | SHA-256 |
|---|---|---|---|
| Severson_NatureEnergy_2019.pdf | https://web.mit.edu/braatzgroup/Severson_NatureEnergy_2019.pdf | 6,060,054 | ad87cfc07a96277c96173f969b2f0b45275bc27c1fc3131b494fe24df7a1879d |
| qt117547r6_noSplash_5559418521be7b1b864a46aa25b416b9.pdf (Attia et al.) | https://escholarship.org/content/qt117547r6/qt117547r6_noSplash_5559418521be7b1b864a46aa25b416b9.pdf | 3,867,928 | d148e885567bb10a21605c8dfe5db1cbb1ec4f129a2d843ef646b85b705ff92a |
| metadata.zip (the 2025 dataset, CC BY-NC 4.0; Part A0.9 only) | https://data.matr.io/11/metadata.zip | 51,659 | eb28c535ea1a899191e8e739dab70ba5ece10a80772c9754a773464fb5740ef4 |

Not taken: the four closed-loop rounds of the second paper (oed_0 to oed_3), which hold no measured lifetime and no temperature; the first paper's low-rate data and structured-data archive; the JSON archives of each batch.

Derived from the structs and committed (CC BY 4.0, with this attribution): `a0/<label>_a0_cells.csv`, `a0/<label>_a0_cycles.csv` and `a0/<label>_a0_structure.json`, written where the files are by `src/stage6/a0_read.py` (labels b1, b2, b3 and b9 for the four files in the order above): each entry's descriptors and, for cycles 1 to 100 only, its per-cycle Tavg, Tmax, Tmin and IR. After Part A0 was committed: `stage6_lifetimes.csv`, from `src/stage6/l_read.py` (D94).
