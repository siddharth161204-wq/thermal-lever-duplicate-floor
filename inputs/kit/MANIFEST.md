# KIT LG INR18650HG2 ageing dataset: what was taken

Source: M. Luh and T. Blank, "Comprehensive battery aging dataset: capacity and impedance fade measurements of a lithium-ion NMC/C-SiO cell [dataset, version 2: result data]", RADAR4KIT, 2024, doi:10.35097/1969 (https://radar.kit.edu/radar/en/dataset/CoNQplSNoVeXExyV). Licence: CC BY 4.0. Version 1 of the dataset, which the data descriptor describes: doi:10.35097/1947. Data descriptor: M. Luh and T. Blank, Scientific Data (2024), doi:10.1038/s41597-024-03831-x.

Archive downloaded by the author on 21 October 2026: `10.35097-1969.tar`, 333,402,624 bytes, MD5 aa60c8fd08755a253ac7479dde41c93d (equal to the value on the RADAR4KIT page), SHA-256 5bb11e1a1d3eb5d8a13585bf32a9cfefcf34b94a488b5e0727994861f98d272e. It is a BagIt bag; the five payload files used were extracted and each matched the bag's own `manifest-md5.txt`:

| File | Bytes | MD5 (bag manifest, verified) |
|---|---|---|
| data/dataset/cfg.zip | 97,252 | b894469e86712b6021bb5364f7296891 |
| data/dataset/cell_eocv2.zip | 93,167,507 | b4f42914c8eaf176fd47e33e6f112d7c |
| data/dataset/Data Structure v05.3.xlsx | 186,499 | 576711692fd3b5aaa6e0bf7f09bff1e9 |
| data/dataset/Cycling Experiment Cell Overview (2024-09-26 15-31).xlsx | 74,470 | 7a5fc1770ee595c6c44a51a0bf344e65 |
| data/dataset/Readme v2.pdf | 1,226,669 | 9c77a65ecdb274ecf2544bcfde462a86 |

Not extracted or opened: cell_eisv2.zip, cell_plsv2.zip, result_plot_examples.7z, "Battery Aging Data Plausibility Check v05.3 (colored).xlsx", WLTP_Class_3b_cycle_power_profile.xlsx, the overview's PDF copy, and the separate log-data record (doi:10.35097/kww7jv8ajuvchcah). The raw files are not redistributed here; `inputs/kit/cyclic_cells.csv` holds the configuration of the 144 cyclic-ageing cells, read from cfg.zip.
