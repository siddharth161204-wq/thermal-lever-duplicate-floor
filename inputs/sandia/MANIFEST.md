# Sandia 18650 cycling data (Preger et al. 2020): what was taken

Source: Battery Archive (www.batteryarchive.org), which holds the Sandia National Laboratories study of commercial 18650 LFP (A123 APR18650M1A), NCA (Panasonic NCR18650B) and NMC (LG 18650HG2) cells. The archive offers no file downloads; on request it shared three zip files through Google Drive on 23 October 2026, which the author downloaded the same day.

| Zip | Bytes | SHA-256 |
|---|---|---|
| SNL NMC.zip | 181,470,424 | d84fd3940eb802a19ae2c28902263e8dc5d4d320b034d1406206e8ba56664137 |
| SNL NCA.zip | 90,993,145 | 39a4cad9b283e94a1b50283f01364ee168b38e1108c72b75b3ca209697890cee |
| SNL LFP.zip | 263,826,451 | 006a335cbcdabc858a85ab0cdbc59a7001150751cf22abe8a7132c85ef63223d |

The zips hold 172 CSV files, a cycle-data file and a time-series file for each of 86 cells (LFP 30, NCA 24, NMC 32); every file's name, size and SHA-256 is in `files_sha256.csv`. The files' own column names, matched by `src/stage2/sandia_io.py`: cycle data Cycle_Index, Start_Time, End_Time, Test_Time (s), Min_Current (A), Max_Current (A), Min_Voltage (V), Max_Voltage (V), Charge_Capacity (Ah), Discharge_Capacity (Ah), Charge_Energy (Wh), Discharge_Energy (Wh); time series Date_Time, Test_Time (s), Cycle_Index, Current (A), Voltage (V), Charge_Capacity (Ah), Discharge_Capacity (Ah), Charge_Energy (Wh), Discharge_Energy (Wh), Environment_Temperature (C), Cell_Temperature (C).

Conditions of use, as the archive's e-mail of 23 October 2026 states them: cite the article that describes the experiments; cite www.batteryarchive.org as the source of the data; send the archive a reference to work that uses the data. The raw files are not redistributed here; per-cell derived values (lifetimes, temperature summaries) are published in `results/` with this citation.

Citation: Y. Preger, H. M. Barkholtz, A. Fresquez, D. L. Campbell, B. W. Juba, J. Romàn-Kustas, S. R. Ferreira and B. Chalamala, Degradation of Commercial Lithium-Ion Cells as a Function of Chemistry and Cycling Conditions, Journal of The Electrochemical Society 167, 120532 (2020), doi:10.1149/1945-7111/abae37 (as printed on the article's IOPscience page; CC BY-NC-ND 4.0). Data source: www.batteryarchive.org.
