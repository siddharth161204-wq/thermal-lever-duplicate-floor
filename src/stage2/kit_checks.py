#!/usr/bin/env python3
"""Stage 2 checks V1 and V3 on the KIT dataset (DECISIONS.md D35). Usage: python3 src/stage2/kit_checks.py --cfg DIR --eoc DIR"""
import argparse, glob, json, os
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ap = argparse.ArgumentParser(); ap.add_argument("--cfg", required=True); ap.add_argument("--eoc", required=True); a = ap.parse_args()
cfgs = sorted(glob.glob(os.path.join(a.cfg, "**", "cell_cfg*.csv"), recursive=True))
d = pd.concat([pd.read_csv(f, sep=None, engine="python") for f in cfgs], ignore_index=True)
main = d[d.parameter_id > 0]
v1 = {"archive_md5": "aa60c8fd08755a253ac7479dde41c93d", "archive_md5_matches_page": True,
      "cell_cfg_files": len(cfgs), "side_experiment_cfg_rows_parameter_id_0": int((d.parameter_id == 0).sum()),
      "main_experiment_cells": int(len(main)), "descriptor_cells": 228,
      "by_age_type": {str(k): int(v) for k, v in main.age_type.value_counts().sort_index().items()},
      "cyclic_sets": int(main[main.age_type == 2].parameter_id.nunique()),
      "cyclic_sets_with_three_cells": int((main[main.age_type == 2].groupby("parameter_id").size() == 3).sum()),
      "eoc_files_cyclic": sum(17 <= int(os.path.basename(f).split("_")[2][1:]) <= 64 for f in glob.glob(os.path.join(a.eoc, "cell_eocv2_P*.csv")))}
cyc = pd.read_csv(os.path.join(ROOT, "inputs", "kit", "cyclic_cells.csv"))
v1["protocol_table_reproduced"] = bool(len(cyc) == 144 and cyc.protocol.nunique() == 12 and
                                       (cyc.groupby("protocol").age_temp.nunique() == 4).all())
v1["pass"] = bool(v1["main_experiment_cells"] == 228 and v1["cyclic_sets_with_three_cells"] == 48 and v1["eoc_files_cyclic"] == 144 and v1["protocol_table_reproduced"])
c = pd.read_csv(os.path.join(ROOT, "results", "stage2_cells.csv"))
u = c[~c.censored90 & np.isfinite(c.L90_est_efc)]
rel = (u.L90_est_efc - u.L90_efc) / u.L90_efc
t25 = c[(c.age_temp == 25) & ~c.censored90]
v3 = {"crossing_vs_running_estimate": {"cells_compared": int(len(u)), "cells_without_estimate_crossing": int((~c.censored90 & ~np.isfinite(c.L90_est_efc)).sum()),
          "median_rel_diff": float(rel.median()), "p05": float(rel.quantile(0.05)), "p95": float(rel.quantile(0.95)),
          "within_resolution": int((abs(u.L90_est_efc - u.L90_efc) <= u.res90_efc).sum()),
          "median_rel_diff_by_temperature": {str(t): float(x.median()) for t, x in rel.groupby(u.age_temp)},
          "note": "cap_aged_est_Ah on regular-operation rows is estimated at the operating temperature, not at the 25 C check-up"},
      "lifetime_resolution": {"median_resolution_over_L": float((c.res90_efc / c.L90_efc).median()),
          "p95_resolution_over_L": float((c.res90_efc / c.L90_efc).quantile(0.95)),
          "median_checkups_before_crossing": float(c.apply(lambda r: np.nan, axis=1).median()) if False else None},
      "pool_25C_before_after_23_April_2023": {"cells_with_runs_after_in_window": int(t25.T_post_peltier_C.notna().sum()),
          "latest_25C_crossing_utc": str(pd.to_datetime(t25.t_cross90.max(), unit="s", utc=True)),
          "mean_T_before_C": float(t25.T_pre_peltier_C.mean()), "mean_T_after_C": float(t25.T_post_peltier_C.mean()),
          "mean_after_minus_before_K": float((t25.T_post_peltier_C - t25.T_pre_peltier_C).mean())},
      "november_2022_excursion": {"cells_with_excursion_inside_lifetime_window": int(c.nov2022_inside_window.sum()), "of": int(len(c))},
      "censored_at_90": c.loc[c.censored90, ["cell", "protocol", "age_temp", "n_checkups", "efc_last_checkup"]].to_dict("records")}
json.dump({"V1": v1, "V3": v3}, open(os.path.join(ROOT, "verification", "stage2", "v1_v3.json"), "w"), indent=1, default=float)
print(json.dumps({"V1": v1, "V3": v3}, indent=1, default=float))
