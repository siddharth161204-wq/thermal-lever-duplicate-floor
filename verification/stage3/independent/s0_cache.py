#!/usr/bin/env python3
"""Stage 3 independent check, step 0: read every raw Sandia time series once and cache the columns this check
uses as float64 arrays, in the order the rows are logged.

Written for the independent V2 check of DECISIONS.md D49; it uses none of the pipeline's functions.
The raw files are untrusted downloaded data: they are only parsed as CSV, never modified, and this script is run
with `python3 -I` from outside the data folders.

Usage: python3 -I s0_cache.py <cache_dir>
"""
import glob, hashlib, json, os, sys
import numpy as np
import pandas as pd

RAW = ["/tmp/sandia/raw_LFP", "/tmp/sandia/raw_NCA", "/tmp/sandia/raw_NMC"]
COLS = {"Test_Time (s)": "t", "Cycle_Index": "cyc", "Current (A)": "I", "Voltage (V)": "V",
        "Charge_Capacity (Ah)": "Qc", "Discharge_Capacity (Ah)": "Qd",
        "Environment_Temperature (C)": "Tenv", "Cell_Temperature (C)": "Tcell"}


def cell_id_from_file(path):
    # SNL_18650_LFP_15C_0-100_0.5-1C_a_timeseries.csv -> SNL_18650_LFP_15C_0-100_0.5/1C_a
    name = os.path.basename(path)[: -len("_timeseries.csv")]
    parts = name.split("_")
    # parts: SNL, 18650, CHEM, 15C, 0-100, 0.5-1C, a
    parts[5] = parts[5].replace("-", "/", 1)
    return "_".join(parts)


def main(cache_dir):
    os.makedirs(cache_dir, exist_ok=True)
    files = sorted(f for r in RAW for f in glob.glob(os.path.join(r, "*", "*_timeseries.csv")))
    index = {}
    for f in files:
        cid = cell_id_from_file(f)
        h = hashlib.sha256()
        with open(f, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        df = pd.read_csv(f, usecols=list(COLS), dtype={c: "float64" for c in COLS}, engine="c")
        arr = {COLS[c]: df[c].to_numpy(dtype=np.float64) for c in COLS}
        out = os.path.join(cache_dir, cid.replace("/", "-") + ".npz")
        np.savez(out, **arr)
        index[cid] = {"file": f, "sha256": h.hexdigest(), "rows": int(len(df)), "cache": out,
                      "nan_counts": {k: int(np.isnan(v).sum()) for k, v in arr.items()}}
        print(cid, len(df), index[cid]["nan_counts"], flush=True)
    with open(os.path.join(cache_dir, "index.json"), "w") as fh:
        json.dump(index, fh, indent=1)


if __name__ == "__main__":
    main(sys.argv[1])
