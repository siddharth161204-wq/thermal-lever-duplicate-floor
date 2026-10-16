#!/usr/bin/env python3
"""Stage 2: reading the Sandia National Laboratories 18650 files of batteryarchive.org (Preger et al. 2020).

Cell IDs, nominal conditions, rated capacities, file discovery, column matching and a chunked reader for the two
file types, the cycle data and the time series. Nothing in this module computes a lifetime or a temperature
statistic; it only reads structure.
"""
import glob, hashlib, os, re

import numpy as np
import pandas as pd

# Rated capacity per chemistry, Ah: the "ah" field of every SNL cell in batteryarchive.org's cycle-test cell list
# (query of the public list, 16 October 2026; DECISIONS.md, Stage 2). The same value for every cell of a chemistry.
RATED_AH = {"LFP": 1.1, "NCA": 3.2, "NMC": 3.0}
CELL_MODEL = {"LFP": "A123 APR18650M1A", "NCA": "Panasonic NCR18650B", "NMC": "LG 18650HG2"}

# SNL_18650_LFP_15C_0-100_0.5/1C_a as the archive lists it; file names may carry "-" in place of "/".
CELL_RE = re.compile(r"SNL_18650_(LFP|NCA|NMC)_(\d+)C_(\d+)-(\d+)_([\d.]+)[/-]([\d.]+)C_([a-z])")


def parse_cell_id(text):
    """Conditions encoded in an archive cell ID (or a file name containing one)."""
    m = CELL_RE.search(text)
    if not m:
        return None
    chem, T, lo, hi, cch, cdis, letter = m.groups()
    cid = f"SNL_18650_{chem}_{T}C_{lo}-{hi}_{cch}/{cdis}C_{letter}"
    return {"cell_id": cid, "key": cid.replace("/", "-"), "chemistry": chem, "ambient_C": float(T),
            "soc_min": float(lo), "soc_max": float(hi), "c_charge": float(cch), "c_discharge": float(cdis),
            "replicate": letter, "rated_Ah": RATED_AH[chem],
            "group": f"{chem}_{T}C_{lo}-{hi}_{cch}/{cdis}C"}


def sha256(path, block=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def norm(name):
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


# Column candidates, matched after removing case, spaces, brackets and units. The archive's standard names first.
COLS = {
    "time": ["Test_Time (s)", "test_time", "testtime", "Test_Time(s)"],
    "datetime": ["Date_Time", "date_time", "datetime"],
    "cycle": ["Cycle_Index", "cycle_index", "cycleindex"],
    "current": ["Current (A)", "current", "i"],
    "voltage": ["Voltage (V)", "voltage", "v"],
    "q_ch": ["Charge_Capacity (Ah)", "charge_capacity", "ah_c"],
    "q_dis": ["Discharge_Capacity (Ah)", "discharge_capacity", "ah_d"],
    "e_ch": ["Charge_Energy (Wh)", "charge_energy", "e_c"],
    "e_dis": ["Discharge_Energy (Wh)", "discharge_energy", "e_d"],
    "t_cell": ["Cell_Temperature (C)", "cell_temperature", "celltemperature"],
    "t_env": ["Environment_Temperature (C)", "environment_temperature", "env_temperature", "environmenttemperature"],
}
CYCLE_COLS = {
    "cycle": ["Cycle_Index", "cycle_index"],
    "start_time": ["Start_Time", "start_time"],
    "end_time": ["End_Time", "end_time"],
    "time": ["Test_Time (s)", "test_time"],
    "i_min": ["Min_Current (A)", "min_current"],
    "i_max": ["Max_Current (A)", "max_current"],
    "v_min": ["Min_Voltage (V)", "min_voltage"],
    "v_max": ["Max_Voltage (V)", "max_voltage"],
    "q_ch": ["Charge_Capacity (Ah)", "charge_capacity"],
    "q_dis": ["Discharge_Capacity (Ah)", "discharge_capacity"],
    "e_ch": ["Charge_Energy (Wh)", "charge_energy"],
    "e_dis": ["Discharge_Energy (Wh)", "discharge_energy"],
}


def match_columns(columns, table=COLS):
    """Map each logical name to the file's column, or None. Exact normalised match first, then a unit-stripped match."""
    nc = {norm(c): c for c in columns}
    stripped = {norm(re.sub(r"\(.*?\)", "", c)): c for c in columns}
    out = {}
    for key, cands in table.items():
        hit = None
        for cand in cands:
            if norm(cand) in nc:
                hit = nc[norm(cand)]
                break
        if hit is None:
            for cand in cands:
                k = norm(re.sub(r"\(.*?\)", "", cand))
                if k in stripped:
                    hit = stripped[k]
                    break
        out[key] = hit
    return out


def header(path):
    return list(pd.read_csv(path, nrows=0).columns)


def discover(root):
    """Pair every SNL cell found under root with its cycle-data and time-series files (CSV, possibly zipped)."""
    files = [p for p in glob.glob(os.path.join(root, "**", "*"), recursive=True) if os.path.isfile(p)]
    cells = {}
    for p in sorted(files):
        meta = parse_cell_id(os.path.basename(p))
        if meta is None:
            continue
        base = os.path.basename(p).lower()
        kind = "timeseries" if ("timeseries" in base or "time_series" in base) else ("cycle" if "cycle" in base else "other")
        c = cells.setdefault(meta["key"], dict(meta, files={}))
        c["files"].setdefault(kind, []).append(p)
    return cells


def iter_timeseries(path, cols, chunksize=500_000):
    """Yield chunks of the time series as dicts of float arrays for the requested logical columns, in file order."""
    m = match_columns(header(path))
    use = {k: m[k] for k in cols if m.get(k)}
    missing = [k for k in cols if not m.get(k)]
    for chunk in pd.read_csv(path, usecols=list(use.values()), chunksize=chunksize, low_memory=False):
        yield {k: pd.to_numeric(chunk[c], errors="coerce").to_numpy(dtype=float) for k, c in use.items()}, missing


def read_cycle_data(path):
    df = pd.read_csv(path, low_memory=False)
    m = match_columns(df.columns, CYCLE_COLS)
    out = pd.DataFrame({k: pd.to_numeric(df[c], errors="coerce") if k not in ("start_time", "end_time") else df[c]
                        for k, c in m.items() if c})
    return out, m
