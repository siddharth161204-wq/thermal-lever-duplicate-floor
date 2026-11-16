#!/usr/bin/env python3
"""Stage 6, Part A0: the structure, the descriptors and the early-cycle temperatures and resistances of one batch file
(DECISIONS.md D84).

For each cell of the struct array `batch` it reads the fields policy_readable, policy, channel_id and barcode, and from
the cell's summary the cycle numbers of its first 101 entries and, for the entries numbered 1 to 100, Tavg, Tmax, Tmin
and IR. It never opens cycle_life, cycles or Vdlin, nor the summary's QDischarge, QCharge or chargetime, and it writes
no array length: a cell's own number of cycles is not read. Runs where the files are, on Python 3.10 with numpy and
h5py only.

Usage: python3 -I a0_read.py <file.mat> <label> <output directory>
Writes <label>_a0_cells.csv, <label>_a0_cycles.csv and <label>_a0_structure.json in the output directory.
"""
import csv, json, os, platform, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import h5py                 # noqa: E402
import numpy as np          # noqa: E402
import matv73 as M          # noqa: E402

DESCRIPTORS = ("policy_readable", "policy", "channel_id", "barcode")
SERIES = ("Tavg", "Tmax", "Tmin", "IR")
NEVER_BATCH = ("cycle_life", "cycles", "Vdlin")
NEVER_SUMMARY = ("QDischarge", "QCharge", "chargetime")
FIRST, LAST = 1, 100        # the cycles read; the measures use 2 to 100 (D86)
PEEK = 101                  # cycle numbers of the first 101 summary entries, enough whether numbering starts at 0, 1 or 2


def main(path, label, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    f = h5py.File(path, "r")
    b = f["batch"]
    batch_keys = sorted(b.keys())
    n = int(b["summary"].shape[0])
    structure = {
        "file": os.path.basename(path), "label": label, "root_keys": sorted(f.keys()), "batch_keys": batch_keys,
        "batch_field_shapes": {k: list(b[k].shape) for k in DESCRIPTORS + ("summary",) if k in batch_keys}, "cells": n,
        "batch_date": M.text(f, f["batch_date"])[0] if "batch_date" in f else None,
        "read_batch_fields": list(DESCRIPTORS) + ["summary"], "read_summary_fields": ["cycle"] + list(SERIES),
        "never_read_batch_fields": list(NEVER_BATCH), "never_read_summary_fields": list(NEVER_SUMMARY),
        "cycles_read": [FIRST, LAST], "cycle_numbers_peeked": PEEK,
        "python": platform.python_version(), "h5py": h5py.__version__, "numpy": np.__version__,
        "summary_keys": {}, "orientation": {}, "descriptor_classes": {}, "notes": []}
    cells, cyc = [], []
    for i in range(n):
        d = {}
        for k in DESCRIPTORS:
            if k in b:
                val, cls = M.text(f, b[k][i, 0])
            else:
                val, cls = None, "absent"
            d[k], d[k + "_class"] = val, cls
            structure["descriptor_classes"].setdefault(k, {}).setdefault(cls, 0)
            structure["descriptor_classes"][k][cls] += 1
        s = M.deref(f, b["summary"][i, 0])
        skeys = tuple(sorted(s.keys()))
        structure["summary_keys"].setdefault("|".join(skeys), 0)
        structure["summary_keys"]["|".join(skeys)] += 1
        c = M.prefix(s["cycle"], PEEK)
        o = M.orientation(s["cycle"])
        structure["orientation"].setdefault("cycle:" + o, 0)
        structure["orientation"]["cycle:" + o] += 1
        if len(c) > 1 and not np.all(np.diff(c) > 0):
            structure["notes"].append(f"{label}c{i}: cycle numbers not increasing in the first entries")
        pos = np.nonzero((c >= FIRST) & (c <= LAST))[0]
        vals = {}
        for k in SERIES:
            if k not in s:
                vals[k] = None
                continue
            o = M.orientation(s[k])
            structure["orientation"].setdefault(k + ":" + o, 0)
            structure["orientation"][k + ":" + o] += 1
            v = M.prefix(s[k], int(pos[-1]) + 1 if len(pos) else 0)
            vals[k] = v[pos[pos < len(v)]] if len(pos) else np.array([])
        got = set(int(round(x)) for x in c[pos])
        window = set(range(2, LAST + 1))
        cells.append({"label": label, "index": i, "key": f"{label}c{i}",
                      "policy_readable": d["policy_readable"], "policy": d["policy"],
                      "channel_id": d["channel_id"], "channel_class": d["channel_id_class"],
                      "barcode": d["barcode"], "barcode_class": d["barcode_class"],
                      "cycle_first": repr(float(c[0])) if len(c) else "",
                      "cycle_is_index": bool(len(c) and np.array_equal(c[:len(pos)], np.arange(1, len(pos) + 1))),
                      "window_cycles": len(window & got), "window_complete": window <= got})
        for j, p in enumerate(pos):
            row = {"label": label, "index": i, "cycle": int(round(c[p]))}
            for k in SERIES:
                v = vals[k]
                row[k] = repr(float(v[j])) if v is not None and j < len(v) else ""
            cyc.append(row)
    f.close()
    with open(os.path.join(out_dir, f"{label}_a0_cells.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(cells[0].keys()))
        w.writeheader()
        w.writerows(cells)
    with open(os.path.join(out_dir, f"{label}_a0_cycles.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["label", "index", "cycle"] + list(SERIES))
        w.writeheader()
        w.writerows(cyc)
    with open(os.path.join(out_dir, f"{label}_a0_structure.json"), "w") as fh:
        json.dump(structure, fh, indent=1)
    print(label, structure["batch_date"], "cells", n, "| batch keys", batch_keys, "| summary keys", structure["summary_keys"],
          "| orientation", structure["orientation"], "| classes", structure["descriptor_classes"],
          "| window complete", sum(r["window_complete"] for r in cells), "of", n, "| notes", structure["notes"][:5])


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
