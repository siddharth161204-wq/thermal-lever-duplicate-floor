#!/usr/bin/env python3
"""Stage 6, the lifetimes of one batch file (DECISIONS.md D94), run only after Part A0 was committed and pushed.

For each cell of the struct array `batch` it reads the field cycle_life and the summary's QDischarge, and writes the
stored cycle life, the number of summary entries, the first cycle whose discharge capacity is below 0.88 Ah (the
authors' MATLAB loader's label, `LoadData.m`), the first such cycle followed by a lower one (the rule of the authors'
`cell_analysis.m` that set the stored field), and the last discharge capacity (the MATLAB loader's 0.885 Ah test).
Runs where the files are, on Python 3.10 with numpy and h5py only.

Usage: python3 -I l_read.py <file.mat> <label> <output directory>
Writes <label>_life.csv in the output directory.
"""
import csv, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import h5py                 # noqa: E402
import numpy as np          # noqa: E402
import matv73 as M          # noqa: E402

THRESHOLD = 0.88


def vector(ds):
    a = np.asarray(ds[()], dtype=float)
    return a.ravel()


def main(path, label, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    f = h5py.File(path, "r")
    b = f["batch"]
    n = int(b["summary"].shape[0])
    rows = []
    for i in range(n):
        cl = M.deref(f, b["cycle_life"][i, 0])
        clv = vector(cl)
        stored = float(clv[0]) if clv.size == 1 and not int(cl.attrs.get("MATLAB_empty", 0)) else float("nan")
        s = M.deref(f, b["summary"][i, 0])
        q = vector(s["QDischarge"])
        below = np.nonzero(q < THRESHOLD)[0]
        first = int(below[0]) + 1 if len(below) else ""
        confirmed = ""
        for k in below:                                  # cell_analysis.m: the next discharge is lower still
            if k + 1 < len(q) and q[k + 1] < q[k]:
                confirmed = int(k) + 1
                break
        rows.append({"label": label, "index": i, "cycle_life_stored": repr(stored), "summary_entries": len(q),
                     "first_below_0p88": first, "first_below_0p88_confirmed": confirmed,
                     "last_QDischarge": repr(float(q[-1])) if len(q) else ""})
    f.close()
    with open(os.path.join(out_dir, f"{label}_life.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(label, n, "cells written")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
