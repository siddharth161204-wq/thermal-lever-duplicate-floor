#!/usr/bin/env python3
"""Stage 6, the lifetimes table (DECISIONS.md D94), built after Part A0 was committed and pushed.

From the outputs of `l_read.py` (inputs/matr/life/<label>_life.csv): each entering cell's registered cycle life L, the
stored cycle_life, as the authors' Python loader and validation code take it, and for a carried cell the stored value of
its 2017-05-12 entry plus the length the loader adds for its continuation (662, 981, 1060, 208 and 482 cycles, the
loader's add_len, checked in V1 against the continuations' numbers of summary entries). Beside it, not used: the MATLAB
loader's label (`LoadData.m`): if the last discharge capacity is below 0.88 Ah, the first cycle below it, otherwise one
more than the number of cycles, over the combined cycles of a carried cell.

Usage: python3 -I src/stage6/l_table.py
Writes inputs/matr/stage6_lifetimes.csv and adds the column cycle_life to results/stage6_cells.csv.
"""
import csv, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s6_common as C                   # noqa: E402
import a0_cells as A                    # noqa: E402

ADD_LEN = {0: 662, 1: 981, 2: 1060, 3: 208, 4: 482}          # the Python loader's add_len, b1c0 to b1c4


def read(label):
    p = C.path(f"inputs/matr/life/{label}_life.csv")
    return {int(r["index"]): r for r in csv.DictReader(open(p))}


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def main():
    rows, sets = A.build(C.path("inputs/matr/a0"))
    life = {lab: read(lab) for lab in A.BATCHES}
    cont = {to[1]: i for i, to in A.BATCHES["b2"]["carried"].items()}          # b1 index -> b2 continuation index
    out = []
    for r in rows:
        if not r["enters"]:
            continue
        e = life[r["label"]][r["index"]]
        stored = num(e["cycle_life_stored"])
        n = int(e["summary_entries"])
        first = int(e["first_below_0p88"]) if e["first_below_0p88"] else None
        rec = {"cell": r["cell"], "batch": r["batch"], "index": r["index"], "authors_key": r["authors_key"],
               "cycle_life_stored": stored, "summary_entries": n, "continuation_entries": "",
               "first_below_0p88_confirmed": e["first_below_0p88_confirmed"], "last_QDischarge": e["last_QDischarge"]}
        # LoadData.m: if the last discharge capacity is below 0.88 Ah, the first cycle below it; else the number of
        # cycles plus one; over the combined cycles of a carried cell
        if r["label"] == "b1" and r["index"] in cont:
            c = life["b2"][cont[r["index"]]]
            m = int(c["summary_entries"])
            rec["continuation_entries"] = m
            rec["cycle_life"] = stored + ADD_LEN[r["index"]]
            if num(c["last_QDischarge"]) < 0.88:
                lab = first if first is not None else n + int(c["first_below_0p88"])
            else:
                lab = n + m + 1
        else:
            rec["cycle_life"] = stored
            lab = first if (num(e["last_QDischarge"]) < 0.88 and first is not None) else n + 1
        rec["matlab_label"] = lab
        out.append(rec)
    cols = ["cell", "batch", "index", "authors_key", "cycle_life", "cycle_life_stored", "summary_entries",
            "continuation_entries", "matlab_label", "first_below_0p88_confirmed", "last_QDischarge"]
    with open(C.path("inputs/matr/stage6_lifetimes.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for x in out:
            w.writerow({k: (repr(v) if isinstance(v, float) else v) for k, v in x.items()})
    L = {x["cell"]: x["cycle_life"] for x in out}
    for r in rows:
        r["cycle_life"] = L.get(r["cell"], "")
    A.write(rows, sets, C.path("results"), extra=("cycle_life",))
    bad = [x["cell"] for x in out if not (isinstance(x["cycle_life"], float) and math.isfinite(x["cycle_life"])
                                          and x["cycle_life"] > 0)]
    print(len(out), "entering cells;", len(bad), "without a finite cycle life:", bad,
          "; MATLAB label differs for", sum(1 for x in out if x["matlab_label"] != x["cycle_life"]))


if __name__ == "__main__":
    main()
