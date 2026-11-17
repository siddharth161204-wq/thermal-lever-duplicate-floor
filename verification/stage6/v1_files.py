#!/usr/bin/env python3
"""V1 of Stage 6 (DECISIONS.md D96 (1)): the files, and the cells, policies and cycle lives of the authors' own code.

(1) Every struct's size and SHA-256 in inputs/matr/MANIFEST.md equal those computed where the files are
(inputs/matr/files.json), each size equals the end of file its HDF5 superblock records, and the extracts were read from
those files. (2) The authors' Python loader (`Load Data.ipynb`, rdbraatz/data-driven-prediction-of-battery-cycle-life-
before-capacity-degradation at 1ef13d2), ported to the committed extracts: its deletions and carry-over give 41, 43 and
40 cells, which are the cells that enter here, with the same policies; its cycle lives, the stored field plus add_len
(662, 981, 1060, 208, 482) for b1c0 to b1c4, equal the registered ones; add_len equals the continuations' numbers of
summary entries; and the MATLAB loader's label (`LoadData.m`) agrees for every one of the 124 cells. (3) The validation
batch: the stored cycle lives of each protocol's five cells equal, as a set, the authors' `Validation/final_results.csv`
(petersen.../battery-parameter-spaces at ec9508f; the repository has no licence file, so the file is not copied here and
is read from the path given).

Usage: python3 -I verification/stage6/v1_files.py [path to final_results.csv]
"""
import csv, hashlib, json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
M = os.path.join(ROOT, "inputs", "matr")
FINAL_SHA = "94359f14e0066e1464e0e06d0ed3c4ceaccc8ea779b2fd492509cea02a22e5bb"
ADD_LEN = [662, 981, 1060, 208, 482]


def manifest():
    rows = {}
    for line in open(os.path.join(M, "MANIFEST.md")):
        m = re.match(r"\| (\d{4}-\d\d-\d\d)_batchdata_updated_struct_errorcorrect\.mat \| (\S+) \| ([\d,]+) \| ([0-9a-f]{64}) \|", line)
        if m:
            rows[m.group(1)] = {"address": m.group(2), "bytes": int(m.group(3).replace(",", "")), "sha256": m.group(4)}
    return rows


def main(final_path):
    out, ok = {"check": "V1 (D96 (1))"}, True
    # (1) the files
    fj = json.load(open(os.path.join(M, "files.json")))["files"]
    man = manifest()
    dates = {"b1": "2017-05-12", "b2": "2017-06-30", "b3": "2018-04-12", "b9": "2019-01-24"}
    files = []
    for lab, d in dates.items():
        f, mrow = fj[lab], man[d]
        st = json.load(open(os.path.join(M, "a0", f"{lab}_a0_structure.json")))
        good = (f["bytes"] == mrow["bytes"] and f["sha256"] == mrow["sha256"] and f["bytes"] == f["hdf5_end_of_file"]
                and st["file"] == f["local_name"])
        ok &= good
        files.append({"batch": d, "bytes": f["bytes"], "sha256": f["sha256"], "manifest_agrees": good,
                      "extract_read_from": st["file"]})
    out["files"] = files
    # (2) the authors' Python loader, ported
    ext = {lab: {int(r["index"]): r for r in csv.DictReader(open(os.path.join(M, "a0", f"{lab}_a0_cells.csv")))}
           for lab in ("b1", "b2", "b3")}
    life = {lab: {int(r["index"]): r for r in csv.DictReader(open(os.path.join(M, "life", f"{lab}_life.csv")))}
            for lab in ("b1", "b2", "b3", "b9")}
    batch = {}
    for lab in ("b1", "b2", "b3"):
        for i, r in ext[lab].items():
            batch[f"{lab}c{i}"] = {"policy": r["policy_readable"], "cycle_life": float(life[lab][i]["cycle_life_stored"]),
                                   "entries": int(life[lab][i]["summary_entries"])}
    for k in ("b1c8", "b1c10", "b1c12", "b1c13", "b1c22"):
        del batch[k]
    cont = {}
    for i, (b2k, b1k) in enumerate(zip(["b2c7", "b2c8", "b2c9", "b2c15", "b2c16"], ["b1c0", "b1c1", "b1c2", "b1c3", "b1c4"])):
        batch[b1k]["cycle_life"] = batch[b1k]["cycle_life"] + ADD_LEN[i]
        cont[b1k] = batch[b2k]["entries"]
        del batch[b2k]
    for k in ("b3c37", "b3c2", "b3c23", "b3c32", "b3c42", "b3c43"):
        del batch[k]
    n = {lab: sum(1 for k in batch if k.startswith(lab)) for lab in ("b1", "b2", "b3")}
    cells = {r["authors_key"]: r for r in csv.DictReader(open(os.path.join(ROOT, "results", "stage6_cells.csv")))
             if r["enters"] == "True" and r["authors_key"]}
    lt = {r["authors_key"]: r for r in csv.DictReader(open(os.path.join(M, "stage6_lifetimes.csv"))) if r["authors_key"]}
    same_cells = set(batch) == set(cells)
    same_policy = all(batch[k]["policy"] == cells[k]["policy"] for k in batch if k in cells)
    same_life = all(batch[k]["cycle_life"] == float(cells[k]["cycle_life"]) == float(lt[k]["cycle_life"]) for k in batch
                    if k in cells)
    add_ok = [cont[f"b1c{i}"] == ADD_LEN[i] for i in range(5)]
    matlab = sum(1 for k in batch if float(lt[k]["matlab_label"]) != batch[k]["cycle_life"])
    good = n == {"b1": 41, "b2": 43, "b3": 40} and same_cells and same_policy and same_life and all(add_ok) and matlab == 0
    ok &= good
    out["authors_python_loader"] = {"cells": n, "same_cells_as_registered": same_cells, "same_policies": same_policy,
                                    "same_cycle_lives": same_life, "continuation_entries": [cont[f"b1c{i}"] for i in range(5)],
                                    "add_len": ADD_LEN, "add_len_equals_continuations": add_ok,
                                    "matlab_label_differs_for": matlab, "pass": bool(good)}
    # (3) the validation batch against the authors' final_results.csv
    if final_path and os.path.exists(final_path):
        raw = open(final_path, "rb").read()
        sha = hashlib.sha256(raw).hexdigest()
        rows = list(csv.reader(raw.decode("utf-8-sig").splitlines()))[1:]
        authors = {tuple(float(x) for x in r[:3]): sorted(float(x) for x in r[3:8]) for r in rows if r}
        b9 = {int(r["index"]): r for r in csv.DictReader(open(os.path.join(M, "a0", "b9_a0_cells.csv")))}
        ours = {}
        for i, r in b9.items():
            cc = tuple(float(x) for x in r["policy_readable"].split("-")[:3])
            ours.setdefault(cc, []).append(float(life["b9"][i]["cycle_life_stored"]))
        ours = {k: sorted(v) for k, v in ours.items()}
        same = ours == authors
        ok &= same and sha == FINAL_SHA
        out["validation_batch"] = {"final_results_sha256": sha, "expected_sha256": FINAL_SHA, "protocols": len(authors),
                                   "cells": sum(len(v) for v in authors.values()), "same_cycle_lives_by_protocol": same}
    else:
        out["validation_batch"] = {"note": "final_results.csv not given; not checked"}
        ok = False
    out["pass"] = bool(ok)
    json.dump(out, open(os.path.join(HERE, "v1.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))
    print("V1 pass:", ok)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
