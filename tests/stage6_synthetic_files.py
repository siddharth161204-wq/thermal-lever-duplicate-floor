#!/usr/bin/env python3
"""Synthetic batch files in the layout of the fast-charging structs (DECISIONS.md D84, D96), for Stage 6's tests and V4.

Writes MATLAB v7.3-style HDF5 files: a root group '#refs#' holding every referenced object, a group 'batch' whose
fields (Vdlin, barcode, channel_id, cycle_life, cycles, policy, policy_readable, summary) are n x 1 datasets of
object references, and a char dataset 'batch_date'. Each summary is a group of column vectors stored with HDF5 shape
(1, k): cycle (1 to k), QDischarge, QCharge, IR, Tmax, Tavg, Tmin and chargetime. A cell's lifetime follows a planted
model, ln L = c_policy + b_batch - s_w (T - T_policy) - s_b T_policy + g ln R + e, so that the readers and Parts B and C
can be checked against known values; its discharge capacity falls below 0.88 Ah first at cycle L.

Usage: python3 -I tests/stage6_synthetic_files.py <output directory> [seed]
"""
import json, math, os, sys

import h5py
import numpy as np

REF = h5py.ref_dtype


def put(refs, name, data, cls, **attrs):
    ds = refs.create_dataset(name, data=data)
    ds.attrs["MATLAB_class"] = np.bytes_(cls)
    for k, v in attrs.items():
        ds.attrs[k] = v
    return ds.ref


def char(refs, name, s):
    a = np.array([ord(c) for c in s], dtype=np.uint16).reshape(len(s), 1)
    return put(refs, name, a, "char", MATLAB_int_decode=np.int32(2))


def column(values):
    return np.asarray(values, dtype=float).reshape(1, -1)


def write_batch(path, date, cells):
    """cells: dicts with policy, policy_readable, channel, barcode, life, T (per cycle), IR (per cycle), swing."""
    with h5py.File(path, "w", userblock_size=512) as f:
        refs = f.create_group("#refs#")
        f.create_group("#subsystem#")
        bd = np.array([ord(c) for c in date], dtype=np.uint16).reshape(len(date), 1)
        d = f.create_dataset("batch_date", data=bd)
        d.attrs["MATLAB_class"] = np.bytes_("char")
        g = f.create_group("batch")
        g.attrs["MATLAB_class"] = np.bytes_("struct")
        fields = {k: [] for k in ("Vdlin", "barcode", "channel_id", "cycle_life", "cycles", "policy", "policy_readable",
                                  "summary")}
        for i, c in enumerate(cells):
            k = len(c["T"])
            fields["policy"].append(char(refs, f"p{i}", c["policy"]))
            fields["policy_readable"].append(char(refs, f"pr{i}", c["policy_readable"]))
            fields["channel_id"].append(char(refs, f"ch{i}", str(c["channel"])))
            fields["barcode"].append(char(refs, f"bc{i}", c["barcode"]))
            fields["cycle_life"].append(put(refs, f"cl{i}", np.array([[float(c["life"])]]), "double"))
            fields["Vdlin"].append(put(refs, f"v{i}", column(np.linspace(3.5, 2.0, 10)), "double"))
            cy = refs.create_group(f"cy{i}")
            cy.attrs["MATLAB_class"] = np.bytes_("struct")
            cy.create_dataset("I", data=np.zeros((1, 1)))
            fields["cycles"].append(cy.ref)
            s = refs.create_group(f"s{i}")
            s.attrs["MATLAB_class"] = np.bytes_("struct")
            n = np.arange(1, k + 1)
            q = 1.07 - 0.19 * (n / c["life"]) ** 3                      # crosses 0.88 between cycles L - 1 and L
            q[n >= c["life"]] = np.minimum(q[n >= c["life"]], 0.8799)
            T = np.asarray(c["T"], dtype=float)
            for name, v in (("cycle", n), ("QDischarge", q), ("QCharge", q + 0.001), ("IR", c["IR"]),
                            ("Tmax", T + c["swing"] / 2), ("Tavg", T), ("Tmin", T - c["swing"] / 2),
                            ("chargetime", np.full(k, 10.0))):
                ds = s.create_dataset(name, data=column(v))
                ds.attrs["MATLAB_class"] = np.bytes_("double")
            fields["summary"].append(s.ref)
        for k, v in fields.items():
            ds = g.create_dataset(k, data=np.array(v, dtype=REF).reshape(-1, 1), dtype=REF)
            ds.attrs["MATLAB_class"] = np.bytes_("cell" if k != "summary" else "struct")
    with open(path, "r+b") as fh:                                        # MATLAB's 116-byte text header
        fh.write(b"MATLAB 7.3 MAT-file, Platform: synthetic, Created on: test HDF5 schema 1.00 .".ljust(116, b" "))


def make(out_dir, seed=20261116, planted=None):
    """Three synthetic batches; returns the planted values and the cells' true temperatures and lifetimes."""
    p = {"s_w": 0.05, "s_b": 0.12, "g": -0.40, "sigma_c": 0.06, "noise_T": 0.40, "coupling_K_per_lnR": 8.0,
         "channel_sd_K": 0.0, "cell_sd_K": 1.2}
    if planted:
        p.update(planted)
    rng = np.random.default_rng(seed)
    os.makedirs(out_dir, exist_ok=True)
    truth = {"planted": p, "cells": []}
    chan = p["channel_sd_K"] * rng.standard_normal(60)                  # a temperature each channel carries
    for bi, (date, n_pol, per) in enumerate((("2030-01-01", 12, 3), ("2030-02-01", 10, 4), ("2030-03-01", 9, 5))):
        cells, ch = [], 1
        for j in range(n_pol):
            T_pol = 31.0 + 6.0 * rng.random()                          # the policy's mean cell temperature
            c_pol = 6.6 + 0.10 * rng.standard_normal()                  # a policy effect apart from heat
            for m in range(per):
                lnR = math.log(0.016) + 0.05 * rng.standard_normal()
                T_cell = (T_pol + chan[ch] + p["cell_sd_K"] * rng.standard_normal() +
                          p["coupling_K_per_lnR"] * (lnR - math.log(0.016)))
                lnL = (c_pol + 0.05 * bi - p["s_w"] * (T_cell - T_pol) - p["s_b"] * (T_pol - 34.0) + p["g"] *
                       (lnR - math.log(0.016)) + p["sigma_c"] * rng.standard_normal())
                life = int(max(120, round(math.exp(lnL))))
                k = life + 5
                T = T_cell + p["noise_T"] * rng.standard_normal(k)
                IR = np.exp(lnR + 0.02 * rng.standard_normal(k))
                IR[rng.random(k) < 0.02] = 0.0                          # an invalid reading now and then
                cells.append({"policy": f"POL{j}", "policy_readable": f"{3 + j % 5}C(50%)-{4 + j % 3}C",
                              "channel": ch, "barcode": f"EL{bi}{ch:04d}", "life": life, "T": T, "IR": IR,
                              "swing": 6.0 + rng.random()})
                truth["cells"].append({"batch": date, "index": len(cells) - 1, "policy": f"POL{j}", "T_cell": T_cell,
                                       "T_policy": T_pol, "lnR": lnR, "life": life, "channel": ch})
                ch += 1
        write_batch(os.path.join(out_dir, f"{date}_batchdata_synthetic.mat"), date, cells)
    with open(os.path.join(out_dir, "truth.json"), "w") as fh:
        json.dump(truth, fh, indent=1)
    return truth


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    t = make(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 20261116)
    print("wrote", len(t["cells"]), "cells in 3 batches to", sys.argv[1])
