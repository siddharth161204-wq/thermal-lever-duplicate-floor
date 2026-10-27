#!/usr/bin/env python3
"""Stage 4, the archived inputs (DECISIONS.md D56 and D57).

Every lever, interval and floor is read at full precision from the archived Stage 1 to 3 files and checked against
the values D56 records at the rounding shown there; the bands and the levers the brief names are built as D57 fixes
them. Nothing here is computed from a Stage 4 run.
"""
import csv, json, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
ROOT = os.path.dirname(SRC)
sys.path.insert(0, SRC)

RATES = [0.25, 0.5, 1.0, 2.0]                 # C/4, C/2, 1C and 2C (D58)
RATE_LABEL = {0.25: "C/4", 0.5: "C/2", 1.0: "1C", 2.0: "2C"}
AMBIENTS_C = [35.0, 25.0, 10.0]
H_LIST = [5.0, 10.0, 25.0]
BAND_OF_AMBIENT = {35.0: "warm", 25.0: "warm", 10.0: "cold"}
BAND_SETPOINTS_C = {"warm": (25.0, 40.0), "cold": (10.0, 25.0)}   # the KIT secants' set points, for the label of D57

R_35_ARCHIVED = 0.01871279570954531           # inputs/lgm50_resistance.json, registered (D56)

# D56 (3): KIT levers in % per K with the 95 % interval, 25 to 40 C and 10 to 25 C, at the rounding of DECISIONS.md
KIT_D56 = {
    "0-100 0.33/0.33C": ((+1.3029, +0.4122, +2.1937), (-4.0180, -4.9813, -3.0547)),
    "0-100 0.33/1.00C": ((-0.6236, -1.5140, +0.2667), (-3.1041, -4.0560, -2.1523)),
    "0-100 1.00/1.00C": ((-2.6850, -3.5876, -1.7825), (-3.3410, -4.2788, -2.4033)),
    "0-100 1.67/1.00C": ((-3.4548, -4.4542, -2.4553), (-3.5800, -4.6130, -2.5471)),
    "10-100 0.33/0.33C": ((+1.3598, +0.4592, +2.2605), (-2.1876, -3.1272, -1.2481)),
    "10-100 0.33/1.00C": ((+0.6724, -0.2212, +1.5659), (-3.6478, -4.5758, -2.7198)),
    "10-100 1.00/1.00C": ((-0.4222, -1.3097, +0.4653), (-4.1691, -5.0866, -3.2517)),
    "10-100 1.67/1.00C": ((-1.9148, -2.8032, -1.0265), (-11.5906, -12.5150, -10.6661)),
    "10-90 0.33/0.33C": ((+1.9097, +1.0043, +2.8152), (-1.7564, -2.7006, -0.8121)),
    "10-90 0.33/1.00C": ((+1.0898, +0.1874, +1.9922), (-2.2142, -3.1580, -1.2704)),
    "10-90 1.00/1.00C": ((-0.1935, -1.0698, +0.6828), (-3.8750, -4.8438, -2.9062)),
    "10-90 1.67/1.00C": ((-0.9952, -1.8798, -0.1106), (-5.2373, -6.1621, -4.3126)),
}
# D56 (4): Sandia, registered 1C series
SANDIA_D56 = {
    ("NMC", "25 to 35 C"): (-2.4432, -4.3314, -0.4448),
    ("NMC", "15 to 25 C"): (-5.6227, -8.4354, -3.1004),
    ("NCA", "15 to 25 C"): (-1.3151, -5.2199, +0.7423),
    ("LFP", "15 to 25 C"): (+3.3718, +2.9491, +3.8270),
}
# D56 (5): floors in %, mean, median, 95th percentile, pairs
FLOORS_D56 = {"Geslin": (2.9910, 1.7457, 8.7176, 45), "KIT": (7.3323, 4.0902, 20.4751, 47)}


def _pct4(x):
    return round(100.0 * float(x), 4)


def _check(name, got, want):
    if any(abs(g - w) > 1e-9 for g, w in zip(got, want)):
        raise SystemExit(f"archived value differs from DECISIONS.md D56: {name}: {got} against {want}")


def read_kit():
    path = os.path.join(ROOT, "results", "stage2_lever.json")
    j = json.load(open(path))
    out = {"warm": [], "cold": []}
    for key, band, secant, idx in (("lever_25_40", "warm", "25 to 40 C", 0), ("lever_10_25", "cold", "10 to 25 C", 1)):
        seen = []
        for r in j[key]:
            if not r["determinable"]:
                raise SystemExit("a KIT lever is not determinable: " + r["protocol"])
            want = KIT_D56[r["protocol"]][idx]
            _check(f"KIT {r['protocol']} {secant}", (_pct4(r["s_per_K"]), _pct4(r["ci95"][0]), _pct4(r["ci95"][1])), want)
            out[band].append({"band": band, "source": "KIT", "label": "KIT " + r["protocol"], "protocol": r["protocol"],
                              "chemistry": "NMC (LG INR18650HG2)", "secant": secant, "s": float(r["s_per_K"]),
                              "lo": float(r["ci95"][0]), "hi": float(r["ci95"][1]),
                              "file": "results/stage2_lever.json", "key": key})
            seen.append(r["protocol"])
        if sorted(seen) != sorted(KIT_D56):
            raise SystemExit("the KIT protocols differ from D56: " + key)
    return out


def read_sandia():
    path = os.path.join(ROOT, "results", "stage3_lever.csv")
    rows = [r for r in csv.DictReader(open(path)) if r["series"] == "registered 1C"]
    got = {}
    for r in rows:
        k = (r["chemistry"], r["secant"])
        if k in SANDIA_D56:
            if r["determinable"] != "True":
                raise SystemExit("a Sandia lever is not determinable: " + str(k))
            s, lo, hi = float(r["s_per_K"]), float(r["ci95_lo"]), float(r["ci95_hi"])
            _check(f"Sandia {k}", (_pct4(s), _pct4(lo), _pct4(hi)), SANDIA_D56[k])
            got[k] = {"source": "Sandia", "chemistry": r["chemistry"], "secant": r["secant"], "s": s, "lo": lo, "hi": hi,
                      "file": "results/stage3_lever.csv", "key": "series registered 1C"}
    if set(got) != set(SANDIA_D56):
        raise SystemExit("Sandia levers missing: " + str(set(SANDIA_D56) - set(got)))
    return got


def read_floors():
    g = json.load(open(os.path.join(ROOT, "inputs", "floor.json")))
    k = json.load(open(os.path.join(ROOT, "results", "stage2_lever.json")))["floor"]
    out = {"Geslin": {"mean": float(g["mean"]), "median": float(g["median"]), "p95": float(g["p95"]), "pairs": int(g["n_pairs"]),
                      "file": "inputs/floor.json"},
           "KIT": {"mean": float(k["mean"]), "median": float(k["median"]), "p95": float(k["p95"]), "pairs": int(k["pairs"]),
                   "file": "results/stage2_lever.json (floor)"}}
    for name, (m, med, q, n) in FLOORS_D56.items():
        f = out[name]
        _check("floor " + name, (_pct4(f["mean"]), _pct4(f["median"]), _pct4(f["p95"]), f["pairs"]), (m, med, q, n))
    return out


def read_r35():
    j = json.load(open(os.path.join(ROOT, "inputs", "lgm50_resistance.json")))
    r = float(j["registered"]["R_DC_ohm"])
    if r != R_35_ARCHIVED:
        raise SystemExit("the archived 35 C resistance differs from D56")
    return r


def bands():
    """D57: the warm band (13 levers) and the cold band (15 levers)."""
    kit = read_kit()
    san = read_sandia()
    warm = list(kit["warm"])
    x = san[("NMC", "25 to 35 C")]
    warm.append({"band": "warm", "label": "Sandia NMC 25 to 35 C", "protocol": "0-100 0.5/1C", **x,
                 "chemistry": "NMC (LG 18650HG2)"})
    cold = list(kit["cold"])
    for chem, name in (("NMC", "NMC (LG 18650HG2)"), ("NCA", "NCA (Panasonic NCR18650B)"), ("LFP", "LFP (A123 APR18650M1A)")):
        x = san[(chem, "15 to 25 C")]
        cold.append({"band": "cold", "label": f"Sandia {chem} 15 to 25 C", "protocol": "0-100 0.5/1C", **x, "chemistry": name})
    for b in (warm, cold):
        for i, lev in enumerate(b):
            lev["index"] = i
    if len(warm) != 13 or len(cold) != 15:
        raise SystemExit("band sizes differ from D57")
    return {"warm": warm, "cold": cold}


def median_lever(levers):
    """The lever at the median of |s| (D57); the bands have odd counts, so it is one measured lever."""
    a = np.array([abs(l["s"]) for l in levers])
    if len(a) % 2 != 1:
        raise SystemExit("median_lever needs an odd count")
    order = np.argsort(a, kind="stable")
    return levers[int(order[len(a) // 2])]


def named_levers(b=None):
    """D57: the median (central) lever of each band and the three levers of C1, with the conditional reading."""
    b = b or bands()
    warm_kit = [l for l in b["warm"] if l["source"] == "KIT"]
    cold_kit = [l for l in b["cold"] if l["source"] == "KIT"]
    one = [l for l in warm_kit if l["protocol"].endswith("1.00/1.00C")]
    dis1 = [l for l in warm_kit if l["protocol"].endswith("/1.00C")]
    if len(one) != 3 or len(dis1) != 9:
        raise SystemExit("protocol counts differ from D57")
    a = np.sort(np.abs([l["s"] for l in cold_kit]))
    if np.any([l["s"] >= 0 for l in cold_kit]):
        raise SystemExit("a KIT 10 to 25 C lever is not negative")
    cold_kit_central = -float(0.5 * (a[5] + a[6]))          # numpy's median of twelve, the mean of the sixth and seventh
    if abs(cold_kit_central + float(np.median(a))) > 1e-15:
        raise SystemExit("central cold KIT value")
    sandia_1c = [l for l in b["warm"] if l["source"] == "Sandia"][0]
    return {
        "median_warm": median_lever(b["warm"]),
        "median_cold": median_lever(b["cold"]),
        "c1": [
            {"name": "central KIT value at 1C", "reading": "the three 1.00/1.00C protocols, median of |s|",
             "conditional": False, **_lever_copy(median_lever(one))},
            {"name": "Sandia 1C value", "reading": "NMC 25 to 35 C, registered 1C series", "conditional": False,
             **_lever_copy(sandia_1c)},
            {"name": "central cold KIT value", "reading": "median of |s| over the twelve KIT 10 to 25 C levers (mean of the sixth and seventh)",
             "conditional": False, "label": "KIT 10 to 25 C median", "source": "KIT", "s": cold_kit_central, "lo": None, "hi": None,
             "file": "results/stage2_lever.json", "key": "lever_10_25"},
            {"name": "KIT value at a 1C discharge", "reading": "the nine protocols with a 1C discharge, median of |s|",
             "conditional": True, **_lever_copy(median_lever(dis1))},
        ],
    }


def _lever_copy(l):
    return {k: l[k] for k in ("label", "source", "s", "lo", "hi", "file", "key")}


def load_all():
    b = bands()
    return {"bands": b, "named": named_levers(b), "floors": read_floors(), "R35": read_r35()}


if __name__ == "__main__":
    d = load_all()
    print(json.dumps({"warm": [(l["label"], l["s"]) for l in d["bands"]["warm"]],
                      "cold": [(l["label"], l["s"]) for l in d["bands"]["cold"]],
                      "median_warm": d["named"]["median_warm"]["label"], "median_cold": d["named"]["median_cold"]["label"],
                      "c1": [(x["name"], x["label"], x["s"]) for x in d["named"]["c1"]], "floors": d["floors"], "R35": d["R35"]},
                     indent=1))
