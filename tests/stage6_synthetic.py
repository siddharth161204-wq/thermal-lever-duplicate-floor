#!/usr/bin/env python3
"""Stage 6's code on synthetic batch files (DECISIONS.md D96, V4), before it reads any real lifetime.

Writes three synthetic batch files with a planted within-set lever, a planted between-policy lever and a planted
resistance effect (`tests/stage6_synthetic_files.py`), runs the readers on them as they run where the real files are,
builds the cell table and runs Parts B and C, and passes where each planted value lies inside its 95 % interval: the
temperature's lever with the resistance held (B2), the resistance's own slope (B2), and the between-policy lever (C1).
It also checks that the A0 reader opens none of the lifetime or capacity fields, by recording every name it reads.

Usage: python3 -I tests/stage6_synthetic.py [scratch directory]
"""
import json, math, os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC6 = os.path.join(ROOT, "src", "stage6")
sys.path.insert(0, SRC6)
sys.path.insert(0, HERE)
import h5py                                     # noqa: E402
import stage6_synthetic_files as SF             # noqa: E402
import a0_cells as A                            # noqa: E402
import a0_analysis as AN                        # noqa: E402
import s6_parts as P                            # noqa: E402

CONFIG = {"s1": {"date": "2030-01-01", "paper": "Severson", "authors_prefix": None, "drop": {}},
          "s2": {"date": "2030-02-01", "paper": "Severson", "authors_prefix": None, "drop": {}},
          "s3": {"date": "2030-03-01", "paper": "Attia", "authors_prefix": None, "drop": {}}}


def guard_a0_reader(path):
    """Run a0_read.main in this process with every h5py name lookup recorded."""
    import a0_read
    seen = []
    orig = h5py.Group.__getitem__

    def rec(self, name):
        if isinstance(name, str):
            seen.append(name.split("/")[-1])
        return orig(self, name)
    h5py.Group.__getitem__ = rec
    try:
        a0_read.main(path, "guard", tempfile.mkdtemp())
    finally:
        h5py.Group.__getitem__ = orig
    bad = sorted(set(seen) & set(a0_read.NEVER_BATCH + a0_read.NEVER_SUMMARY))
    return bad, sorted(set(seen))


def scenario(scratch, planted, seed):
    os.makedirs(scratch, exist_ok=True)
    data, ext = os.path.join(scratch, "files"), os.path.join(scratch, "extract")
    truth = SF.make(data, seed, planted)
    p = truth["planted"]
    for lab, cfg in CONFIG.items():
        f = os.path.join(data, f"{cfg['date']}_batchdata_synthetic.mat")
        for script in ("a0_read.py", "l_read.py"):
            subprocess.run([sys.executable, "-I", os.path.join(SRC6, script), f, lab, ext], check=True,
                           stdout=subprocess.DEVNULL)
    bad, names = guard_a0_reader(os.path.join(data, "2030-01-01_batchdata_synthetic.mat"))
    rows, sets = A.build(ext, CONFIG)
    life = {}
    import csv
    for lab, cfg in CONFIG.items():
        for r in csv.DictReader(open(os.path.join(ext, f"{lab}_life.csv"))):
            life[f"{cfg['date']}:{r['index']}"] = float(r["cycle_life_stored"])
    a0 = {"reliability": AN.reliability(rows), "channel": AN.channel(rows)}
    res = P.run(rows, life, a0)
    reg = res["registered"]["pooled"]
    return p, bad, reg, a0


def inside(x, v):
    return x["ci95"][0] <= v <= x["ci95"][1]


def main(scratch):
    p, bad, reg, a0 = scenario(os.path.join(scratch, "confounded"), None, 20261116)
    checks = {
        "A0 reader opened no lifetime or capacity field": (not bad, bad),
        "B2: temperature's lever with R held covers the planted s_w":
            (reg["B2"]["s_w_R_held_per_K"]["ci95"][0] <= p["s_w"] <= reg["B2"]["s_w_R_held_per_K"]["ci95"][1],
             reg["B2"]["s_w_R_held_per_K"]),
        "B2: resistance's slope covers the planted g":
            (reg["B2"]["lnL_per_lnR"]["ci95"][0] <= p["g"] <= reg["B2"]["lnL_per_lnR"]["ci95"][1], reg["B2"]["lnL_per_lnR"]),
        "C1: between-policy lever covers the planted s_b":
            (reg["C1"]["s_b_per_K"]["ci95"][0] <= p["s_b"] <= reg["C1"]["s_b_per_K"]["ci95"][1], reg["C1"]["s_b_per_K"]),
        "B4: residual sd covers the planted cell noise":
            (reg["B4"]["sd_cell"]["ci95"][0] <= p["sigma_c"] * 1.0 <= reg["B4"]["sd_cell"]["ci95"][1] * 1.2,
             reg["B4"]["sd_cell"]),
    }
    # a second case: no link between a cell's temperature and its resistance, so B1 alone recovers the lever
    p2, _, reg2, _ = scenario(os.path.join(scratch, "independent"), {"coupling_K_per_lnR": 0.0}, 20261118)
    checks["B1: s_w covers the planted within-set lever when T does not follow R"] = (
        inside(reg2["B1"]["s_w_per_K"], p2["s_w"]), reg2["B1"]["s_w_per_K"])
    checks["B1 pairs form covers it too"] = (inside(reg2["B1"]["pairs_form_s_per_K"], p2["s_w"]), reg2["B1"]["pairs_form_s_per_K"])
    # a third: each channel carries a temperature into every batch; A0.7 finds it and B5 recovers the lever
    p3, _, reg3, a03 = scenario(os.path.join(scratch, "channel"), {"coupling_K_per_lnR": 0.0, "channel_sd_K": 1.5,
                                                                    "cell_sd_K": 0.4}, 20261119)
    checks["A0.7: the planted channel temperature is found"] = (a03["channel"]["channel_carries_temperature"],
                                                                {k: v for k, v in a03["channel"].items() if k != "instrument"})
    checks["B5: the instrumented lever covers the planted lever"] = (
        "B5" in reg3 and inside(reg3["B5"]["s_iv_per_K"], p3["s_w"]), reg3.get("B5"))
    ok = all(v[0] for v in checks.values())
    for k, v in checks.items():
        print(("pass " if v[0] else "FAIL ") + k, json.dumps(v[1])[:300])
    print("B1 s_w (biased by the resistance it carries, as planted):", reg["B1"]["s_w_per_K"]["value"],
          "| sign test", reg["B3"], "| channel", {k: v for k, v in a0["channel"].items() if k != "instrument"})
    print("synthetic test pass:", ok)
    return ok


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1] if len(sys.argv) > 1 else tempfile.mkdtemp()) else 1)
