#!/usr/bin/env python3
"""V1 of Stage 5 (DECISIONS.md D77): the Part A inputs are the archived ones.

(1) The SHA-256 of every file of D71 equals that of the same path at tag v1.3.0 (main at 40445cf), with its hash at the
release that first archived it reported beside. (2) Every lifetime, temperature difference, lever and interval in
results/stage5_pairs.csv equals the archived value exactly; KIT dT_K equals T_C(warm) - T_C(cool) of
results/stage2_cells.csv within 1e-12 K, L_warm and L_cool equal the cells' L90_efc, and no scored pair holds a
censored cell; the Sandia orientation, dT and d equal warm, cool, dT_K and obs_ln_ratio of
results/stage3_pairs_oriented_measure_a.csv within 1e-12.

Usage: python3 verification/stage5/v1_inputs.py
"""
import csv, hashlib, json, math, os, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
FILES = ["results/stage2_warmer_pairs.csv", "results/stage2_cells.csv", "results/stage2_lever.json",
         "results/stage3_pairs.csv", "results/stage3_pairs_oriented_measure_a.csv", "results/stage3_lever.csv",
         "inputs/floor.json", "results/stage4_experiment.json", "results/stage4_paper_numbers.json",
         "src/stage4/common.py", "src/stage4/experiment.py"]
TAGS = ["v1.0.0", "v1.1.0", "v1.2.0", "v1.3.0"]


def sha(b):
    return hashlib.sha256(b).hexdigest()


def at(tag, rel):
    r = subprocess.run(["git", "-C", ROOT, "show", f"{tag}:{rel}"], capture_output=True)
    return r.stdout if r.returncode == 0 else None


def f(x):
    return float(x) if x not in ("", None) else None


def main():
    out = {"check": "V1 (D77)", "files": [], "values": {}}
    ok = True
    for rel in FILES:
        now = sha(open(os.path.join(ROOT, rel), "rb").read())
        v13 = at("v1.3.0", rel)
        first = next((t for t in TAGS if at(t, rel) is not None), None)
        h13 = sha(v13) if v13 is not None else None
        hfirst = sha(at(first, rel)) if first else None
        same = now == h13
        ok &= same
        out["files"].append({"file": rel, "sha256": now, "sha256_v1.3.0": h13, "equal_to_v1.3.0": same,
                             "first_release": first, "sha256_first_release": hfirst, "equal_to_first_release": now == hfirst})
    P = list(csv.DictReader(open(os.path.join(ROOT, "results", "stage5_pairs.csv"))))
    # KIT
    wp = list(csv.DictReader(open(os.path.join(ROOT, "results", "stage2_warmer_pairs.csv"))))
    cells = {r["cell"]: r for r in csv.DictReader(open(os.path.join(ROOT, "results", "stage2_cells.csv")))}
    lev = json.load(open(os.path.join(ROOT, "results", "stage2_lever.json")))
    L = {("warm", r["protocol"]): r for r in lev["lever_25_40"]}
    L.update({("cold", r["protocol"]): r for r in lev["lever_10_25"]})
    kp = [r for r in P if r["dataset"] == "KIT"]
    bad = []
    max_dT_vs_cells = 0.0
    if len(kp) != len(wp):
        bad.append("KIT pair count")
    for a, b in zip(kp, wp):
        if (a["warm"], a["cool"]) != (b["warm"], b["cool"]) or f(a["dT_K"]) != f(b["dT_K"]) or \
                f(a["L_warm_efc"]) != f(b["L_warm"]) or f(a["L_cool_efc"]) != f(b["L_cool"]):
            bad.append(("KIT pair", a["warm"], a["cool"]))
        cw, cc = cells[a["warm"]], cells[a["cool"]]
        max_dT_vs_cells = max(max_dT_vs_cells, abs(f(a["dT_K"]) - (f(cw["T_C"]) - f(cc["T_C"]))))
        if f(a["L_warm_efc"]) != f(cw["L90_efc"]) or f(a["L_cool_efc"]) != f(cc["L90_efc"]):
            bad.append(("KIT L90", a["warm"], a["cool"]))
        if a["scored"] == "True":
            if cw["censored90"] == "True" or cc["censored90"] == "True":
                bad.append(("KIT censored", a["warm"], a["cool"]))
            l = L[(a["band"], a["protocol"])]
            if f(a["s_per_K"]) != l["s_per_K"] or f(a["s_lo_per_K"]) != l["ci95"][0] or f(a["s_hi_per_K"]) != l["ci95"][1]:
                bad.append(("KIT lever", a["warm"], a["cool"]))
    ok_dT = max_dT_vs_cells <= 1e-12
    out["values"]["KIT"] = {"pairs": len(kp), "mismatches": [str(x) for x in bad], "max_abs_dT_vs_cells_K": max_dT_vs_cells,
                            "dT_vs_cells_within_1e-12": ok_dT}
    ok &= not bad and ok_dT
    # Sandia
    sp = {(r["group"], r["cell_a"], r["cell_b"]): r for r in csv.DictReader(open(os.path.join(ROOT, "results", "stage3_pairs.csv")))}
    orient = {(r["warm"], r["cool"]): r for r in csv.DictReader(open(os.path.join(ROOT, "results", "stage3_pairs_oriented_measure_a.csv")))}
    sl = {(r["chemistry"], r["secant"]): r for r in csv.DictReader(open(os.path.join(ROOT, "results", "stage3_lever.csv")))
          if r["series"] == "registered 1C"}
    snl = [r for r in P if r["dataset"] == "Sandia"]
    bad = []
    max_dT, max_d = 0.0, 0.0
    oriented = 0
    for a in snl:
        if a["dT_K"] == "":
            continue
        oriented += 1
        key = next(k for k in sp if k[0] == a["set"] and {k[1], k[2]} == {a["warm"], a["cool"]})
        r = sp[key]
        La, Lb, ra, rb = f(r["L_a_efc"]), f(r["L_b_efc"]), f(r["rise_a_K"]), f(r["rise_b_K"])
        if a["warm"] == r["cell_a"]:
            want = (La, Lb, ra - rb)
        else:
            want = (Lb, La, rb - ra)
        if (f(a["L_warm_efc"]), f(a["L_cool_efc"]), f(a["dT_K"])) != want:
            bad.append(("Sandia values", a["warm"], a["cool"]))
        o = orient.get((a["warm"], a["cool"]))
        if o is None:
            bad.append(("Sandia orientation", a["warm"], a["cool"]))
            continue
        max_dT = max(max_dT, abs(f(a["dT_K"]) - f(o["dT_K"])))
        max_d = max(max_d, abs(f(a["d"]) - f(o["obs_ln_ratio"])))
        if a["scored"] == "True":
            chem = a["chemistry"]
            sec = "25 to 35 C" if a["band"] == "warm" else "15 to 25 C"
            l = sl[(chem, sec)]
            if f(a["s_per_K"]) != f(l["s_per_K"]) or f(a["s_lo_per_K"]) != f(l["ci95_lo"]) or f(a["s_hi_per_K"]) != f(l["ci95_hi"]):
                bad.append(("Sandia lever", a["warm"], a["cool"]))
            if r["censored_a"] == "True" or r["censored_b"] == "True":
                bad.append(("Sandia censored", a["warm"], a["cool"]))
    ok_s = not bad and max_dT <= 1e-12 and max_d <= 1e-12 and oriented == len(orient)
    out["values"]["Sandia"] = {"pairs": len(snl), "oriented_by_measure_a": oriented, "in_stage3_oriented_file": len(orient),
                               "mismatches": [str(x) for x in bad], "max_abs_dT_vs_stage3_K": max_dT,
                               "max_abs_d_vs_stage3_obs_ln_ratio": max_d}
    ok &= ok_s
    out["pass"] = bool(ok)
    json.dump(out, open(os.path.join(HERE, "v1.json"), "w"), indent=1)
    for x in out["files"]:
        print(x["file"], x["equal_to_v1.3.0"], x["first_release"], x["equal_to_first_release"])
    print(json.dumps(out["values"], indent=1))
    print("V1 pass:", out["pass"])


if __name__ == "__main__":
    main()
