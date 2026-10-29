#!/usr/bin/env python3
"""V3 of Stage 5 (DECISIONS.md D77): two pairs' p and d by hand.

Plain Python only (no numpy and nothing of src/stage5): the first KIT pair at 25 C in the order of
results/stage2_warmer_pairs.csv (P041-1 and P041-2), and the first Sandia NMC pair at 25 C in the order of
results/stage3_pairs.csv (0-100 %, 0.5/1C, cells a and b). From the archived values: the orientation, dT, the band
(the narrower interval at 25 C), the lever, p = -s dT and d = ln(L_warm / L_cool), each written out; compared with
results/stage5_pairs.csv, passing at agreement to 1e-12 relative.

Usage: python3 verification/stage5/v3_hand.py
"""
import csv, json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def kit_case():
    row = next(r for r in csv.DictReader(open(os.path.join(ROOT, "results", "stage2_warmer_pairs.csv"))) if r["age_temp"] == "25")
    lev = json.load(open(os.path.join(ROOT, "results", "stage2_lever.json")))
    warm = next(r for r in lev["lever_25_40"] if r["protocol"] == row["protocol"])
    cold = next(r for r in lev["lever_10_25"] if r["protocol"] == row["protocol"])
    ww = warm["ci95"][1] - warm["ci95"][0]
    cw = cold["ci95"][1] - cold["ci95"][0]
    band, l = ("warm", warm) if ww < cw else ("cold", cold)
    dT = float(row["dT_K"])
    Lw, Lc = float(row["L_warm"]), float(row["L_cool"])
    s = l["s_per_K"]
    p = -s * dT
    d = math.log(Lw) - math.log(Lc)
    steps = [
        f"pair: set P{int(row['parameter_id']):03d} at {row['age_temp']} C, protocol {row['protocol']}; warm {row['warm']}, cool {row['cool']} (the archive's orientation by T, D29, D32)",
        f"dT = {dT!r} K (dT_K of results/stage2_warmer_pairs.csv)",
        f"interval widths: 25 to 40 C {warm['ci95'][1]!r} - {warm['ci95'][0]!r} = {ww!r}; 10 to 25 C {cold['ci95'][1]!r} - {cold['ci95'][0]!r} = {cw!r}; band {band}",
        f"s = {s!r} per K",
        f"p = -s dT = -({s!r}) x {dT!r} = {p!r}",
        f"d = ln {Lw!r} - ln {Lc!r} = {math.log(Lw)!r} - {math.log(Lc)!r} = {d!r}",
    ]
    return {"dataset": "KIT", "warm": row["warm"], "cool": row["cool"], "band": band, "dT_K": dT, "s_per_K": s, "p": p, "d": d,
            "steps": steps}


def sandia_case():
    # D77 names the pair: NMC, 25 C, 0-100 %, 0.5/1C, cells a and b, the first NMC pair at 25 C that is scored (the
    # pairs of the 0.5/0.5C groups before it in the file have no measure (a))
    row = next(r for r in csv.DictReader(open(os.path.join(ROOT, "results", "stage3_pairs.csv")))
               if r["group"] == "NMC_25C_0-100_0.5/1C" and r["cell_a"].endswith("_a") and r["cell_b"].endswith("_b"))
    lev = {r["secant"]: r for r in csv.DictReader(open(os.path.join(ROOT, "results", "stage3_lever.csv")))
           if r["series"] == "registered 1C" and r["chemistry"] == "NMC"}
    ra, rb = float(row["rise_a_K"]), float(row["rise_b_K"])
    if ra > rb:
        w, c, Lw, Lc, dT = row["cell_a"], row["cell_b"], float(row["L_a_efc"]), float(row["L_b_efc"]), ra - rb
    else:
        w, c, Lw, Lc, dT = row["cell_b"], row["cell_a"], float(row["L_b_efc"]), float(row["L_a_efc"]), rb - ra
    warm, cold = lev["25 to 35 C"], lev["15 to 25 C"]
    ww = float(warm["ci95_hi"]) - float(warm["ci95_lo"])
    cw = float(cold["ci95_hi"]) - float(cold["ci95_lo"])
    band, l = ("warm", warm) if ww < cw else ("cold", cold)
    s = float(l["s_per_K"])
    p = -s * dT
    d = math.log(Lw) - math.log(Lc)
    steps = [
        f"pair: group {row['group']}, cells {row['cell_a'][-1]} and {row['cell_b'][-1]}",
        f"measure (a): {ra!r} K and {rb!r} K, so the warm cell is {w[-1]} and dT = {max(ra, rb)!r} - {min(ra, rb)!r} = {dT!r} K",
        f"interval widths: 25 to 35 C {float(warm['ci95_hi'])!r} - {float(warm['ci95_lo'])!r} = {ww!r}; 15 to 25 C {float(cold['ci95_hi'])!r} - {float(cold['ci95_lo'])!r} = {cw!r}; band {band}",
        f"s = {s!r} per K (NMC, {l['secant']}, registered 1C)",
        f"p = -s dT = -({s!r}) x {dT!r} = {p!r}",
        f"d = ln {Lw!r} - ln {Lc!r} = {math.log(Lw)!r} - {math.log(Lc)!r} = {d!r}",
    ]
    return {"dataset": "Sandia", "warm": w, "cool": c, "band": band, "dT_K": dT, "s_per_K": s, "p": p, "d": d, "steps": steps}


def main():
    P = list(csv.DictReader(open(os.path.join(ROOT, "results", "stage5_pairs.csv"))))
    out = {"check": "V3 (D77): two pairs by hand", "cases": []}
    ok = True
    for case in (kit_case(), sandia_case()):
        row = next(r for r in P if r["dataset"] == case["dataset"] and r["warm"] == case["warm"] and r["cool"] == case["cool"])
        diffs = {k: rel(case[k], float(row[k])) for k in ("dT_K", "s_per_K", "p", "d")}
        same_band = case["band"] == row["band"]
        good = max(diffs.values()) <= 1e-12 and same_band and row["scored"] == "True"
        ok &= good
        case.update({"pipeline": {k: float(row[k]) for k in ("dT_K", "s_per_K", "p", "d")}, "pipeline_band": row["band"],
                     "rel_diff": diffs, "pass": good})
        out["cases"].append(case)
        print("\n".join(case["steps"]))
        print("relative differences:", diffs, "band", same_band, "pass", good, "\n")
    out["pass"] = bool(ok)
    json.dump(out, open(os.path.join(HERE, "v3.json"), "w"), indent=1)
    print("V3 pass:", ok)


if __name__ == "__main__":
    main()
