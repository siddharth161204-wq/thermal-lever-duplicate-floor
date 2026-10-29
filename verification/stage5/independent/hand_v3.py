#!/usr/bin/env python3
"""Independent Stage 5 check, V3 (D77) by hand, in plain Python (csv, json, math; no numpy, nothing of src/).

The two pairs D77 names: the first KIT pair at 25 C in the order of results/stage2_warmer_pairs.csv (P041-1 and
P041-2), and the Sandia NMC pair at 25 C, 0-100 %, 0.5/1C, cells a and b, of results/stage3_pairs.csv.
Every step is written out; the result is compared with this checker's own out/pairs.csv to 1e-12 relative.

Output: out/v3_hand.json.
Usage: python3 -B verification/stage5/independent/hand_v3.py
"""
import csv, json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))


def rd(rel):
    return os.path.join(ROOT, rel)


steps = []


def say(s):
    steps.append(s)
    print(s)


# ---------------- KIT
pairs = list(csv.DictReader(open(rd("results/stage2_warmer_pairs.csv"))))
first25 = None
for i, r in enumerate(pairs):
    if r["age_temp"] == "25":
        first25 = (i, r)
        break
i, r = first25
say(f"KIT: the first 25 C pair is row {i + 1} of results/stage2_warmer_pairs.csv: set {r['parameter_id']}, protocol "
    f"{r['protocol']}, warm {r['warm']}, cool {r['cool']}")
cells = {c["cell"]: c for c in csv.DictReader(open(rd("results/stage2_cells.csv")))}
Tw, Tc = float(cells[r["warm"]]["T_C"]), float(cells[r["cool"]]["T_C"])
say(f"  T_C {r['warm']} = {Tw!r} C, T_C {r['cool']} = {Tc!r} C; the warm cell is the one with the higher T_C: "
    f"{'ok' if Tw > Tc else 'NOT ok'}")
dT_cells = Tw - Tc
dT = float(r["dT_K"])
say(f"  dT = {Tw!r} - {Tc!r} = {dT_cells!r} K; archived dT_K = {dT!r} K; difference {abs(dT - dT_cells)!r}")
Lw, Lc = float(cells[r["warm"]]["L90_efc"]), float(cells[r["cool"]]["L90_efc"])
say(f"  L90 warm = {Lw!r} EFC, L90 cool = {Lc!r} EFC (archived pair file: {r['L_warm']}, {r['L_cool']}); censored90: "
    f"{cells[r['warm']]['censored90']}, {cells[r['cool']]['censored90']}")
lev = json.load(open(rd("results/stage2_lever.json")))
w = [x for x in lev["lever_25_40"] if x["protocol"] == r["protocol"]][0]
c = [x for x in lev["lever_10_25"] if x["protocol"] == r["protocol"]][0]
ww, wc = w["ci95"][1] - w["ci95"][0], c["ci95"][1] - c["ci95"][0]
band = "warm" if ww < wc else "cold"
say(f"  25 C band: the 25 to 40 C interval [{w['ci95'][0]!r}, {w['ci95'][1]!r}] is {ww!r} wide, the 10 to 25 C "
    f"interval [{c['ci95'][0]!r}, {c['ci95'][1]!r}] is {wc!r} wide; the narrower is {band}")
s = w["s_per_K"] if band == "warm" else c["s_per_K"]
say(f"  lever: {('25 to 40 C' if band == 'warm' else '10 to 25 C')} {r['protocol']}, s = {s!r} per K")
p = -s * dT
d = math.log(Lw / Lc)
say(f"  p = -s dT = -({s!r}) x {dT!r} = {p!r}")
say(f"  d = ln(L_warm / L_cool) = ln({Lw!r} / {Lc!r}) = ln({Lw / Lc!r}) = {d!r}")
res = max(float(cells[r["warm"]]["res90_efc"]), float(cells[r["cool"]]["res90_efc"]))
tie = abs(Lw - Lc) < res
say(f"  tie: |L_warm - L_cool| = {abs(Lw - Lc)!r} against the larger res90_efc {res!r}: {'tie' if tie else 'untied'}")
kit = {"set": "P%03d" % int(r["parameter_id"]), "warm": r["warm"], "cool": r["cool"], "dT_K": dT, "band": band,
       "s_per_K": s, "p": p, "d": d, "tie": tie, "L_warm": Lw, "L_cool": Lc}

# ---------------- Sandia
sp = list(csv.DictReader(open(rd("results/stage3_pairs.csv"))))
literal_first = [(j, x) for j, x in enumerate(sp) if x["chemistry"] == "NMC" and x["ambient_C"] == "25.0"][0]
say(f"Sandia: the first NMC pair at 25 C in file order is row {literal_first[0] + 1}, group {literal_first[1]['group']} "
    f"(rise_a_K '{literal_first[1]['rise_a_K']}', rise_b_K '{literal_first[1]['rise_b_K']}'): it has no measure (a); "
    f"D77 names the 0-100 %, 0.5/1C pair of cells a and b, used here")
j, x = [(j, x) for j, x in enumerate(sp) if x["group"] == "NMC_25C_0-100_0.5/1C" and x["cell_a"].endswith("_a")
        and x["cell_b"].endswith("_b")][0]
say(f"  row {j + 1}: group {x['group']}, cells {x['cell_a']} and {x['cell_b']}")
ra, rb = float(x["rise_a_K"]), float(x["rise_b_K"])
La, Lb = float(x["L_a_efc"]), float(x["L_b_efc"])
if ra > rb:
    warm, cool, Lw, Lc, dT = x["cell_a"], x["cell_b"], La, Lb, ra - rb
else:
    warm, cool, Lw, Lc, dT = x["cell_b"], x["cell_a"], Lb, La, rb - ra
say(f"  measure (a): rise_a = {ra!r} K, rise_b = {rb!r} K; warm = {warm}; dT = {dT!r} K")
say(f"  lifetimes: L_a = {La!r}, L_b = {Lb!r} EFC (censored {x['censored_a']}, {x['censored_b']})")
lv = {(q["chemistry"], q["secant"]): q for q in csv.DictReader(open(rd("results/stage3_lever.csv"))) if q["series"] == "registered 1C"}
hi = lv[("NMC", "25 to 35 C")]
lo = lv[("NMC", "15 to 25 C")]
wh = float(hi["ci95_hi"]) - float(hi["ci95_lo"])
wl = float(lo["ci95_hi"]) - float(lo["ci95_lo"])
band = "warm" if wh < wl else "cold"
say(f"  25 C band for NMC: 25 to 35 C interval {wh!r} wide, 15 to 25 C interval {wl!r} wide; the narrower is {band}")
s = float(hi["s_per_K"]) if band == "warm" else float(lo["s_per_K"])
say(f"  lever: NMC {'25 to 35 C' if band == 'warm' else '15 to 25 C'} (registered 1C), s = {s!r} per K")
p = -s * dT
d = math.log(Lw / Lc)
say(f"  p = -s dT = -({s!r}) x {dT!r} = {p!r}")
say(f"  d = ln(L_warm / L_cool) = ln({Lw!r} / {Lc!r}) = {d!r}")
say(f"  tie flag (below_resolution, archived): {x['below_resolution']} (|dL| = {abs(La - Lb)!r}, resolution {x['resolution_efc']})")
san = {"set": x["group"], "warm": warm, "cool": cool, "dT_K": dT, "band": band, "s_per_K": s, "p": p, "d": d,
       "tie": x["below_resolution"] == "True", "L_warm": Lw, "L_cool": Lc}

# ---------------- against this checker's pair table
mine = list(csv.DictReader(open(os.path.join(HERE, "out", "pairs.csv"))))
cmp = {}
for name, h in (("KIT", kit), ("Sandia", san)):
    m = [y for y in mine if y["set"] == h["set"] and y["warm"] == h["warm"] and y["cool"] == h["cool"]]
    if len(m) != 1:
        raise SystemExit("pair not found in out/pairs.csv")
    m = m[0]
    worst = 0.0
    for k in ("dT_K", "s_per_K", "p", "d", "L_warm", "L_cool"):
        a, b = float(m[k]), h[k]
        worst = max(worst, abs(a - b) / max(1e-300, abs(b)))
    cmp[name] = {"max_rel_diff": worst, "band_equal": m["band"] == h["band"], "tie_equal": (m["tie"] == "True") == h["tie"],
                 "pass_1e-12": worst <= 1e-12 and m["band"] == h["band"]}
say(f"against out/pairs.csv: {cmp}")
json.dump({"check": "V3 (D77) by hand, plain Python", "kit_pair": kit, "sandia_pair": san, "against_own_pairs_csv": cmp,
           "literal_first_nmc_25C_pair": literal_first[1]["group"] + " (no measure (a); not scored)", "steps": steps},
          open(os.path.join(HERE, "out", "v3_hand.json"), "w"), indent=1)
