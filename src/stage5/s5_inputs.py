#!/usr/bin/env python3
"""Stage 5, Part A: the pairs read from the archives, with their levers, p and d (DECISIONS.md D71 to D73).

Reads only the archived files named in D71; every lever, floor, shape and median is checked against D56, D57 and the
Stage 4 files before it is used. Writes results/stage5_pairs.csv: one row for each of the 142 KIT and 86 Sandia pairs,
scored or not, with the reason.

Usage: python3 src/stage5/s5_inputs.py
"""
import csv, json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
ROOT = os.path.dirname(SRC)
SRC4 = os.path.join(SRC, "stage4")
if SRC4 not in sys.path:
    sys.path.insert(0, SRC4)
import common as CM4                    # noqa: E402  src/stage4/common.py, unchanged
import experiment as EX4                # noqa: E402  src/stage4/experiment.py, unchanged

SEED = 20261029

# D71: every file Stage 5 reads (V1 hashes each against tag v1.3.0)
FILES = [
    "results/stage2_warmer_pairs.csv",
    "results/stage2_cells.csv",
    "results/stage2_lever.json",
    "results/stage3_pairs.csv",
    "results/stage3_pairs_oriented_measure_a.csv",
    "results/stage3_lever.csv",
    "inputs/floor.json",
    "results/stage4_experiment.json",
    "results/stage4_paper_numbers.json",
    "src/stage4/common.py",
    "src/stage4/experiment.py",
]

KIT_CHEM = "NMC (LG INR18650HG2)"
WARM_SECANT_KIT, COLD_SECANT_KIT = "25 to 40 C", "10 to 25 C"
WARM_SECANT_SNL, COLD_SECANT_SNL = "25 to 35 C", "15 to 25 C"
NCA_15C_1C_GROUP = "NCA_15C_0-100_0.5/1C"      # D74's bracket: cell b's sensor did not follow the cell (D55)


def path(rel):
    return os.path.join(ROOT, rel)


# ---------------------------------------------------------------------------------------------------------------- levers
def kit_levers():
    """Per protocol, the warm (25 to 40 C) and cold (10 to 25 C) lever with its interval, as D56 (3) lists them."""
    b = CM4.read_kit()                  # stops if any value differs from D56 at its rounding
    out = {}
    for band in ("warm", "cold"):
        for l in b[band]:
            out.setdefault(l["protocol"], {})[band] = {"s": l["s"], "lo": l["lo"], "hi": l["hi"],
                                                       "secant": l["secant"], "file": l["file"], "key": l["key"]}
    if len(out) != 12 or any(set(v) != {"warm", "cold"} for v in out.values()):
        raise SystemExit("KIT levers: twelve protocols with both secants expected")
    return out


def sandia_levers():
    """The registered 1C series of results/stage3_lever.csv by chemistry and secant, determinable or not."""
    CM4.read_sandia()                   # checks the four levers D56 (4) lists
    out = {}
    for r in csv.DictReader(open(path("results/stage3_lever.csv"))):
        if r["series"] != "registered 1C":
            continue
        det = r["determinable"] == "True"
        out[(r["chemistry"], r["secant"])] = {
            "determinable": det,
            "s": float(r["s_per_K"]) if det else None,
            "lo": float(r["ci95_lo"]) if det else None,
            "hi": float(r["ci95_hi"]) if det else None,
            "secant": r["secant"], "file": "results/stage3_lever.csv", "key": "series registered 1C"}
    if len(out) != 6:
        raise SystemExit("Sandia levers: six chemistry-secant rows expected")
    nca = out[("NCA", WARM_SECANT_SNL)]
    if not nca["determinable"] or round(100.0 * nca["s"], 2) != -2.04:     # D57 names it at -2.04 % per K
        raise SystemExit("Sandia NCA 25 to 35 C lever differs from D57")
    if out[("LFP", WARM_SECANT_SNL)]["determinable"]:
        raise SystemExit("LFP 25 to 35 C is not determinable in the archive (D52)")
    return out


def width(l):
    return l["hi"] - l["lo"]


def kit_band(temp, lev):
    """D72: 40 C warm, 10 C cold, 0 C none, 25 C the narrower interval of the protocol's two."""
    if temp == 40:
        return "warm", "40 C set: 25 to 40 C lever"
    if temp == 10:
        return "cold", "10 C set: 10 to 25 C lever"
    if temp == 25:
        ww, cw = width(lev["warm"]), width(lev["cold"])
        b = "warm" if ww < cw else "cold"
        return b, f"25 C set: narrower interval ({100 * ww:.4f} against {100 * cw:.4f} % per K wide)"
    if temp == 0:
        return None, "0 C set: no measured lever reaches 0 C"
    raise SystemExit(f"unexpected KIT temperature {temp}")


def sandia_band(chem, ambient, lev):
    """D72: 15 C the 15 to 25 C lever, 35 C the 25 to 35 C lever, 25 C the narrower determinable interval."""
    warm, cold = lev[(chem, WARM_SECANT_SNL)], lev[(chem, COLD_SECANT_SNL)]
    if ambient == 15.0:
        return ("cold", cold, "15 C set: 15 to 25 C lever") if cold["determinable"] else (None, None, "15 C set: 15 to 25 C lever not determinable")
    if ambient == 35.0:
        return ("warm", warm, "35 C set: 25 to 35 C lever") if warm["determinable"] else (None, None, "35 C set: 25 to 35 C lever not determinable")
    if ambient == 25.0:
        if warm["determinable"] and cold["determinable"]:
            ww, cw = width(warm), width(cold)
            if ww < cw:
                return "warm", warm, f"25 C set: narrower interval ({100 * ww:.4f} against {100 * cw:.4f} % per K wide)"
            return "cold", cold, f"25 C set: narrower interval ({100 * cw:.4f} against {100 * ww:.4f} % per K wide)"
        if cold["determinable"]:
            return "cold", cold, "25 C set: the 25 to 35 C lever is not determinable, the 15 to 25 C one is used"
        if warm["determinable"]:
            return "warm", warm, "25 C set: the 15 to 25 C lever is not determinable, the 25 to 35 C one is used"
        return None, None, "25 C set: no determinable lever"
    raise SystemExit(f"unexpected Sandia ambient {ambient}")


# ----------------------------------------------------------------------------------------------------------------- pairs
def _b(x):
    return str(x).strip() == "True"


def kit_pairs():
    lev = kit_levers()
    cells = {r["cell"]: r for r in csv.DictReader(open(path("results/stage2_cells.csv")))}
    n_cens_cells = sum(_b(r["censored90"]) for r in cells.values())
    rows = []
    for r in csv.DictReader(open(path("results/stage2_warmer_pairs.csv"))):
        temp = int(r["age_temp"])
        w, c = r["warm"], r["cool"]
        cw, cc = cells[w], cells[c]
        if _b(cw["censored90"]) or _b(cc["censored90"]):
            raise SystemExit("a KIT pair in the archive holds a censored cell")
        dT = float(r["dT_K"])
        Lw, Lc = float(r["L_warm"]), float(r["L_cool"])
        if not dT > 0:
            raise SystemExit("KIT pair with dT not above zero")
        band, reason = kit_band(temp, lev[r["protocol"]])
        res_w, res_c = float(cw["res90_efc"]), float(cc["res90_efc"])
        row = {"dataset": "KIT", "set": "P%03d" % int(r["parameter_id"]), "temp_C": float(temp),
               "protocol": r["protocol"], "chemistry": KIT_CHEM, "band": band or "none", "band_reason": reason,
               "warm": w, "cool": c, "dT_K": dT, "L_warm_efc": Lw, "L_cool_efc": Lc,
               "resolution_efc": max(res_w, res_c), "tie": abs(Lw - Lc) < max(res_w, res_c),
               "d": math.log(Lw / Lc)}
        if band:
            l = lev[r["protocol"]][band]
            row.update({"lever": f"KIT {r['protocol']} {l['secant']}", "s_per_K": l["s"], "s_lo_per_K": l["lo"],
                        "s_hi_per_K": l["hi"], "p": -l["s"] * dT, "scored": True})
        else:
            row.update({"lever": "", "s_per_K": None, "s_lo_per_K": None, "s_hi_per_K": None, "p": None, "scored": False})
        rows.append(row)
    n_sets_total = len({r["parameter_id"] for r in csv.DictReader(open(path("results/stage2_cells.csv")))})
    info = {"pairs_in_file": len(rows), "censored_cells_at_90pct": n_cens_cells,
            "pairs_not_in_file_for_censoring": 3 * n_sets_total - len(rows), "sets": n_sets_total}
    return rows, info


def sandia_pairs():
    lev = sandia_levers()
    rows = []
    for r in csv.DictReader(open(path("results/stage3_pairs.csv"))):
        chem, amb = r["chemistry"], float(r["ambient_C"])
        La, Lb = float(r["L_a_efc"]), float(r["L_b_efc"])
        has_a = r["rise_a_K"] not in ("", "nan") and r["rise_b_K"] not in ("", "nan")
        cens = _b(r["censored_a"]) or _b(r["censored_b"])
        tie = _b(r["below_resolution"])
        base = {"dataset": "Sandia", "set": r["group"], "temp_C": amb,
                "protocol": f"{r['window']} {float(r['c_charge']):g}/{float(r['c_discharge']):g}C", "chemistry": chem,
                "resolution_efc": float(r["resolution_efc"]), "tie": tie}
        if not has_a:
            rows.append({**base, "band": "none", "band_reason": "no measure (a) for both cells (not oriented: warm and cool are cells a and b)", "warm": r["cell_a"],
                         "cool": r["cell_b"], "dT_K": None, "L_warm_efc": La, "L_cool_efc": Lb, "d": None, "lever": "",
                         "s_per_K": None, "s_lo_per_K": None, "s_hi_per_K": None, "p": None, "scored": False})
            continue
        ra, rb = float(r["rise_a_K"]), float(r["rise_b_K"])
        if ra == rb:
            raise SystemExit("a Sandia pair with equal measure (a) cannot be oriented")
        if ra > rb:
            w, c, Lw, Lc, dT = r["cell_a"], r["cell_b"], La, Lb, ra - rb
        else:
            w, c, Lw, Lc, dT = r["cell_b"], r["cell_a"], Lb, La, rb - ra
        band, l, reason = sandia_band(chem, amb, lev)
        row = {**base, "band": band or "none", "band_reason": reason, "warm": w, "cool": c, "dT_K": dT,
               "L_warm_efc": Lw, "L_cool_efc": Lc, "d": math.log(Lw / Lc)}
        if band and not cens:
            row.update({"lever": f"Sandia {chem} {l['secant']}", "s_per_K": l["s"], "s_lo_per_K": l["lo"],
                        "s_hi_per_K": l["hi"], "p": -l["s"] * dT, "scored": True})
        else:
            if cens:
                row["band_reason"] += "; censored lifetime"
            row.update({"lever": "", "s_per_K": None, "s_lo_per_K": None, "s_hi_per_K": None, "p": None, "scored": False})
        # measure (b), Stage 3's registered secondary (D74's unscored line)
        qa, qb = r["rest_a_C"], r["rest_b_C"]
        if qa not in ("", "nan") and qb not in ("", "nan") and float(qa) != float(qb):
            qa, qb = float(qa), float(qb)
            if qa > qb:
                wb, cb, Lwb, Lcb, dTb = r["cell_a"], r["cell_b"], La, Lb, qa - qb
            else:
                wb, cb, Lwb, Lcb, dTb = r["cell_b"], r["cell_a"], Lb, La, qb - qa
            row.update({"warm_b": wb, "cool_b": cb, "dT_b_K": dTb, "d_b": math.log(Lwb / Lcb),
                        "p_b": (-row["s_per_K"] * dTb) if row["scored"] else None})
        rows.append(row)
    return rows


# ------------------------------------------------------------------------------------------------- noise and the medians
def noise():
    """The floors (D56 (5)) and Stage 4's registered shape (D61), checked against results/stage4_experiment.json."""
    floors = CM4.read_floors()
    ex = json.load(open(path("results/stage4_experiment.json")))["shapes"]
    out = {}
    for name, f in floors.items():
        ln = EX4.lognormal_shape(f["mean"], f["p95"])
        if ln["sd_e"] != ex[name]["lognormal"]["sd_e"] or ln["sigma"] != ex[name]["lognormal"]["sigma"] or ln["mu"] != ex[name]["lognormal"]["mu"]:
            raise SystemExit(f"the {name} shape differs from results/stage4_experiment.json")
        sd_norm = f["mean"] * math.sqrt(math.pi / 2.0)
        if abs(sd_norm - ex[name]["normal"]["sd_e"]) > 1e-15:
            raise SystemExit(f"the {name} normal shape differs from results/stage4_experiment.json")
        out[name] = {"mean": f["mean"], "p95": f["p95"], "median": f["median"], "pairs": f["pairs"], "file": f["file"],
                     "mu": ln["mu"], "sigma": ln["sigma"], "sd_e": ln["sd_e"], "sd_e_normal_shape": sd_norm}
    if out["Geslin"]["sd_e"] != 0.04454298465298393 or out["KIT"]["sd_e"] != 0.10412996462328801:
        raise SystemExit("sd(e) differs from D71 (5)")
    return out


def medians():
    """D57's band medians: labels from results/stage4_experiment.json, values from results/stage2_lever.json."""
    named = json.load(open(path("results/stage4_experiment.json")))["named"]
    pn = json.load(open(path("results/stage4_paper_numbers.json")))["median_levers_pct_per_K"]
    lev = kit_levers()
    out = {}
    for band, key, want in (("warm", "median_warm", 1.3029), ("cold", "median_cold", -3.5800)):
        label = named[key]
        if not label.startswith("KIT "):
            raise SystemExit("band median is not a KIT lever")
        prot = label[4:]
        s = lev[prot][band]["s"]
        if pn[band][0] != label or abs(100.0 * s - pn[band][1]) > 1e-12 or round(100.0 * s, 4) != want:
            raise SystemExit(f"band median {band} differs from D57 or results/stage4_paper_numbers.json")
        out[band] = {"label": label, "protocol": prot, "secant": lev[prot][band]["secant"], "s_per_K": s,
                     "abs_s_per_K": abs(s), "file": "results/stage4_experiment.json (named) and results/stage2_lever.json"}
    return out


# ------------------------------------------------------------------------------------------------------------------ main
COLUMNS = ["dataset", "set", "temp_C", "protocol", "chemistry", "band", "band_reason", "scored", "warm", "cool", "dT_K",
           "L_warm_efc", "L_cool_efc", "d", "lever", "s_per_K", "s_lo_per_K", "s_hi_per_K", "p", "resolution_efc", "tie",
           "warm_b", "cool_b", "dT_b_K", "d_b", "p_b"]


def _cell(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "True" if v else "False"
    if isinstance(v, float):
        return repr(v)
    return str(v)


def load_all():
    kit, kit_info = kit_pairs()
    snl = sandia_pairs()
    return {"pairs": kit + snl, "kit_info": kit_info, "noise": noise(), "medians": medians()}


def write_pairs(pairs, out):
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(COLUMNS)
        for r in pairs:
            w.writerow([_cell(r.get(k)) for k in COLUMNS])


if __name__ == "__main__":
    d = load_all()
    write_pairs(d["pairs"], path("results/stage5_pairs.csv"))
    P = d["pairs"]
    for ds in ("KIT", "Sandia"):
        s = [r for r in P if r["dataset"] == ds]
        sc = [r for r in s if r["scored"]]
        print(ds, "pairs", len(s), "scored", len(sc), "sets", len({r["set"] for r in sc}),
              "warm", sum(r["band"] == "warm" for r in sc), "cold", sum(r["band"] == "cold" for r in sc),
              "ties among scored", sum(r["tie"] for r in sc))
    print(json.dumps(d["kit_info"]))
