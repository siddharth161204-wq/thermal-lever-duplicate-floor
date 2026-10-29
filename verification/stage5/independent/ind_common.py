#!/usr/bin/env python3
"""Independent Stage 5 check: shared readers, the pair table and the noise shapes.

Written by the independent checker from DECISIONS.md D71 to D78, the brief and the archived inputs that D71 names.
It imports nothing from src/ (no Stage 1 to 5 module) and reads no Stage 5 pipeline file.

Conventions (D73): s in natural-log units per kelvin, positive when the warmer cell reaches the threshold sooner;
dT = T_warm - T_cool > 0 by the orientation the measured temperatures set; d = ln(L_warm / L_cool); p = -s dT.
"""
import csv, hashlib, json, math, os, subprocess

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = os.path.join(HERE, "out")

BASE_COMMIT = "da613a1488838fd8a73032882d4f9fc5d6a0750d"
ARCHIVE_TAG = "v1.3.0"          # main at 40445cf, the archive Stage 5 starts from (D77)

# the files D71 names (the inputs this check may read)
FILES_D71 = [
    "results/stage2_warmer_pairs.csv",
    "results/stage2_cells.csv",
    "results/stage2_lever.json",
    "results/stage3_pairs.csv",
    "results/stage3_pairs_oriented_measure_a.csv",
    "results/stage3_lever.csv",
    "inputs/floor.json",
    "results/stage4_experiment.json",
    "results/stage4_paper_numbers.json",
]

# Values written in the task and in D56, D57 and D71, used only to stop if an archived value differs.
SD_E_TASK = {"Geslin": 0.04454298465298393, "KIT": 0.10412996462328801}
S_MED_TASK = {"warm": 0.013029423647732621, "cold": -0.03580047101563345}
D56_KIT = {   # % per K at 4 decimals: (s, lo, hi); 25 to 40 C then 10 to 25 C
    "warm": {
        "0-100 0.33/0.33C": (1.3029, 0.4122, 2.1937), "0-100 0.33/1.00C": (-0.6236, -1.5140, 0.2667),
        "0-100 1.00/1.00C": (-2.6850, -3.5876, -1.7825), "0-100 1.67/1.00C": (-3.4548, -4.4542, -2.4553),
        "10-100 0.33/0.33C": (1.3598, 0.4592, 2.2605), "10-100 0.33/1.00C": (0.6724, -0.2212, 1.5659),
        "10-100 1.00/1.00C": (-0.4222, -1.3097, 0.4653), "10-100 1.67/1.00C": (-1.9148, -2.8032, -1.0265),
        "10-90 0.33/0.33C": (1.9097, 1.0043, 2.8152), "10-90 0.33/1.00C": (1.0898, 0.1874, 1.9922),
        "10-90 1.00/1.00C": (-0.1935, -1.0698, 0.6828), "10-90 1.67/1.00C": (-0.9952, -1.8798, -0.1106)},
    "cold": {
        "0-100 0.33/0.33C": (-4.0180, -4.9813, -3.0547), "0-100 0.33/1.00C": (-3.1041, -4.0560, -2.1523),
        "0-100 1.00/1.00C": (-3.3410, -4.2788, -2.4033), "0-100 1.67/1.00C": (-3.5800, -4.6130, -2.5471),
        "10-100 0.33/0.33C": (-2.1876, -3.1272, -1.2481), "10-100 0.33/1.00C": (-3.6478, -4.5758, -2.7198),
        "10-100 1.00/1.00C": (-4.1691, -5.0866, -3.2517), "10-100 1.67/1.00C": (-11.5906, -12.5150, -10.6661),
        "10-90 0.33/0.33C": (-1.7564, -2.7006, -0.8121), "10-90 0.33/1.00C": (-2.2142, -3.1580, -1.2704),
        "10-90 1.00/1.00C": (-3.8750, -4.8438, -2.9062), "10-90 1.67/1.00C": (-5.2373, -6.1621, -4.3126)},
}
D56_SANDIA = {("NMC", "25 to 35 C"): (-2.4432, -4.3314, -0.4448), ("NMC", "15 to 25 C"): (-5.6227, -8.4354, -3.1004),
              ("NCA", "15 to 25 C"): (-1.3151, -5.2199, 0.7423), ("LFP", "15 to 25 C"): (3.3718, 2.9491, 3.8270)}
D56_FLOORS = {"Geslin": (2.9910, 1.7457, 8.7176, 45), "KIT": (7.3323, 4.0902, 20.4751, 47)}

Z95 = float(stats.norm.ppf(0.95))
SAN_BAND_SECANT = {"warm": "25 to 35 C", "cold": "15 to 25 C"}
KIT_BAND_LEVER = {"warm": "lever_25_40", "cold": "lever_10_25"}


def path(rel):
    return os.path.join(ROOT, rel)


def sha256_file(fn):
    h = hashlib.sha256()
    with open(fn, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_at(rev, rel):
    """SHA-256 of a file as stored at a release that precedes the base commit (None if absent there)."""
    try:
        blob = subprocess.run(["git", "-C", ROOT, "show", f"{rev}:{rel}"], capture_output=True, check=True).stdout
    except subprocess.CalledProcessError:
        return None
    return hashlib.sha256(blob).hexdigest()


def stop(msg):
    raise SystemExit("STOP: " + msg)


def check4(name, got, want):
    """got in natural-log units per K (or fraction); want in % at 4 decimals."""
    if round(100.0 * got, 4) != round(want, 4):
        stop(f"{name}: archived {100.0 * got!r} % does not round to {want}")


# ------------------------------------------------------------------ readers
def read_kit_levers():
    j = json.load(open(path("results/stage2_lever.json")))
    out = {}
    for band, key in KIT_BAND_LEVER.items():
        out[band] = {}
        for r in j[key]:
            s, lo, hi = float(r["s_per_K"]), float(r["ci95"][0]), float(r["ci95"][1])
            out[band][r["protocol"]] = {"determinable": bool(r["determinable"]), "s": s, "lo": lo, "hi": hi,
                                        "secant": "25 to 40 C" if band == "warm" else "10 to 25 C"}
            want = D56_KIT[band][r["protocol"]]
            check4(f"KIT {band} {r['protocol']} s", s, want[0])
            check4(f"KIT {band} {r['protocol']} lo", lo, want[1])
            check4(f"KIT {band} {r['protocol']} hi", hi, want[2])
        if len(out[band]) != 12:
            stop("KIT lever count")
    return out


def read_sandia_levers():
    out = {}
    for r in csv.DictReader(open(path("results/stage3_lever.csv"))):
        if r["series"] != "registered 1C":
            continue
        det = r["determinable"] == "True"
        rec = {"determinable": det, "s": float(r["s_per_K"]) if det else None,
               "lo": float(r["ci95_lo"]) if det else None, "hi": float(r["ci95_hi"]) if det else None}
        out[(r["chemistry"], r["secant"])] = rec
        if (r["chemistry"], r["secant"]) in D56_SANDIA:
            want = D56_SANDIA[(r["chemistry"], r["secant"])]
            check4(f"Sandia {r['chemistry']} {r['secant']} s", rec["s"], want[0])
            check4(f"Sandia {r['chemistry']} {r['secant']} lo", rec["lo"], want[1])
            check4(f"Sandia {r['chemistry']} {r['secant']} hi", rec["hi"], want[2])
    if out[("LFP", "25 to 35 C")]["determinable"]:
        stop("LFP 25 to 35 C should not be determinable (D52)")
    return out


def read_floors():
    g = json.load(open(path("inputs/floor.json")))
    k = json.load(open(path("results/stage2_lever.json")))["floor"]
    out = {"Geslin": {"mean": float(g["mean"]), "median": float(g["median"]), "p95": float(g["p95"]), "pairs": int(g["n_pairs"])},
           "KIT": {"mean": float(k["mean"]), "median": float(k["median"]), "p95": float(k["p95"]), "pairs": int(k["pairs"])}}
    for name, (m, med, q, n) in D56_FLOORS.items():
        f = out[name]
        check4(name + " floor mean", f["mean"], m)
        check4(name + " floor median", f["median"], med)
        check4(name + " floor p95", f["p95"], q)
        if f["pairs"] != n:
            stop(name + " floor pairs")
    return out


def lognormal_shape(m, q):
    """D61, rebuilt here: |e| log-normal with mean m and 95th percentile q. sigma^2/2 - z sigma + ln(q/m) = 0, the
    smaller root; mu = ln m - sigma^2/2; sd(e) = sqrt(E e^2) = exp(mu + sigma^2)."""
    L = math.log(q / m)
    disc = Z95 * Z95 - 2.0 * L
    if disc < 0:
        stop("no log-normal with this mean and 95th percentile")
    r1, r2 = Z95 - math.sqrt(disc), Z95 + math.sqrt(disc)
    sigma = r1
    mu = math.log(m) - 0.5 * sigma * sigma
    return {"sigma": sigma, "mu": mu, "roots": [r1, r2], "sd_e": math.exp(mu + sigma * sigma),
            "median": math.exp(mu), "mean_check": math.exp(mu + 0.5 * sigma * sigma), "p95_check": math.exp(mu + Z95 * sigma),
            "sd_e_normal_shape": m * math.sqrt(math.pi / 2.0)}


def shapes():
    fl = read_floors()
    out = {}
    for name, f in fl.items():
        sh = lognormal_shape(f["mean"], f["p95"])
        if abs(sh["sd_e"] - SD_E_TASK[name]) > 1e-15 * SD_E_TASK[name] * 10:
            stop(f"{name} sd(e) rebuilt {sh['sd_e']!r} differs from {SD_E_TASK[name]!r}")
        sh["floor"] = f
        out[name] = sh
    # against Stage 4's archived shapes (D71 (5))
    j = json.load(open(path("results/stage4_experiment.json")))["shapes"]
    for name in out:
        a = j[name]["lognormal"]
        for key in ("sigma", "mu", "sd_e"):
            if abs(a[key] - out[name][key]) > 1e-14 * max(1.0, abs(a[key])):
                stop(f"{name} {key}: rebuilt {out[name][key]!r} against Stage 4's {a[key]!r}")
        if abs(j[name]["normal"]["sd_e"] - out[name]["sd_e_normal_shape"]) > 1e-15:
            stop(f"{name} normal-shape sd(e)")
    return out


def band_medians():
    """D71 (6): labels from stage4_experiment.json (named), values from stage2_lever.json, checked against
    stage4_paper_numbers.json (median_levers_pct_per_K)."""
    named = json.load(open(path("results/stage4_experiment.json")))["named"]
    paper = json.load(open(path("results/stage4_paper_numbers.json")))["median_levers_pct_per_K"]
    kit = read_kit_levers()
    out = {}
    for band in ("warm", "cold"):
        label = named["median_" + band]
        if not label.startswith("KIT "):
            stop("median label " + label)
        proto = label[4:]
        s = kit[band][proto]["s"]
        if paper[band][0] != label or abs(paper[band][1] - 100.0 * s) > 1e-12:
            stop(f"band median {band}: {label} {s} against paper numbers {paper[band]}")
        if s != S_MED_TASK[band]:
            stop(f"band median {band} {s!r} differs from the task's {S_MED_TASK[band]!r}")
        out[band] = {"label": label, "s": s}
    return out


# ------------------------------------------------------------------ the pairs (D72, D73)
def _f(x):
    return float(x) if x not in ("", None) else None


def kit_band_of_25C(kit):
    """The 25 C sets take the band whose archived interval is narrower, protocol by protocol (D72)."""
    out = {}
    for proto in kit["warm"]:
        w = kit["warm"][proto]["hi"] - kit["warm"][proto]["lo"]
        c = kit["cold"][proto]["hi"] - kit["cold"][proto]["lo"]
        out[proto] = {"band": "warm" if w < c else "cold", "width_warm": w, "width_cold": c}
        if w == c:
            stop("equal widths at 25 C for " + proto)
    return out


def sandia_band_of_25C(san):
    out = {}
    for chem in ("LFP", "NCA", "NMC"):
        def width(sec):
            x = san[(chem, sec)]
            return (x["hi"] - x["lo"]) if x["determinable"] else math.inf
        w, c = width("25 to 35 C"), width("15 to 25 C")
        out[chem] = {"band": "warm" if w < c else "cold", "width_warm": w, "width_cold": c}
    return out


def build_pairs():
    kit = read_kit_levers()
    san = read_sandia_levers()
    k25 = kit_band_of_25C(kit)
    s25 = sandia_band_of_25C(san)
    cells = {r["cell"]: r for r in csv.DictReader(open(path("results/stage2_cells.csv")))}
    rows = []
    for r in csv.DictReader(open(path("results/stage2_warmer_pairs.csv"))):
        w, c = cells[r["warm"]], cells[r["cool"]]
        T = int(r["age_temp"])
        proto = r["protocol"]
        dT = float(r["dT_K"])
        Lw, Lc = float(r["L_warm"]), float(r["L_cool"])
        if Lw != float(w["L90_efc"]) or Lc != float(c["L90_efc"]):
            stop("KIT lifetime mismatch " + r["warm"] + " " + r["cool"])
        dT_cells = float(w["T_C"]) - float(c["T_C"])
        if abs(dT - dT_cells) > 1e-12:
            stop("KIT dT mismatch")
        if not dT > 0:
            stop("KIT dT not positive")
        cens = (w["censored90"] == "True") or (c["censored90"] == "True")
        res_w, res_c = float(w["res90_efc"]), float(c["res90_efc"])
        res = max(res_w, res_c)
        tie = abs(Lw - Lc) < res
        if T == 0:
            band, reason = "none", "0 C set: no secant of results/stage2_lever.json reaches 0 C (D72)"
        elif T == 10:
            band, reason = "cold", ""
        elif T == 25:
            band, reason = k25[proto]["band"], ""
        elif T == 40:
            band, reason = "warm", ""
        else:
            stop("KIT temperature " + str(T))
        lev = kit[band][proto] if band != "none" else None
        scored = band != "none" and lev["determinable"] and not cens
        s = lev["s"] if lev else None
        d = float(np.log(Lw / Lc))
        rows.append({
            "dataset": "KIT", "set": "P%03d" % int(r["parameter_id"]), "temperature_C": float(T), "protocol_or_chemistry": proto,
            "window": proto.split()[0], "rates": proto.split()[1],
            "band": band if scored else "none", "set_band": band, "reason_not_scored": reason if not scored else "",
            "warm": r["warm"], "cool": r["cool"], "dT_K": dT, "L_warm": Lw, "L_cool": Lc, "d": d,
            "lever": (lev["secant"] + " " + proto) if lev else "", "s_per_K": s,
            "s_ci95_lo": lev["lo"] if lev else None, "s_ci95_hi": lev["hi"] if lev else None,
            "p": (-s * dT) if (lev and scored) else None,
            "res_warm_efc": res_w, "res_cool_efc": res_c, "resolution_efc": res, "tie": bool(tie),
            "censored": bool(cens), "scored": bool(scored),
            "warm_b": "", "cool_b": "", "dT_b_K": None, "d_b": None, "p_b": None,
            "warm_cell_T_C": float(w["T_C"]), "cool_cell_T_C": float(c["T_C"]),
        })
    for r in csv.DictReader(open(path("results/stage3_pairs.csv"))):
        chem, amb = r["chemistry"], float(r["ambient_C"])
        La, Lb = _f(r["L_a_efc"]), _f(r["L_b_efc"])
        ra, rb = _f(r["rise_a_K"]), _f(r["rise_b_K"])
        qa, qb = _f(r["rest_a_C"]), _f(r["rest_b_C"])
        cens = (r["censored_a"] == "True") or (r["censored_b"] == "True")
        if amb == 15.0:
            band = "cold"
        elif amb == 35.0:
            band = "warm"
        elif amb == 25.0:
            band = s25[chem]["band"]
        else:
            stop("Sandia ambient")
        lev = san[(chem, SAN_BAND_SECANT[band])]
        has_a = ra is not None and rb is not None
        reason = ""
        if not has_a:
            reason = "measure (a) missing for a cell (0.5C discharge group, D50 (5), or LFP cell c at 35 C, D52)"
        elif not lev["determinable"]:
            reason = f"{chem} {SAN_BAND_SECANT[band]} lever not determinable (D52)"
        elif cens:
            reason = "censored lifetime"
        scored = has_a and lev["determinable"] and not cens
        if has_a:
            if ra == rb:
                stop("tied measure (a) in " + r["group"])
            if ra > rb:
                warm, cool, Lw, Lc, dT = r["cell_a"], r["cell_b"], La, Lb, ra - rb
            else:
                warm, cool, Lw, Lc, dT = r["cell_b"], r["cell_a"], Lb, La, rb - ra
            d = float(np.log(Lw / Lc))
        else:
            warm = cool = ""
            Lw = Lc = dT = d = None
        has_b = qa is not None and qb is not None
        if has_b:
            if qa == qb:
                stop("tied measure (b) in " + r["group"])
            if qa > qb:
                wb, cb, Lwb, Lcb, dTb = r["cell_a"], r["cell_b"], La, Lb, qa - qb
            else:
                wb, cb, Lwb, Lcb, dTb = r["cell_b"], r["cell_a"], Lb, La, qb - qa
            db = float(np.log(Lwb / Lcb))
        else:
            wb = cb = ""
            dTb = db = None
        s = lev["s"] if lev["determinable"] else None
        rows.append({
            "dataset": "Sandia", "set": r["group"], "temperature_C": amb, "protocol_or_chemistry": chem,
            "window": r["window"], "rates": f"{r['c_charge']}/{r['c_discharge']}C",
            "band": band if scored else "none", "set_band": band, "reason_not_scored": reason,
            "warm": warm, "cool": cool, "dT_K": dT, "L_warm": Lw, "L_cool": Lc, "d": d,
            "lever": f"{chem} {SAN_BAND_SECANT[band]} (registered 1C)", "s_per_K": s,
            "s_ci95_lo": lev["lo"], "s_ci95_hi": lev["hi"],
            "p": (-s * dT) if scored else None,
            "res_warm_efc": None, "res_cool_efc": None, "resolution_efc": float(r["resolution_efc"]),
            "tie": r["below_resolution"] == "True", "censored": bool(cens), "scored": bool(scored),
            "warm_b": wb, "cool_b": cb, "dT_b_K": dTb, "d_b": db,
            "p_b": (-s * dTb) if (scored and has_b) else None,
            "warm_cell_T_C": None, "cool_cell_T_C": None,
            "cells_in_group": int(r["cells_in_group"]), "cell_a": r["cell_a"], "cell_b": r["cell_b"],
            "L_a_efc": La, "L_b_efc": Lb, "rise_a_K": ra, "rise_b_K": rb, "rest_a_C": qa, "rest_b_C": qb,
        })
    return rows, {"kit_25C": k25, "sandia_25C": s25}


def v1_oriented_check(rows):
    """D71 (3): the oriented measure-(a) file, read by V1 only: orientation, dT and d against it within 1e-12."""
    ori = {r["group"] + "|" + "|".join(sorted([r["warm"], r["cool"]])): r
           for r in csv.DictReader(open(path("results/stage3_pairs_oriented_measure_a.csv")))}
    worst = 0.0
    n = 0
    for x in rows:
        if x["dataset"] != "Sandia" or x["dT_K"] is None:
            continue
        key = x["set"] + "|" + "|".join(sorted([x["warm"], x["cool"]]))
        o = ori[key]
        if o["warm"] != x["warm"] or o["cool"] != x["cool"]:
            stop("Sandia orientation differs from the oriented file at " + key)
        worst = max(worst, abs(float(o["dT_K"]) - x["dT_K"]), abs(float(o["obs_ln_ratio"]) - x["d"]))
        n += 1
    if worst > 1e-12:
        stop(f"oriented file differs by {worst}")
    return {"pairs_compared": n, "oriented_file_rows": len(ori), "max_abs_diff_dT_and_d": worst}


def select(rows, dataset=None, band=None, scored=True, temperature=None, untied=False, exclude_set=None):
    out = []
    for x in rows:
        if dataset and x["dataset"] != dataset:
            continue
        if scored and not x["scored"]:
            continue
        if band and x["band"] != band:
            continue
        if temperature is not None and x["temperature_C"] != temperature:
            continue
        if untied and x["tie"]:
            continue
        if exclude_set and x["set"] in exclude_set:
            continue
        out.append(x)
    return out


def write_json(name, obj):
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, name), "w") as f:
        json.dump(obj, f, indent=1, default=_json_default)


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def write_csv(name, rows, keys=None):
    os.makedirs(OUT, exist_ok=True)
    if keys is None:
        keys = []
        for r in rows:
            for k in r:
                if k not in keys:
                    keys.append(k)
    with open(os.path.join(OUT, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else (repr(r[k]) if isinstance(r.get(k), float) else r.get(k))) for k in keys})


def rng_for(*key):
    """The checker's own generator: SFC64 from a SeedSequence of [BASE_SEED, *key] (not the pipeline's PCG64
    default_rng with [20261029, k], nor V2's Mersenne Twister with 20261030)."""
    return np.random.Generator(np.random.SFC64(np.random.SeedSequence([BASE_SEED, *key])))


BASE_SEED = 5_550_129
