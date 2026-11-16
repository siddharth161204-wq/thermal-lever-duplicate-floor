#!/usr/bin/env python3
"""Stage 6, Parts B and C on the lifetimes (DECISIONS.md D93 to D95), run only after Part A0 was committed and pushed.

From the cell table of Part A0 (`a0_cells.build`), each entering cell's cycle life (D94) and Part A0's reliability and
channel results (`results/stage6_a0.json`): B1, the within-set lever s_w, with the pairs form beside it; B2, with
ln R held, over the cells with a valid resistance; B3, the sign test; B4, the residual spread beside the floors; B5,
the instrumented slope, only if A0.7 found that the channel carries a temperature; C1, the between-policy lever s_b;
C2, s_b - s_w; C3, the share of the between-policy variance. Every interval is the percentile interval of 10,000 draws
of sets (`s6_lever`, D93 (6)), seeds [20261116, k] with k as below; the brackets and unscored lines of D93 (7) are the
same analysis with one change each.

Usage: python3 -I src/stage6/s6_parts.py [extract directory] [lifetimes CSV] [results directory]
Writes stage6_within.json and stage6_between.json in the results directory.
"""
import csv, json, math, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s6_common as C                   # noqa: E402
import s6_lever as Lv                   # noqa: E402
import a0_cells as A                    # noqa: E402

SCOPES = {"pooled": ("Severson", "Attia"), "Severson": ("Severson",), "Attia": ("Attia",)}
SEED_K = {"pooled": 1, "Severson": 2, "Attia": 3}
BRACKETS = {"T_max": ("T_max", False), "T_2_50": ("T_2_50", False), "T_min": ("T_min", False),
            "without_flagged": ("T_reg", True)}
WITHIN_KEYS = ("B1", "B2", "B3", "B4", "B5")
BETWEEN_KEYS = ("C1", "C2", "C3")


def read_life(path):
    """cell -> registered cycle life (D94 (2)), from inputs/matr/stage6_lifetimes.csv."""
    out = {}
    for r in csv.DictReader(open(path)):
        v = r.get("cycle_life", "")
        out[r["cell"]] = float(v) if v not in ("", "nan", "None") else float("nan")
    return out


def cells_for(rows, life, papers, tkey="T_reg", unflagged=False, rkey=None, instrument=None):
    """Entering cells of the papers with a finite measure and lifetime (and R or instrument if asked)."""
    out = []
    for r in rows:
        if not r["enters"] or r["paper"] not in papers or (unflagged and r["contact_flag"]):
            continue
        L, T = life.get(r["cell"], float("nan")), r[tkey]
        if not (math.isfinite(L) and L > 0 and math.isfinite(T)):
            continue
        c = {"cell": r["cell"], "set": r["set"], "batch": r["batch"], "paper": r["paper"], "T": T, "y": math.log(L)}
        if rkey:
            if not math.isfinite(r[rkey]):
                continue
            c["r"] = r[rkey]
        if instrument is not None:
            z = instrument.get(r["cell"], float("nan"))
            if z is None or not math.isfinite(z):
                continue
            c["z"] = z
        out.append(c)
    return out


def summ(point, draws, key):
    ci, mc, used = Lv.percentile_ci(draws[key])
    return {"value": float(point[key][0]), "ci95": ci, "mc_sd_of_ends": mc, "draws_used": used}


def core(cells, k):
    """B1, B3, B4, C1 to C3 from one set of draws."""
    sets = Lv.set_table(cells, "T", "y")
    multi = [s for s in sets if s["n"] >= 2]
    if not multi:
        return {"note": "no set of two or more cells"}
    o = Lv.ones(sets)
    W = Lv.draw_weights(sets, C.N_BOOT, np.random.default_rng([C.SEED, k]))
    pw, dw = Lv.within(sets, o), Lv.within(sets, W)
    pb, db = Lv.between(sets, o, pw["s_w"]), Lv.between(sets, W, dw["s_w"])
    sd2 = Lv.percentile_ci(dw["sd_cell"])
    out = {"seed": [C.SEED, k], "draws": C.N_BOOT, "cells_within": int(sum(s["n"] for s in multi)),
           "sets_within": len(multi), "cells_all": int(sum(s["n"] for s in sets)), "sets_all": len(sets),
           "batches": sorted(set(s["batch"] for s in sets)),
           "B1": {"s_w_per_K": summ(pw, dw, "s_w"), "pairs_form_s_per_K": summ(pw, dw, "s_pairs"),
                  "pairs": int(sum(s["n"] * (s["n"] - 1) // 2 for s in multi))},
           "B3": Lv.sign_test(sets),
           "B4": {"sd_cell": summ(pw, dw, "sd_cell"), "df": float(pw["df"][0]),
                  "sqrt2_sd_cell": float(math.sqrt(2.0) * pw["sd_cell"][0]),
                  "sqrt2_sd_cell_ci95": [math.sqrt(2.0) * x for x in sd2[0]]},
           "C1": {"s_b_per_K": summ(pb, db, "s_b")},
           "C2": {"s_b_minus_s_w_per_K": summ(pb, db, "s_b_minus_s_w")},
           "C3": {"share": summ(pb, db, "share"), "share_residual_form": summ(pb, db, "share_residual_form")}}
    return out


def resistance_part(cells, k):
    """B2 over the cells given (all with a valid R), with s_w over the same cells beside it."""
    sets = Lv.set_table(cells, "T", "y", "r")
    multi = [s for s in sets if s["n"] >= 2]
    if not multi:
        return {"note": "no set of two or more cells with a valid resistance"}
    o = Lv.ones(sets)
    W = Lv.draw_weights(sets, C.N_BOOT, np.random.default_rng([C.SEED, k]))
    pr, dr = Lv.with_resistance(sets, o), Lv.with_resistance(sets, W)
    return {"seed": [C.SEED, k], "cells": int(sum(s["n"] for s in multi)), "sets": len(multi),
            "s_w_R_held_per_K": summ(pr, dr, "s_w_R_held"), "lnL_per_lnR": summ(pr, dr, "lnL_per_lnR"),
            "s_w_same_cells_per_K": summ(pr, dr, "s_w_same_cells"),
            "lnR_per_K_within": summ(pr, dr, "lnR_per_K_within")}


def instrument_part(cells, k):
    sets = Lv.set_table(cells, "T", "y", None, "z")
    multi = [s for s in sets if s["n"] >= 2]
    if not multi:
        return {"note": "no set of two or more cells with an instrument"}
    o = Lv.ones(sets)
    W = Lv.draw_weights(sets, C.N_BOOT, np.random.default_rng([C.SEED, k]))
    pi, di = Lv.instrumented(sets, o), Lv.instrumented(sets, W)
    F = float(pi["first_stage_F"][0])
    return {"seed": [C.SEED, k], "cells": int(sum(s["n"] for s in multi)), "sets": len(multi),
            "s_iv_per_K": summ(pi, di, "s_iv"), "first_stage_F": F, "first_stage_slope": float(pi["first_stage_slope"][0]),
            "registered": bool(F > 10.0), "status": "registered (F > 10)" if F > 10.0 else "exploratory (F <= 10)"}


def analysis(rows, life, papers, k, tkey="T_reg", unflagged=False, rkey="lnR", ch=None, rel=None):
    out = core(cells_for(rows, life, papers, tkey, unflagged), k)
    if "note" in out:
        return out
    cr = cells_for(rows, life, papers, tkey, unflagged, rkey=rkey)
    out["B2"] = resistance_part(cr, k + 100) if cr else {"note": "no cell with a valid resistance (D90)"}
    if ch is not None and ch.get("channel_carries_temperature"):
        out["B5"] = instrument_part(cells_for(rows, life, papers, tkey, unflagged, instrument=ch["instrument"]), k + 200)
    if rel:
        rb = {"r_oe": rel, "s_w_over_r": out["B1"]["s_w_per_K"]["value"] / rel,
              "s_b_over_r": out["C1"]["s_b_per_K"]["value"] / rel}
        if "s_w_R_held_per_K" in out["B2"]:
            rb["s_w_R_held_over_r"] = out["B2"]["s_w_R_held_per_K"]["value"] / rel
        out["reliability_bracket"] = rb
    return out


def multi_only_between(rows, life, papers, k):
    cells = cells_for(rows, life, papers)
    size = {}
    for c in cells:
        size[c["set"]] = size.get(c["set"], 0) + 1
    keep = [c for c in cells if size[c["set"]] >= 2]
    r = core(keep, k)
    return {kk: r[kk] for kk in ("seed", "sets_all", "cells_all", "C1", "C2", "C3") if kk in r}


def run(rows, life, a0):
    rel = {s: a0["reliability"][s]["odd_even"]["r"] for s in SCOPES}
    ch = a0["channel"]
    res = {"registered": {}, "brackets": {}, "unscored": {}}
    for s, papers in SCOPES.items():
        res["registered"][s] = analysis(rows, life, papers, SEED_K[s], ch=ch, rel=rel[s])
        for j, (name, (tkey, unf)) in enumerate(BRACKETS.items()):
            res["brackets"].setdefault(name, {})[s] = analysis(rows, life, papers, 10 * (j + 1) + SEED_K[s], tkey, unf)
        u = res["unscored"]
        u.setdefault("T_rise", {})[s] = analysis(rows, life, papers, 60 + SEED_K[s], "T_rise")
        cr = cells_for(rows, life, papers, rkey="lnR_2_10")
        u.setdefault("R_cycles_2_10", {})[s] = resistance_part(cr, 170 + SEED_K[s]) if cr else {"note": "no valid R"}
        u.setdefault("C1_sets_of_two_or_more", {})[s] = multi_only_between(rows, life, papers, 80 + SEED_K[s])
    return res


def split(res):
    """The within-set and between-policy halves of the results, for the two output files."""
    def pick(d, keys):
        if not isinstance(d, dict):
            return d
        keep = {"seed", "draws", "batches", "note", "reliability_bracket", "cells_within", "sets_within", "cells_all",
                "sets_all", "cells", "sets"}
        return {k: v for k, v in d.items() if k in keys or k in keep}
    within, between = {}, {}
    for part in ("registered", "brackets", "unscored"):
        src = res[part]
        if part == "registered":
            within[part] = {s: pick(v, WITHIN_KEYS) for s, v in src.items()}
            between[part] = {s: pick(v, BETWEEN_KEYS) for s, v in src.items()}
        else:
            within[part] = {b: {s: pick(v, WITHIN_KEYS) for s, v in d.items()} for b, d in src.items()
                            if b != "C1_sets_of_two_or_more"}
            between[part] = {b: {s: pick(v, BETWEEN_KEYS) for s, v in d.items()} for b, d in src.items()
                             if b not in ("R_cycles_2_10",)}
    return within, between


if __name__ == "__main__":
    ext = sys.argv[1] if len(sys.argv) > 1 else C.path("inputs/matr/a0")
    lif = sys.argv[2] if len(sys.argv) > 2 else C.path("inputs/matr/stage6_lifetimes.csv")
    out = sys.argv[3] if len(sys.argv) > 3 else C.path("results")
    rows, sets = A.build(ext)
    a0 = json.load(open(os.path.join(out, "stage6_a0.json")))
    res = run(rows, read_life(lif), a0)
    w, b = split(res)
    json.dump(w, open(os.path.join(out, "stage6_within.json"), "w"), indent=1)
    json.dump(b, open(os.path.join(out, "stage6_between.json"), "w"), indent=1)
    for s in SCOPES:
        r = res["registered"][s]
        if "B1" in r:
            print(s, "s_w", r["B1"]["s_w_per_K"]["value"], r["B1"]["s_w_per_K"]["ci95"], "s_b",
                  r["C1"]["s_b_per_K"]["value"], r["C1"]["s_b_per_K"]["ci95"])
