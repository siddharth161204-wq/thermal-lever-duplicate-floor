#!/usr/bin/env python3
"""V3 of Stage 4 (DECISIONS.md D62): the implied divergences of two cases recomputed by hand.

Plain Python only (no numpy and no Stage 4 function): the two cases fixed in D62, 35 C at 1C and h = 10 and 10 C at
2C and h = 5. Each case's spread is read from results/stage4_module.csv; the levers are read from the archived files
directly (results/stage2_lever.json and results/stage3_lever.csv). Every lever's |s| x spread is formed, the band's
values sorted, the median taken as the middle value and the 95th percentile by linear interpolation at position
0.95 (n - 1), each written out; compared with results/stage4_divergence.csv, passing at agreement to 1e-12.

Usage: python3 verification/stage4/v3_hand.py
"""
import csv, json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CASES = [(35.0, 1.0, 10.0), (10.0, 2.0, 5.0)]


def levers(band):
    kit = json.load(open(os.path.join(ROOT, "results", "stage2_lever.json")))
    key = "lever_25_40" if band == "warm" else "lever_10_25"
    out = [("KIT " + r["protocol"], r["s_per_K"]) for r in kit[key]]
    for r in csv.DictReader(open(os.path.join(ROOT, "results", "stage3_lever.csv"))):
        if r["series"] != "registered 1C":
            continue
        if band == "warm" and r["chemistry"] == "NMC" and r["secant"] == "25 to 35 C":
            out.append(("Sandia NMC 25 to 35 C", float(r["s_per_K"])))
        if band == "cold" and r["secant"] == "15 to 25 C":
            out.append((f"Sandia {r['chemistry']} 15 to 25 C", float(r["s_per_K"])))
    return out


def main():
    mod = list(csv.DictReader(open(os.path.join(ROOT, "results", "stage4_module.csv"))))
    div = list(csv.DictReader(open(os.path.join(ROOT, "results", "stage4_divergence.csv"))))
    out = {"check": "V3 (D62): two cases by hand", "cases": []}
    ok_all = True
    for T, c, h in CASES:
        m = [r for r in mod if float(r["ambient_C"]) == T and float(r["crate"]) == c and float(r["h"]) == h]
        assert len(m) == 1
        spread = float(m[0]["spread_K"])
        band = "warm" if T >= 25.0 else "cold"
        lev = levers(band)
        vals = [(name, abs(s) * spread) for name, s in lev]
        xs = sorted(v for _, v in vals)
        n = len(xs)
        med = xs[n // 2] if n % 2 else 0.5 * (xs[n // 2 - 1] + xs[n // 2])
        pos = 0.95 * (n - 1)
        lo = int(math.floor(pos))
        frac = pos - lo
        p95 = xs[lo] + frac * (xs[lo + 1] - xs[lo])
        rows = [r for r in div if float(r["ambient_C"]) == T and float(r["crate"]) == c and float(r["h"]) == h]
        piped = {r["lever"]: float(r["divergence"]) for r in rows}
        diffs = [abs(piped[name] - v) / v for name, v in vals]
        d_med = abs(piped["band median"] - med) / med
        d_p95 = abs(piped["band 95th percentile"] - p95) / p95
        ok = max(diffs) <= 1e-12 and d_med <= 1e-12 and d_p95 <= 1e-12 and len(vals) == (13 if band == "warm" else 15)
        ok_all &= ok
        out["cases"].append({
            "case": f"{T:g} C, {c:g}C, h = {h:g}", "band": band, "levers": n, "spread_K": spread,
            "by_lever": [{"lever": name, "abs_s_per_K": abs(dict(lev)[name]), "divergence": v} for name, v in vals],
            "sorted": xs,
            "median_written_out": f"middle of {n} sorted values, position {n // 2 + 1}: {med!r}",
            "p95_written_out": f"position 0.95 x ({n} - 1) = {pos:.2f}: {xs[lo]!r} + {frac:.2f} x ({xs[lo + 1]!r} - {xs[lo]!r}) = {p95!r}",
            "median": med, "p95": p95, "max_rel_diff_levers": max(diffs), "rel_diff_median": d_med, "rel_diff_p95": d_p95, "pass": ok})
        print(out["cases"][-1]["case"], "median", med, "p95", p95, "pass", ok)
    out["pass"] = ok_all
    json.dump(out, open(os.path.join(HERE, "v3.json"), "w"), indent=1)
    print("V3 pass:", ok_all)


if __name__ == "__main__":
    main()
