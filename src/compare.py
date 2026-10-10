#!/usr/bin/env python3
"""Step 2(c), Step 3 and the scoring: the implied lifetime divergence, where a module or a cell first offers the floor,
the dataset's own conditions, the secondary results of DECISIONS.md D12 and D16 to D19, and the scorecard of
predictions.md. Reads results/lever.json and results/resolved_cell.json; writes results/comparison.json,
results/implied_divergence.csv and results/prediction_scorecard.json.

Usage: python3 src/compare.py --results results
"""
import argparse, csv, json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import thermal_lever as TL

RATE_LABEL = {1 / 16: "C/16", 1 / 10: "C/10", 1 / 5: "C/5", 1 / 2: "C/2", 1.0: "1C", 2.0: "2C"}
DATASET_RATES = [1 / 16, 1 / 10, 1 / 5, 1 / 2]
LAWS = list(TL.LAWS)


def lab(c):
    return min(RATE_LABEL.items(), key=lambda kv: abs(kv[0] - c))[1]


def law_life(law, T_K, c):
    try:
        L = TL.LAWS[law](T_K, c)
        return L if np.isfinite(L) else float("nan")
    except (ValueError, KeyError):
        return float("nan")


def exact_divergence(law, T_cells_C, c):
    """Law evaluated at each cell's temperature: pair difference hottest/coldest, worst-of-8 shortfall (D18)."""
    L = np.array([law_life(law, t + 273.15, c) for t in T_cells_C])
    if not np.all(np.isfinite(L)):
        return float("nan"), float("nan")
    i_hot, i_cold = int(np.argmax(T_cells_C)), int(np.argmin(T_cells_C))
    pair = abs(L[i_hot] - L[i_cold]) / (0.5 * (L[i_hot] + L[i_cold]))
    short = (L.mean() - L.min()) / L.mean()
    return float(pair), float(short)


def main(res_dir):
    lev = json.load(open(os.path.join(res_dir, "lever.json")))
    cel = json.load(open(os.path.join(res_dir, "resolved_cell.json")))
    floors = lev["floors"]
    s = {law: {round(e["crate"], 6): e["s_per_K"] for e in lev["levers"][law]} for law in LAWS}
    dtf = {law: {round(e["crate"], 6): e["dT_floor_K"] for e in lev["levers"][law]} for law in LAWS}
    mod = {(round(r["crate"], 6), r["h"]): r for r in lev["module"]}
    cell = {(round(r["crate"], 6), r["h"]): r for r in cel["rows"] if r["base_h"] is None}
    cell_plate = {round(r["crate"], 6): r for r in cel["rows"] if r["base_h"] is not None}
    plate = {(round(r["crate"], 6), r["h"]): r for r in lev["module_coldplate"]}

    rows = []
    for law in LAWS:
        for c in TL.CRATES:
            ck = round(c, 6)
            for h in TL.H_LIST:
                m, rc = mod[(ck, h)], cell[(ck, h)]
                sv = s[law][ck]
                pair, short = exact_divergence(law, m["T_cells_C"], c)
                row = {"law": law, "rate": lab(c), "crate": c, "h": h, "q_cell_W": m["q_cell_W"], "s_per_K": sv,
                       "dT_module_K": m["dT_module_K"], "dT_module_hot_minus_mean_K": m["dT_hot_minus_mean_K"],
                       "dT_cell_K": rc["dT_cell_K"], "dT_cell_hot_minus_mean_K": rc["dT_hot_minus_mean_K"],
                       "implied_module": abs(sv) * m["dT_module_K"], "implied_module_hot_minus_mean": abs(sv) * m["dT_hot_minus_mean_K"],
                       "implied_cell": abs(sv) * rc["dT_cell_K"],
                       "exact_pair_hottest_coldest": pair, "exact_worst_of_8_shortfall": short,
                       "T_module_max_C": m["T_max_C"], "T_cell_max_C": rc["T_max_C"]}
                for k, f in floors.items():
                    row[f"module_exceeds_{k}"] = m["dT_module_K"] > abs(dtf[law][ck][k])
                    row[f"cell_exceeds_{k}"] = rc["dT_cell_K"] > abs(dtf[law][ck][k])
                rows.append(row)
    with open(os.path.join(res_dir, "implied_divergence.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    # Step 3: lowest rate on the grid at which the spread exceeds |dT_floor(f)|
    def lowest(law, h, k, which):
        for c in TL.CRATES:
            r = [x for x in rows if x["law"] == law and abs(x["crate"] - c) < 1e-9 and x["h"] == h][0]
            if r[f"{which}_exceeds_{k}"]:
                return lab(c)
        return "none up to 2C"
    step3 = {law: {str(h): {k: {"module": lowest(law, h, k, "module"), "cell": lowest(law, h, k, "cell")} for k in floors}
                   for h in TL.H_LIST} for law in LAWS}

    # the dataset's own conditions: average rates C/16 to C/2, h = 5 to 25, the brief's heat
    ds = [x for x in rows if x["crate"] <= 0.5 + 1e-9]
    dataset = {law: {"max_implied_module": max(x["implied_module"] for x in ds if x["law"] == law),
                     "max_implied_cell": max(x["implied_cell"] for x in ds if x["law"] == law),
                     "max_dT_module_K": max(x["dT_module_K"] for x in ds if x["law"] == law)} for law in LAWS}

    # D19: the cycle-averaged heat of the dataset's protocol (CC charge at C/2 included, CC discharge at c)
    i1 = TL.CAP_AH
    cyc = []
    for c in DATASET_RATES:
        q = TL.R_DC * i1 ** 2 * (0.5 + c) / (2.0 + 1.0 / c)
        for h in TL.H_LIST:
            T, bal, _ = TL.module_steady(q, h)
            d = float(T.max() - T.min())
            cyc.append({"rate": lab(c), "h": h, "q_cycle_avg_W": q, "ratio_to_brief_heat": q / TL.heat_per_cell(c),
                        "dT_module_K": d, **{f"implied_{law}": abs(s[law][round(c, 6)]) * d for law in LAWS}})

    # secondary: Law 1 per printed charge rate, and one-sided secants at the dataset's C/2 charge
    l1 = {}
    for cr in (0.2, 0.4, 0.5, 0.6):
        lnL = {t: TL.ln_L_kucinskis(t, cr) for t in (25.0, 35.0, 45.0)}
        sc = -(lnL[45.0] - lnL[25.0]) / 20.0
        l1[str(cr)] = {"L_25_35_45": [math.exp(lnL[t]) for t in (25.0, 35.0, 45.0)], "s_centred": sc,
                       "s_25_to_35": -(lnL[35.0] - lnL[25.0]) / 10.0, "s_35_to_45": -(lnL[45.0] - lnL[35.0]) / 10.0,
                       "dT_floor_K": {k: f / sc for k, f in floors.items()}}
    # secondary: assumed noise on the printed values (D16), Monte Carlo
    rng = np.random.default_rng(20261010)
    sig = 0.0265 / math.sqrt(2.0)
    n = 20000
    lnL = {(cr, t): TL.ln_L_kucinskis(t, cr) for cr in (0.4, 0.6) for t in (25.0, 45.0)}
    draws = {key: v + rng.normal(0.0, sig, n) for key, v in lnL.items()}
    s_mc = -(0.5 * (draws[(0.4, 45.0)] + draws[(0.6, 45.0)]) - 0.5 * (draws[(0.4, 25.0)] + draws[(0.6, 25.0)])) / 20.0
    l1_mc = {"sigma_lnL_per_printed_value": sig, "draws": n, "s_2.5_50_97.5": [float(np.percentile(s_mc, q)) for q in (2.5, 50, 97.5)],
             "fraction_s_positive": float((s_mc > 0).mean()),
             "abs_dT_floor_p95_2.5_50_97.5_K": [float(np.percentile(np.abs(floors["p95_8.7pct"] / s_mc), q)) for q in (2.5, 50, 97.5)]}
    # secondary: Law 2 by term, and with the text values
    T35 = TL.T_AMB
    l2 = {}
    for c in TL.CRATES:
        a, b, parts, _ = TL.schimpe_cycle_coeffs(T35, c)
        terms = {}
        for name, tt in (("calendar only", ("cal",)), ("high-temperature cycling only", ("high",)),
                         ("low-temperature cycling only", ("low",)), ("cycling terms only", ("high", "low", "hs"))):
            terms[name] = TL.sensitivity(lambda T, cc: TL.lifetime_law2(T, cc, terms=tt), c)
        pt = TL.schimpe_params("text")
        l2[lab(c)] = {"share_of_loss_at_35C": {k: v / a for k, v in parts.items()}, "s_by_term": terms,
                      "s_text_values": TL.sensitivity(lambda T, cc: TL.lifetime_law2(T, cc, p=pt), c),
                      "L_cycles_35C": TL.lifetime_law2(T35, c)}
    # secondary: module sensitivity at C/2, h = 10 (strip conductance, radiation as received, resistance range)
    q2 = TL.heat_per_cell(0.5)
    sens = {}
    for g in (0.0, TL.G_STRIP, 2 * TL.G_STRIP):
        T, _, _ = TL.module_steady(q2, 10.0, g_strip=g)
        sens[f"strip_{1e3 * g:.2f}_mW_per_K"] = float(T.max() - T.min())
    sens["radiation_as_received"] = float(_as_received_spread(q2, 10.0))
    res_json = json.load(open(os.path.join(TL.INPUTS, "lgm50_resistance.json")))
    for pul in res_json["pulses"]:
        if abs(pul["c_rate"] - 0.5) < 1e-9:
            T, _, _ = TL.module_steady((0.5 * TL.CAP_AH) ** 2 * pul["R_ohm"], 10.0)
            sens[f"R_{pul['pulse_s']}s_pulse_{1e3 * pul['R_ohm']:.2f}_mOhm"] = float(T.max() - T.min())

    # scoring (predictions.md)
    score = {}
    p1 = {law: {lab(c): dtf[law][round(c, 6)]["p95_8.7pct"] for c in DATASET_RATES} for law in LAWS}
    score["P1"] = {"statement": "Delta T_floor(8.7 %) at 35 C lies between 1 and 3 K for every law used",
                   "measured_K": p1, "central_K": {"law_1_kucinskis2022": -8.8, "law_2_schimpe2018": 1.7},
                   "held_per_law": {law: all(1.0 <= v <= 3.0 for v in p1[law].values()) for law in LAWS}}
    score["P1"]["held"] = all(score["P1"]["held_per_law"].values())
    v2 = mod[(0.5, 10.0)]["dT_module_K"]
    score["P2"] = {"statement": "At C/2 with h = 10, Delta T_module is below 1 K", "measured_K": v2, "central_K": 0.39,
                   "by_h_K": {str(h): mod[(0.5, h)]["dT_module_K"] for h in TL.H_LIST}, "held": v2 < 1.0}
    v3 = {str(h): cell[(0.5, h)]["dT_cell_K"] for h in TL.H_LIST}
    score["P3"] = {"statement": "At C/2 the within-cell spread Delta T_cell is below 0.5 K (scored at h = 10, must hold at 5 and 25)",
                   "measured_K": v3, "central_K": 0.12, "held": all(v < 0.5 for v in v3.values())}
    grid = [lab(c) for c in TL.CRATES]
    p4 = {}
    for law in LAWS:
        r10, r5 = step3[law]["10.0"]["p95_8.7pct"]["module"], step3[law]["5.0"]["p95_8.7pct"]["module"]
        first = r10 == "none up to 2C" or grid.index(r10) >= grid.index("1C")
        second = (r10 != "none up to 2C" and r5 != "none up to 2C" and grid.index(r5) == grid.index(r10) - 1)
        p4[law] = {"h10": r10, "h5": r5, "first_clause": first, "second_clause": second}
    score["P4"] = {"statement": "The lowest rate at which Delta T_module exceeds Delta T_floor(8.7 %) with h = 10 is 1C or higher, and with h = 5 one step lower",
                   "measured": p4, "central": {"law_1_kucinskis2022": {"h10": "none up to 2C", "h5": "none up to 2C"},
                                               "law_2_schimpe2018": {"h10": "2C", "h5": "2C"}},
                   "held": all(v["first_clause"] and v["second_clause"] for v in p4.values())}
    ratio = {lab(c): max(abs(p1[l][lab(c)]) for l in LAWS) / min(abs(p1[l][lab(c)]) for l in LAWS) for c in DATASET_RATES}
    score["P5"] = {"statement": "The two laws' Delta T_floor(8.7 %) differ by less than a factor of 2", "measured_ratio": ratio,
                   "central_ratio": 5.2, "held": all(v < 2.0 for v in ratio.values())}

    out = {"floors": floors, "s_per_K": {law: {lab(c): s[law][round(c, 6)] for c in TL.CRATES} for law in LAWS},
           "dT_floor_K": {law: {lab(c): dtf[law][round(c, 6)] for c in TL.CRATES} for law in LAWS},
           "step3_lowest_rate": step3, "dataset_conditions_brief_heat": dataset, "dataset_conditions_cycle_averaged_heat_D19": cyc,
           "coldplate": {"module": {f"{lab(c)} h{h:g}": plate[(round(c, 6), h)]["dT_module_K"] for c in TL.CRATES for h in TL.H_LIST},
                         "cell_h10": {lab(c): cell_plate[round(c, 6)]["dT_cell_K"] for c in TL.CRATES}},
           "law1_by_charge_rate": l1, "law1_assumed_noise_mc": l1_mc, "law2_by_term": l2, "module_sensitivity_C2_h10_K": sens}
    with open(os.path.join(res_dir, "comparison.json"), "w") as f:
        json.dump(out, f, indent=1)
    with open(os.path.join(res_dir, "prediction_scorecard.json"), "w") as f:
        json.dump(score, f, indent=1)
    return out, score


def _as_received_spread(q, h):
    """The received script's radiation form (eps sigma F A_face exchange, room radiation over the exposed area), for
    the record only (D10)."""
    T_amb = TL.T_AMB
    n = 8
    a_conv, _ = TL._areas(n, TL.NEIGHBOURS, False)
    a_rad = a_conv.copy()
    g_pair = TL.G_STRIP + TL.G_GAP
    r_pair = TL.EPS * TL.SIGMA * TL.F_VIEW * TL.A_FACE
    T = np.full(n, T_amb + 1.0)
    for _ in range(100):
        F = q - h * a_conv * (T - T_amb) - TL.EPS * TL.SIGMA * a_rad * (T ** 4 - T_amb ** 4)
        J = np.diag(-h * a_conv - 4 * TL.EPS * TL.SIGMA * a_rad * T ** 3)
        for a, b in TL.NEIGHBOURS:
            flow = g_pair * (T[a] - T[b]) + r_pair * (T[a] ** 4 - T[b] ** 4)
            F[a] -= flow; F[b] += flow
            da = g_pair + 4 * r_pair * T[a] ** 3; db = g_pair + 4 * r_pair * T[b] ** 3
            J[a, a] -= da; J[a, b] += db; J[b, a] += da; J[b, b] -= db
        dT = np.linalg.solve(J, -F); T = T + dT
        if np.abs(dT).max() < 1e-11:
            break
    return T.max() - T.min()


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--results", default="results"); a = ap.parse_args()
    out, score = main(a.results)
    print(json.dumps(score, indent=1))
