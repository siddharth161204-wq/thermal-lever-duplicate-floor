"""Compare this check's outputs (out/*.json) with the published results and write the tables of REPORT.md.

Relative difference = (repository value - my value) / |my value|.
Verdicts: "agrees" when |rel| < 0.1 % (deterministic quantities; most agree to 1e-9 or better);
"small difference explained" when the difference is within 1 % (3 % for resolved-cell quantities) and its cause is
identified (named in the note column); "disagrees" otherwise, and for any differing discrete outcome
(a position, a lowest rate, a prediction score).
Usage: python3 make_report.py   (needs out/law1.json, law2.json, module.json, resolved_cell.json, scoring.json)
       python3 make_report.py --fill < template.md > filled.md   (fills {{name}} placeholders from out/named_tables.json)
"""
import csv
import json
import math
import os

from common import FLOOR_KEYS, H_LIST, OUT, RATES, RATE_NAMES, REPO, load_json, repo_json, repo_verif

ROWS = []          # every comparison: section, quantity, repo, mine, rel, verdict, note
FILL_MODE = False  # set by "--fill": read a template on stdin, print it filled, write nothing else to stdout


def rel(repo, mine):
    if isinstance(repo, (str, bool)) or isinstance(mine, (str, bool)) or repo is None or mine is None:
        return None
    if mine == 0:
        return 0.0 if repo == 0 else float("inf")
    return (repo - mine) / abs(mine)


def add(section, quantity, repo, mine, tol_small=0.01, note_small="", note=""):
    r = rel(repo, mine)
    if r is None:
        verdict = "agrees" if repo == mine else "disagrees"
    elif isinstance(repo, float) and math.isnan(repo):
        verdict = "agrees" if (isinstance(mine, float) and math.isnan(mine)) else "disagrees"
    elif abs(r) < 1e-3:
        verdict = "agrees"
    elif abs(r) <= tol_small and note_small:
        verdict = "small difference explained"
        note = (note + "; " if note else "") + note_small
    else:
        verdict = "disagrees"
    ROWS.append({"section": section, "quantity": quantity, "repo": repo, "mine": mine, "rel": r,
                 "verdict": verdict, "note": note})
    return verdict


def add_maxrel(section, quantity, max_rel, note=""):
    """A row summarising many values by their largest absolute relative difference."""
    verdict = "agrees" if max_rel < 1e-3 else "disagrees"
    ROWS.append({"section": section, "quantity": quantity, "repo": "41 values", "mine": "41 values", "rel": max_rel,
                 "verdict": verdict, "note": note})


# published mesh geometry (mesh_cell2.py, refine 1), used only to interpret reported node positions
_RINGS_MM = [0.0, 0.65, 1.3, 2.0, 2.65, 3.3, 4.0, 4.7, 5.4, 6.1, 6.8, 7.5, 8.05, 8.55, 9.0, 9.4, 9.73, 10.0, 10.2,
             10.35, 10.46, 10.55]


def _layers_mm(Lz=70.1, first=0.2, growth=1.3, dmax=2.5):
    half = Lz / 2
    z = [0.0]
    d = first
    while z[-1] + d < half - 0.5 * min(d, dmax):
        z.append(z[-1] + d)
        d = min(d * growth, dmax)
    return z + [half]                        # distances from the nearer end


_LAYERS_MM = _layers_mm()


def snap(r_mm, z_mm, Lz=70.1):
    """Nearest node of the published mesh, as (ring radius, distance from the nearer end)."""
    e = min(z_mm, Lz - z_mm)
    rr = min(_RINGS_MM, key=lambda x: abs(x - r_mm))
    ee = min(_LAYERS_MM, key=lambda x: abs(x - e))
    return round(rr, 3), round(ee, 3)


def add_position(section, quantity, rep_rz, my_rz):
    a = snap(*rep_rz)
    b = snap(*my_rz)
    verdict = "agrees" if a == b else "disagrees"
    ROWS.append({"section": section, "quantity": quantity,
                 "repo": f"({rep_rz[0]:.2f}, {rep_rz[1]:.2f})", "mine": f"({my_rz[0]:.2f}, {my_rz[1]:.2f})", "rel": None,
                 "verdict": verdict,
                 "note": f"nearest published-mesh node to my point: r = {b[0]} mm, {b[1]} mm from the nearer end; reported node: r = {a[0]} mm, {a[1]} mm from the nearer end"})


def fmt(v):
    if isinstance(v, bool) or v is None:
        return str(v)
    if isinstance(v, str):
        return v
    if isinstance(v, float) and math.isnan(v):
        return "nan"
    if v == 0:
        return "0"
    a = abs(v)
    if a >= 1e4 or a < 1e-4:
        return f"{v:.6e}"
    return f"{v:.9g}"


def frel(r):
    if r is None:
        return "n/a"
    if math.isinf(r):
        return "inf"
    if r == 0:
        return "0"
    return f"{r:+.2e}"


def table(section_prefix, title=None):
    rows = [r for r in ROWS if r["section"] == section_prefix]
    lines = []
    if title:
        lines.append(f"**{title}**\n")
    lines.append("| Quantity | Repository | This check | Rel. diff. | Verdict | Note |")
    lines.append("|---|---|---|---|---|---|")
    for r in rows:
        lines.append(f"| {r['quantity']} | {fmt(r['repo'])} | {fmt(r['mine'])} | {frel(r['rel'])} | {r['verdict']} | {r['note']} |")
    return "\n".join(lines) + "\n"


def summary_line(section_prefix):
    rows = [r for r in ROWS if r["section"] == section_prefix]
    n = len(rows)
    cnt = {}
    for r in rows:
        cnt[r["verdict"]] = cnt.get(r["verdict"], 0) + 1
    finite = [abs(r["rel"]) for r in rows if r["rel"] is not None and not math.isinf(r["rel"])]
    mx = max(finite) if finite else 0.0
    return n, cnt, mx


def main():
    L1 = load_json(os.path.join(OUT, "law1.json"))
    L2 = load_json(os.path.join(OUT, "law2.json"))
    MOD = load_json(os.path.join(OUT, "module.json"))
    RC = load_json(os.path.join(OUT, "resolved_cell.json"))
    SC = load_json(os.path.join(OUT, "scoring.json"))
    lever = repo_json("lever.json")
    comp = repo_json("comparison.json")
    card = repo_json("prediction_scorecard.json")
    rcell = repo_json("resolved_cell.json")
    v12 = repo_verif("v1_v2.json")
    v3 = repo_verif("v3_resolved_cell.json")
    tabs = {}

    # ------------------------------------------------------------------ Law 1
    sec = "law1"
    pl1 = {round(x["crate"], 6): x for x in lever["levers"]["law_1_kucinskis2022"]}
    pc = comp["law1_by_charge_rate"]["0.5"]
    for i, T in enumerate((25, 35, 45)):
        add(sec, f"L({T} C) at C/2 charge, cycles", pc["L_25_35_45"][i], L1["C2"][f"L_{T}"])
    add(sec, "L(35 C) in levers (every rate)", pl1[0.5]["L_at_35C"], L1["C2"]["L_35"])
    add(sec, "s(35 C), per K (centred secant)", pl1[0.5]["s_per_K"], L1["C2"]["s_centred_secant"])
    same = all(abs(pl1[k]["s_per_K"] - pl1[0.5]["s_per_K"]) == 0 for k in pl1)
    add(sec, "same s at all six rates (no discharge-rate dependence, D3)", same, True)
    for k in FLOOR_KEYS:
        add(sec, f"Delta T_floor = f/s, {k}, K", pl1[0.5]["dT_floor_K"][k], L1["C2"]["dT_floor_K"][k])
    for k in FLOOR_KEYS:
        add(sec, f"exact pair 2 atanh(f/2)/s, {k}, K", pl1[0.5]["dT_floor_exact_pair_K"][k], L1["C2"]["dT_floor_exact_pair_K"][k])
    for r in ("0.2", "0.4", "0.6"):
        a, b = comp["law1_by_charge_rate"][r], L1["by_charge_rate"][r]
        for i, T in enumerate((25, 35, 45)):
            add(sec, f"charge {r}C: L({T} C)", a["L_25_35_45"][i], b["L_25_35_45"][i])
        add(sec, f"charge {r}C: s centred, per K", a["s_centred"], b["s_centred"])
        add(sec, f"charge {r}C: s 25 to 35 C, per K", a["s_25_to_35"], b["s_25_to_35"])
        add(sec, f"charge {r}C: s 35 to 45 C, per K", a["s_35_to_45"], b["s_35_to_45"])
        for k in FLOOR_KEYS:
            add(sec, f"charge {r}C: Delta T_floor {k}, K", a["dT_floor_K"][k], b["dT_floor_K"][k])
    curve_rep = lever["lifetime_curves_at_C2"]["law_1_kucinskis2022"]
    mx = max(abs(rel(a["L"], b["L"])) for a, b in zip(curve_rep, L1["curve_C2"]))
    add_maxrel(sec, "lifetime curve at C/2, 15 to 55 C (41 points): largest abs. rel. diff.", mx,
               note="every point compared")
    tabs["law1"] = table(sec)

    sec = "law1mc"
    mc_r = comp["law1_assumed_noise_mc"]
    mc_m = L1["assumed_noise_mc"]
    add(sec, "sigma of ln L per printed value", mc_r["sigma_lnL_per_printed_value"], mc_m["sigma_lnL_per_printed_value"])
    big = mc_m["draws"]["2000000"]
    small = mc_m["draws"]["20000"]
    note_mc = "Monte Carlo sampling noise (20,000 draws in the repository, 2,000,000 here)"
    for i, p in enumerate(("2.5", "50", "97.5")):
        add(sec, f"s percentile {p}, per K", mc_r["s_2.5_50_97.5"][i], big["s_2.5_50_97.5"][i], note_small=note_mc)
    for i, p in enumerate(("2.5", "50", "97.5")):
        add(sec, f"abs Delta T_floor(8.7 %) percentile {p}, K", mc_r["abs_dT_floor_p95_2.5_50_97.5_K"][i],
            big["abs_dT_floor_p95_2.5_50_97.5_K"][i], note_small=note_mc)
    add(sec, "fraction of draws with s > 0", mc_r["fraction_s_positive"], big["fraction_s_positive"])
    add(sec, "D16 hand figure: standard error of the secant, per K", 0.0013, mc_m["se_s_analytic_four_values"],
        note="my propagation over the four printed values that enter the C/2 secant: sigma/20; "
             f"MC std {big['s_std']:.6f}; the D16 form sqrt(2) x 0.0187/20 = {mc_m['se_s_single_printed_rate_two_values']:.5f} is the error of a single printed rate's secant")
    tabs["law1mc"] = table(sec)

    # ------------------------------------------------------------------ Law 2
    sec = "law2"
    pl2 = {round(x["crate"], 6): x for x in lever["levers"]["law_2_schimpe2018"]}
    for c, n in zip(RATES, RATE_NAMES):
        a, b = pl2[round(c, 6)], L2["rates"][n]
        add(sec, f"{n}: L(35 C), cycles to 10 % loss", a["L_at_35C"], b["L_35C"])
        add(sec, f"{n}: s(35 C), per K (centred difference +-0.5 K)", a["s_per_K"], b["s"],
            note=f"+-0.25 K gives {b['s_pm0.25K']:.7f}, exact derivative {b['s_pm0.001K']:.7f}")
        for k in FLOOR_KEYS:
            add(sec, f"{n}: Delta T_floor {k}, K", a["dT_floor_K"][k], b["dT_floor_K"][k])
        for k in FLOOR_KEYS:
            add(sec, f"{n}: exact pair {k}, K", a["dT_floor_exact_pair_K"][k], b["dT_floor_exact_pair_K"][k])
    tabs["law2"] = table(sec)

    sec = "law2terms"
    for n in RATE_NAMES:
        a, b = comp["law2_by_term"][n], L2["rates"][n]
        for k_r, k_m in (("cal", "cal"), ("high", "high"), ("low", "low")):
            add(sec, f"{n}: share of loss at L(35 C), {k_r}", a["share_of_loss_at_35C"][k_r], b["share_of_loss_at_35C"][k_m],
                note=("mine over all four terms; the repository's shares sum to 1 over the three square-root terms "
                      f"(over three, mine is {b['share_over_three_sqrt_terms'][k_m]:.10f}); fourth term "
                      f"{b['share_of_loss_at_35C']['high_soc']:.1e}, not reported") if k_r == "cal" else "")
        for k in ("calendar only", "high-temperature cycling only", "low-temperature cycling only", "cycling terms only"):
            add(sec, f"{n}: s by term, {k}", a["s_by_term"][k], b["s_by_term"][k])
        add(sec, f"{n}: s with the text values (secondary)", a["s_text_values"], b["s_text_values"])
        add(sec, f"{n}: L(35 C) in law2_by_term", a["L_cycles_35C"], b["L_35C"])
    curve_rep = lever["lifetime_curves_at_C2"]["law_2_schimpe2018"]
    mx = max(abs(rel(a["L"], b["L"])) for a, b in zip(curve_rep, L2["curve_C2"]))
    add_maxrel(sec, "lifetime curve at C/2, 15 to 55 C (41 points): largest abs. rel. diff.", mx,
               note=f"includes the fourth term: without it L(15 C) would be {L2['curve_C2'][0]['L_without_high_soc']:.1f} cycles "
                    f"(published {curve_rep[0]['L']:.1f}, mine {L2['curve_C2'][0]['L']:.1f})")
    tabs["law2terms"] = table(sec)

    sec = "law2v1"
    for T in ("25", "45"):
        a = v12["V1_law2_schimpe2018"][f"{T}C"]
        b = L2["V1"][f"{T}.0"]
        add(sec, f"V1 {T} C, 2800 FEC at 1C: loss, % (full Eq. 23, four terms)", a["model_loss_pct"], b["loss_pct_full_eq23"],
            note=f"three square-root terms only: {b['loss_pct_three_sqrt_terms']:.6f} %; fourth term {b['parts_pct']['high_soc']:.4f} points")
        for k_r, k_m in (("cal", "cal"), ("high", "high"), ("low", "low")):
            add(sec, f"V1 {T} C: {k_r} part, %", a["parts_pct"][k_r], b["parts_pct"][k_m])
        add(sec, f"V1 {T} C: text values, three terms, %", a["model_loss_pct_text_values"], b["loss_pct_three_terms_text_values"])
        lo_, hi_ = a["accept_range_pct"]
        add(sec, f"V1 {T} C: pass (accept {lo_:.0f} to {hi_:.0f} %)", a["pass"], lo_ <= b["loss_pct_full_eq23"] <= hi_)
    tabs["law2v1"] = table(sec)

    # ------------------------------------------------------------------ module
    sec = "module"
    mm = {(r["crate"], r["h"]): r for r in MOD["module"]}
    for a in lever["module"]:
        b = mm[(a["crate"], a["h"])]
        n = RATE_NAMES[RATES.index(a["crate"])]
        tag = f"{n}, h={a['h']:.0f}"
        add(sec, f"{tag}: q per cell, W", a["q_cell_W"], b["q_cell_W"])
        add(sec, f"{tag}: Delta T_module, K", a["dT_module_K"], b["dT_module_K"])
        add(sec, f"{tag}: T_max, C", a["T_max_C"], b["T_max_C"])
        add(sec, f"{tag}: T_min, C", a["T_min_C"], b["T_min_C"])
        add(sec, f"{tag}: T_mean, C", a["T_mean_C"], b["T_mean_C"])
        add(sec, f"{tag}: hottest minus mean, K", a["dT_hot_minus_mean_K"], b["dT_hot_minus_mean_K"])
        add(sec, f"{tag}: hottest position", a["hot_position"], b["hot_position"])
        add(sec, f"{tag}: coldest position", a["cold_position"], b["cold_position"])
        add(sec, f"{tag}: energy balance closes (|residual| < 1e-12 W)", abs(a["balance_W"]) < 1e-12, abs(b["balance_W"]) < 1e-12)
    sec = "coldplate"
    mc = {(r["crate"], r["h"]): r for r in MOD["module_coldplate"]}
    for a in lever["module_coldplate"]:
        b = mc[(a["crate"], a["h"])]
        n = RATE_NAMES[RATES.index(a["crate"])]
        tag = f"{n}, h={a['h']:.0f}"
        add(sec, f"{tag}: Delta T_module, K", a["dT_module_K"], b["dT_module_K"])
        add(sec, f"{tag}: T_max, C", a["T_max_C"], b["T_max_C"])
        add(sec, f"{tag}: T_mean, C", a["T_mean_C"], b["T_mean_C"])
        add(sec, f"{tag}: hottest minus mean, K", a["dT_hot_minus_mean_K"], b["dT_hot_minus_mean_K"])
        add(sec, f"{tag}: hottest / coldest position", f"{a['hot_position']} / {a['cold_position']}",
            f"{b['hot_position']} / {b['cold_position']}")
        add(sec, f"{tag}: energy balance closes (|residual| < 1e-12 W)", abs(a["balance_W"]) < 1e-12, abs(b["balance_W"]) < 1e-12)
    sec = "modulemisc"
    for h in H_LIST:
        a = [x for x in v12["V2_module"]["b_strip_1e3_uniformity"] if x["h"] == h][0]
        b = [x for x in MOD["V2_strip_1e3"] if x["h"] == h][0]
        add(sec, f"V2 (b): strip 1e3 W/K, 2C, h={h:.0f}: max minus min, K", a["max_minus_min_K"], b["max_minus_min_K"])
    for k, v in comp["module_sensitivity_C2_h10_K"].items():
        add(sec, f"sensitivity C/2, h=10: {k}, K", v, MOD["sensitivity_C2_h10"][k],
            note=("the received form: pair exchange eps sigma F (A_lat/4), room radiation over the convective area; "
                  f"with eps^2 instead of eps it would be {MOD['sensitivity_C2_h10']['radiation_quarter_area_eps2_variant']:.6f} K")
            if k == "radiation_as_received" else "")
    for c, n in zip(RATES, RATE_NAMES):
        for h in H_LIST:
            a = comp["coldplate"]["module"][f"{n} h{int(h)}"]
            add(sec, f"comparison.json coldplate.module {n} h{int(h)}, K", a, mc[(c, h)]["dT_module_K"])
    tabs["module"] = table("module")
    tabs["coldplate"] = table("coldplate")
    tabs["modulemisc"] = table("modulemisc")

    # two-node sanity check: the published eight-node spread against my exact corner/inner reduction
    sec = "twonode"
    central = {5.0: 0.37, 10.0: 0.39, 25.0: 0.35}
    for h in H_LIST:
        tn = MOD["two_node_C2"][str(h)]
        add(sec, f"C/2, h={h:.0f}: published eight-node Delta T_module vs my two-node (corner/inner) reduction, K",
            mm[(0.5, h)]["dT_module_K"] if False else [x for x in lever["module"] if x["crate"] == 0.5 and x["h"] == h][0]["dT_module_K"],
            tn["nonlinear_two_node"]["spread"], note="symmetry makes the reduction exact")
    tabs["twonode"] = table(sec)
    info = ["| h, W/(m2 K) | predictions.md hand estimate, K | my two-node, linearised at 35 C, K | my two-node, linearised at the module mean, K | exact, K |",
            "|---|---|---|---|---|"]
    for h in H_LIST:
        tn = MOD["two_node_C2"][str(h)]
        info.append(f"| {h:.0f} | {central[h]:.2f} | {tn['linearised_at_35C']['spread']:.4f} | {tn['linearised_at_mean']['spread']:.4f} | {tn['nonlinear_two_node']['spread']:.4f} |")
    tabs["twonode_info"] = "\n".join(info) + "\n"

    # ------------------------------------------------------------------ resolved cell
    sec = "cellconv"
    disc = "discretisation of the published mesh (Section 4, item 3)"
    conv = RC["convergence"]
    pub = {(r["crate"], r["h"], r["base_h"]): r for r in rcell["rows"]}
    conv_lines = ["| Case | Quantity | Nr=20 | Nr=40 | Nr=80 | Nr=160 | Nr=320 | Observed order | Extrapolated | Repository (refine 1) | Rel. diff. (repo vs extrapolated) |",
                  "|---|---|---|---|---|---|---|---|---|---|---|"]
    for key, cc in conv.items():
        c, h, bh = cc["crate"], cc["h"], cc["base_h"]
        n = RATE_NAMES[RATES.index(c)]
        if bh is not None and c not in (0.5, 2.0):
            continue
        p = pub[(c, h, bh)]
        for q, label in (("dT_cell_K", "hottest minus coldest, K"), ("dT_hot_minus_mean_K", "hottest minus volume mean, K")):
            seq = [s[q] for s in cc["sequence"]]
            ri = cc["richardson"][q]
            conv_lines.append(f"| {n}, h={h:.0f}{', base 600' if bh else ''} | {label} | " + " | ".join(f"{v:.7f}" for v in seq) +
                              f" | {('n/a' if math.isnan(ri['order']) else format(ri['order'], '.2f'))} | {ri['extrapolated']:.7f} | {p[q]:.7f} | {frel(rel(p[q], ri['extrapolated']))} |")
    tabs["cellconv"] = "\n".join(conv_lines) + "\n"

    for key, cc in conv.items():
        c, h, bh = cc["crate"], cc["h"], cc["base_h"]
        n = RATE_NAMES[RATES.index(c)]
        p = pub[(c, h, bh)]
        ri = cc["richardson"]
        tag = f"{n}, h={h:.0f}{', base 600 W/(m2 K)' if bh else ''}"
        add(sec, f"{tag}: hottest minus coldest, K", p["dT_cell_K"], ri["dT_cell_K"]["extrapolated"], tol_small=0.03, note_small=disc)
        add(sec, f"{tag}: hottest minus volume mean, K", p["dT_hot_minus_mean_K"], ri["dT_hot_minus_mean_K"]["extrapolated"], tol_small=0.03, note_small=disc)
    for case in v3["cases"]:
        h = case["h"]
        cc = conv[f"c=0.5_h={h}_base=None"]
        add(sec, f"V3 (d) C/2, h={h:.0f}: refine 2 hottest minus coldest, K", case["d_dT_cell_refine2_K"],
            cc["richardson"]["dT_cell_K"]["extrapolated"], tol_small=0.03, note_small="discretisation (refine 2 is closer to the converged value)")
    tabs["cellconv_cmp"] = table(sec)

    sec = "cellall"
    myall = {(r["crate"], r["h"]): r for r in RC["all_cases"]}
    for p in rcell["rows"]:
        c, h, bh = p["crate"], p["h"], p["base_h"]
        n = RATE_NAMES[RATES.index(c)]
        if bh is None:
            b = myall[(c, h)]
        else:
            seqs = conv[f"c={c}_h={h}_base={bh}"]["sequence"]
            b = [s for s in seqs if s["nr"] == 160][0]
        tag = f"{n}, h={h:.0f}{', base' if bh else ''}"
        add(sec, f"{tag}: hottest minus coldest, K", p["dT_cell_K"], b["dT_cell_K"], tol_small=0.03, note_small=disc)
        add(sec, f"{tag}: hottest minus mean, K", p["dT_hot_minus_mean_K"], b["dT_hot_minus_mean_K"], tol_small=0.03, note_small=disc)
        add(sec, f"{tag}: T_max rise, K", p["T_max_C"] - 35.0, b["T_max_C"] - 35.0, tol_small=0.03, note_small=disc)
        add(sec, f"{tag}: volume-mean rise, K", p["T_vol_mean_C"] - 35.0, b["T_vol_mean_C"] - 35.0, tol_small=0.03, note_small=disc)
        add(sec, f"{tag}: surface-mean rise, K", p["T_surface_mean_C"] - 35.0, b["T_surface_mean_C"] - 35.0, tol_small=0.03, note_small=disc)
        add(sec, f"{tag}: lumped-cell rise, K", p["T_lumped_C"] - 35.0, b["T_lumped_C"] - 35.0)
        add_position(sec, f"{tag}: hottest point (r, z) mm", (p["hot"]["r_mm"], p["hot"]["z_mm"]), (b["hot_r_mm"], b["hot_z_mm"]))
        add_position(sec, f"{tag}: coldest point (r, z) mm", (p["cold"]["r_mm"], p["cold"]["z_mm"]), (b["cold_r_mm"], b["cold_z_mm"]))
    tabs["cellall"] = table(sec)

    sec = "cellv3"
    for case in v3["cases"]:
        h = case["h"]
        cc = conv[f"c=0.5_h={h}_base=None"]
        fin = [s for s in cc["sequence"] if s["nr"] == 320][0]
        add(sec, f"V3 (b) C/2, h={h:.0f}: (surface mean - lumped)/rise", case["b_surface_mean_minus_lumped_over_rise"],
            fin["V3b_surface_mean_minus_lumped_over_rise"], tol_small=0.03, note_small="second-order radiation effect of a size set by the surface temperature variation (Section 3.4)", note="")
        add(sec, f"V3 (c) C/2, h={h:.0f}: (volume mean - lumped)/rise", case["c_volume_mean_minus_lumped_over_rise"],
            fin["V3c_volume_mean_minus_lumped_over_rise"], tol_small=0.03, note_small=disc)
        add(sec, f"V3 (c) C/2, h={h:.0f}: brief's literal 1 % test passes", case["c_brief_literal_pass_1pct"],
            abs(fin["V3c_volume_mean_minus_lumped_over_rise"]) < 0.01)
    e = v3["e_analytic_radial"]
    an = RC["analytic_infinite_cylinder"][-1]
    add(sec, "V3 (e) analytic centre-minus-surface (infinite cylinder, C/2 heat), K", e["center_to_surface_analytic_K"],
        an["centre_minus_surface_exact"])
    tabs["cellv3"] = table(sec)

    sec = "cellbase"
    for c, n in zip(RATES, RATE_NAMES):
        cc = conv[f"c={c}_h=10.0_base=600.0"]
        add(sec, f"base-cooled cell, h=10, {n}: hottest minus coldest, K (comparison.json)", comp["coldplate"]["cell_h10"][n],
            cc["richardson"]["dT_cell_K"]["extrapolated"], tol_small=0.03, note_small=disc)
    tabs["cellbase"] = table(sec)

    # ------------------------------------------------------------------ comparison and scoring
    sec = "implied"
    rep_rows = list(csv.DictReader(open(os.path.join(REPO, "results", "implied_divergence.csv"))))
    mine_rows = {(r["law"], r["rate"], r["h"]): r for r in SC["rows"]}
    for a in rep_rows:
        b = mine_rows[(a["law"], a["rate"], float(a["h"]))]
        tag = f"{'Law 1' if a['law'].startswith('law_1') else 'Law 2'} {a['rate']} h={float(a['h']):.0f}"
        add(sec, f"{tag}: implied module", float(a["implied_module"]), b["implied_module"])
        add(sec, f"{tag}: implied cell", float(a["implied_cell"]), b["implied_cell"], tol_small=0.03, note_small=disc)
        ep = float(a["exact_pair_hottest_coldest"])
        note_rng = "" if b["cells_inside_law_range"] else "cells outside the law's range (D18 says report only inside)"
        add(sec, f"{tag}: implied module, hottest minus mean", float(a["implied_module_hot_minus_mean"]), b["implied_module_hot_minus_mean"])
        add(sec, f"{tag}: dT_module, K", float(a["dT_module_K"]), b["dT_module_K"])
        add(sec, f"{tag}: dT_cell, K", float(a["dT_cell_K"]), b["dT_cell_K"], tol_small=0.03, note_small=disc)
        add(sec, f"{tag}: s, per K", float(a["s_per_K"]), b["s_per_K"])
        add(sec, f"{tag}: exact pair (D18)", ep, b["exact_pair_hottest_coldest"] if (b["cells_inside_law_range"] or not math.isnan(ep)) else float("nan"),
            note=note_rng)
        ew = float(a["exact_worst_of_8_shortfall"])
        add(sec, f"{tag}: exact worst-of-8 shortfall (D18)", ew,
            b["exact_worst_of_8_shortfall"] if (b["cells_inside_law_range"] or not math.isnan(ew)) else float("nan"), note=note_rng)
        for k in FLOOR_KEYS:
            add(sec, f"{tag}: module exceeds {k}", a[f"module_exceeds_{k}"] == "True", b[f"module_exceeds_{k}"])
            add(sec, f"{tag}: cell exceeds {k}", a[f"cell_exceeds_{k}"] == "True", b[f"cell_exceeds_{k}"])
    tabs["implied_summary"] = None

    sec = "step3"
    for law in ("law_1_kucinskis2022", "law_2_schimpe2018"):
        for h in H_LIST:
            for k in FLOOR_KEYS:
                for w in ("module", "cell"):
                    add(sec, f"{'Law 1' if law.startswith('law_1') else 'Law 2'}, h={h:.0f}, {k}, {w}",
                        comp["step3_lowest_rate"][law][str(h)][k][w], SC["step3_lowest_rate"][law][str(h)][k][w])
    for law in ("law_1_kucinskis2022", "law_2_schimpe2018"):
        for q in ("max_implied_module", "max_implied_cell", "max_dT_module_K"):
            add(sec, f"dataset conditions (C/16 to C/2, every h), {('Law 1' if law.startswith('law_1') else 'Law 2')}: {q}",
                comp["dataset_conditions_brief_heat"][law][q], SC["dataset_conditions_brief_heat"][law][q],
                tol_small=0.03, note_small=disc if "cell" in q else "")
    tabs["step3"] = table(sec)

    sec = "d19"
    md = {(r["rate"], r["h"]): r for r in SC["D19"]}
    for a in comp["dataset_conditions_cycle_averaged_heat_D19"]:
        b = md[(a["rate"], a["h"])]
        tag = f"{a['rate']}, h={a['h']:.0f}"
        add(sec, f"{tag}: cycle-averaged heat, W", a["q_cycle_avg_W"], b["q_cycle_avg_W"])
        add(sec, f"{tag}: ratio to the brief's heat", a["ratio_to_brief_heat"], b["ratio_to_brief_heat"])
        add(sec, f"{tag}: Delta T_module, K", a["dT_module_K"], b["dT_module_K"])
        add(sec, f"{tag}: implied, Law 1", a["implied_law_1_kucinskis2022"], b["implied_law_1_kucinskis2022"])
        add(sec, f"{tag}: implied, Law 2", a["implied_law_2_schimpe2018"], b["implied_law_2_schimpe2018"])
    tabs["d19"] = table(sec)

    sec = "scores"
    S = SC["scorecard"]
    for law in ("law_1_kucinskis2022", "law_2_schimpe2018"):
        for n in RATE_NAMES[:4]:
            add(sec, f"P1 value {('Law 1' if law.startswith('law_1') else 'Law 2')} {n}, K", card["P1"]["measured_K"][law][n], S["P1"]["values_K"][law][n])
        add(sec, f"P1 held for {('Law 1' if law.startswith('law_1') else 'Law 2')}", card["P1"]["held_per_law"][law], S["P1"]["held_per_law"][law])
    add(sec, "P1 held", card["P1"]["held"], S["P1"]["held"])
    add(sec, "P2 value (C/2, h=10), K", card["P2"]["measured_K"], S["P2"]["value_K"])
    add(sec, "P2 held", card["P2"]["held"], S["P2"]["held"])
    for h in H_LIST:
        add(sec, f"P3 value C/2, h={h:.0f}, K", card["P3"]["measured_K"][str(h)], S["P3"]["values_K"][str(h)], tol_small=0.03, note_small=disc)
    add(sec, "P3 held", card["P3"]["held"], S["P3"]["held"])
    for law in ("law_1_kucinskis2022", "law_2_schimpe2018"):
        lab = "Law 1" if law.startswith("law_1") else "Law 2"
        for q in ("h10", "h5", "first_clause", "second_clause"):
            add(sec, f"P4 {lab} {q}", card["P4"]["measured"][law][q], S["P4"][law][q])
    add(sec, "P4 held", card["P4"]["held"], S["P4"]["held"])
    for n in RATE_NAMES[:4]:
        add(sec, f"P5 ratio {n}", card["P5"]["measured_ratio"][n], S["P5"]["ratio"][n])
    add(sec, "P5 held", card["P5"]["held"], S["P5"]["held"])
    tabs["scores"] = table(sec)

    # ------------------------------------------------------------------ compact body tables
    body = {}
    t = ["| Rate | L(35 C) repo | L(35 C) mine | s repo, per K | s mine, per K | Delta T_floor 3.0 % repo / mine, K | 4.1 % repo / mine, K | 8.7 % repo / mine, K | largest abs. rel. diff. |",
         "|---|---|---|---|---|---|---|---|---|"]
    for c, n in zip(RATES, RATE_NAMES):
        a, b = pl2[round(c, 6)], L2["rates"][n]
        rels = [abs(rel(a["L_at_35C"], b["L_35C"])), abs(rel(a["s_per_K"], b["s"]))] + \
               [abs(rel(a["dT_floor_K"][k], b["dT_floor_K"][k])) for k in FLOOR_KEYS]
        t.append(f"| {n} | {a['L_at_35C']:.4f} | {b['L_35C']:.4f} | {a['s_per_K']:.7f} | {b['s']:.7f} | "
                 + " | ".join(f"{a['dT_floor_K'][k]:.5f} / {b['dT_floor_K'][k]:.5f}" for k in FLOOR_KEYS) + f" | {max(rels):.1e} |")
    body["law2_compact"] = "\n".join(t) + "\n"
    a, b = comp["law2_by_term"]["C/2"], L2["rates"]["C/2"]
    t = ["| Term (C/2 discharge, at L(35 C)) | Share, repo | Share, mine (over three square-root terms) | Share, mine (over all four terms) | s of the term alone, repo / mine, per K |",
         "|---|---|---|---|---|"]
    names = {"cal": "calendar only", "high": "high-temperature cycling only", "low": "low-temperature cycling only"}
    for k in ("cal", "high", "low"):
        t.append(f"| {names[k]} | {a['share_of_loss_at_35C'][k]:.6f} | {b['share_over_three_sqrt_terms'][k]:.6f} | {b['share_of_loss_at_35C'][k]:.6f} | "
                 f"{a['s_by_term'][names[k]]:.6f} / {b['s_by_term'][names[k]]:.6f} |")
    t.append(f"| low-temperature high-SOC (fourth) | not reported | n/a | {b['share_of_loss_at_35C']['high_soc']:.1e} | n/a |")
    body["law2_terms_compact"] = "\n".join(t) + "\n"

    def module_compact(rep_list, mine_map):
        t = ["| Rate | h | q per cell, W | Delta T repo, K | Delta T mine, K | rel. diff. | T_max repo / mine, C | T_mean repo / mine, C | hottest / coldest (same in both) |",
             "|---|---|---|---|---|---|---|---|---|"]
        for a in rep_list:
            b = mine_map[(a["crate"], a["h"])]
            n = RATE_NAMES[RATES.index(a["crate"])]
            same = (a["hot_position"] == b["hot_position"]) and (a["cold_position"] == b["cold_position"])
            t.append(f"| {n} | {a['h']:.0f} | {a['q_cell_W']:.6g} | {a['dT_module_K']:.9g} | {b['dT_module_K']:.9g} | {frel(rel(a['dT_module_K'], b['dT_module_K']))} | "
                     f"{a['T_max_C']:.6f} / {b['T_max_C']:.6f} | {a['T_mean_C']:.6f} / {b['T_mean_C']:.6f} | "
                     f"{a['hot_position']} / {a['cold_position']}{'' if same else ' (DIFFERENT)'} |")
        return "\n".join(t) + "\n"
    body["module_compact"] = module_compact(lever["module"], mm)
    body["coldplate_compact"] = module_compact(lever["module_coldplate"], mc)

    t = ["| Case | Delta T_cell repo, K | mine (Nr=160), K | rel. diff. | hottest minus mean repo, K | mine, K | rel. diff. | volume-mean rise repo / mine, K | surface-mean rise repo / mine, K |",
         "|---|---|---|---|---|---|---|---|---|"]
    for p in rcell["rows"]:
        c, h, bh = p["crate"], p["h"], p["base_h"]
        n = RATE_NAMES[RATES.index(c)]
        if bh is None:
            b = myall[(c, h)]
        else:
            b = [s for s in conv[f"c={c}_h={h}_base={bh}"]["sequence"] if s["nr"] == 160][0]
        t.append(f"| {n}, h={h:.0f}{', base 600' if bh else ''} | {p['dT_cell_K']:.6g} | {b['dT_cell_K']:.6g} | {frel(rel(p['dT_cell_K'], b['dT_cell_K']))} | "
                 f"{p['dT_hot_minus_mean_K']:.6g} | {b['dT_hot_minus_mean_K']:.6g} | {frel(rel(p['dT_hot_minus_mean_K'], b['dT_hot_minus_mean_K']))} | "
                 f"{p['T_vol_mean_C'] - 35:.6g} / {b['T_vol_mean_C'] - 35:.6g} | {p['T_surface_mean_C'] - 35:.6g} / {b['T_surface_mean_C'] - 35:.6g} |")
    body["cell_compact"] = "\n".join(t) + "\n"

    t = ["| Law | Rate | s, per K | Delta T_module, K | implied module repo / mine | Delta T_cell repo / mine, K | implied cell repo / mine |",
         "|---|---|---|---|---|---|---|"]
    for a in rep_rows:
        if float(a["h"]) != 10.0:
            continue
        b = mine_rows[(a["law"], a["rate"], 10.0)]
        t.append(f"| {'Law 1' if a['law'].startswith('law_1') else 'Law 2'} | {a['rate']} | {float(a['s_per_K']):.6f} | {float(a['dT_module_K']):.6g} | "
                 f"{float(a['implied_module']):.6g} / {b['implied_module']:.6g} | {float(a['dT_cell_K']):.6g} / {b['dT_cell_K']:.6g} | "
                 f"{float(a['implied_cell']):.6g} / {b['implied_cell']:.6g} |")
    body["implied_h10"] = "\n".join(t) + "\n"

    t = ["| Law | h | floor | module: repo / mine | cell: repo / mine |", "|---|---|---|---|---|"]
    for law in ("law_1_kucinskis2022", "law_2_schimpe2018"):
        for h in H_LIST:
            for k in FLOOR_KEYS:
                ra = comp["step3_lowest_rate"][law][str(h)][k]
                rb = SC["step3_lowest_rate"][law][str(h)][k]
                t.append(f"| {'Law 1' if law.startswith('law_1') else 'Law 2'} | {h:.0f} | {k} | {ra['module']} / {rb['module']} | {ra['cell']} / {rb['cell']} |")
    body["step3_compact"] = "\n".join(t) + "\n"

    t = ["| Rate | h | heat repo / mine, W | ratio to the brief's heat | Delta T_module repo / mine, K | implied Law 1 repo / mine | implied Law 2 repo / mine |",
         "|---|---|---|---|---|---|---|"]
    for a in comp["dataset_conditions_cycle_averaged_heat_D19"]:
        b = md[(a["rate"], a["h"])]
        t.append(f"| {a['rate']} | {a['h']:.0f} | {a['q_cycle_avg_W']:.6g} / {b['q_cycle_avg_W']:.6g} | {a['ratio_to_brief_heat']:.4f} / {b['ratio_to_brief_heat']:.4f} | "
                 f"{a['dT_module_K']:.6g} / {b['dT_module_K']:.6g} | {a['implied_law_1_kucinskis2022']:.5g} / {b['implied_law_1_kucinskis2022']:.5g} | "
                 f"{a['implied_law_2_schimpe2018']:.5g} / {b['implied_law_2_schimpe2018']:.5g} |")
    body["d19_compact"] = "\n".join(t) + "\n"

    t = ["| Prediction | Repository value | My value | Repository score | My score |", "|---|---|---|---|---|"]
    t.append(f"| P1: Delta T_floor(8.7 %) in [1, 3] K for every law | Law 1 {card['P1']['measured_K']['law_1_kucinskis2022']['C/2']:.4f} K; Law 2 {card['P1']['measured_K']['law_2_schimpe2018']['C/16']:.4f} to {card['P1']['measured_K']['law_2_schimpe2018']['C/2']:.4f} K | "
             f"Law 1 {S['P1']['values_K']['law_1_kucinskis2022']['C/2']:.4f} K; Law 2 {S['P1']['values_K']['law_2_schimpe2018']['C/16']:.4f} to {S['P1']['values_K']['law_2_schimpe2018']['C/2']:.4f} K | "
             f"{'held' if card['P1']['held'] else 'failed'} | {'held' if S['P1']['held'] else 'failed'} |")
    t.append(f"| P2: Delta T_module(C/2, h=10) below 1 K | {card['P2']['measured_K']:.6f} K | {S['P2']['value_K']:.6f} K | {'held' if card['P2']['held'] else 'failed'} | {'held' if S['P2']['held'] else 'failed'} |")
    t.append(f"| P3: Delta T_cell(C/2) below 0.5 K at h = 5, 10, 25 | " + ", ".join(f"{card['P3']['measured_K'][str(h)]:.5f}" for h in H_LIST) + " K | "
             + ", ".join(f"{S['P3']['values_K'][str(h)]:.5f}" for h in H_LIST) + f" K | {'held' if card['P3']['held'] else 'failed'} | {'held' if S['P3']['held'] else 'failed'} |")
    t.append(f"| P4: lowest rate (8.7 % floor) at h=10 is 1C or higher, at h=5 one step lower | Law 1: {card['P4']['measured']['law_1_kucinskis2022']['h10']} / {card['P4']['measured']['law_1_kucinskis2022']['h5']}; Law 2: {card['P4']['measured']['law_2_schimpe2018']['h10']} / {card['P4']['measured']['law_2_schimpe2018']['h5']} | "
             f"Law 1: {S['P4']['law_1_kucinskis2022']['h10']} / {S['P4']['law_1_kucinskis2022']['h5']}; Law 2: {S['P4']['law_2_schimpe2018']['h10']} / {S['P4']['law_2_schimpe2018']['h5']} | "
             f"{'held' if card['P4']['held'] else 'failed'} (first clause held, second failed, both laws) | {'held' if S['P4']['held'] else 'failed'} (same clauses) |")
    t.append(f"| P5: the laws' Delta T_floor(8.7 %) within a factor of 2 | ratio {card['P5']['measured_ratio']['C/16']:.4f} to {card['P5']['measured_ratio']['C/2']:.4f} | ratio {S['P5']['ratio']['C/16']:.4f} to {S['P5']['ratio']['C/2']:.4f} | "
             f"{'held' if card['P5']['held'] else 'failed'} | {'held' if S['P5']['held'] else 'failed'} |")
    body["scores_compact"] = "\n".join(t) + "\n"
    with open(os.path.join(OUT, "body_tables.json"), "w") as f:
        json.dump(body, f, indent=1)
    keep = ("L(25 C)", "L(35 C) at", "L(45 C)", "s(35 C)", "same s", "Delta T_floor = f/s", "exact pair 2 atanh",
            "lifetime curve")
    sel = [r for r in ROWS if r["section"] == "law1" and (r["quantity"].startswith(keep) or
           ("s centred" in r["quantity"]) or ("Delta T_floor p95" in r["quantity"]))]
    t = ["| Quantity | Repository | This check | Rel. diff. | Verdict |", "|---|---|---|---|---|"]
    for r in sel:
        t.append(f"| {r['quantity']} | {fmt(r['repo'])} | {fmt(r['mine'])} | {frel(r['rel'])} | {r['verdict']} |")
    body["law1_compact"] = "\n".join(t) + "\n"
    sel = [r for r in ROWS if r["section"] == "modulemisc" and not r["quantity"].startswith("comparison.json coldplate")]
    t = ["| Quantity | Repository | This check | Rel. diff. | Verdict |", "|---|---|---|---|---|"]
    for r in sel:
        t.append(f"| {r['quantity']} | {fmt(r['repo'])} | {fmt(r['mine'])} | {frel(r['rel'])} | {r['verdict']} |")
    body["modulemisc_short"] = "\n".join(t) + "\n"
    tabs.update({"BODY_" + k: v for k, v in body.items()})

    # ------------------------------------------------------------------ write
    with open(os.path.join(OUT, "comparison_all.json"), "w") as f:
        json.dump(ROWS, f, indent=1, default=str)
    with open(os.path.join(OUT, "tables.md"), "w") as f:
        for k, v in tabs.items():
            if v is None:
                continue
            f.write(f"\n<!-- {k} -->\n{v}\n")
    # verdict summary by section
    secs = []
    for r in ROWS:
        if r["section"] not in secs:
            secs.append(r["section"])
    summ = ["| Section | Quantities compared | agrees | small difference explained | disagrees | Largest abs. rel. diff. |",
            "|---|---|---|---|---|---|"]
    for s_ in secs:
        n, cnt, mx = summary_line(s_)
        summ.append(f"| {s_} | {n} | {cnt.get('agrees', 0)} | {cnt.get('small difference explained', 0)} | {cnt.get('disagrees', 0)} | {mx:.2e} |")
    with open(os.path.join(OUT, "summary.md"), "w") as f:
        f.write("\n".join(summ) + "\n")
    # descriptive summary for the report
    names = {"law1": "Law 1 lever (D3, D17)", "law1mc": "Law 1 assumed-noise Monte Carlo and D16 hand figure",
             "law2": "Law 2 lever (D4, D17)", "law2terms": "Law 2 term shares, term levers, text values, lifetime curve",
             "law2v1": "Law 2 V1 (D15)", "module": "2x4 module in air (D9 to D11)", "coldplate": "1x8 module on a cold plate (D12)",
             "modulemisc": "V2 uniformity, module sensitivities, comparison.json cold-plate copy",
             "twonode": "two-node reduction at C/2", "cellconv": "resolved cell, converged values (C/2, 2C, base-cooled)",
             "cellall": "resolved cell, all 24 published rows", "cellv3": "resolved cell, V3 (b), (c), (e)",
             "cellbase": "base-cooled cell, comparison.json", "implied": "implied divergence (implied_divergence.csv)",
             "step3": "Step 3 lowest rates and dataset conditions", "d19": "D19 cycle-averaged heat",
             "scores": "prediction scorecard"}
    summ2 = ["| Part | Quantities compared | agrees | small difference explained | disagrees | Largest abs. rel. diff. |",
             "|---|---|---|---|---|---|"]
    for s_ in secs:
        n, cnt, mx = summary_line(s_)
        summ2.append(f"| {names.get(s_, s_)} | {n} | {cnt.get('agrees', 0)} | {cnt.get('small difference explained', 0)} | "
                     f"{cnt.get('disagrees', 0)} | {mx:.2e} |")
    tot = {"n": len(ROWS)}
    for r in ROWS:
        tot[r["verdict"]] = tot.get(r["verdict"], 0) + 1
    summ2.append(f"| **all** | **{tot['n']}** | **{tot.get('agrees', 0)}** | **{tot.get('small difference explained', 0)}** | "
                 f"**{tot.get('disagrees', 0)}** | |")
    tabs["summary"] = "\n".join(summ2) + "\n"
    ds_rows = [r for r in ROWS if r["section"] == "step3" and r["quantity"].startswith("dataset conditions")]
    t = ["| Quantity | Repository | This check | Rel. diff. | Verdict | Note |", "|---|---|---|---|---|---|"]
    for r in ds_rows:
        t.append(f"| {r['quantity']} | {fmt(r['repo'])} | {fmt(r['mine'])} | {frel(r['rel'])} | {r['verdict']} | {r['note']} |")
    tabs["dataset_conditions"] = "\n".join(t) + "\n"
    with open(os.path.join(OUT, "named_tables.json"), "w") as f:
        json.dump({k: v for k, v in tabs.items() if v is not None}, f, indent=1)
    if FILL_MODE:
        # fill a report template read from stdin ({{name}} placeholders) and print it; nothing is written
        import re
        import sys
        txt = sys.stdin.read()
        for key in re.findall(r"\{\{([A-Za-z0-9_]+)\}\}", txt):
            if key not in tabs or tabs[key] is None:
                raise KeyError(f"no table named {key}")
            txt = txt.replace("{{" + key + "}}", tabs[key].rstrip("\n"))
        bad = [ch for ch in ("\u2014", "\u2013", "\u2212") if ch in txt]
        if bad:
            raise ValueError(f"dash characters found in the report: {bad}")
        sys.stdout.write(txt)
        return
    print("\n".join(summ))
    bad = [r for r in ROWS if r["verdict"] == "disagrees"]
    print("\nDISAGREES:")
    for r in bad:
        print(r["section"], "|", r["quantity"], "|", r["repo"], "|", r["mine"], "|", r["rel"], "|", r["note"])
    smallx = [r for r in ROWS if r["verdict"] == "small difference explained"]
    print("\nSMALL (max per section):")
    by = {}
    for r in smallx:
        if r["section"] not in by or abs(r["rel"]) > abs(by[r["section"]]["rel"]):
            by[r["section"]] = r
    for s_, r in by.items():
        print(s_, "|", r["quantity"], "|", r["repo"], "|", r["mine"], "|", f"{r['rel']:+.3e}")


if __name__ == "__main__":
    import sys
    FILL_MODE = "--fill" in sys.argv[1:]
    main()
