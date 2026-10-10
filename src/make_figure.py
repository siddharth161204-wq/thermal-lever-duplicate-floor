#!/usr/bin/env python3
"""The figure: (a) Delta T_floor against the floor f for each law, with the three floors marked and the uncertainty
bands; (b) Delta T_module (lumped, solid) and Delta T_cell (resolved, dashed) against C-rate for h = 5, 10, 25, with
Delta T_floor(3.0 %) and Delta T_floor(8.7 %) of each law as horizontal lines; (c) the implied lifetime divergence
s x Delta T against rate, with the floor lines.

Usage: python3 src/make_figure.py --results results
"""
import argparse, json, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import thermal_lever as TL

INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a897f", "#e1e0d9"
LAW_COLOR = {"law_1_kucinskis2022": "#eb6834", "law_2_schimpe2018": "#2a78d6"}      # categorical slots 2 and 1
LAW_NAME = {"law_1_kucinskis2022": "Law 1: NCA/Si-graphite 21700",
            "law_2_schimpe2018": "Law 2: LiFePO4/graphite 26650"}
H_STYLE = {5.0: ("o", INK), 10.0: ("s", INK2), 25.0: ("^", MUTED)}
RATES = TL.CRATES
TICKS = ["C/16", "C/10", "C/5", "C/2", "1C", "2C"]

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.0, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "axes.spines.top": False, "axes.spines.right": False})


def main(res_dir):
    lev = json.load(open(os.path.join(res_dir, "lever.json")))
    cmp_ = json.load(open(os.path.join(res_dir, "comparison.json")))
    cel = json.load(open(os.path.join(res_dir, "resolved_cell.json")))
    floors = lev["floors"]
    s = {law: {round(e["crate"], 6): e["s_per_K"] for e in lev["levers"][law]} for law in TL.LAWS}
    fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.55), dpi=220)
    fig.subplots_adjust(left=0.065, right=0.99, bottom=0.17, top=0.86, wspace=0.34)

    # (a) |Delta T_floor| against f
    a = ax[0]
    f = np.linspace(0.002, 0.15, 300)
    s1 = s["law_1_kucinskis2022"][0.5]
    s2 = {c: s["law_2_schimpe2018"][round(c, 6)] for c in (1 / 16, 1 / 10, 1 / 5, 1 / 2)}
    mc = cmp_["law1_assumed_noise_mc"]["s_2.5_50_97.5"]
    a.axvspan(100 * 0.0565, 100 * 0.1394, color=GRID, alpha=0.6, lw=0)
    lo1 = np.abs(f / mc[0]); hi1 = np.where(mc[2] < 0, np.abs(f / mc[2]), 1e3)
    a.fill_between(100 * f, lo1, np.minimum(hi1, 1e3), color=LAW_COLOR["law_1_kucinskis2022"], alpha=0.16, lw=0)
    a.plot(100 * f, np.abs(f / s1), color=LAW_COLOR["law_1_kucinskis2022"], lw=1.4)
    a.fill_between(100 * f, f / max(s2.values()), f / min(s2.values()), color=LAW_COLOR["law_2_schimpe2018"], alpha=0.25, lw=0)
    a.plot(100 * f, f / s2[0.5], color=LAW_COLOR["law_2_schimpe2018"], lw=1.4)
    for k, lab in (("mean_3.0pct", "3.0 %"), ("worst_of_8_4.1pct", "4.1 %"), ("p95_8.7pct", "8.7 %")):
        x = 100 * floors[k]
        a.axvline(x, color=INK2, lw=0.6, ls=":")
        a.text(x, 33, lab, rotation=90, va="top", ha="right", color=INK2, fontsize=6)
        for law, sv in (("law_1_kucinskis2022", s1), ("law_2_schimpe2018", s2[0.5])):
            a.plot([x], [abs(floors[k] / sv)], marker="o", ms=3.4, color=LAW_COLOR[law], mec="white", mew=0.6, zorder=5)
    a.set_yscale("log"); a.set_ylim(0.1, 40); a.set_xlim(0, 15)
    a.set_xlabel("floor f, lifetime difference (%)"); a.set_ylabel("|Delta T_floor| (K)")
    a.text(14.6, abs(0.145 / s1) * 0.82, "Law 1", color=INK, ha="right", va="top", fontsize=6.3)
    a.text(14.6, 0.145 / s2[0.5] * 0.8, "Law 2", color=INK, ha="right", va="top", fontsize=6.3)
    a.set_title("(a) temperature difference the floor needs", loc="left", fontsize=7.2, color=INK)
    a.grid(axis="y", color=GRID, lw=0.5)

    # (b) spreads against rate
    b = ax[1]
    x = np.array(RATES)
    mod = {(round(r["crate"], 6), r["h"]): r for r in lev["module"]}
    cell = {(round(r["crate"], 6), r["h"]): r for r in cel["rows"] if r["base_h"] is None}
    b.axvspan(1 / 16 * 0.85, 0.5 * 1.15, color=GRID, alpha=0.45, lw=0)
    for h, (mk, col) in H_STYLE.items():
        b.plot(x, [mod[(round(c, 6), h)]["dT_module_K"] for c in RATES], color=col, lw=1.2, marker=mk, ms=3.2)
        b.plot(x, [cell[(round(c, 6), h)]["dT_cell_K"] for c in RATES], color=col, lw=1.0, ls="--", marker=mk, ms=2.6, mfc="white")
    for law in TL.LAWS:
        for k, ls in (("mean_3.0pct", (0, (1, 1.2))), ("p95_8.7pct", "-")):
            v = abs(floors[k] / s[law][0.5])
            b.axhline(v, color=LAW_COLOR[law], lw=0.9, ls=ls)
            b.text(2.6, v, ("L1 " if law.startswith("law_1") else "L2 ") + ("3.0 %" if k.startswith("mean") else "8.7 %"),
                   color=INK2, fontsize=5.6, va="center", ha="left")
    b.set_xscale("log"); b.set_yscale("log")
    b.set_xticks(x); b.set_xticklabels(TICKS); b.minorticks_off()
    b.set_xlim(1 / 16 * 0.8, 2.0 * 1.25); b.set_ylim(1e-3, 30)
    b.set_xlabel("C-rate (average)"); b.set_ylabel("spread (K)")
    b.set_title("(b) spread a module or a cell offers", loc="left", fontsize=7.2, color=INK)
    b.grid(axis="y", color=GRID, lw=0.5)
    hl = [Line2D([], [], color=INK2, lw=1.2, label="module, lumped"), Line2D([], [], color=INK2, lw=1.0, ls="--", label="cell, resolved")]
    hl += [Line2D([], [], color=col, marker=mk, lw=0, ms=3.2, label=f"h = {h:g}") for h, (mk, col) in H_STYLE.items()]
    b.legend(handles=hl, loc="upper left", fontsize=5.6, frameon=False, handlelength=1.8, borderaxespad=0.2)

    # (c) implied divergence at h = 10 (lines) with the h = 5 to 25 range (band)
    cax = ax[2]
    cax.axvspan(1 / 16 * 0.85, 0.5 * 1.15, color=GRID, alpha=0.45, lw=0)
    for law in TL.LAWS:
        sv = np.array([abs(s[law][round(c, 6)]) for c in RATES])
        for kind, ls in (("module", "-"), ("cell", "--")):
            src = mod if kind == "module" else cell
            key = "dT_module_K" if kind == "module" else "dT_cell_K"
            vals = np.array([[src[(round(c, 6), h)][key] for c in RATES] for h in TL.H_LIST])
            cax.fill_between(x, 100 * sv * vals.min(axis=0), 100 * sv * vals.max(axis=0), color=LAW_COLOR[law], alpha=0.12, lw=0)
            cax.plot(x, 100 * sv * vals[1], color=LAW_COLOR[law], lw=1.2 if kind == "module" else 1.0, ls=ls)
    for k, lab in (("mean_3.0pct", "3.0 %"), ("worst_of_8_4.1pct", "4.1 %"), ("p95_8.7pct", "8.7 %")):
        cax.axhline(100 * floors[k], color=INK2, lw=0.6, ls=":")
        cax.text(1 / 16 * 0.84, 100 * floors[k] * 1.08, lab, color=INK2, fontsize=5.6, va="bottom", ha="left")
    cax.set_xscale("log"); cax.set_yscale("log")
    cax.set_xticks(x); cax.set_xticklabels(TICKS); cax.minorticks_off()
    cax.set_xlim(1 / 16 * 0.8, 2.0 * 1.25); cax.set_ylim(1e-3, 100)
    cax.set_xlabel("C-rate (average)"); cax.set_ylabel("implied lifetime divergence (%)")
    cax.set_title("(c) s x spread against the floor", loc="left", fontsize=7.2, color=INK)
    cax.grid(axis="y", color=GRID, lw=0.5)
    hl = [Line2D([], [], color=LAW_COLOR[l], lw=1.2, label=LAW_NAME[l]) for l in TL.LAWS]
    hl += [Line2D([], [], color=INK2, lw=1.2, label="module (h = 10; band h = 5 to 25)"), Line2D([], [], color=INK2, lw=1.0, ls="--", label="within cell")]
    cax.legend(handles=hl, loc="lower right", fontsize=5.4, frameon=False, handlelength=1.8, borderaxespad=0.2)

    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(res_dir, f"thermal_lever_figure.{ext}"))
    print("figure written")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--results", default="results"); a = ap.parse_args()
    main(a.results)
