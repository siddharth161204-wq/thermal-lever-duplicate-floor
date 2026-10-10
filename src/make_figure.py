#!/usr/bin/env python3
"""The figure: (a) Delta T_floor against the floor f for each law, with the three floors marked and the uncertainty
bands; (b) Delta T_module (lumped, solid) and Delta T_cell (resolved, dashed) against C-rate for h = 5, 10, 25, with
Delta T_floor(3.0 %) and Delta T_floor(8.7 %) of each law as horizontal lines; (c) the implied lifetime divergence
s x Delta T against rate, with the floor lines.

Usage: python3 src/make_figure.py --results results
"""
import argparse, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import thermal_lever as TL

INK, INK2, MUTED, GRID, SHADE = "#0b0b0b", "#52514e", "#8a897f", "#e1e0d9", "#f0efec"
L1, L2 = "law_1_kucinskis2022", "law_2_schimpe2018"
COL = {L1: "#eb6834", L2: "#2a78d6"}                       # categorical slots 2 and 1 of the reference palette
H_STYLE = {5.0: ("o", INK), 10.0: ("s", INK2), 25.0: ("^", MUTED)}
RATES = TL.CRATES
TICKS = ["C/16", "C/10", "C/5", "C/2", "1C", "2C"]

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.6, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "xtick.labelsize": 6.0, "ytick.labelsize": 6.0,
                     "axes.spines.top": False, "axes.spines.right": False})


def rate_axis(ax, ylim):
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xticks(RATES); ax.set_xticklabels(TICKS); ax.minorticks_off()
    ax.set_xlim(1 / 16 * 0.78, 2.0 * 1.3); ax.set_ylim(*ylim)
    ax.axvspan(1 / 16 * 0.86, 0.5 * 1.16, color=SHADE, lw=0, zorder=0)
    ax.text(np.sqrt(0.2 * 0.5), ylim[0] * 1.35, "the dataset's rates", ha="center", va="bottom", fontsize=5.4, color=INK2)
    ax.grid(axis="y", color=GRID, lw=0.5)
    ax.set_xlabel("average C-rate")


def main(res_dir):
    lev = json.load(open(os.path.join(res_dir, "lever.json")))
    cmp_ = json.load(open(os.path.join(res_dir, "comparison.json")))
    cel = json.load(open(os.path.join(res_dir, "resolved_cell.json")))
    fl = lev["floors"]
    s = {law: {round(e["crate"], 6): e["s_per_K"] for e in lev["levers"][law]} for law in TL.LAWS}
    mod = {(round(r["crate"], 6), r["h"]): r for r in lev["module"]}
    cell = {(round(r["crate"], 6), r["h"]): r for r in cel["rows"] if r["base_h"] is None}

    fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.75), dpi=220)
    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.28, top=0.9, wspace=0.42)

    # (a) |Delta T_floor| against f
    a = ax[0]
    f = np.linspace(0.002, 0.15, 400)
    s1 = s[L1][0.5]
    s2 = {c: s[L2][round(c, 6)] for c in (1 / 16, 1 / 10, 1 / 5, 1 / 2)}
    mc = cmp_["law1_assumed_noise_mc"]["s_2.5_50_97.5"]
    a.axvspan(100 * fl_bs(lev)[0], 100 * fl_bs(lev)[1], color=SHADE, lw=0, zorder=0)
    a.text(11.85, 0.13, "95th percentile,\nbootstrap 5.7 to 13.9 %", fontsize=5.4, color=INK2, ha="center", va="bottom")
    a.fill_between(100 * f, np.abs(f / mc[0]), np.abs(f / mc[2]), color=COL[L1], alpha=0.18, lw=0)
    a.plot(100 * f, np.abs(f / s1), color=COL[L1], lw=1.4)
    a.fill_between(100 * f, f / max(s2.values()), f / min(s2.values()), color=COL[L2], alpha=0.3, lw=0)
    a.plot(100 * f, f / s2[0.5], color=COL[L2], lw=1.4)
    for k, lab in (("mean_3.0pct", "3.0 %"), ("worst_of_8_4.1pct", "4.1 %"), ("p95_8.7pct", "8.7 %")):
        x = 100 * fl[k]
        a.axvline(x, color=INK2, lw=0.6, ls=":", zorder=1)
        a.text(x - 0.15, 34, lab, rotation=90, va="top", ha="right", color=INK2, fontsize=5.6)
        for law, sv in ((L1, s1), (L2, s2[0.5])):
            a.plot([x], [abs(fl[k] / sv)], marker="o", ms=3.6, color=COL[law], mec="white", mew=0.6, zorder=5)
    a.text(14.7, abs(0.147 / s1) * 1.12, "Law 1", color=INK, ha="right", va="bottom", fontsize=6.2)
    a.text(14.7, 0.147 / s2[0.5] * 1.12, "Law 2", color=INK, ha="right", va="bottom", fontsize=6.2)
    a.set_yscale("log"); a.set_ylim(0.1, 40); a.set_xlim(0, 15)
    a.set_xlabel("floor f, lifetime difference (%)"); a.set_ylabel("|ΔT$_{floor}$| (K)")
    a.set_title("(a) the temperature difference the floor needs", loc="left", fontsize=6.9, color=INK)
    a.grid(axis="y", color=GRID, lw=0.5)

    # (b) spreads against rate, with the floor temperature differences
    b = ax[1]
    rate_axis(b, (1e-3, 40))
    for h, (mk, col) in H_STYLE.items():
        b.plot(RATES, [mod[(round(c, 6), h)]["dT_module_K"] for c in RATES], color=col, lw=1.1, marker=mk, ms=3.0, zorder=3)
        b.plot(RATES, [cell[(round(c, 6), h)]["dT_cell_K"] for c in RATES], color=col, lw=0.9, ls="--", marker=mk, ms=2.6,
               mfc="white", zorder=3)
    for law, tag in ((L1, "Law 1"), (L2, "Law 2")):
        for k, ls, lab in (("mean_3.0pct", (0, (1, 1.3)), "3.0 %"), ("p95_8.7pct", "-", "8.7 %")):
            v = abs(fl[k] / s[law][0.5])
            b.axhline(v, color=COL[law], lw=0.9, ls=ls, zorder=2)
            below = (law == L2)
            b.text(1 / 16 * 0.82, v * (0.88 if below else 1.12), f"{tag}, {lab}", color=INK2, fontsize=5.4,
                   va="top" if below else "bottom", ha="left")
    b.set_ylabel("temperature difference (K)")
    b.set_title("(b) what a module and a cell offer", loc="left", fontsize=6.9, color=INK)

    # (c) implied divergence at h = 10 (lines) with the h = 5 to 25 range (band)
    c = ax[2]
    rate_axis(c, (1e-3, 100))
    for law in TL.LAWS:
        sv = np.array([abs(s[law][round(r, 6)]) for r in RATES])
        for kind, ls, lw in (("module", "-", 1.3), ("cell", "--", 1.0)):
            src, key = (mod, "dT_module_K") if kind == "module" else (cell, "dT_cell_K")
            vals = np.array([[src[(round(r, 6), h)][key] for r in RATES] for h in TL.H_LIST])
            c.fill_between(RATES, 100 * sv * vals.min(axis=0), 100 * sv * vals.max(axis=0), color=COL[law], alpha=0.14, lw=0)
            c.plot(RATES, 100 * sv * vals[1], color=COL[law], lw=lw, ls=ls, zorder=3)
    for k, lab, va, fac in (("mean_3.0pct", "3.0 % mean", "top", 0.9), ("worst_of_8_4.1pct", "4.1 % worst of 8", "bottom", 1.05),
                            ("p95_8.7pct", "8.7 % 95th percentile", "bottom", 1.08)):
        c.axhline(100 * fl[k], color=INK2, lw=0.6, ls=":", zorder=2)
        c.text(1 / 16 * 0.82, 100 * fl[k] * fac, lab, color=INK2, fontsize=5.4, va=va, ha="left")
    c.set_ylabel("implied lifetime divergence (%)")
    c.set_title("(c) s × temperature difference", loc="left", fontsize=6.9, color=INK)

    handles = [Line2D([], [], color=COL[L1], lw=1.4, label="Law 1: NCA+LNO / Si-graphite 21700, C/2 charge"),
               Line2D([], [], color=COL[L2], lw=1.4, label="Law 2: LiFePO4 / graphite 26650, full model"),
               Line2D([], [], color=INK2, lw=1.1, label="eight-cell module, lumped"),
               Line2D([], [], color=INK2, lw=0.9, ls="--", label="one cell, resolved")]
    handles += [Line2D([], [], color=col, marker=mk, lw=0, ms=3.2, label=f"h = {h:g} W/(m\u00b2 K)") for h, (mk, col) in H_STYLE.items()]
    handles += [Patch(color=COL[L1], alpha=0.25, lw=0, label="bands: (a) Law 1 assumed noise, Law 2 C/16 to C/2; (c) h = 5 to 25")]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=5.5, frameon=False, handlelength=2.0,
               columnspacing=1.2, bbox_to_anchor=(0.5, 0.0))
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(res_dir, f"thermal_lever_figure.{ext}"))
    print("figure written")


def fl_bs(lev):
    with open(os.path.join(TL.INPUTS, "floor.json")) as f:
        return json.load(f)["p95_bootstrap_95"]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--results", default="results"); a = ap.parse_args()
    main(a.results)
