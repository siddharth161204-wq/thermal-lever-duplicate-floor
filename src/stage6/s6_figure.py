#!/usr/bin/env python3
"""Stage 6 figure, three panels (DECISIONS.md D95 (4)).

(a) Within each set of two or more cells, every cell's deviation of ln L from its set's mean against its deviation of
T, the two papers in two colours, the cells the contact rule flags as rings; the fitted within-set line (slope -s_w)
inside its 95 % interval, and Sandia's levers for this cell (Stage 3) as thin lines through the origin.
(b) Every set's mean of ln L against its mean of T, each as a deviation from its batch's weighted means, sized by its
number of cells and coloured by batch; the fitted between-policy line (slope -s_b) and, dashed, the within-set lever.
(c) The smallest |s| detectable at 80 % power: Stage 6 (each paper, pooled, and pooled without the flagged cells)
beside Stage 5 (KIT, Sandia, pooled), under the KIT floor (filled) and the Geslin floor (open), with Sandia's levers
and this stage's |s_w| marked.

Usage: python3 -I src/stage6/s6_figure.py [results directory]
Writes stage6_figure.pdf and stage6_figure.png in the results directory.
"""
import csv, json, math, os, sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s6_common as C                   # noqa: E402

INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a897f", "#e1e0d9"
PAPER_COL = {"Severson": "#2a78d6", "Attia": "#eb6834"}
BATCH_COL = {"2017-05-12": "#2a78d6", "2017-06-30": "#7fb2ec", "2018-04-12": "#1d4f91", "2019-01-24": "#eb6834"}
SANDIA = [("registered 1C", "15 to 25 C", "Sandia 15 to 25 C, 1C"), ("secondary 2C", "15 to 25 C", "Sandia 15 to 25 C, 2C"),
          ("secondary 2C", "25 to 35 C", "Sandia 25 to 35 C, 2C")]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.4, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "xtick.labelsize": 5.8, "ytick.labelsize": 5.8,
                     "axes.spines.top": False, "axes.spines.right": False})


def cells(res_dir):
    out = []
    for r in csv.DictReader(open(os.path.join(res_dir, "stage6_cells.csv"))):
        if r["enters"] != "True" or r.get("cycle_life", "") in ("", "nan"):
            continue
        out.append({"set": r["set"], "batch": r["batch"], "paper": r["paper"], "T": float(r["T_reg"]),
                    "y": math.log(float(r["cycle_life"])), "flag": r["contact_flag"] == "True"})
    return out


def panel_a(ax, cl, W, lev):
    by = {}
    for c in cl:
        by.setdefault(c["set"], []).append(c)
    xs, ys = [], []
    for m in by.values():
        if len(m) < 2:
            continue
        Tm, ym = np.mean([c["T"] for c in m]), np.mean([c["y"] for c in m])
        for c in m:
            c["dx"], c["dy"] = c["T"] - Tm, c["y"] - ym
            xs.append(c["dx"]); ys.append(c["dy"])
    span = 1.08 * max(abs(min(xs)), abs(max(xs)))
    x = np.linspace(-span, span, 50)
    sw = W["B1"]["s_w_per_K"]
    lo, hi = sw["ci95"]
    ax.fill_between(x, np.minimum(-lo * x, -hi * x), np.maximum(-lo * x, -hi * x), color=GRID, lw=0, zorder=1)
    ax.axhline(0, color=GRID, lw=0.5, zorder=0)
    ax.axvline(0, color=GRID, lw=0.5, zorder=0)
    for (series, sec, lab), ls in zip(SANDIA, ((0, (1, 1.2)), (0, (3, 1.5)), (0, (5, 1.5, 1, 1.5)))):
        ax.plot(x, -lev[(series, sec)][0] * x, color=MUTED, lw=0.7, ls=ls, zorder=2)
    ax.plot(x, -sw["value"] * x, color=INK, lw=1.1, zorder=4)
    for c in [c for c in cl if "dx" in c]:
        col = PAPER_COL[c["paper"]]
        if c["flag"]:
            ax.plot(c["dx"], c["dy"], ls="none", marker="o", ms=2.8, mfc="white", mec=col, mew=0.6, zorder=5)
        else:
            ax.plot(c["dx"], c["dy"], ls="none", marker="o", ms=3.0, mfc=col, mec="white", mew=0.4, zorder=6)
    ax.set_xlim(-span, span)
    ym = 1.1 * max(abs(min(ys)), abs(max(ys)))
    ax.set_ylim(-ym, ym)
    ax.set_xlabel("T minus the set's mean (K)")
    ax.set_ylabel("ln L minus the set's mean")
    ax.set_title("a  Within a set", loc="left", fontsize=6.8, fontweight="bold")
    h = [Line2D([], [], ls="none", marker="o", ms=3, mfc=PAPER_COL["Severson"], mec="white", label="Severson et al. cells"),
         Line2D([], [], ls="none", marker="o", ms=3, mfc=PAPER_COL["Attia"], mec="white", label="Attia et al. cells"),
         Line2D([], [], ls="none", marker="o", ms=3, mfc="white", mec=INK2, label="flagged by the contact rule"),
         Line2D([], [], color=INK, lw=1.1, label=f"fit, s_w = {100 * sw['value']:+.2f} % per K"),
         Line2D([], [], color=MUTED, lw=0.7, ls=(0, (1, 1.2)), label="Sandia's levers, same cell")]
    ax.legend(handles=h, frameon=False, fontsize=5.2, loc="lower left", handlelength=1.6)


def panel_b(ax, cl, Bt, W):
    by = {}
    for c in cl:
        by.setdefault(c["set"], []).append(c)
    sets = [{"batch": m[0]["batch"], "n": len(m), "T": np.mean([c["T"] for c in m]), "y": np.mean([c["y"] for c in m])}
            for m in by.values()]
    for b in sorted(set(s["batch"] for s in sets)):
        ss = [s for s in sets if s["batch"] == b]
        n = np.array([s["n"] for s in ss], dtype=float)
        Tm = (n * np.array([s["T"] for s in ss])).sum() / n.sum()
        ym = (n * np.array([s["y"] for s in ss])).sum() / n.sum()
        for s in ss:
            s["dx"], s["dy"] = s["T"] - Tm, s["y"] - ym
    xs = np.array([s["dx"] for s in sets])
    span = 1.08 * np.abs(xs).max()
    x = np.linspace(-span, span, 50)
    ax.axhline(0, color=GRID, lw=0.5, zorder=0)
    ax.axvline(0, color=GRID, lw=0.5, zorder=0)
    sb = Bt["C1"]["s_b_per_K"]
    lo, hi = sb["ci95"]
    ax.fill_between(x, np.minimum(-lo * x, -hi * x), np.maximum(-lo * x, -hi * x), color=GRID, lw=0, zorder=1)
    ax.plot(x, -sb["value"] * x, color=INK, lw=1.1, zorder=4)
    ax.plot(x, -W["B1"]["s_w_per_K"]["value"] * x, color=INK2, lw=0.8, ls=(0, (3, 2)), zorder=3)
    for s in sets:
        ax.scatter(s["dx"], s["dy"], s=6 + 5 * s["n"], c=BATCH_COL.get(s["batch"], INK2), edgecolors="white",
                   linewidths=0.4, zorder=5)
    ax.set_xlim(-span, span)
    ax.set_xlabel("set's mean T minus its batch's (K)")
    ax.set_ylabel("set's mean ln L minus its batch's")
    ax.set_title("b  Between policies", loc="left", fontsize=6.8, fontweight="bold")
    h = [Line2D([], [], ls="none", marker="o", ms=3.4, mfc=c, mec="white", label=b) for b, c in BATCH_COL.items()
         if any(s["batch"] == b for s in sets)]
    h += [Line2D([], [], color=INK, lw=1.1, label=f"fit, s_b = {100 * sb['value']:+.2f} % per K"),
          Line2D([], [], color=INK2, lw=0.8, ls=(0, (3, 2)), label="the within-set lever")]
    ax.legend(handles=h, frameon=False, fontsize=5.2, loc="lower left", handlelength=1.6)


def panel_c(ax, P, s5, lev, sw):
    rows = []
    for case, lab in (("KIT", "Stage 5, KIT replicates"), ("Sandia", "Stage 5, Sandia replicates"),
                      ("pooled", "Stage 5, pooled")):
        rows.append((lab, s5[(case, "KIT")]["detectable_s_per_K"], s5[(case, "Geslin")]["detectable_s_per_K"]))
    get = lambda scope, floor: next(p for p in P if p["scope"] == scope and p["floor"] == floor)
    for scope, lab in (("Severson", "Stage 6, Severson et al."), ("Attia", "Stage 6, Attia et al."),
                       ("pooled", "Stage 6, pooled")):
        rows.append((lab, get(scope, "KIT")["detectable_s_per_K"], get(scope, "Geslin")["detectable_s_per_K"]))
    wf = get("pooled", "KIT").get("without_flagged")
    if wf:
        rows.append(("Stage 6, pooled, without flagged cells", wf["detectable_s_per_K"],
                     get("pooled", "Geslin")["without_flagged"]["detectable_s_per_K"]))
    y = np.arange(len(rows))[::-1]
    for (lab, k, g), yy in zip(rows, y):
        ax.plot([100 * g, 100 * k], [yy, yy], color=GRID, lw=1.2, zorder=1)
        ax.plot(100 * k, yy, ls="none", marker="o", ms=3.6, mfc=INK, mec=INK, zorder=3)
        ax.plot(100 * g, yy, ls="none", marker="o", ms=3.6, mfc="white", mec=INK, mew=0.7, zorder=3)
    for (series, sec, lab), ls in zip(SANDIA, ((0, (1, 1.2)), (0, (3, 1.5)), (0, (5, 1.5, 1, 1.5)))):
        ax.axvline(100 * lev[(series, sec)][0], color=MUTED, lw=0.7, ls=ls, zorder=0)
    ax.axvline(100 * abs(sw), color=PAPER_COL["Attia"], lw=0.8, zorder=0)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=5.4)
    ax.set_xlabel("smallest |s| detectable at 80 % power (% per K)")
    ax.set_title("c  What the data could detect", loc="left", fontsize=6.8, fontweight="bold")
    h = [Line2D([], [], ls="none", marker="o", ms=3.4, mfc=INK, mec=INK, label="KIT floor"),
         Line2D([], [], ls="none", marker="o", ms=3.4, mfc="white", mec=INK, label="Geslin floor"),
         Line2D([], [], color=MUTED, lw=0.7, ls=(0, (1, 1.2)), label="Sandia's levers"),
         Line2D([], [], color=PAPER_COL["Attia"], lw=0.8, label="|s_w| measured here")]
    ax.legend(handles=h, frameon=False, fontsize=5.2, loc="upper right", handlelength=1.6)


def main(res_dir):
    arch = C.archive()
    W = json.load(open(os.path.join(res_dir, "stage6_within.json")))["registered"]["pooled"]
    Bt = json.load(open(os.path.join(res_dir, "stage6_between.json")))["registered"]["pooled"]
    P = json.load(open(os.path.join(res_dir, "stage6_power.json")))["cases"]
    cl = cells(res_dir)
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.55), gridspec_kw={"width_ratios": [1, 1, 1.12]})
    panel_a(axs[0], cl, W, arch["levers"])
    panel_b(axs[1], cl, Bt, W)
    panel_c(axs[2], P, arch["stage5"], arch["levers"], W["B1"]["s_w_per_K"]["value"])
    fig.tight_layout(w_pad=1.2)
    fig.savefig(os.path.join(res_dir, "stage6_figure.pdf"))
    fig.savefig(os.path.join(res_dir, "stage6_figure.png"), dpi=300)
    print("wrote", os.path.join(res_dir, "stage6_figure.pdf"))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else C.path("results"))
