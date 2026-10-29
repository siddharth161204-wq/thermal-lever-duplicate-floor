#!/usr/bin/env python3
"""Stage 5 figure, two panels (DECISIONS.md D76). (a) d against p for the 106 scored KIT pairs, warm and cold band,
with the line beta = 1 and the fitted line d = beta p inside its 95 % interval; the pairs tied within their check
spacing drawn as rings. (b) The calculated detectable |s| against the number of pairs at the observed dT spread of
KIT and of Sandia, under the KIT and Geslin floors, each dataset's available pairs marked at its calculated (filled)
and simulated (ring) value, the band medians as horizontal marks and the 24 KIT between-set |s| as ticks at the edge.

Usage: python3 src/stage5/s5_figure.py
"""
import csv, json, os, sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s5_inputs as A                   # noqa: E402

INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a897f", "#e1e0d9"
WARM, COLD, AQUA = "#eb6834", "#2a78d6", "#1baf7a"
DS_COL = {"KIT": INK, "Sandia": AQUA}
FL_LS = {"KIT": "-", "Geslin": (0, (4, 1.6))}

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.4, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
                     "xtick.labelsize": 5.8, "ytick.labelsize": 5.8, "axes.spines.top": False, "axes.spines.right": False})


def f(x):
    return float(x) if x not in ("", None) else None


def panel_a(ax, pairs, beta):
    K = beta["dataset"]["KIT"]
    b, (lo, hi) = K["beta"], K["beta_boot"]["ci95"]
    kit = [r for r in pairs if r["dataset"] == "KIT" and r["scored"] == "True"]
    p = np.array([f(r["p"]) for r in kit])
    x0, x1 = p.min() - 0.012, p.max() + 0.008
    x = np.linspace(x0, x1, 200)
    ax.fill_between(x, np.minimum(lo * x, hi * x), np.maximum(lo * x, hi * x), color=GRID, lw=0, zorder=1)
    ax.axhline(0, color=GRID, lw=0.5, zorder=0)
    ax.axvline(0, color=GRID, lw=0.5, zorder=0)
    ax.plot(x, x, color=INK2, lw=0.9, ls=(0, (3, 2)), zorder=3)
    ax.plot(x, b * x, color=INK, lw=1.1, zorder=4)
    for band, col in (("warm", WARM), ("cold", COLD)):
        rr = [r for r in kit if r["band"] == band]
        tied = [r for r in rr if r["tie"] == "True"]
        free = [r for r in rr if r["tie"] != "True"]
        ax.plot([f(r["p"]) for r in tied], [f(r["d"]) for r in tied], ls="none", marker="o", ms=2.9, mfc="white", mec=col,
                mew=0.7, zorder=5)
        ax.plot([f(r["p"]) for r in free], [f(r["d"]) for r in free], ls="none", marker="o", ms=3.3, mfc=col, mec="white",
                mew=0.5, zorder=6)
    d = np.array([f(r["d"]) for r in kit])
    ym = 1.1 * np.abs(d).max()
    ax.set_xlim(x0, x1)
    ax.set_ylim(-ym, ym)
    ax.set_xlabel("predicted by the between-set lever, p = -s dT")
    ax.set_ylabel("observed within the set, d = ln(L_warm / L_cool)")
    nw = sum(r["band"] == "warm" for r in kit)
    nc = sum(r["band"] == "cold" for r in kit)
    nt = sum(r["tie"] == "True" for r in kit)
    hs = [Line2D([], [], ls="none", marker="o", ms=3.3, mfc=WARM, mec="white", label=f"warm band ({nw} pairs)"),
          Line2D([], [], ls="none", marker="o", ms=3.3, mfc=COLD, mec="white", label=f"cold band ({nc} pairs)"),
          Line2D([], [], ls="none", marker="o", ms=2.9, mfc="white", mec=INK2, mew=0.7, label=f"ring: tied ({nt} of {len(kit)})"),
          Line2D([], [], color=INK2, lw=0.9, ls=(0, (3, 2)), label="beta = 1"),
          Line2D([], [], color=INK, lw=1.1, label=f"fitted, beta = {b:+.2f}"),
          Patch(fc=GRID, ec="none", label=f"95 %: {lo:+.2f} to {hi:+.2f}")]
    ax.legend(handles=hs, frameon=False, fontsize=4.9, loc="lower right", handlelength=1.8, handletextpad=0.4,
              borderaxespad=0.3)
    ax.set_title("(a) KIT replicate pairs: observed against predicted", loc="left", fontsize=7.0, color=INK)


def panel_b(ax, power, levers):
    cur = power["curves"]
    for ds in ("KIT", "Sandia"):
        for fl in ("KIT", "Geslin"):
            cc = [c for c in cur if c["dataset"] == ds and c["floor"] == fl]
            ax.plot([c["pairs"] for c in cc], [100 * c["detectable_s_per_K"] for c in cc], color=DS_COL[ds], ls=FL_LS[fl],
                    lw=1.1 if ds == "KIT" else 1.0, zorder=3)
            c1 = [x for x in power["C1"] if x["registered"] and x["case"] == ds and x["floor"] == fl][0]
            ax.plot([c1["pairs"]], [100 * c1["detectable_s_per_K"]], ls="none", marker="o", ms=3.4, mfc=DS_COL[ds],
                    mec="white", mew=0.5, zorder=6)
            ax.plot([c1["pairs"]], [100 * c1["sim_lognormal_detectable_s_per_K"]], ls="none", marker="o", ms=5.0,
                    mfc="none", mec=DS_COL[ds], mew=0.6, zorder=5)
    med = power["band_medians"]
    for band, col in (("warm", WARM), ("cold", COLD)):
        ax.axhline(100 * med[band]["abs_s_per_K"], color=col, lw=0.9, zorder=2)
    y0, y1 = 0.3, 100.0
    for band, col in (("warm", WARM), ("cold", COLD)):
        for s in levers[band]:
            if y0 <= 100 * abs(s) <= y1:
                ax.plot([1.0, 1.035], [100 * abs(s)] * 2, color=col, lw=0.8, transform=ax.get_yaxis_transform(),
                        clip_on=False, zorder=7)
    k = [x for x in power["C1"] if x["registered"] and x["case"] == "KIT" and x["floor"] == "KIT"][0]
    s_ = [x for x in power["C1"] if x["registered"] and x["case"] == "Sandia" and x["floor"] == "KIT"][0]
    ax.annotate(f"KIT pairs, dT {k['dT_rms_K']:.2f} K rms", (k["pairs"], 100 * k["detectable_s_per_K"]), xytext=(5, 6),
                textcoords="offset points", fontsize=4.9, color=INK)
    s_g = [x for x in power["C1"] if x["registered"] and x["case"] == "Sandia" and x["floor"] == "Geslin"][0]
    ax.annotate(f"Sandia pairs, dT {s_['dT_rms_K']:.2f} K rms", (s_g["pairs"], 100 * s_g["detectable_s_per_K"]),
                xytext=(-7, -9), textcoords="offset points", fontsize=4.9, color=INK, ha="right")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(3, 3000)
    ax.set_xticks([3, 10, 30, 100, 300, 1000, 3000])
    ax.set_xticklabels(["3", "10", "30", "100", "300", "1000", "3000"])
    ax.set_ylim(y0, y1)
    ax.set_yticks([0.3, 1, 3, 10, 30, 100])
    ax.set_yticklabels(["0.3", "1", "3", "10", "30", "100"])
    ax.minorticks_off()
    ax.set_xlabel("replicate pairs at the observed dT spread")
    ax.set_ylabel("smallest |s| detectable at 80 % power (% per K)")
    hs = [Line2D([], [], color=INK2, lw=1.1, label="against the KIT floor"),
          Line2D([], [], color=INK2, lw=1.1, ls=FL_LS["Geslin"], label="against the Geslin floor"),
          Line2D([], [], ls="none", marker="o", ms=3.4, mfc=INK2, mec="white", label="the pairs available: calculated"),
          Line2D([], [], ls="none", marker="o", ms=5.0, mfc="none", mec=INK2, mew=0.6, label="simulated, log-normal floor"),
          Line2D([], [], color=WARM, lw=0.9, label=f"warm band median, {100 * med['warm']['abs_s_per_K']:.2f} % per K"),
          Line2D([], [], color=COLD, lw=0.9, label=f"cold band median, {100 * med['cold']['abs_s_per_K']:.2f} % per K"),
          Line2D([], [], color=INK2, lw=0.8, label="ticks at the edge: KIT between-set |s|, by band")]
    ax.legend(handles=hs, frameon=False, fontsize=4.9, loc="upper right", bbox_to_anchor=(0.97, 1.0), handlelength=2.0,
              handletextpad=0.4)
    ax.set_title("(b) what the replicates could detect", loc="left", fontsize=7.0, color=INK)


def main():
    R = os.path.join(A.ROOT, "results")
    pairs = list(csv.DictReader(open(os.path.join(R, "stage5_pairs.csv"))))
    beta = json.load(open(os.path.join(R, "stage5_beta.json")))
    power = json.load(open(os.path.join(R, "stage5_power.json")))
    lev = A.kit_levers()
    levers = {b: [lev[p][b]["s"] for p in lev] for b in ("warm", "cold")}
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.75), dpi=220)
    fig.subplots_adjust(left=0.07, right=0.965, bottom=0.16, top=0.9, wspace=0.25)
    panel_a(ax[0], pairs, beta)
    panel_b(ax[1], power, levers)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(R, f"stage5_figure.{ext}"))
    print("written", os.path.join(R, "stage5_figure.png"))


if __name__ == "__main__":
    main()
