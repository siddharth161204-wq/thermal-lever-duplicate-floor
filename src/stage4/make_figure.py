#!/usr/bin/env python3
"""Stage 4 figure, two panels. (a) The median implied divergence |s| x spread over the band's levers against the
average rate, at 35, 25 and 10 C ambient and h = 5, 10 and 25 W/(m2 K), with the band's 95th percentile at h = 10
and the four floor lines. (b) The duplicate pairs needed to detect the thermal share against dT, in the warm (35 C)
and cold (10 C) chambers under their central levers, against the Geslin and KIT floors in the registered log-normal
shape, with the cases of C2 marked.

Usage: python3 src/stage4/make_figure.py --results results
"""
import argparse, csv, json, os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a897f", "#e1e0d9"
AMB_COL = {35.0: "#eb6834", 25.0: "#1baf7a", 10.0: "#2a78d6"}
H_STYLE = {5.0: (0, (4, 1.6)), 10.0: "-", 25.0: (0, (1, 1.4))}
H_MARK = {5.0: "o", 10.0: "s", 25.0: "^"}
RATES = [0.25, 0.5, 1.0, 2.0]
RATE_LAB = ["C/4", "C/2", "1C", "2C"]

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.4, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
                     "xtick.labelsize": 5.8, "ytick.labelsize": 5.8, "axes.spines.top": False, "axes.spines.right": False})


def panel_a(ax, mp):
    cases = mp["cases"]
    fl = mp["floors"]
    for T in (35.0, 25.0, 10.0):
        for h in (5.0, 10.0, 25.0):
            cc = sorted([c for c in cases if c["ambient_C"] == T and c["h"] == h], key=lambda c: c["crate"])
            y = [100 * c["median"] for c in cc]
            ax.plot(RATES, y, ls=H_STYLE[h], color=AMB_COL[T], lw=1.1 if h == 10.0 else 0.9, zorder=4,
                    marker=H_MARK[h], ms=2.6, mfc=AMB_COL[T] if h == 10.0 else "white", mec=AMB_COL[T], mew=0.6)
        cc = sorted([c for c in cases if c["ambient_C"] == T and c["h"] == 10.0], key=lambda c: c["crate"])
        ax.plot(RATES, [100 * c["p95"] for c in cc], color=AMB_COL[T], lw=0.5, alpha=0.55, zorder=3)
    lines = [("Geslin", "mean", "Geslin mean 3.0 %", -1), ("Geslin", "p95", "Geslin 95th percentile 8.7 %", +1),
             ("KIT", "mean", "KIT mean 7.3 %", -1), ("KIT", "p95", "KIT 95th percentile 20.5 %", +1)]
    for name, key, lab, side in lines:
        v = 100 * fl[name][key]
        ax.axhline(v, color=INK2 if name == "Geslin" else MUTED, lw=0.6, ls=(0, (3, 2)), zorder=2)
        ax.annotate(lab, (0.235, v), xytext=(0, 1.6 * side), textcoords="offset points", fontsize=4.8, color=INK2,
                    va="bottom" if side > 0 else "top", ha="left", zorder=6,
                    bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none", alpha=0.9))
    ax.set_xscale("log", base=2)
    ax.set_xticks(RATES)
    ax.set_xticklabels(RATE_LAB)
    ax.minorticks_off()
    ax.set_xlim(0.22, 2.3)
    ax.set_yscale("log")
    ax.set_ylim(0.08, 200)
    ax.set_yticks([0.1, 1, 10, 100])
    ax.set_yticklabels(["0.1", "1", "10", "100"])
    ax.set_xlabel("average rate")
    ax.set_ylabel("implied lifetime divergence, |s| x spread (%)")
    hs = [Line2D([], [], color=AMB_COL[T], lw=1.1, label=f"{T:g} C ambient") for T in (35.0, 25.0, 10.0)]
    hs += [Line2D([], [], color=INK2, lw=0.5, alpha=0.6, label="95th percentile (h = 10)")]
    hs += [Line2D([], [], color=INK2, ls=H_STYLE[h], marker=H_MARK[h], ms=2.6, mfc=INK2 if h == 10.0 else "white",
                  mec=INK2, mew=0.6, lw=0.9, label=f"h = {h:g} W/(m2 K)") for h in (5.0, 10.0, 25.0)]
    ax.legend(handles=hs, frameon=False, fontsize=4.9, loc="upper left", bbox_to_anchor=(0.0, 1.0), ncol=2,
              columnspacing=0.8, handlelength=2.0, handletextpad=0.4)
    ax.set_title("(a) the measured levers in the module: median over the band", loc="left", fontsize=7.0, color=INK)


def panel_b(ax, curves, exp, named):
    col = {"warm": AMB_COL[35.0], "cold": AMB_COL[10.0]}
    for ch in ("warm", "cold"):
        for fl, ls in (("Geslin", "-"), ("KIT", (0, (4, 1.6)))):
            cc = [c for c in curves if c["chamber"] == ch and c["floor"] == fl]
            x = [float(c["dT_K"]) for c in cc]
            y = [int(c["n_pairs"]) for c in cc]
            ax.plot(x, y, color=col[ch], ls=ls, lw=1.1, drawstyle="steps-post", zorder=3)
            pts = [e for e in exp if e["part"] == "C2" and e["central"] == "True" and e["chamber"].startswith(ch) and e["floor"] == fl]
            ax.plot([float(e["dT_K"]) for e in pts], [int(e["n_pairs"]) for e in pts], ls="none", marker="o", ms=3.4,
                    mfc=col[ch] if fl == "Geslin" else "white", mec=col[ch], mew=0.7, zorder=5)
            for e in pts:
                if fl == "Geslin":
                    ax.annotate(f"{int(e['n_pairs'])}", (float(e["dT_K"]), int(e["n_pairs"])),
                                xytext=(3, 4) if ch == "warm" else (3, -8), textcoords="offset points", fontsize=5.2, color=INK)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.5, 10.5)
    ax.set_xticks([0.5, 1, 2, 5, 10])
    ax.set_xticklabels(["0.5", "1", "2", "5", "10"])
    ax.minorticks_off()
    ax.set_ylim(1.5, 9000)
    ax.set_yticks([2, 10, 100, 1000])
    ax.set_yticklabels(["2", "10", "100", "1000"])
    ax.set_xlabel("temperature difference held within each pair, dT (K)")
    ax.set_ylabel("duplicate pairs for 80 % power")
    w, c = named["median_warm"], named["median_cold"]
    hs = [Line2D([], [], color=col["warm"], lw=1.1, label=f"warm chamber, 35 C: {w['label']}, {100 * abs(w['s']):.2f} % per K"),
          Line2D([], [], color=col["cold"], lw=1.1, label=f"cold chamber, 10 C: {c['label']}, {100 * abs(c['s']):.2f} % per K"),
          Line2D([], [], color=INK2, lw=1.1, label="against the Geslin floor"),
          Line2D([], [], color=INK2, lw=1.1, ls=(0, (4, 1.6)), label="against the KIT floor")]
    ax.legend(handles=hs, frameon=False, fontsize=4.9, loc="upper right", handlelength=2.2, handletextpad=0.4)
    ax.set_title("(b) the experiment: pairs needed, paired t-test at 5 %", loc="left", fontsize=7.0, color=INK)


def main(res):
    mp = json.load(open(os.path.join(res, "stage4_map.json")))
    curves = list(csv.DictReader(open(os.path.join(res, "stage4_experiment_curves.csv"))))
    exp = list(csv.DictReader(open(os.path.join(res, "stage4_experiment.csv"))))
    ej = json.load(open(os.path.join(res, "stage4_experiment.json")))
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import common as CM
    named = CM.named_levers()
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.75), dpi=220)
    fig.subplots_adjust(left=0.065, right=0.985, bottom=0.16, top=0.9, wspace=0.22)
    panel_a(ax[0], mp)
    panel_b(ax[1], curves, exp, named)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(res, f"stage4_figure.{ext}"))
    print("written", os.path.join(res, "stage4_figure.png"), ej["cases"], "experiment rows")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    main(ap.parse_args().results)
