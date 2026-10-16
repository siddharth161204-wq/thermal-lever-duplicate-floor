#!/usr/bin/env python3
"""Stage 2 figure, three panels. (a) ln lifetime against ambient temperature on the registered 1C temperature series,
each chemistry relative to its own 25 C mean, with Law 1 (at NCA's threshold) and Law 2 (at LFP's threshold) over the
same temperatures, likewise relative to 25 C. (b) The replicate floor per chemistry: every pair's relative lifetime
difference, its mean and 95th percentile, beside the Geslin floor. (c) Observed against predicted ln lifetime ratio
(warmer cell over cooler cell by self-heating) per pair, filled where the signs agree, open where they disagree,
grey where the lifetime difference is below its resolution.

Usage: python3 src/stage2/make_figure.py --results results
"""
import argparse, json, os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a897f", "#e1e0d9"
CHEM_COL = {"LFP": "#2a78d6", "NCA": "#eb6834", "NMC": "#1baf7a"}     # reference palette slots 1 to 3 (validated)
LAW_COL = {"law_1_kucinskis2022": "#eb6834", "law_2_schimpe2018": "#2a78d6"}   # as in the Stage 1 figure
CHEMS = ["LFP", "NCA", "NMC"]

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.6, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "xtick.labelsize": 6.0, "ytick.labelsize": 6.0,
                     "axes.spines.top": False, "axes.spines.right": False})


def panel_a(ax, lev):
    reg = lev["registered"]
    for chem in CHEMS:
        e = reg.get(chem, {})
        if not e.get("determinable"):
            continue
        m25 = e["mean_lnL"]["25C"]
        Ts = [T for T in (15, 25, 35) if f"{T}C" in e["mean_lnL"]]
        for T in Ts:
            Ls = e["L_efc"][f"{T}C"]
            ax.scatter([T + {"LFP": -0.6, "NCA": 0.0, "NMC": 0.6}[chem]] * len(Ls), np.log(Ls) - m25, s=7,
                       facecolor="none", edgecolor=CHEM_COL[chem], lw=0.6, zorder=3)
        ys = [e["mean_lnL"][f"{T}C"] - m25 for T in Ts]
        ax.plot(Ts, ys, "-o", color=CHEM_COL[chem], lw=1.4, ms=3.6, zorder=4)
        ax.annotate(f"{chem}  s = {100 * e['s_25_35_per_K']:.1f} %/K", (Ts[-1], ys[-1]), xytext=(4, 0),
                    textcoords="offset points", fontsize=5.6, color=INK, va="center")
    for law, src in (("law_1_kucinskis2022", "NCA"), ("law_2_schimpe2018", "LFP")):
        curves = reg.get(src, {}).get("law_curves") or []
        if not curves:
            continue
        key = "ln_L_law1" if law.startswith("law_1") else "ln_L_law2"
        T = np.array([c["T_C"] for c in curves])
        y = np.array([c[key] for c in curves])
        y25 = np.interp(25.0, T, y)
        ax.plot(T, y - y25, "--", color=LAW_COL[law], lw=1.0, zorder=2)
    ax.axhline(0, color=GRID, lw=0.6, zorder=0)
    ax.set_xlim(12, 44)
    ax.set_xticks([15, 25, 35])
    ax.set_xlabel("ambient temperature (C)")
    ax.set_ylabel("ln L minus its mean at 25 C")
    ax.legend(handles=[Line2D([], [], color=INK2, lw=1.4, marker="o", ms=3, label="measured, 1C series (cells open)"),
                       Line2D([], [], color=LAW_COL["law_1_kucinskis2022"], lw=1.0, ls="--", label="Law 1 (NCA threshold)"),
                       Line2D([], [], color=LAW_COL["law_2_schimpe2018"], lw=1.0, ls="--", label="Law 2 (LFP threshold)")],
              frameon=False, fontsize=5.4, loc="lower left")
    ax.set_title("(a) the lever", loc="left", fontsize=7.2, color=INK)


def panel_b(ax, pairs, floor):
    rng = np.random.default_rng(7)
    counts = {}
    for i, chem in enumerate(CHEMS):
        p = pairs[(pairs.chemistry == chem) & pairs.rel_diff.notna()]
        if p.empty:
            continue
        y = 100 * p.rel_diff.to_numpy()
        x = i + rng.uniform(-0.18, 0.18, y.size)
        ax.scatter(x, y, s=8, facecolor=CHEM_COL[chem], edgecolor="white", lw=0.4, zorder=3)
        ax.plot([i - 0.28, i + 0.28], [y.mean()] * 2, color=INK, lw=1.2, zorder=4)
        ax.scatter([i + 0.34], [np.percentile(y, 95)], marker="<", s=12, color=INK, zorder=4)
        counts[i] = len(y)
    for val, lab in ((100 * floor["mean"], "Geslin mean 3.0 %"), (100 * floor["p95"], "Geslin 95th percentile 8.7 %")):
        ax.axhline(val, color=MUTED, lw=0.8, ls=":", zorder=1)
        ax.text(2.45, val, lab, fontsize=5.2, color=INK2, va="bottom", ha="right")
    ax.set_xticks(range(3))
    ax.set_xticklabels([f"{c}\n{counts.get(i, 0)} pairs" for i, c in enumerate(CHEMS)])
    ax.set_xlim(-0.6, 2.6)
    top = max(10.0, 1.15 * 100 * float(np.nanmax(pairs.rel_diff))) if pairs.rel_diff.notna().any() else 10.0
    ax.set_ylim(0, top)
    ax.set_ylabel("pair lifetime difference (%)")
    ax.legend(handles=[Line2D([], [], color=INK, lw=1.2, label="mean"),
                       Line2D([], [], color=INK, marker="<", ls="none", ms=3.5, label="95th percentile")],
              frameon=False, fontsize=5.4, loc="upper left", bbox_to_anchor=(0.0, 0.84))
    ax.set_title("(b) the replicate floor", loc="left", fontsize=7.2, color=INK)


def panel_c(ax, o, share):
    o = o[o.obs_ln_ratio.notna()]
    for _, r in o.iterrows():
        x, y = r.pred_ln_ratio, r.obs_ln_ratio
        mk = "o" if r.rate_group == "1C or below" else "^"
        if r.outcome == "tie":
            ax.scatter(x, y, marker=mk, s=10, facecolor="none", edgecolor=MUTED, lw=0.6, zorder=2)
        elif np.sign(x) == np.sign(y):
            ax.scatter(x, y, marker=mk, s=12, facecolor=CHEM_COL[r.chemistry], edgecolor="white", lw=0.4, zorder=3)
        else:
            ax.scatter(x, y, marker=mk, s=12, facecolor="none", edgecolor=CHEM_COL[r.chemistry], lw=0.8, zorder=3)
    lim = max(1e-4, np.nanmax(np.abs(np.r_[o.pred_ln_ratio, o.obs_ln_ratio]))) * 1.08 if len(o) else 0.1
    ax.plot([-lim, lim], [-lim, lim], color=INK2, lw=0.7, ls="--", zorder=1)
    ax.axhline(0, color=GRID, lw=0.6, zorder=0)
    ax.axvline(0, color=GRID, lw=0.6, zorder=0)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    a = share["all"]
    p = a["sign_test_p_one_sided"]
    ptxt = "p < 0.001" if p < 0.001 else f"p = {p:.3f}"
    txt = (f"warmer shorter in {a['warmer_shorter']} of {a['warmer_shorter'] + a['warmer_longer']} pairs (ties {a['ties']})\n"
           f"sign test {ptxt}; slope {a.get('slope', float('nan')):.2f}")
    ax.text(0.02, 0.98, txt, transform=ax.transAxes, fontsize=5.3, color=INK, va="top")
    ax.set_xlabel("predicted ln(L warmer / L cooler) = -s dT")
    ax.set_ylabel("observed ln(L warmer / L cooler)")
    hs = [Line2D([], [], marker="o", ls="none", ms=3.6, mfc=CHEM_COL[c], mec="white", label=c) for c in CHEMS]
    hs += [Line2D([], [], marker="o", ls="none", ms=3.6, mfc="none", mec=INK2, label="signs disagree"),
           Line2D([], [], marker="^", ls="none", ms=3.6, mfc=INK2, mec="white", label="2C and 3C")]
    ax.legend(handles=hs, frameon=False, fontsize=5.0, loc="lower right", handletextpad=0.2)
    ax.set_title("(c) the thermal share", loc="left", fontsize=7.2, color=INK)


def main(res):
    lev = json.load(open(os.path.join(res, "stage2_lever.json")))
    pairs = pd.read_csv(os.path.join(res, "stage2_pairs.csv"))
    o = pd.read_csv(os.path.join(res, "stage2_pairs_oriented_measure_a.csv"))
    share = json.load(open(os.path.join(res, "stage2_thermal_share.json")))["measure_a"]
    floor = json.load(open(os.path.join(ROOT, "inputs", "floor.json")))
    fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.75), dpi=220)
    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.17, top=0.9, wspace=0.36)
    panel_a(ax[0], lev)
    panel_b(ax[1], pairs, floor)
    panel_c(ax[2], o, share)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(res, f"stage2_figure.{ext}"))
    print("written", os.path.join(res, "stage2_figure.png"))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    main(ap.parse_args().results)
