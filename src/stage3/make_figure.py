#!/usr/bin/env python3
"""Stage 3 figure, three panels. (a) ln lifetime against the logged cycling temperature on the registered 1C series,
each chemistry relative to its own 25 C mean, with Law 1 (at the NCA and NMC threshold) and Law 2 (at LFP's) relative
to 25 C, and the KIT NMC protocols at 0.33/1C and 1/1C. (b) The replicate floor per chemistry, every pair (filled when
its difference exceeds its check spacing), with its mean and 95th percentile, beside the KIT floor and the Geslin
floor. (c) Observed against predicted ln lifetime ratio (warmer over cooler by measure (a)) per pair, filled where
the signs agree.

Usage: python3 src/stage3/make_figure.py --results results
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
CHEM_COL = {"LFP": "#2a78d6", "NCA": "#eb6834", "NMC": "#1baf7a"}
LAW_COL = {"law1": "#7a4fb5", "law2": "#52514e"}
CHEMS = ["LFP", "NCA", "NMC"]

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.4, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "xtick.labelsize": 5.8, "ytick.labelsize": 5.8,
                     "axes.spines.top": False, "axes.spines.right": False})


def panel_a(ax, lev, kit):
    reg = lev["registered_1C"]
    for chem in CHEMS:
        e = reg[chem]
        m25 = e["mean_lnL"]["25C"]
        xs, ys = [], []
        for T in (15, 25, 35):
            cells = e["cells"][f"{T}C"]
            Tc = np.array([c["T_cyc_C"] if c["T_cyc_C"] is not None else np.nan for c in cells], float)
            y = np.log([c["L_efc"] for c in cells]) - m25
            ok = np.isfinite(Tc)
            ax.scatter(Tc[ok], y[ok], s=7, facecolor="none", edgecolor=CHEM_COL[chem], lw=0.6, zorder=3)
            if (~ok).any():                                   # a cell without a cycling temperature, at its setpoint
                ax.scatter([float(T)] * int((~ok).sum()), y[~ok], s=12, marker="x", color=CHEM_COL[chem], lw=0.7, zorder=3)
            xs.append(float(np.nanmean(Tc)) if ok.any() else float(T))
            ys.append(float(np.mean(y)))
        ax.plot(xs[:2], ys[:2], "-o", color=CHEM_COL[chem], lw=1.3, ms=3.2, zorder=4)
        det = e["s_25_35"] is not None
        ax.plot(xs[1:], ys[1:], "-" if det else ":", color=CHEM_COL[chem], lw=1.3, zorder=4)
        ax.plot(xs[2:], ys[2:], "o", color=CHEM_COL[chem], ms=3.2, mfc=CHEM_COL[chem] if det else "white", zorder=5)
        s = e["s_25_35"]["s"] * 100 if det else None
        lab = f"{chem}  {s:+.1f} %/K" if det else f"{chem}\nnot determinable"
        dx, dy, ha = {"LFP": (4, -2, "left"), "NCA": (3, -5, "left"), "NMC": (3, 5, "left")}[chem]
        ax.annotate(lab, (xs[2], ys[2]), xytext=(dx, dy), textcoords="offset points", fontsize=5.2, color=INK, va="center", ha=ha)
    for key, src, col, name in (("ln_L_law1", "NCA", LAW_COL["law1"], "Law 1, 90 %"), ("ln_L_law2", "LFP", LAW_COL["law2"], "Law 2, 95 %")):
        cur = [c for c in reg[src]["law_curves"] if c[key] is not None and np.isfinite(c[key])]
        T = np.array([c["T_C"] for c in cur]); y = np.array([c[key] for c in cur])
        ax.plot(T, y - np.interp(25.0, T, y), ls=(0, (4, 2)), color=col, lw=1.0, alpha=0.9, zorder=6)
    for k, p in enumerate(("0-100 0.33/1.00C", "0-100 1.00/1.00C")):
        hi = next(r for r in kit if r["protocol"] == p and r["secant"].startswith("25"))
        lo = next(r for r in kit if r["protocol"] == p and r["secant"].startswith("10"))
        T = [lo["T_lo_C"], hi["T_lo_C"], hi["T_hi_C"]]
        y = [lo["s_per_K"] * (lo["T_hi_C"] - lo["T_lo_C"]), 0.0, -hi["s_per_K"] * (hi["T_hi_C"] - hi["T_lo_C"])]
        ax.plot(T, y, ":", color=CHEM_COL["NMC"], lw=0.7, zorder=1)
        ax.scatter(T, y, marker="D", s=9, facecolor="white", edgecolor=CHEM_COL["NMC"], lw=0.6, zorder=2)
    ax.axhline(0, color=GRID, lw=0.6, zorder=0)
    ax.set_xlim(9, 47)
    ax.set_xlabel("logged cycling temperature (C)")
    ax.set_ylabel("ln L minus its mean at 25 C")
    ax.legend(handles=[Line2D([], [], color=INK2, lw=1.3, marker="o", ms=3, label="Sandia 1C series, group means"),
                       Line2D([], [], color=INK2, marker="x", ls="none", ms=3.5, label="cell without temperature"),
                       Line2D([], [], color=CHEM_COL["NMC"], marker="D", ls=":", ms=2.8, mfc="white", lw=0.7, label="KIT HG2, 0.33/1C and 1/1C"),
                       Line2D([], [], color=LAW_COL["law1"], ls="--", lw=0.9, label="Law 1 (90 %)"),
                       Line2D([], [], color=LAW_COL["law2"], ls="--", lw=0.9, label="Law 2 (95 %)")],
              frameon=False, fontsize=4.9, loc="lower left", handlelength=1.8)
    ax.set_title("(a) the lever", loc="left", fontsize=7.0, color=INK)


def panel_b(ax, pairs, kit_pairs, geslin, kit_floor):
    rng = np.random.default_rng(7)
    labels = []
    cols = CHEMS + ["KIT", "Geslin"]
    for i, name in enumerate(cols):
        if name in CHEMS:
            p = pairs[(pairs.chemistry == name) & pairs.rel_diff.notna()]
            y = 100 * p.rel_diff.to_numpy()
            above = ~p.below_resolution.astype(bool).to_numpy()
            x = i + rng.uniform(-0.2, 0.2, y.size)
            ax.scatter(x[~above], y[~above], s=7, facecolor="none", edgecolor=CHEM_COL[name], lw=0.6, zorder=3)
            ax.scatter(x[above], y[above], s=8, facecolor=CHEM_COL[name], edgecolor="white", lw=0.3, zorder=3)
            mean, p95, n = y.mean(), np.percentile(y, 95), len(y)
        elif name == "KIT":
            y = 100 * kit_pairs.rel_diff.to_numpy()
            ax.scatter(i + rng.uniform(-0.2, 0.2, y.size), y, s=6, facecolor=MUTED, edgecolor="white", lw=0.3, zorder=3)
            mean, p95, n = 100 * kit_floor["mean"], 100 * kit_floor["p95"], kit_floor["pairs"]
        else:
            mean, p95, n = 100 * geslin["mean"], 100 * geslin["p95"], geslin["n_pairs"]
        ax.plot([i - 0.28, i + 0.28], [mean] * 2, color=INK, lw=1.2, zorder=4)
        ax.scatter([i + 0.34], [p95], marker="<", s=11, color=INK, zorder=4)
        labels.append(f"{name}\n{n} pairs")
    ax.set_yscale("log")
    ax.set_ylim(0.1, 2000)
    ax.set_yticks([0.1, 1, 10, 100])
    ax.set_yticklabels(["0.1", "1", "10", "100"])
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(labels, fontsize=5.4)
    ax.set_xlim(-0.6, len(cols) - 0.4)
    ax.set_ylabel("pair lifetime difference (%)")
    ax.legend(handles=[Line2D([], [], color=INK, lw=1.2, label="mean"),
                       Line2D([], [], color=INK, marker="<", ls="none", ms=3.3, label="95th percentile"),
                       Line2D([], [], marker="o", ls="none", ms=3, mfc="none", mec=INK2, label="below its check spacing")],
              frameon=False, fontsize=4.9, loc="upper center", ncol=3, columnspacing=0.8, handletextpad=0.3)
    ax.set_title("(b) the replicate floor", loc="left", fontsize=7.0, color=INK)


def panel_c(ax, o, share):
    o = o[o.obs_ln_ratio.notna() & o.pred_ln_ratio.notna()]
    for _, r in o.iterrows():
        x, y = r.pred_ln_ratio, r.obs_ln_ratio
        mk = "o" if r.rate_group == "1C or below" else "^"
        agree = np.sign(x) == np.sign(y)
        ax.scatter(x, y, marker=mk, s=11, facecolor=CHEM_COL[r.chemistry] if agree else "none",
                   edgecolor=CHEM_COL[r.chemistry] if not agree else "white", lw=0.7 if not agree else 0.3, zorder=3)
    xl = max(1e-3, float(np.nanmax(np.abs(o.pred_ln_ratio)))) * 1.15
    yl = max(1e-3, float(np.nanmax(np.abs(o.obs_ln_ratio)))) * 1.1
    ax.plot([-xl, xl], [-xl, xl], color=INK2, lw=0.7, ls="--", zorder=1)
    ax.axhline(0, color=GRID, lw=0.6, zorder=0)
    ax.axvline(0, color=GRID, lw=0.6, zorder=0)
    ax.set_xlim(-xl, xl)
    ax.set_ylim(-yl, 1.45 * yl)
    a = share["all"]
    agree = int((np.sign(o.pred_ln_ratio) == np.sign(o.obs_ln_ratio)).sum())
    txt = (f"signs agree in {agree} of {len(o)} pairs with a prediction\n"
           f"{a['ties']} of {a['pairs']} pairs differ by less than their check spacing")
    ax.text(0.02, 0.98, txt, transform=ax.transAxes, fontsize=5.0, color=INK, va="top")
    ax.set_xlabel("predicted ln(L warmer / L cooler) = -s dT")
    ax.set_ylabel("observed ln(L warmer / L cooler)")
    hs = [Line2D([], [], marker="o", ls="none", ms=3.4, mfc=CHEM_COL[c], mec="white", label=c) for c in ("NCA", "NMC")]
    hs += [Line2D([], [], marker="o", ls="none", ms=3.4, mfc="none", mec=INK2, label="signs disagree"),
           Line2D([], [], marker="^", ls="none", ms=3.4, mfc=INK2, mec="white", label="2C and 3C")]
    ax.legend(handles=hs, frameon=False, fontsize=4.9, loc="lower left", handletextpad=0.2)
    ax.set_title("(c) the thermal share", loc="left", fontsize=7.0, color=INK)


def main(res):
    lev = json.load(open(os.path.join(res, "stage3_lever.json")))
    pairs = pd.read_csv(os.path.join(res, "stage3_pairs.csv"))
    o = pd.read_csv(os.path.join(res, "stage3_pairs_oriented_measure_a.csv"))
    ts = json.load(open(os.path.join(res, "stage3_thermal_share.json")))
    kit_pairs = pd.read_csv(os.path.join(res, "stage2_pairs.csv"))
    fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.75), dpi=220)
    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.2, top=0.9, wspace=0.34)
    panel_a(ax[0], lev, lev["kit_levers"])
    panel_b(ax[1], pairs, kit_pairs, ts["geslin_floor"], ts["kit_floor"])
    panel_c(ax[2], o, ts["measure_a"])
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(res, f"stage3_figure.{ext}"))
    print("written", os.path.join(res, "stage3_figure.png"))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    main(ap.parse_args().results)
