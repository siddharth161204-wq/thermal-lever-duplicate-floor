#!/usr/bin/env python3
"""Stage 2 figure (KIT dataset): (a) the 25 to 40 C lever s per protocol with its 95 % interval (filled) and the 10 to
25 C lever (open), the two Stage 1 laws marked; (b) the replicate floor, every pre-specified pair's lifetime difference
at 90 % by temperature, its mean and 95th percentile, beside the Geslin floor. Usage: python3 src/stage2/make_figure_kit.py"""
import json, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a897f", "#e1e0d9"
WIN_COL = {"0-100": "#2a78d6", "10-100": "#eb6834", "10-90": "#1baf7a"}     # reference palette slots 1 to 3
TEMP_COL = {0: "#2a78d6", 10: "#1baf7a", 25: "#eda100", 40: "#eb6834"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.6, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.6, "xtick.labelsize": 5.8,
                     "ytick.labelsize": 6.0, "axes.spines.top": False, "axes.spines.right": False})
r = json.load(open(os.path.join(ROOT, "results", "stage2_lever.json")))
fl = json.load(open(os.path.join(ROOT, "inputs", "floor.json")))
pairs = pd.read_csv(os.path.join(ROOT, "results", "stage2_pairs.csv"))
fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.9), dpi=220, gridspec_kw={"width_ratios": [1.9, 1]})
fig.subplots_adjust(left=0.075, right=0.985, bottom=0.29, top=0.9, wspace=0.28)
a = ax[0]
hi = {x["protocol"]: x for x in r["lever_25_40"]}
lo = {x["protocol"]: x for x in r["lever_10_25"]}
order = [f"{w} {c}C" for w in ("0-100", "10-100", "10-90") for c in ("0.33/0.33", "0.33/1.00", "1.00/1.00", "1.67/1.00")]
for i, p in enumerate(order):
    col = WIN_COL[p.split()[0]]
    x = hi[p]
    if x["determinable"]:
        a.errorbar(i - 0.12, 100 * x["s_per_K"], yerr=[[100 * (x["s_per_K"] - x["ci95"][0])], [100 * (x["ci95"][1] - x["s_per_K"])]],
                   fmt="o", ms=3.6, color=col, ecolor=col, elinewidth=1.0, capsize=0)
    y = lo[p]
    if y["determinable"]:
        a.plot(i + 0.14, 100 * y["s_per_K"], "o", ms=3.2, mfc="white", mec=col, mew=0.9)
for v, lab, ls in ((5.2, "Law 2, 5.2 %/K", "--"), (-0.99, "Law 1, -0.99 %/K", ":")):
    a.axhline(v, color=INK2, lw=0.8, ls=ls)
    a.text(-0.5, v, lab, fontsize=5.3, color=INK2, ha="left", va="bottom")
a.axhline(0, color=GRID, lw=0.6)
a.set_xticks(range(12))
a.set_xticklabels([p.split()[1] for p in order], rotation=50, ha="right")
for k, w in enumerate(("0-100 %", "10-100 %", "10-90 %")):
    a.text(1.5 + 4 * k, -0.36, w, transform=a.get_xaxis_transform(), ha="center", fontsize=6.0, color=INK)
a.set_xlim(-0.6, 11.6)
a.set_ylabel("lever s (% per K)")
a.legend(handles=[Line2D([], [], marker="o", ls="none", color=INK2, ms=3.6, label="25 to 40 C, with 95 % interval"),
                  Line2D([], [], marker="o", ls="none", mfc="white", mec=INK2, ms=3.2, label="10 to 25 C")],
         frameon=False, fontsize=5.4, loc="lower left")
a.set_title("(a) the lever by protocol (charge/discharge rate)", loc="left", fontsize=7.2, color=INK)
b = ax[1]
rng = np.random.default_rng(7)
y = 100 * pairs.rel_diff.to_numpy()
b.scatter(rng.uniform(-0.2, 0.2, y.size), y, s=9, c=[TEMP_COL[t] for t in pairs.age_temp], edgecolor="white", lw=0.4, zorder=3)
f = r["floor"]
b.plot([-0.3, 0.3], [100 * f["mean"]] * 2, color=INK, lw=1.2)
b.scatter([0.36], [100 * f["p95"]], marker="<", s=12, color=INK)
b.scatter([1.0], [100 * fl["mean"]], marker="_", s=120, color=INK, lw=1.2)
b.scatter([1.06], [100 * fl["p95"]], marker="<", s=12, color=INK)
b.set_xticks([0, 1]); b.set_xticklabels([f"KIT HG2\n{f['pairs']} pairs", "Geslin\n45 pairs"])
b.set_xlim(-0.6, 1.5); b.set_ylim(0, max(25, 1.1 * y.max()))
b.set_ylabel("pair lifetime difference at 90 % (%)")
b.legend(handles=[Line2D([], [], color=INK, lw=1.2, label="mean"), Line2D([], [], color=INK, marker="<", ls="none", ms=3.5, label="95th percentile")]
         + [Line2D([], [], marker="o", ls="none", mfc=c, mec="white", ms=3.6, label=f"{t} C") for t, c in TEMP_COL.items()],
         frameon=False, fontsize=5.2, loc="upper right")
b.set_title("(b) the replicate floor", loc="left", fontsize=7.2, color=INK)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(ROOT, "results", f"stage2_figure.{ext}"))
print("written results/stage2_figure.png")
