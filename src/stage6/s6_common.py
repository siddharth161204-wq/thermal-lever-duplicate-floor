"""Stage 6, what every part shares (DECISIONS.md D82 to D96): the archived values, the fixed-effects fit, the power.

The archived values are read at full precision from the Stage 3 to 5 files and checked against D82; the code stops if
any differs. Stage 4's noncentral t power (`src/stage4/experiment.py`, unchanged) is imported, as Stage 5 did.
"""
import csv, json, math, os, sys

import numpy as np
from scipy import optimize

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
ROOT = os.path.dirname(SRC)
SRC4 = os.path.join(SRC, "stage4")
if SRC4 not in sys.path:
    sys.path.insert(0, SRC4)
import common as CM4                    # noqa: E402  src/stage4/common.py, unchanged
import experiment as EX4                # noqa: E402  src/stage4/experiment.py, unchanged

SEED = 20261116
N_BOOT = 10000
ALPHA, POWER = EX4.ALPHA, EX4.POWER
CHAMBER_C = 30.0                        # the chamber's set point, both papers (D85)

# D82 (2), at full precision
LEVERS_D82 = {
    ("registered 1C", "15 to 25 C"): (0.03371837126587428, 0.02949074115460169, 0.038270002135397496),
    ("secondary 2C", "15 to 25 C"): (0.02049870998889622, 0.011837002110908613, 0.029076510430014992),
    ("secondary 2C", "25 to 35 C"): (0.11005822444213981, 0.08302214841830174, 0.13729950535845456),
}
SD_E_D82 = {"Geslin": 0.04454298465298393, "KIT": 0.10412996462328801}
STAGE5_DETECTABLE_D82 = {("KIT", "KIT"): 0.06559944091124832, ("KIT", "Geslin"): 0.028061037956985384}


def path(rel):
    return os.path.join(ROOT, rel)


def archive():
    """Sandia's LFP levers, the floors with Stage 4's registered shape, and Stage 5's detectable |s|, checked."""
    lev = {}
    for r in csv.DictReader(open(path("results/stage3_lever.csv"))):
        if r["chemistry"] == "LFP" and r["determinable"] == "True":
            lev[(r["series"], r["secant"])] = (float(r["s_per_K"]), float(r["ci95_lo"]), float(r["ci95_hi"]),
                                               int(float(r["n_lo"])), int(float(r["n_hi"])))
    for k, v in LEVERS_D82.items():
        if k not in lev or lev[k][:3] != v:
            raise SystemExit(f"Stage 3 LFP lever {k} differs from D82")
    floors = CM4.read_floors()
    ex = json.load(open(path("results/stage4_experiment.json")))["shapes"]
    shapes = {}
    for name in ("Geslin", "KIT"):
        ln = EX4.lognormal_shape(floors[name]["mean"], floors[name]["p95"])
        if ln["sd_e"] != ex[name]["lognormal"]["sd_e"] or ln["sd_e"] != SD_E_D82[name]:
            raise SystemExit(f"the {name} floor's sd(e) differs from D82 or results/stage4_experiment.json")
        shapes[name] = {"kind": "lognormal", "mu": ln["mu"], "sigma": ln["sigma"], "sd_e": ln["sd_e"],
                        "sd_cell": ln["sd_e"] / math.sqrt(2.0), "mean_abs_e": floors[name]["mean"],
                        "p95_abs_e": floors[name]["p95"], "pairs": floors[name]["pairs"],
                        "sd_e_normal_shape": ex[name]["normal"]["sd_e"]}
    s5 = {}
    for x in json.load(open(path("results/stage5_paper_numbers.json")))["C1"]:
        if x.get("registered"):
            s5[(x["case"], x["floor"])] = x
    for k, v in STAGE5_DETECTABLE_D82.items():
        if s5[k]["detectable_s_per_K"] != v:
            raise SystemExit(f"Stage 5 detectable |s| {k} differs from D82")
    return {"levers": lev, "shapes": shapes, "stage5": s5}


# --------------------------------------------------------------------------------------------- the fixed-effects fit
def demean(x, groups):
    """x minus its group mean, for numpy arrays x (1-d or 2-d by rows) and integer group labels."""
    x = np.asarray(x, dtype=float)
    g = np.asarray(groups)
    out = x.copy()
    for u in np.unique(g):
        m = g == u
        out[m] = x[m] - x[m].mean(axis=0)
    return out


def fe_fit(y, X, groups):
    """Least squares of y on the columns of X with one constant per group. Returns the coefficients, the residuals,
    the residual degrees of freedom N - J - p and the residual standard deviation."""
    yd = demean(y, groups)
    Xd = demean(np.asarray(X, dtype=float).reshape(len(yd), -1), groups)
    coef, *_ = np.linalg.lstsq(Xd, yd, rcond=None)
    res = yd - Xd @ coef
    J = len(np.unique(groups))
    df = len(yd) - J - Xd.shape[1]
    sd = math.sqrt(float(res @ res) / df) if df > 0 else float("nan")
    return coef, res, df, sd


# ------------------------------------------------------------------------------------------------------------ power
def fe_power(s, Sxx, N, J, sd_cell):
    """Two-sided 5 % power of the fixed-effects slope's t-test (N cells, J sets, one slope) at lever |s| per K:
    the noncentral t with N - J - 1 degrees of freedom and noncentrality |s| sqrt(Sxx) / sd_cell, Sxx the within-set
    sum of squared temperature deviations. Stage 4's power_t(n, delta, sd) is that t with n - 1 degrees of freedom
    and noncentrality sqrt(n) delta / sd, so it is called with n = N - J and delta = |s| sqrt(Sxx / (N - J))."""
    n = N - J
    return EX4.power_t(n, abs(s) * math.sqrt(Sxx / n), sd_cell)


def fe_detectable(Sxx, N, J, sd_cell):
    """The |s| per K at which fe_power is 0.80 (Brent's method, 1e-12 relative)."""
    hi = 1e-4
    while fe_power(hi, Sxx, N, J, sd_cell) < POWER:
        hi *= 2.0
    return optimize.brentq(lambda s: fe_power(s, Sxx, N, J, sd_cell) - POWER, 0.0, hi, rtol=1e-12, xtol=1e-15)


def pairs_detectable(n_pairs, dT_rms, sd_e):
    """Stage 5's calculation (D75): power_t(n, |s| dT_rms, sd(e)) at 0.80, the pairs taken as independent."""
    f = lambda s: EX4.power_t(n_pairs, s * dT_rms, sd_e) - POWER
    hi = 1e-4
    while f(hi) < 0:
        hi *= 2.0
    return optimize.brentq(f, 0.0, hi, rtol=1e-12, xtol=1e-15)


def shape_draws(shape, size, rng):
    """Draws of a pair difference e from Stage 4's registered shape (log-normal magnitude, random sign), or normal."""
    if shape["kind"] == "lognormal":
        mag = np.exp(shape["mu"] + shape["sigma"] * rng.standard_normal(size))
        return mag * np.where(rng.random(size) < 0.5, -1.0, 1.0)
    return shape["sd_e"] * rng.standard_normal(size)
