"""Stage 6, the estimators of Parts B and C and their bootstrap (DECISIONS.md D90 to D93).

Everything is computed from per-set sums, so that a bootstrap draw of sets is a weighted sum: for set j with n_j cells,
temperature T (the measure in use), y = ln L and, where needed, r = ln R and an instrument z, the within-set sums of
squares and products Sxx, Sxy, Syy, Srr, Sxr, Sry, Szz, Szx, Szy (deviations from the set's own means), and the set's
means. A draw takes, within each stratum (a batch, and within it the sets of one cell apart from the sets of two or
more), as many sets as the stratum has, with replacement; the weight of a set is the number of times it is drawn.
Signs: s = -d ln L / dT (D82 (3)), so s_w = -(the within-set slope) and s_b = -(the between-set slope).
"""
import math

import numpy as np
from scipy import stats


def set_table(cells, tkey="T", ykey="y", rkey=None, zkey=None):
    """Per-set sums from a list of cell dicts carrying 'set', 'batch', 'paper' and the named numeric keys."""
    by = {}
    for c in cells:
        by.setdefault(c["set"], []).append(c)
    out = []
    for sid, m in sorted(by.items()):
        T = np.array([c[tkey] for c in m], dtype=float)
        y = np.array([c[ykey] for c in m], dtype=float) if ykey else np.zeros(len(m))
        d = {"set": sid, "batch": m[0]["batch"], "paper": m[0]["paper"], "n": len(m), "Tbar": T.mean(), "ybar": y.mean(),
             "Sxx": float(((T - T.mean()) ** 2).sum()), "Sxy": float(((T - T.mean()) * (y - y.mean())).sum()),
             "Syy": float(((y - y.mean()) ** 2).sum()), "cells": [c["cell"] for c in m],
             "T": T.tolist(), "y": y.tolist()}
        if rkey:
            r = np.array([c[rkey] for c in m], dtype=float)
            d.update({"Srr": float(((r - r.mean()) ** 2).sum()), "Sxr": float(((T - T.mean()) * (r - r.mean())).sum()),
                      "Sry": float(((r - r.mean()) * (y - y.mean())).sum())})
        if zkey:
            z = np.array([c[zkey] for c in m], dtype=float)
            d.update({"Szz": float(((z - z.mean()) ** 2).sum()), "Szx": float(((z - z.mean()) * (T - T.mean())).sum()),
                      "Szy": float(((z - z.mean()) * (y - y.mean())).sum())})
        out.append(d)
    return out


def col(sets, key):
    return np.array([s[key] for s in sets], dtype=float)


def strata(sets):
    """Indices of the sets by stratum: batch, and within it one-cell sets apart from the others."""
    st = {}
    for i, s in enumerate(sets):
        st.setdefault((s["batch"], s["n"] >= 2), []).append(i)
    return [np.array(v) for _, v in sorted(st.items())]


def draw_weights(sets, B, rng):
    """B x J matrix of how many times each set is drawn."""
    J = len(sets)
    W = np.zeros((B, J))
    for idx in strata(sets):
        m = len(idx)
        picks = idx[rng.integers(0, m, size=(B, m))]
        np.add.at(W, (np.repeat(np.arange(B), m), picks.ravel()), 1.0)
    return W


def ones(sets):
    return np.ones((1, len(sets)))


# ------------------------------------------------------------------------------------------- the statistics, vectorized
def within(sets, W):
    """B1 and B4 over the sets of two or more: s_w, the pairs form's lever, the residual sd and its df."""
    multi = col(sets, "n") >= 2
    Wm = W * multi
    n = col(sets, "n")
    Sxx, Sxy, Syy = Wm @ col(sets, "Sxx"), Wm @ col(sets, "Sxy"), Wm @ col(sets, "Syy")
    N, J = Wm @ n, Wm.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        b = Sxy / Sxx
        ssr = Syy - Sxy * Sxy / Sxx
        df = N - J - 1
        sd = np.sqrt(ssr / df)
        pairs = (Wm @ (n * col(sets, "Sxy"))) / (Wm @ (n * col(sets, "Sxx")))
    return {"s_w": -b, "s_pairs": -pairs, "sd_cell": sd, "df": df, "N": N, "J": J, "Sxx": Sxx}


def with_resistance(sets, W):
    """B2: ln L on T and ln R with one constant per set: the temperature's lever with R held, and d ln L / d ln R."""
    multi = col(sets, "n") >= 2
    Wm = W * multi
    Sxx, Srr, Sxr = Wm @ col(sets, "Sxx"), Wm @ col(sets, "Srr"), Wm @ col(sets, "Sxr")
    Sxy, Sry, Syy = Wm @ col(sets, "Sxy"), Wm @ col(sets, "Sry"), Wm @ col(sets, "Syy")
    with np.errstate(divide="ignore", invalid="ignore"):
        det = Sxx * Srr - Sxr * Sxr
        bT = (Srr * Sxy - Sxr * Sry) / det
        bR = (Sxx * Sry - Sxr * Sxy) / det
        bT_alone = Sxy / Sxx
        r_on_T = Sxr / Sxx
    return {"s_w_R_held": -bT, "lnL_per_lnR": bR, "s_w_same_cells": -bT_alone, "lnR_per_K_within": r_on_T}


def instrumented(sets, W):
    """B5: the slope with the channel's temperature from other batches as the instrument, and the first stage's F."""
    multi = col(sets, "n") >= 2
    Wm = W * multi
    n = col(sets, "n")
    Szz, Szx, Szy, Sxx = Wm @ col(sets, "Szz"), Wm @ col(sets, "Szx"), Wm @ col(sets, "Szy"), Wm @ col(sets, "Sxx")
    N, J = Wm @ n, Wm.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        b_iv = Szy / Szx
        g = Szx / Szz
        ss1 = Sxx - Szx * Szx / Szz
        F = g * g * Szz / (ss1 / (N - J - 1))
    return {"s_iv": -b_iv, "first_stage_slope": g, "first_stage_F": F}


def between(sets, W, s_w=None):
    """C1 and C3: the set means of ln L on the set means of T, one constant per batch, weighted by the set's cells."""
    n = col(sets, "n")
    Tb, yb = col(sets, "Tbar"), col(sets, "ybar")
    batches = sorted(set(s["batch"] for s in sets))
    CTT = CTy = Cyy = 0.0
    for b in batches:
        mb = np.array([s["batch"] == b for s in sets], dtype=float)
        Wn = W * (n * mb)
        SW, ST, Sy = Wn.sum(axis=1), Wn @ Tb, Wn @ yb
        STT, STy, Syy = Wn @ (Tb * Tb), Wn @ (Tb * yb), Wn @ (yb * yb)
        with np.errstate(divide="ignore", invalid="ignore"):
            CTT = CTT + np.where(SW > 0, STT - ST * ST / SW, 0.0)
            CTy = CTy + np.where(SW > 0, STy - ST * Sy / SW, 0.0)
            Cyy = Cyy + np.where(SW > 0, Syy - Sy * Sy / SW, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = {"s_b": -(CTy / CTT), "between_var_lnL_sum": Cyy, "between_var_T_sum": CTT}
        if s_w is not None:
            out["share"] = s_w * s_w * CTT / Cyy
            out["share_residual_form"] = -(2.0 * s_w * CTy + s_w * s_w * CTT) / Cyy
            out["s_b_minus_s_w"] = out["s_b"] - s_w
    return out


def percentile_ci(vals):
    v = np.asarray(vals, dtype=float)
    v = v[np.isfinite(v)]
    lo, hi = np.percentile(v, [2.5, 97.5])
    B = len(v)
    srt = np.sort(v)
    mc = []
    for q in (0.025, 0.975):                             # D77's Monte Carlo spread of one percentile end
        k = q * (B - 1)
        h = math.sqrt(B * q * (1 - q))
        a, b = int(max(0, math.floor(k - h))), int(min(B - 1, math.ceil(k + h)))
        mc.append(float((srt[b] - srt[a]) / 2.0))
    return [float(lo), float(hi)], mc, int(len(v))


# ------------------------------------------------------------------------------------------------------- B3, exactly
def sign_test(sets):
    """Within each set of two or more: is the warmest cell the first to die? Sets whose highest temperature or shortest
    life is shared are left out and counted. Exact two-sided test of the count against the Poisson-binomial with
    chance 1 / n per set (scipy's binomtest rule: outcomes no more likely than the one observed)."""
    hits, ps, ties = 0, [], 0
    for s in sets:
        if s["n"] < 2:
            continue
        T, y = np.array(s["T"]), np.array(s["y"])
        if (T == T.max()).sum() > 1 or (y == y.min()).sum() > 1:
            ties += 1
            continue
        hits += int(np.argmax(T) == np.argmin(y))
        ps.append(1.0 / s["n"])
    pmf = np.array([1.0])
    for p in ps:
        pmf = np.convolve(pmf, [1.0 - p, p])
    obs = pmf[hits]
    p_two = float(pmf[pmf <= obs * (1.0 + 1e-7)].sum())
    return {"sets": len(ps), "warmest_first": hits, "share": hits / len(ps) if ps else float("nan"),
            "chance_share": float(np.mean(ps)) if ps else float("nan"), "expected": float(np.sum(ps)),
            "p_two_sided": min(1.0, p_two), "left_out_ties": ties}
