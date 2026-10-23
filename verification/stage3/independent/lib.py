"""Stage 3 independent check: the row, step, ampere-hour and temperature primitives, written from the text of
DECISIONS.md D40, D41, D43 and D50 without the pipeline's code (src/stage2, src/stage3 were not opened).

Every function works on the cache written by s0_cache.py (the raw columns as float64, in logged order).
"""
import json, os
import numpy as np

CACHE = os.environ.get("S3V2_CACHE", "/tmp/claude-0/-home-claude/6d63e4c8-32e9-55cd-b414-99eb5512ec94/scratchpad/cache")
RATED = {"LFP": 1.1, "NCA": 3.2, "NMC": 3.0}            # Ah, D39
CUTOFF = {"LFP": 2.0, "NCA": 2.5, "NMC": 2.0}           # V, D41
THRESHOLD = {"LFP": 0.95, "NCA": 0.90, "NMC": 0.90}     # D42
REST, CHARGE, DISCHARGE = 0, 1, -1


def index():
    with open(os.path.join(CACHE, "index.json")) as fh:
        return json.load(fh)


def parse_id(cid):
    # SNL_18650_LFP_15C_0-100_0.5/1C_a
    p = cid.split("_")
    chem, amb, win, rates, rep = p[2], p[3], p[4], p[5], p[6]
    c_ch, c_dis = rates.rstrip("C").split("/")
    lo, hi = win.split("-")
    return {"cell_id": cid, "chemistry": chem, "ambient_C": float(amb.rstrip("C")), "soc_min": float(lo),
            "soc_max": float(hi), "window": win, "c_charge": float(c_ch), "c_discharge": float(c_dis),
            "replicate": rep, "rated_Ah": RATED[chem], "group": f"{chem}_{amb}_{win}_{rates}"}


def load(cid, idx=None):
    idx = idx or index()
    d = np.load(idx[cid]["cache"])
    return {k: d[k] for k in d.files}


def row_structure(cyc):
    """Late records (cycle index exactly one below the running maximum) and rows more than one cycle below it,
    on the rows as logged (D50 (1))."""
    rm = np.maximum.accumulate(cyc)
    late = cyc == rm - 1
    other = cyc < rm - 1
    return rm, late, other


def ordered_rows(a):
    """D50 (1): drop the rows more than one cycle below the running maximum, then a stable sort by cycle index."""
    _, late, other = row_structure(a["cyc"])
    keep = ~other
    o = np.argsort(a["cyc"][keep], kind="stable")
    out = {k: v[keep][o] for k, v in a.items()}
    out["_logged_row"] = np.flatnonzero(keep)[o]
    return out, int(late.sum()), int(other.sum())


def labels(I, rated_Ah):
    """D40: rest when |I| < 0.005C, otherwise charge or discharge by the sign of the current."""
    return np.where(np.abs(I) < 0.005 * rated_Ah, REST, np.where(I > 0, CHARGE, DISCHARGE)).astype(np.int8)


def steps(lab):
    """Maximal runs of one label in the ordered rows: start and end row of each step (inclusive)."""
    b = np.flatnonzero(np.diff(lab) != 0) + 1
    starts = np.r_[0, b]
    ends = np.r_[b - 1, len(lab) - 1]
    return starts, ends


def step_ah_per_sample(Q, starts, ends):
    """D50 (2): sum of the counter's increases over the transitions from the sample before the step to its last
    sample; a decrease counts as a restart after which the new value counts in full. For a step that starts at the
    first row there is no sample before it, and its first value counts in full."""
    dq = np.diff(Q)
    inc = np.where(dq >= 0, dq, Q[1:])          # contribution of transition j-1 -> j, stored at j-1
    # each step sums its own transitions (rows max(s,1)..e, i.e. inc[max(s,1)-1 : e]); summed per step, not through
    # a global running sum, so that no rounding accumulates across the file
    lo = np.maximum(starts, 1) - 1
    hi = ends                                   # exclusive end in inc
    out = np.zeros(len(starts))
    nonempty = hi > lo
    if nonempty.any() and len(inc):
        seg = np.add.reduceat(inc, lo[nonempty]) if len(lo[nonempty]) else np.zeros(0)
        # reduceat sums inc[lo[i]:lo[i+1]]; the segments are contiguous because steps tile the rows
        out[nonempty] = seg
    out = np.where(starts == 0, out + Q[0], out)
    return out


def step_ah_per_sample_loop(Q, s, e):
    """The same rule written as a plain loop, used to cross-check the vectorised form on the V2 cells."""
    total = 0.0 if s > 0 else Q[0]
    for j in range(max(s, 1), e + 1):
        d = Q[j] - Q[j - 1]
        total += d if d >= 0 else Q[j]
    return total


def counter_ah_registered(Q, starts, ends):
    """The registered rule of D40 as D50 (2) describes process_cells.counter_ah: the last value minus the value
    before the step, or the last alone after a restart (last below before)."""
    before = np.where(starts > 0, Q[np.maximum(starts - 1, 0)], 0.0)
    last = Q[ends]
    return np.where(last >= before, last - before, last)


def step_temperature(t, T, starts, ends, zero_dt="sample_mean"):
    """D50 (3): a step's mean temperature is the trapezoidal mean over its consecutive pairs of samples that both
    hold a value; the mean of its valid samples when it has no such pair; missing when none holds a value.
    Also returns the last valid value, the number of valid samples and the valid-pair duration.

    zero_dt: what to do when a step has valid pairs but their durations sum to zero (two records logged at the same
    instant). "sample_mean" (primary reading) takes the mean of the valid samples, the limit of the trapezoidal mean
    as the duration goes to zero; "nan" leaves the step mean missing (0/0)."""
    n = len(t)
    sid = np.repeat(np.arange(len(starts)), ends - starts + 1)
    valid = ~np.isnan(T)
    pair_ok = valid[:-1] & valid[1:] & (sid[:-1] == sid[1:])
    dt = np.diff(t)
    num = np.where(pair_ok, 0.5 * (np.nan_to_num(T[:-1]) + np.nan_to_num(T[1:])) * dt, 0.0)
    den = np.where(pair_ok, dt, 0.0)
    npair = np.bincount(sid[:-1], weights=pair_ok.astype(float), minlength=len(starts))
    snum = np.bincount(sid[:-1], weights=num, minlength=len(starts))
    sden = np.bincount(sid[:-1], weights=den, minlength=len(starts))
    nval = np.bincount(sid, weights=valid.astype(float), minlength=len(starts))
    sval = np.bincount(sid, weights=np.nan_to_num(T), minlength=len(starts))
    with np.errstate(invalid="ignore", divide="ignore"):
        trap = snum / sden
        smean = sval / nval
    zero = (npair > 0) & (sden == 0)
    if zero_dt == "sample_mean":
        trap = np.where(zero, smean, trap)
    mean = np.where(npair > 0, trap, np.where(nval > 0, smean, np.nan))
    # last valid value of each step
    pos = np.where(valid, np.arange(n), -1)
    lastpos = np.maximum.reduceat(pos, starts) if len(starts) else np.array([], int)
    lastpos = np.where(lastpos >= starts, lastpos, -1)
    last = np.where(lastpos >= 0, T[np.maximum(lastpos, 0)], np.nan)
    return {"mean": mean, "npair": npair, "nvalid": nval, "pair_dt": sden, "last": last,
            "zero_den": zero}


def step_temperature_loop(t, T, s, e):
    """Plain-loop version of step_temperature for one step (cross-check on the V2 cells)."""
    num = den = 0.0
    npair = 0
    for j in range(s, e):
        if not (np.isnan(T[j]) or np.isnan(T[j + 1])):
            num += 0.5 * (T[j] + T[j + 1]) * (t[j + 1] - t[j])
            den += t[j + 1] - t[j]
            npair += 1
    v = [x for x in T[s:e + 1] if not np.isnan(x)]
    if npair and den != 0:
        return num / den
    return float(np.mean(v)) if v else float("nan")
