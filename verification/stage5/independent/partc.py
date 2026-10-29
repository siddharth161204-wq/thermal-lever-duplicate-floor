#!/usr/bin/env python3
"""Independent Stage 5 check, Part C (D75): could the data have shown it.

The test: slope b of d on dT through the origin, t = b sqrt(sum dT^2) / s_r, s_r^2 = sum (d - b dT)^2 / (n - 1),
two-sided 5 % against t with n - 1 df.
C1 calculation: noncentral t, df n - 1, nc = |s| sqrt(sum dT^2) / sd(e); the detectable |s| is where power = 0.80.
  Implemented here with scipy.stats.nct directly (both tails), cross-checked by quadrature over the chi distribution.
C1 simulation: 10,000 synthetic datasets with the observed dT and e from Stage 4's registered shape (log-normal |e|
  with the floor's mean and 95th percentile, smaller root, random sign); the residuals do not depend on s, so one set
  of draws serves every s; the simulated detectable |s| by bisection to 1e-6 relative (and the exact step location).
Beside: normal e of the same sd(e); the normal-shape calculation (sd = m sqrt(pi/2)); per band; KIT over 142 pairs;
  the shared-cell factor F and its line (sd(e) sqrt(F)).
C2: per band, power at the band's median between-set |s|, the pairs needed (calculation, simulation with dT drawn
  with replacement from the band, Stage 4's bisection between 2 and 4 n_calc) and the shared-cell number.

Outputs: out/C.json, out/C1.csv, out/C2.csv.
Usage: PYTHONDONTWRITEBYTECODE=1 python3 verification/stage5/independent/partc.py
"""
import collections, json, math, os, sys, time

import numpy as np
from scipy import integrate, optimize, stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ind_common as C   # noqa: E402

ALPHA, POWER, N_SIM = 0.05, 0.80, 10_000
NAN_LOWER = {"count": 0}


# ------------------------------------------------------------------ noncentral t power (own implementation)
def tcrit(df):
    return float(stats.t.ppf(1.0 - ALPHA / 2.0, df))


def power_quad(df, nc):
    """P(|T| > tc), T = (Z + nc) / W, W = sqrt(V / df), V ~ chi2(df), by quadrature over the quantile of V:
    E[ Phi(nc - tc W) + Phi(-tc W - nc) ]."""
    tc = tcrit(df)

    def f(u):
        v = stats.chi2.ppf(u, df)
        w = math.sqrt(v / df)
        return stats.norm.cdf(nc - tc * w) + stats.norm.cdf(-tc * w - nc)
    val, _ = integrate.quad(f, 0.0, 1.0, limit=400, epsabs=1e-14, epsrel=1e-13)
    return val


def power_nct(df, nc):
    tc = tcrit(df)
    up = float(stats.nct.sf(tc, df, nc))
    lo = float(stats.nct.cdf(-tc, df, nc))
    if not math.isfinite(lo):
        NAN_LOWER["count"] += 1
        # the lower tail is below P(Z + nc < 0) = Phi(-nc); take it from the quadrature of the lower term alone
        lo = float(integrate.quad(lambda u: stats.norm.cdf(-tc * math.sqrt(stats.chi2.ppf(u, df) / df) - nc),
                                  0.0, 1.0, limit=400, epsabs=1e-16)[0])
    if not math.isfinite(up):
        C.stop(f"upper tail undefined at df={df}, nc={nc}")
    return up + lo


def nc_for_power(df, target=POWER):
    return optimize.brentq(lambda x: power_nct(df, x) - target, 0.0, 60.0, xtol=1e-15, rtol=1e-14, maxiter=500)


def detectable_s(n, ss, sd):
    """|s| at power 0.80 for n pairs with sum dT^2 = ss and noise sd."""
    nc = nc_for_power(n - 1)
    return nc * sd / math.sqrt(ss), nc


def n_required(effect, sd):
    """Smallest n >= 2 with power >= 0.80, power at n = noncentral t, df n - 1, nc sqrt(n) effect / sd (Stage 4's
    definition, own search: start at the normal approximation and walk)."""
    pw = lambda n: power_nct(n - 1, math.sqrt(n) * effect / sd)
    if pw(2) >= POWER:
        return 2
    z = stats.norm.ppf(1 - ALPHA / 2) + stats.norm.ppf(POWER)
    n = max(2, int((z * sd / effect) ** 2))
    if pw(n) >= POWER:
        while n > 2 and pw(n - 1) >= POWER:
            n -= 1
    else:
        while pw(n) < POWER:
            n += 1
    if not (pw(n) >= POWER and (n == 2 or pw(n - 1) < POWER)):
        C.stop("n_required search")
    return n


# ------------------------------------------------------------------ simulation
def draw_e(rng, shape, size, kind):
    if kind == "lognormal":
        mag = np.exp(shape["mu"] + shape["sigma"] * rng.standard_normal(size))
        sign = np.where(rng.random(size) < 0.5, -1.0, 1.0)
        return sign * mag
    return shape["sd_e"] * rng.standard_normal(size)


class FixedDesign:
    """10,000 datasets at the observed dT; the slope's noise part and the residual sd per dataset."""

    def __init__(self, dT, shape, kind, rng, n_sim=N_SIM):
        self.dT = np.asarray(dT, float)
        self.n = self.dT.size
        self.ss = float(np.sum(self.dT ** 2))
        e = draw_e(rng, shape, (n_sim, self.n), kind)
        self.be = e @ self.dT / self.ss
        r = e - self.be[:, None] * self.dT[None, :]
        self.sr = np.sqrt(np.sum(r * r, axis=1) / (self.n - 1))
        self.tc = tcrit(self.n - 1)
        self.k = self.tc * self.sr / math.sqrt(self.ss)
        self.n_sim = n_sim

    def power(self, s):
        t = (s + self.be) * math.sqrt(self.ss) / self.sr
        return float(np.mean(np.abs(t) > self.tc))

    def wrong_sign_share(self, s):
        t = (s + self.be) * math.sqrt(self.ss) / self.sr
        return float(np.mean(t < -self.tc))

    def detectable_bisection(self, s_hint):
        lo, hi = 0.0, 2.0 * s_hint
        if self.power(lo) >= POWER:
            C.stop("power at zero")
        while self.power(hi) < POWER:
            lo, hi = hi, 2.0 * hi
        while hi - lo > 1e-6 * hi:
            mid = 0.5 * (lo + hi)
            if self.power(mid) >= POWER:
                hi = mid
            else:
                lo = mid
        return hi

    def detectable_exact(self):
        """The step: smallest s >= 0 (as a limit from above) with power >= 0.80; U_i is where draw i starts rejecting
        upward, V_i where it stops rejecting downward."""
        U = np.sort(self.k - self.be)
        V = -self.k - self.be
        Vs = np.sort(V)
        need = math.ceil(POWER * self.n_sim - 1e-9)
        for j in range(len(U)):
            s = U[j]
            if s < 0:
                continue
            up = j + 1                                  # draws with U < s + 0
            down = len(Vs) - np.searchsorted(Vs, s, side="right")
            if up + down >= need:
                return float(s)
        return None


def sim_power_at_n(n, s, dT_pool, shape, kind, rng, n_sim=N_SIM, chunk=4_000_000):
    tc = tcrit(n - 1)
    rej = done = 0
    per = max(1, chunk // n)
    pool = np.asarray(dT_pool, float)
    while done < n_sim:
        k = min(per, n_sim - done)
        dT = pool[rng.integers(0, pool.size, size=(k, n))]
        e = draw_e(rng, shape, (k, n), kind)
        d = s * dT + e
        ss = np.sum(dT * dT, axis=1)
        b = np.sum(dT * d, axis=1) / ss
        r = d - b[:, None] * dT
        sr = np.sqrt(np.sum(r * r, axis=1) / (n - 1))
        t = b * np.sqrt(ss) / sr
        rej += int(np.sum(np.abs(t) > tc))
        done += k
    return rej / n_sim


def n_simulated(s, dT_pool, shape, kind, n_calc, case_key):
    tried = {}

    def ok(n):
        if n not in tried:
            tried[n] = sim_power_at_n(n, s, dT_pool, shape, kind, C.rng_for(*case_key, n))
        return tried[n] >= POWER

    if ok(2):
        return 2, tried
    lo, hi = 2, max(4 * n_calc, 3)
    while not ok(hi):
        lo, hi = hi, hi * 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ok(mid):
            hi = mid
        else:
            lo = mid
    return hi, tried


# ------------------------------------------------------------------ design quantities
def shared_cell_F(pairs):
    w = collections.defaultdict(float)
    ss = 0.0
    for x in pairs:
        w[(x["dataset"], x["warm"])] += x["dT_K"]
        w[(x["dataset"], x["cool"])] -= x["dT_K"]
        ss += x["dT_K"] ** 2
    return sum(v * v for v in w.values()) / (2.0 * ss), len(w)


def main():
    t0 = time.time()
    rows, _ = C.build_pairs()
    sh = C.shapes()
    med = C.band_medians()
    sel = collections.OrderedDict()
    sel["KIT"] = (C.select(rows, "KIT"), True)
    sel["Sandia"] = (C.select(rows, "Sandia"), True)
    sel["pooled"] = (C.select(rows), True)
    for ds in ("KIT", "Sandia"):
        for b in ("warm", "cold"):
            sel[f"{ds} {b}"] = (C.select(rows, ds, b), False)
    for b in ("warm", "cold"):
        sel[f"pooled {b}"] = (C.select(rows, None, b), False)
    sel["KIT all 142"] = (C.select(rows, "KIT", scored=False), False)
    design = {}
    for name, (ps, reg) in sel.items():
        dT = np.array([x["dT_K"] for x in ps])
        F, ncell = shared_cell_F(ps)
        design[name] = {"n_pairs": len(ps), "n_sets": len({(x['dataset'], x['set']) for x in ps}), "n_cells": ncell,
                        "sum_dT2": float(np.sum(dT ** 2)), "dT_rms_K": float(math.sqrt(np.mean(dT ** 2))),
                        "dT_median_K": float(np.median(dT)), "dT_mean_K": float(np.mean(dT)), "F": F,
                        "registered": reg, "dT": dT}
    print("design", {k: (v["n_pairs"], round(v["dT_rms_K"], 4), round(v["F"], 4)) for k, v in design.items()}, flush=True)

    # ---------------- C1
    c1 = []
    sims = {}
    for ci, (name, dz) in enumerate(design.items()):
        n, ss = dz["n_pairs"], dz["sum_dT2"]
        for fi, (fname, shp) in enumerate(sh.items()):
            s_calc, nc = detectable_s(n, ss, shp["sd_e"])
            s_quad_check = power_quad(n - 1, s_calc * math.sqrt(ss) / shp["sd_e"])
            s_norm_shape, _ = detectable_s(n, ss, shp["sd_e_normal_shape"])
            s_shared = s_calc * math.sqrt(dz["F"])
            key_ln, key_nm = (100 + ci, fi, 1), (100 + ci, fi, 2)
            fd_ln = FixedDesign(dz["dT"], shp, "lognormal", C.rng_for(*key_ln))
            fd_nm = FixedDesign(dz["dT"], shp, "normal", C.rng_for(*key_nm))
            s_sim = fd_ln.detectable_bisection(s_calc)
            s_sim_exact = fd_ln.detectable_exact()
            s_simn = fd_nm.detectable_bisection(s_calc)
            s_simn_exact = fd_nm.detectable_exact()
            sims[(name, fname)] = (fd_ln, fd_nm)
            rel = (s_sim - s_calc) / s_calc
            reln = (s_simn - s_calc) / s_calc
            c1.append({"case": name, "registered": dz["registered"], "floor": fname, "n_pairs": n, "n_sets": dz["n_sets"],
                       "sum_dT2_K2": ss, "dT_rms_K": dz["dT_rms_K"], "dT_median_K": dz["dT_median_K"],
                       "sd_e": shp["sd_e"], "nc_at_power_0.8": nc,
                       "detectable_s_calc_pct_per_K": 100 * s_calc,
                       "power_quadrature_at_calc": s_quad_check,
                       "detectable_s_sim_lognormal_pct_per_K": 100 * s_sim,
                       "detectable_s_sim_lognormal_exact_step_pct_per_K": 100 * s_sim_exact,
                       "v4_rel_diff": rel, "v4_pass": abs(rel) <= 0.05,
                       "sim_lognormal_power_at_calc": fd_ln.power(s_calc),
                       "detectable_s_sim_normal_same_sd_pct_per_K": 100 * s_simn,
                       "detectable_s_sim_normal_exact_step_pct_per_K": 100 * s_simn_exact,
                       "normal_check_rel_diff": reln, "normal_check_pass": abs(reln) <= 0.05,
                       "sd_e_normal_shape": shp["sd_e_normal_shape"],
                       "detectable_s_calc_normal_shape_pct_per_K": 100 * s_norm_shape,
                       "F_shared_cell": dz["F"], "n_cells": dz["n_cells"],
                       "detectable_s_calc_shared_cell_pct_per_K": 100 * s_shared,
                       "sim_size_lognormal_at_zero": fd_ln.power(0.0), "sim_size_normal_at_zero": fd_nm.power(0.0),
                       "seed_lognormal": [C.BASE_SEED, *key_ln], "seed_normal": [C.BASE_SEED, *key_nm]})
            print(f"C1 {name:14s} {fname:6s} n={n:3d} calc={100*s_calc:.5f} sim={100*s_sim:.5f} ({rel:+.4f}) "
                  f"simN={100*s_simn:.5f} ({reln:+.4f}) normal-shape={100*s_norm_shape:.5f} F={dz['F']:.4f} "
                  f"shared={100*s_shared:.5f} quad={s_quad_check:.12f}", flush=True)

    # ---------------- C2
    c2 = []
    for bi, b in enumerate(("warm", "cold")):
        smed = abs(med[b]["s"])
        for di, ds in enumerate(("KIT", "Sandia", "pooled")):
            name = f"{ds} {b}"
            dz = design[name]
            n, ss, rms = dz["n_pairs"], dz["sum_dT2"], dz["dT_rms_K"]
            for fi, (fname, shp) in enumerate(sh.items()):
                sd = shp["sd_e"]
                nc = smed * math.sqrt(ss) / sd
                pw_calc = power_nct(n - 1, nc)
                pw_quad = power_quad(n - 1, nc)
                fd_ln, fd_nm = sims[(name, fname)]
                pw_sim = fd_ln.power(smed)
                pw_sim_n = fd_nm.power(smed)
                pw_shared = power_nct(n - 1, nc / math.sqrt(dz["F"]))
                n_calc = n_required(smed * rms, sd)
                n_shared = n_required(smed * rms, sd * math.sqrt(dz["F"]))
                n_norm_shape = n_required(smed * rms, shp["sd_e_normal_shape"])
                key = (200 + 10 * bi + di, fi, 1)
                t1 = time.time()
                n_sim, tried = n_simulated(smed, dz["dT"], shp, "lognormal", n_calc, key)
                named = (ds == "KIT" and fname == "KIT")
                c2.append({"band": b, "dataset": ds, "floor": fname, "named_in_D75": named,
                           "band_median_label": med[b]["label"], "band_median_abs_s_per_K": smed,
                           "n_pairs_available": n, "dT_rms_K": rms, "sum_dT2_K2": ss, "sd_e": sd, "nc": nc,
                           "power_calc": pw_calc, "power_quadrature": pw_quad, "power_sim_lognormal": pw_sim,
                           "power_sim_normal_same_sd": pw_sim_n, "power_calc_shared_cell": pw_shared,
                           "F_shared_cell": dz["F"],
                           "n_required_calc": n_calc, "n_required_sim_lognormal": n_sim,
                           "n_sim_rel_diff": (n_sim - n_calc) / n_calc,
                           "sim_power_at_n_sim": tried[n_sim], "sim_power_at_n_calc": tried.get(n_calc),
                           "n_tried": sorted(tried), "n_required_shared_cell": n_shared,
                           "n_required_calc_normal_shape": n_norm_shape,
                           "seed_n_sim": [C.BASE_SEED, *key, "n"]})
                print(f"C2 {name:12s} {fname:6s} power calc={pw_calc:.5f} quad={pw_quad:.5f} sim={pw_sim:.4f} "
                      f"n_calc={n_calc} n_sim={n_sim} n_shared={n_shared} ({time.time()-t1:.0f} s)", flush=True)

    doc = {"method": "D75 as read by the independent checker", "alpha": ALPHA, "power": POWER, "n_sim": N_SIM,
           "base_seed": C.BASE_SEED, "generator": "numpy.random.Generator(SFC64(SeedSequence([BASE_SEED, ...])))",
           "nct_lower_tail_nan_count": NAN_LOWER["count"],
           "design": {k: {kk: vv for kk, vv in v.items() if kk != "dT"} for k, v in design.items()},
           "shapes": {k: {kk: v[kk] for kk in ("sigma", "mu", "sd_e", "sd_e_normal_shape")} for k, v in sh.items()},
           "band_medians": med, "C1": c1, "C2": c2,
           "v4": {"cases": len(c1), "pass": sum(r["v4_pass"] for r in c1),
                  "max_abs_rel_diff": max(abs(r["v4_rel_diff"]) for r in c1),
                  "normal_check_pass": sum(r["normal_check_pass"] for r in c1),
                  "normal_check_max_abs_rel_diff": max(abs(r["normal_check_rel_diff"]) for r in c1)}}
    C.write_json("C.json", doc)
    C.write_csv("C1.csv", c1)
    C.write_csv("C2.csv", [{k: (json.dumps(v) if isinstance(v, list) else v) for k, v in r.items()} for r in c2])
    print(json.dumps(doc["v4"]), "nan lower tails:", NAN_LOWER["count"], f"{time.time()-t0:.0f} s")


if __name__ == "__main__":
    main()
