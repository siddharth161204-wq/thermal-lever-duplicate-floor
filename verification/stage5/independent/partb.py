#!/usr/bin/env python3
"""Independent Stage 5 check, Part B (D74): beta, the within-set lever s_w, the sign test and the residual noise,
with the brackets. Bootstrap over sets (KIT parameter sets, Sandia groups), 10,000 draws, percentile 95 % interval
(numpy's linear interpolation), the pooled draws stratified by dataset. The checker's own generator (SFC64, seeds
[BASE_SEED, k]); B4's interval comes from B1's draws, with beta refitted in each draw.

Outputs: out/B.json, out/B.csv.
Usage: PYTHONDONTWRITEBYTECODE=1 python3 verification/stage5/independent/partb.py
"""
import collections, json, math, os, sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ind_common as C   # noqa: E402

N_DRAWS = 10_000
Q = (2.5, 97.5)


def as_xy(pairs, measure="a", key=None):
    """(set key, dT, d, p) per pair; measure (b) for the Sandia bracket."""
    out = []
    for x in pairs:
        if measure == "a":
            dT, d, p = x["dT_K"], x["d"], x["p"]
        else:
            dT, d, p = x["dT_b_K"], x["d_b"], x["p_b"]
        out.append(((x["dataset"], x["set"]) if key is None else key(x), dT, d, p))
    return out


def stratum(xy):
    """Per-set sufficient statistics, sets in first-appearance order."""
    order, acc = [], {}
    for k, dT, d, p in xy:
        if k not in acc:
            order.append(k)
            acc[k] = {"n": 0, "pd": 0.0, "pp": 0.0, "dd": 0.0, "dTd": 0.0, "dTdT": 0.0, "has_p": True}
        a = acc[k]
        a["n"] += 1
        a["dd"] += d * d
        a["dTd"] += dT * d
        a["dTdT"] += dT * dT
        if p is None:
            a["has_p"] = False
        else:
            a["pd"] += p * d
            a["pp"] += p * p
    S = {f: np.array([acc[k][f] for k in order], dtype=float) for f in ("n", "pd", "pp", "dd", "dTd", "dTdT")}
    S["keys"] = order
    S["has_p"] = all(acc[k]["has_p"] for k in order)
    return S


def point(xy):
    dT = np.array([x[1] for x in xy])
    d = np.array([x[2] for x in xy])
    out = {"n_pairs": len(xy), "n_sets": len({x[0] for x in xy})}
    out["s_w_per_K"] = float(-np.sum(dT * d) / np.sum(dT * dT)) if len(xy) else None
    if len(xy) and all(x[3] is not None for x in xy):
        p = np.array([x[3] for x in xy])
        beta = float(np.sum(p * d) / np.sum(p * p))
        r = d - beta * p
        out["beta"] = beta
        out["resid_sd"] = float(math.sqrt(np.sum(r * r) / (len(xy) - 1))) if len(xy) > 1 else None
    else:
        out["beta"] = out["resid_sd"] = None
    return out


def bootstrap(strata, seed_key, n_draws=N_DRAWS):
    rng = C.rng_for(*seed_key)
    tot = {f: np.zeros(n_draws) for f in ("n", "pd", "pp", "dd", "dTd", "dTdT")}
    for S in strata:
        m = len(S["n"])
        idx = rng.integers(0, m, size=(n_draws, m))
        for f in tot:
            tot[f] += S[f][idx].sum(axis=1)
    out = {"draws": n_draws, "seed": [C.BASE_SEED, *seed_key], "generator": "numpy Generator(SFC64(SeedSequence(seed)))",
           "strata_sets": [len(S["n"]) for S in strata]}
    ok_sw = tot["dTdT"] > 0
    sw = -tot["dTd"][ok_sw] / tot["dTdT"][ok_sw]
    out["s_w_left_out"] = int((~ok_sw).sum())
    out["s_w_ci95"] = [float(v) for v in np.percentile(sw, Q)] if sw.size else None
    if all(S["has_p"] for S in strata):
        ok = tot["pp"] > 0
        beta = tot["pd"][ok] / tot["pp"][ok]
        out["beta_left_out"] = int((~ok).sum())
        out["beta_ci95"] = [float(v) for v in np.percentile(beta, Q)]
        n = tot["n"][ok]
        ss = tot["dd"][ok] - 2.0 * beta * tot["pd"][ok] + beta * beta * tot["pp"][ok]
        ok2 = n > 1
        sd = np.sqrt(np.maximum(ss[ok2], 0.0) / (n[ok2] - 1.0))
        out["resid_sd_left_out"] = int((~ok2).sum()) + out["beta_left_out"]
        out["resid_sd_ci95"] = [float(v) for v in np.percentile(sd, Q)]
        out["beta_draws_mean"] = float(beta.mean())
        out["beta_draws_sd"] = float(beta.std(ddof=1))
    return out


def sign_test(xy):
    agree = left = 0
    n = 0
    for _, dT, d, p in xy:
        if d == 0 or p == 0 or p is None:
            left += 1
            continue
        n += 1
        agree += int(np.sign(d) == np.sign(p))
    pv = float(stats.binomtest(agree, n, 0.5, alternative="two-sided").pvalue) if n else None
    return {"pairs_tested": n, "agree": agree, "left_out": left, "share": agree / n if n else None, "p_two_sided": pv}


def analysis(name, strata_pairs, seed_key, measure="a", key=None, registered=True, note=""):
    xys = [as_xy(ps, measure, key) for ps in strata_pairs]
    allxy = [x for xy in xys for x in xy]
    res = {"name": name, "registered": registered, "note": note, "measure": measure}
    res.update(point(allxy))
    res["strata_pairs"] = [len(xy) for xy in xys]
    if len(allxy):
        res["bootstrap"] = bootstrap([stratum(xy) for xy in xys if len(xy)], seed_key)
    if all(x[3] is not None for x in allxy) and len(allxy):
        res["sign"] = sign_test(allxy)
    return res


def floor_verdict(res, shapes, against):
    lo, hi = res["bootstrap"]["resid_sd_ci95"]
    out = {}
    for f in against:
        sd = shapes[f]["sd_e"]
        out[f] = {"sd_e": sd, "inside_interval": bool(lo <= sd <= hi),
                  "point_over_sd_e": res["resid_sd"] / sd}
    return out


def main():
    rows, _ = C.build_pairs()
    sh = C.shapes()
    med = C.band_medians()
    K = C.select(rows, "KIT")
    S = C.select(rows, "Sandia")
    Kw, Kc = C.select(rows, "KIT", "warm"), C.select(rows, "KIT", "cold")
    Sw, Sc = C.select(rows, "Sandia", "warm"), C.select(rows, "Sandia", "cold")
    K0 = C.select(rows, "KIT", scored=False, temperature=0.0)
    Ku = C.select(rows, "KIT", untied=True)
    Kuw = C.select(rows, "KIT", "warm", untied=True)
    nca15 = "NCA_15C_0-100_0.5/1C"
    S_noNCA = C.select(rows, "Sandia", exclude_set={nca15})
    if len(K) != 106 or len(S) != 56 or len(Kw) != 70 or len(Kc) != 36 or len(Sw) != 36 or len(Sc) != 20 or len(K0) != 36:
        C.stop("pair counts differ from D72")

    A = collections.OrderedDict()
    # registered B1 (with B4 and, as extra, s_w over all pairs, from the same draws)
    A["KIT"] = analysis("B1/B4 KIT", [K], (1,))
    A["Sandia"] = analysis("B1/B4 Sandia", [S], (2,))
    A["pooled"] = analysis("B1/B4 pooled (stratified)", [K, S], (3,))
    # registered B2 per band (beta and sd per band from the same draws are extras)
    A["KIT warm"] = analysis("B2 KIT warm", [Kw], (4,))
    A["KIT cold"] = analysis("B2 KIT cold", [Kc], (5,))
    A["Sandia warm"] = analysis("B2 Sandia warm", [Sw], (6,))
    A["Sandia cold"] = analysis("B2 Sandia cold", [Sc], (7,))
    A["pooled warm"] = analysis("B2 pooled warm (stratified)", [Kw, Sw], (8,))
    A["pooled cold"] = analysis("B2 pooled cold (stratified)", [Kc, Sc], (9,))
    A["KIT 0 C"] = analysis("B2 KIT 0 C sets (unscored, no between-set lever)", [K0], (10,), registered=False)
    # brackets (unscored)
    A["untied KIT"] = analysis("bracket: untied pairs, KIT", [Ku], (11,), registered=False,
                               note="bootstrap over the sets holding an untied pair")
    A["untied KIT warm"] = analysis("bracket: untied pairs, KIT warm", [Kuw], (12,), registered=False)
    A["untied pooled"] = analysis("bracket: untied pairs, pooled (Sandia has none)", [Ku, []], (13,), registered=False)
    A["untied pooled warm"] = analysis("bracket: untied pairs, pooled warm", [Kuw, []], (14,), registered=False)
    A["Sandia measure b"] = analysis("bracket: Sandia with measure (b)", [S], (15,), measure="b", registered=False)
    A["Sandia measure b warm"] = analysis("bracket: Sandia with measure (b), warm", [Sw], (151,), measure="b", registered=False)
    A["Sandia measure b cold"] = analysis("bracket: Sandia with measure (b), cold", [Sc], (152,), measure="b", registered=False)
    A["Sandia without NCA 15C 1C"] = analysis("bracket: Sandia without the NCA 15 C 1C pair", [S_noNCA], (16,), registered=False)
    # checker's sensitivities, not in D74 (readings): untied over all sets; protocol clusters
    A["sens untied KIT, all 36 sets drawn"] = analysis_untied_all_sets(rows, (17,))
    proto_key = lambda x: (x["dataset"], x["protocol_or_chemistry"] if x["dataset"] == "KIT" else
                           f"{x['protocol_or_chemistry']} {x['window']} {x['rates']}")
    A["sens KIT, protocol clusters"] = analysis("sensitivity: KIT resampling the 12 protocols (each with its sets)", [K], (20,),
                                                key=proto_key, registered=False)
    A["sens Sandia, protocol clusters"] = analysis("sensitivity: Sandia resampling chemistry x window x rates", [S], (21,),
                                                   key=proto_key, registered=False)
    A["sens pooled, protocol clusters"] = analysis("sensitivity: pooled, protocol clusters", [K, S], (22,), key=proto_key,
                                                   registered=False)
    A["sens KIT warm, protocol clusters"] = analysis("sensitivity: KIT warm, 12 protocols", [Kw], (23,), key=proto_key,
                                                     registered=False)
    A["sens KIT cold, protocol clusters"] = analysis("sensitivity: KIT cold, 12 protocols", [Kc], (24,), key=proto_key,
                                                     registered=False)

    # B4 floor verdicts (registered: KIT vs KIT floor, Geslin beside; Sandia and pooled vs both)
    for name in ("KIT", "Sandia", "pooled", "untied KIT", "untied pooled"):
        A[name]["floor_verdicts"] = floor_verdict(A[name], sh, ("KIT", "Geslin"))
    for name in ("KIT warm", "KIT cold", "Sandia warm", "Sandia cold", "pooled warm", "pooled cold"):
        A[name]["floor_verdicts_extra"] = floor_verdict(A[name], sh, ("KIT", "Geslin"))

    # convert s_w to % per K for reporting, beside the band's median lever
    for name, r in A.items():
        if r.get("s_w_per_K") is not None:
            r["s_w_pct_per_K"] = 100.0 * r["s_w_per_K"]
            ci = r["bootstrap"]["s_w_ci95"]
            r["s_w_ci95_pct_per_K"] = [100.0 * ci[0], 100.0 * ci[1]]
        if r.get("beta") is not None and "beta_ci95" in r.get("bootstrap", {}):
            lo, hi = r["bootstrap"]["beta_ci95"]
            r["beta_ci95"] = [lo, hi]
            r["beta_ci95_width"] = hi - lo
            r["resid_sd_ci95"] = r["bootstrap"]["resid_sd_ci95"]
        for b in ("warm", "cold"):
            if name.endswith(b):
                r["band_median_lever_pct_per_K"] = 100.0 * med[b]["s"]
                r["band_median_lever_label"] = med[b]["label"]

    doc = {"method": "D74 as read by the independent checker", "n_draws": N_DRAWS, "percentiles": Q,
           "base_seed": C.BASE_SEED, "generator": "numpy.random.Generator(SFC64(SeedSequence([BASE_SEED, k])))",
           "sd_e": {k: v["sd_e"] for k, v in sh.items()}, "band_medians": med, "analyses": A}
    C.write_json("B.json", doc)
    flat = []
    for name, r in A.items():
        b = r.get("bootstrap", {})
        sg = r.get("sign", {})
        fv = r.get("floor_verdicts", r.get("floor_verdicts_extra", {}))
        flat.append({"analysis": name, "registered": r["registered"], "measure": r["measure"], "n_pairs": r["n_pairs"],
                     "n_sets": r["n_sets"], "beta": r.get("beta"),
                     "beta_lo": (b.get("beta_ci95") or [None, None])[0], "beta_hi": (b.get("beta_ci95") or [None, None])[1],
                     "s_w_pct_per_K": r.get("s_w_pct_per_K"),
                     "s_w_lo": (r.get("s_w_ci95_pct_per_K") or [None, None])[0],
                     "s_w_hi": (r.get("s_w_ci95_pct_per_K") or [None, None])[1],
                     "band_median_pct_per_K": r.get("band_median_lever_pct_per_K"),
                     "sign_agree": sg.get("agree"), "sign_n": sg.get("pairs_tested"), "sign_left_out": sg.get("left_out"),
                     "sign_share": sg.get("share"), "sign_p_two_sided": sg.get("p_two_sided"),
                     "resid_sd": r.get("resid_sd"), "resid_sd_lo": (b.get("resid_sd_ci95") or [None, None])[0],
                     "resid_sd_hi": (b.get("resid_sd_ci95") or [None, None])[1],
                     "KIT_floor_sd_e_inside": (fv.get("KIT") or {}).get("inside_interval"),
                     "Geslin_floor_sd_e_inside": (fv.get("Geslin") or {}).get("inside_interval"),
                     "seed": json.dumps(b.get("seed")), "beta_left_out": b.get("beta_left_out"),
                     "s_w_left_out": b.get("s_w_left_out")})
    C.write_csv("B.csv", flat)
    for f in flat:
        print({k: (round(v, 6) if isinstance(v, float) else v) for k, v in f.items()
               if k in ("analysis", "n_pairs", "n_sets", "beta", "beta_lo", "beta_hi", "s_w_pct_per_K", "s_w_lo", "s_w_hi",
                        "sign_agree", "sign_n", "sign_p_two_sided", "resid_sd", "resid_sd_lo", "resid_sd_hi",
                        "KIT_floor_sd_e_inside", "Geslin_floor_sd_e_inside")})


def analysis_untied_all_sets(rows, seed_key):
    """Reading not taken (sensitivity): the untied bracket drawing all 36 KIT sets, sets without an untied pair adding
    nothing to a draw."""
    K = C.select(rows, "KIT")
    keys = []
    for x in K:
        if x["set"] not in keys:
            keys.append(x["set"])
    acc = {k: {"n": 0, "pd": 0.0, "pp": 0.0, "dd": 0.0, "dTd": 0.0, "dTdT": 0.0} for k in keys}
    for x in K:
        if x["tie"]:
            continue
        a = acc[x["set"]]
        a["n"] += 1
        a["pd"] += x["p"] * x["d"]
        a["pp"] += x["p"] ** 2
        a["dd"] += x["d"] ** 2
        a["dTd"] += x["dT_K"] * x["d"]
        a["dTdT"] += x["dT_K"] ** 2
    S = {f: np.array([acc[k][f] for k in keys]) for f in ("n", "pd", "pp", "dd", "dTd", "dTdT")}
    S["has_p"] = True
    res = {"name": "sensitivity: untied KIT pairs, drawing all 36 sets", "registered": False, "measure": "a",
           "note": "draws with no untied pair have sum p^2 = 0 and are left out and counted"}
    res.update(point(as_xy([x for x in K if not x["tie"]])))
    res["bootstrap"] = bootstrap([S], seed_key)
    return res


if __name__ == "__main__":
    main()
