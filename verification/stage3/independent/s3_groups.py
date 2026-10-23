#!/usr/bin/env python3
"""Stage 3 independent check, step 3: levers (D44), laws over the same secants, the floor (D45), the thermal share
(D46), the module numbers (D47) and the verdicts P1 to P10 (D48), from the per-cell table of s1_cells.py.

Only Stage 1's laws are imported (src/thermal_lever.py, which the brief allows); nothing from src/stage2 or src/stage3.
Usage: python3 -I s3_groups.py <cells_full.json> <out_json>
"""
import csv, itertools, json, math, os, sys
import numpy as np
from scipy import stats

REPO = "/home/claude/thermal-lever-duplicate-floor"
sys.path.insert(0, os.path.join(REPO, "src"))
import thermal_lever as tl  # Stage 1 laws, allowed by the brief

SEED = 20261024
NBOOT = 10000
CHEMS = ["LFP", "NCA", "NMC"]
THR = {"LFP": 0.95, "NCA": 0.90, "NMC": 0.90}
OFFER_K = None  # read from results/module_spread.csv below
GESLIN = json.load(open(os.path.join(REPO, "inputs", "floor.json")))


def module_offer():
    rows = list(csv.DictReader(open(os.path.join(REPO, "results", "module_spread.csv"))))
    vals = [float(r["dT_module_K"]) for r in rows if float(r["crate"]) <= 0.5]
    return max(vals)


def law1_lnL(T_C, soh):
    return tl.ln_L_kucinskis(T_C, 0.5, soh)


def law2_lnL(T_C, c_dis, soh):
    a, b = tl.schimpe_cycle_coeffs(T_C + 273.15, c_dis, 0.5, None, ("cal", "high", "low", "hs"))[:2]
    loss = 1.0 - soh
    u = loss / a if b <= 0 else (-a + math.sqrt(a * a + 4.0 * b * loss)) / (2.0 * b)
    return math.log(u * u)


def secant(f_lo, f_hi, T_lo, T_hi):
    return -(f_hi - f_lo) / (T_hi - T_lo)


def nanmean_or_nan(x):
    x = np.asarray(x, float)
    return float(np.mean(x)) if len(x) else float("nan")


def main(cells_path, out_path):
    cells = json.load(open(cells_path))
    cells.sort(key=lambda r: r["cell_id"])
    by_id = {c["cell_id"]: c for c in cells}
    offer = module_offer()
    out = {"module_offer_K": offer}

    # ------------------------------------------------------------------ the lever (D44)
    def series(chem, cdis, amb):
        return [c for c in cells if c["chemistry"] == chem and c["window"] == "0-100" and c["c_charge"] == 0.5
                and c["c_discharge"] == cdis and c["ambient_C"] == amb]

    rng = np.random.default_rng(SEED)
    rng_alt = np.random.default_rng(SEED)        # alternative stream: the LFP 1C 25-35 secant draws as well
    levers = {}
    stream_log = []
    for cdis in (1.0, 2.0):
        for chem in CHEMS:
            for lo, hi in ((25.0, 35.0), (15.0, 25.0)):
                g_lo, g_hi = series(chem, cdis, lo), series(chem, cdis, hi)
                key = f"{chem}_{int(cdis)}C_{int(lo)}-{int(hi)}"
                cens = any(c["censored"] for c in g_lo + g_hi)
                lnL_lo = np.array([math.log(c["L_efc"]) for c in g_lo]) if not cens else None
                lnL_hi = np.array([math.log(c["L_efc"]) for c in g_hi]) if not cens else None
                Tl = np.array([c["T_cyc_C"] for c in g_lo])
                Th = np.array([c["T_cyc_C"] for c in g_hi])
                rec = {"n_lo": len(g_lo), "n_hi": len(g_hi), "cells_lo": [c["cell_id"] for c in g_lo],
                       "cells_hi": [c["cell_id"] for c in g_hi], "any_censored": cens,
                       "determinable_D44": (not cens) and len(g_lo) > 0 and len(g_hi) > 0}
                if not rec["determinable_D44"]:
                    levers[key] = rec
                    stream_log.append((key, "no draws (censored)"))
                    continue
                T_lo, T_hi = float(np.mean(Tl)), float(np.mean(Th))     # NaN if any cell lacks T_cyc (strict D43)
                m_lo, m_hi = float(np.mean(lnL_lo)), float(np.mean(lnL_hi))
                rec.update({"T_lo_C": T_lo, "T_hi_C": T_hi, "dT_K": T_hi - T_lo, "mean_lnL_lo": m_lo, "mean_lnL_hi": m_hi,
                            "s_per_K": secant(m_lo, m_hi, T_lo, T_hi) if not (math.isnan(T_lo) or math.isnan(T_hi)) else float("nan"),
                            "s_setpoint_per_K": -(m_hi - m_lo) / (hi - lo)})
                temps_defined = not (math.isnan(T_lo) or math.isnan(T_hi))
                rec["group_temperatures_defined"] = temps_defined
                # alternative readings of D43 where a cell has no T_cyc
                Tl_ok, Th_ok = Tl[~np.isnan(Tl)], Th[~np.isnan(Th)]
                T_lo_avail, T_hi_avail = float(np.mean(Tl_ok)), float(np.mean(Th_ok))
                rec["alt_T_over_cells_with_one"] = {"T_lo_C": T_lo_avail, "T_hi_C": T_hi_avail,
                                                    "s_per_K": secant(m_lo, m_hi, T_lo_avail, T_hi_avail)}
                keep_lo, keep_hi = ~np.isnan(Tl), ~np.isnan(Th)
                rec["alt_drop_cells_without_T_cyc"] = {
                    "s_per_K": secant(float(np.mean(lnL_lo[keep_lo])), float(np.mean(lnL_hi[keep_hi])), T_lo_avail, T_hi_avail),
                    "n_lo": int(keep_lo.sum()), "n_hi": int(keep_hi.sum())}
                # laws over the same secant
                soh = THR[chem]
                for name, f in (("law1", lambda T: law1_lnL(T, soh)), ("law2", lambda T: law2_lnL(T, cdis, soh))):
                    if temps_defined:
                        a_, b_ = f(T_lo), f(T_hi)
                        rec[f"{name}_s_per_K"] = secant(a_, b_, T_lo, T_hi) if not (math.isnan(a_) or math.isnan(b_)) else float("nan")
                    else:
                        rec[f"{name}_s_per_K"] = float("nan")
                    a2, b2 = f(T_lo_avail), f(T_hi_avail)
                    rec["alt_T_over_cells_with_one"][f"{name}_s_per_K"] = secant(a2, b2, T_lo_avail, T_hi_avail) if not (math.isnan(a2) or math.isnan(b2)) else float("nan")
                    a3, b3 = f(lo), f(hi)
                    rec[f"{name}_s_setpoint_per_K"] = secant(a3, b3, lo, hi) if not (math.isnan(a3) or math.isnan(b3)) else float("nan")
                # bootstrap
                def boot(gen):
                    sb = np.empty(NBOOT)
                    ssp = np.empty(NBOOT)
                    for it in range(NBOOT):
                        il = gen.integers(0, len(g_lo), len(g_lo))
                        ih = gen.integers(0, len(g_hi), len(g_hi))
                        ml, mh = lnL_lo[il].mean(), lnL_hi[ih].mean()
                        sb[it] = -(mh - ml) / (Th[ih].mean() - Tl[il].mean())
                        ssp[it] = -(mh - ml) / (hi - lo)
                    return sb, ssp
                if temps_defined:
                    sb, ssp = boot(rng)
                    rec["ci95"] = [float(np.percentile(sb, 2.5)), float(np.percentile(sb, 97.5))]
                    rec["ci95_setpoint"] = [float(np.percentile(ssp, 2.5)), float(np.percentile(ssp, 97.5))]
                    rec["ci95_abs_s"] = [float(np.percentile(np.abs(sb), 2.5)), float(np.percentile(np.abs(sb), 97.5))]
                    stream_log.append((key, "drew"))
                    sb2, ssp2 = boot(rng_alt)
                    rec["ci95_alt_stream"] = [float(np.percentile(sb2, 2.5)), float(np.percentile(sb2, 97.5))]
                else:
                    stream_log.append((key, "no draws (group temperature undefined)"))
                    sb2, ssp2 = boot(rng_alt)     # the alternative stream draws here
                    rec["ci95_setpoint_alt_stream"] = [float(np.percentile(ssp2, 2.5)), float(np.percentile(ssp2, 97.5))]
                    # interval of the alternative reading (T over cells that have one), from a separate generator
                    g3 = np.random.default_rng(SEED)
                    sb3 = np.empty(NBOOT)
                    Tl3 = np.where(np.isnan(Tl), np.nan, Tl)
                    Th3 = np.where(np.isnan(Th), np.nan, Th)
                    for it in range(NBOOT):
                        il = g3.integers(0, len(g_lo), len(g_lo))
                        ih = g3.integers(0, len(g_hi), len(g_hi))
                        sb3[it] = -(lnL_hi[ih].mean() - lnL_lo[il].mean()) / (np.nanmean(Th3[ih]) - np.nanmean(Tl3[il]))
                    rec["alt_T_over_cells_with_one"]["ci95_own_seed"] = [float(np.nanpercentile(sb3, 2.5)), float(np.nanpercentile(sb3, 97.5))]
                    rec["alt_T_over_cells_with_one"]["n_boot_nan"] = int(np.isnan(sb3).sum())
                # implied activation energy (reported)
                if temps_defined:
                    Tm = 0.5 * (T_lo + T_hi) + 273.15
                    rec["Ea_kJ_per_mol"] = rec["s_per_K"] * 8.314462618 * Tm * Tm / 1000.0
                levers[key] = rec
    out["levers"] = levers
    out["rng_stream"] = stream_log

    # KIT values beside NMC (D44), read unchanged
    kit = {}
    for r in csv.DictReader(open(os.path.join(REPO, "results", "stage2_lever.csv"))):
        if r["protocol"] in ("0-100 0.33/1.00C", "0-100 1.00/1.00C"):
            kit[f"{r['protocol']} {r['T_lo_set']}-{r['T_hi_set']}"] = float(r["s_per_K"])
    out["kit_beside_nmc"] = kit
    kitfloor = json.load(open(os.path.join(REPO, "results", "stage2_lever.json")))["floor"]
    out["kit_floor"] = {k: kitfloor[k] for k in ("pairs", "mean", "median", "p95")}

    # ------------------------------------------------------------------ pairs
    groups = {}
    for c in cells:
        groups.setdefault(c["group"], []).append(c)
    pairs = []
    for g, cs in groups.items():
        cs = sorted(cs, key=lambda r: r["cell_id"])
        for x, y in itertools.combinations(cs, 2):
            p = {"group": g, "chemistry": x["chemistry"], "c_discharge": x["c_discharge"], "a": x["cell_id"], "b": y["cell_id"],
                 "censored_any": bool(x["censored"] or y["censored"])}
            if not p["censored_any"]:
                La, Lb = x["L_efc"], y["L_efc"]
                p["rel_diff"] = abs(La - Lb) / (0.5 * (La + Lb))
                p["abs_diff_efc"] = abs(La - Lb)
                p["resolution_efc"] = max(x["resolution_efc"], y["resolution_efc"])
                p["below_resolution"] = bool(abs(La - Lb) < p["resolution_efc"])
            pairs.append(p)
    out["pairs_all"] = pairs

    # ------------------------------------------------------------------ the floor (D45)
    floor = {}
    for chem in CHEMS:
        pc = [p for p in pairs if p["chemistry"] == chem]
        obs = [p for p in pc if not p["censored_any"]]
        d = np.array([p["rel_diff"] for p in obs])
        gset = sorted({p["group"] for p in pc})
        floor[chem] = {"pairs_total": len(pc), "groups_with_pairs": len(gset), "pairs_observed": len(obs),
                       "groups_observed": len({p["group"] for p in obs}),
                       "mean": float(d.mean()) if len(d) else float("nan"),
                       "median": float(np.median(d)) if len(d) else float("nan"),
                       "p95": float(np.percentile(d, 95)) if len(d) else float("nan"),
                       "pairs_below_resolution": int(sum(p["below_resolution"] for p in obs)),
                       "pairs_censored": len(pc) - len(obs)}
        # sensitivity, unregistered: without the three NMC 20-80 % cells whose crossing rests on a 1.8 Ah cycling discharge
    out["floor"] = floor
    out["geslin_floor"] = {"mean": GESLIN["mean"], "p95": GESLIN["p95"], "median": GESLIN["median"], "n_pairs": GESLIN["n_pairs"]}

    # ------------------------------------------------------------------ thermal share (D46)
    s_reg = {chem: levers[f"{chem}_1C_25-35"].get("s_per_K", float("nan")) for chem in CHEMS}
    out["s_registered"] = s_reg

    def orient(measure_key, s_map, include_ties_in_slope=True):
        res = []
        for p in pairs:
            if p["c_discharge"] == 0.5:
                continue
            x, y = by_id[p["a"]], by_id[p["b"]]
            mx, my = x[measure_key], y[measure_key]
            q = {"group": p["group"], "chemistry": p["chemistry"], "c_discharge": p["c_discharge"], "a": p["a"], "b": p["b"]}
            if mx is None or my is None or math.isnan(mx) or math.isnan(my):
                q["status"] = "no_measure"
                res.append(q)
                continue
            w, c = (x, y) if mx > my else (y, x)
            q["warmer"], q["cooler"] = w["cell_id"], c["cell_id"]
            q["dT"] = abs(mx - my)
            s = s_map[p["chemistry"]]
            q["predicted"] = -s * q["dT"] if not math.isnan(s) else float("nan")
            if not w["censored"] and not c["censored"]:
                q["observed"] = math.log(w["L_efc"] / c["L_efc"])
                q["resolution"] = max(w["resolution_efc"], c["resolution_efc"])
                q["tie"] = bool(abs(w["L_efc"] - c["L_efc"]) < q["resolution"])
                q["enters_sign_test"] = not q["tie"]
                q["warmer_shorter"] = bool(w["L_efc"] < c["L_efc"])
                q["status"] = "both_observed"
            elif w["censored"] and c["censored"]:
                q["status"] = "both_censored"
                q["enters_sign_test"] = False
            else:
                cen, unc = (w, c) if w["censored"] else (c, w)
                q["status"] = "one_censored"
                q["enters_sign_test"] = bool(cen["censor_efc"] > unc["L_efc"])
                q["warmer_shorter"] = bool(c["censored"])   # the censored cell lives longer
            res.append(q)
        return res

    def sign(sub):
        ent = [q for q in sub if q.get("enters_sign_test")]
        n = len(ent)
        k = sum(q["warmer_shorter"] for q in ent)
        if n == 0:
            return {"n_entering": 0, "k_warmer_shorter": 0, "share": float("nan"), "p_one_sided": float("nan"),
                    "p_two_sided": float("nan"), "n_ties": sum(1 for q in sub if q.get("tie")), "n_pairs": len(sub)}
        return {"n_entering": n, "k_warmer_shorter": int(k), "share": k / n,
                "p_one_sided": float(stats.binomtest(k, n, 0.5, alternative="greater").pvalue),
                "p_two_sided": float(stats.binomtest(k, n, 0.5, alternative="two-sided").pvalue),
                "n_ties": sum(1 for q in sub if q.get("tie")), "n_pairs": len(sub)}

    def slope_items(sub, include_ties=True):
        return [q for q in sub if q.get("status") == "both_observed" and not math.isnan(q["predicted"])
                and (include_ties or not q["tie"])]

    def slope(items):
        P = np.array([q["predicted"] for q in items])
        O = np.array([q["observed"] for q in items])
        den = float(np.sum(P * P))
        return (float(np.sum(P * O)) / den if den > 0 else float("nan"),
                float(math.sqrt(np.mean(P * P)) / math.sqrt(np.mean(O * O))) if len(items) and np.mean(O * O) > 0 else float("nan"))

    def boot_slope(items, gen):
        gl = sorted({q["group"] for q in items})
        byg = {g: [q for q in items if q["group"] == g] for g in gl}
        P = {g: np.array([q["predicted"] for q in byg[g]]) for g in gl}
        O = {g: np.array([q["observed"] for q in byg[g]]) for g in gl}
        PP = np.array([np.sum(P[g] * P[g]) for g in gl])
        PO = np.array([np.sum(P[g] * O[g]) for g in gl])
        sb = np.empty(NBOOT)
        for it in range(NBOOT):
            ix = gen.integers(0, len(gl), len(gl))
            den = PP[ix].sum()
            sb[it] = PO[ix].sum() / den if den > 0 else np.nan
        return [float(np.nanpercentile(sb, 2.5)), float(np.nanpercentile(sb, 97.5))], len(gl), int(np.isnan(sb).sum())

    def one_per_group(sub):
        seen = {}
        for q in sorted(sub, key=lambda q: (q["group"], q["a"], q["b"])):
            if q["status"] == "no_measure":
                continue
            seen.setdefault(q["group"], q)
        return list(seen.values())

    share = {}
    oriented = {}
    for meas, key in (("a", "measure_a_K"), ("b", "measure_b_C")):
        ori = orient(key, s_reg)
        oriented[meas] = ori
        subsets = [("all", ori),
                   ("rate_1C_or_below", [q for q in ori if q["c_discharge"] <= 1.0]),
                   ("rate_2C_and_3C", [q for q in ori if q["c_discharge"] >= 2.0])]
        subsets += [(f"chem_{ch}", [q for q in ori if q["chemistry"] == ch]) for ch in CHEMS]
        subsets += [("one_pair_per_group", one_per_group(ori))]
        res = {}
        for name, sub in subsets:
            rr = sign(sub)
            rr["n_no_measure"] = sum(1 for q in sub if q["status"] == "no_measure")
            items = slope_items(sub, True)
            sl, rms = slope(items) if items else (float("nan"), float("nan"))
            rr.update({"slope": sl, "rms_pred_over_obs": rms, "n_slope_pairs": len(items)})
            items_nt = slope_items(sub, False)
            if items_nt:
                rr["slope_untied_only"], rr["rms_untied_only"] = slope(items_nt)
            if items:
                ci, ng, nnan = boot_slope(items, rng)
                rr.update({"slope_ci95": ci, "slope_boot_groups": ng, "slope_boot_nan": nnan})
                ci2, _, _ = boot_slope(items, rng_alt)
                rr["slope_ci95_alt_stream"] = ci2
            res[name] = rr
        share[meas] = res
    out["thermal_share"] = share
    out["oriented_pairs"] = oriented

    # ------------------------------------------------------------------ P6: self-heating differences at 1C or below
    p6 = [q for q in oriented["a"] if q["c_discharge"] <= 1.0 and q["status"] != "no_measure"]
    dts = np.array([q["dT"] for q in p6])
    out["P6"] = {"n_pairs": len(p6), "n_below_1K": int(np.sum(dts < 1.0)), "share_below_1K": float(np.mean(dts < 1.0)),
                 "median_dT_K": float(np.median(dts)), "max_dT_K": float(dts.max()),
                 "n_pairs_without_measure": sum(1 for q in oriented["a"] if q["c_discharge"] <= 1.0 and q["status"] == "no_measure")}

    # ------------------------------------------------------------------ module (D47)
    sN = levers["NCA_1C_25-35"]
    mod = {"s_NCA": sN["s_per_K"], "offer_K": offer, "divergence": abs(sN["s_per_K"]) * offer,
           "divergence_ci95_from_s_ci": sorted([abs(sN["ci95"][0]) * offer, abs(sN["ci95"][1]) * offer]),
           "divergence_ci95_from_abs_s": [sN["ci95_abs_s"][0] * offer, sN["ci95_abs_s"][1] * offer]}
    for name, sub in (("NCA_pairs", [q for q in oriented["a"] if q["chemistry"] == "NCA" and q["status"] != "no_measure"]),
                      ("all_pairs", [q for q in oriented["a"] if q["status"] != "no_measure"])):
        v = np.array([abs(sN["s_per_K"]) * q["dT"] for q in sub])
        mod[name] = {"n": len(v), "mean": float(v.mean()), "median": float(np.median(v)), "max": float(v.max()),
                     "p95": float(np.percentile(v, 95)),
                     "n_above_geslin_mean": int(np.sum(v > GESLIN["mean"])), "n_above_geslin_p95": int(np.sum(v > GESLIN["p95"])),
                     "n_above_half_geslin_mean": int(np.sum(v > 0.5 * GESLIN["mean"]))}
    out["module"] = mod

    # ------------------------------------------------------------------ verdicts (D48)
    v = {}
    L = levers
    sNMC, sNMC_lo = L["NMC_1C_25-35"]["s_per_K"], L["NMC_1C_15-25"]["s_per_K"]
    sNCA, sLFP = L["NCA_1C_25-35"]["s_per_K"], L["LFP_1C_25-35"]["s_per_K"]
    det = lambda x: not math.isnan(x)
    v["P1"] = ("holds" if -0.035 <= sNMC <= 0.005 else "did not hold") if det(sNMC) else "not scored"
    v["P2"] = (("holds" if (sNMC_lo < 0 and abs(sNMC_lo) > abs(sNMC)) else "did not hold")
               if det(sNMC) and det(sNMC_lo) else "not scored")
    v["P3"] = ("holds" if 0.0 <= sNCA <= 0.05 else "did not hold") if det(sNCA) else "not scored"
    law2_lfp = L["LFP_1C_25-35"]["law2_s_per_K"]
    if det(sLFP) and det(sNCA) and det(sNMC) and det(law2_lfp):
        ratio = sLFP / law2_lfp
        v["P4"] = "holds" if (sLFP > sNCA and sLFP > sNMC and np.sign(sLFP) == np.sign(law2_lfp) and 1 / 3 <= ratio <= 3) else "did not hold"
    else:
        v["P4"] = "not scored"
    ok5 = [floor[ch]["pairs_observed"] >= 5 for ch in CHEMS]
    v["P5"] = "holds" if all((0.015 <= floor[ch]["mean"] <= 0.10) for ch in CHEMS if floor[ch]["pairs_observed"] >= 5) else "did not hold"
    v["P6"] = "holds" if out["P6"]["share_below_1K"] >= 0.8 else "did not hold"
    sa = share["a"]["all"]
    v["P7"] = ("holds" if 0.35 <= sa["share"] <= 0.65 else "did not hold") if sa["n_entering"] > 0 else "not scored"
    r1, r2 = share["a"]["rate_1C_or_below"], share["a"]["rate_2C_and_3C"]
    if r1["n_entering"] >= 5 and r2["n_entering"] >= 5:
        v["P8"] = "holds" if r2["share"] > r1["share"] else "did not hold"
    else:
        v["P8"] = "not scored"
    v["P9"] = ("holds" if abs(sNCA) * offer < 0.015 else "did not hold") if det(sNCA) else "not scored"
    lows = [L[f"{ch}_1C_15-25"]["s_per_K"] for ch in CHEMS]
    v["P10"] = ("holds" if all(x < 0 for x in lows) else "did not hold") if all(det(x) for x in lows) else "not scored"
    out["verdicts"] = v

    # P4 under the alternative readings (point (i) of the brief)
    alt = {}
    lf = L["LFP_1C_25-35"]
    for name, sl, law2 in (("T_cyc_over_cells_with_one", lf["alt_T_over_cells_with_one"]["s_per_K"], lf["alt_T_over_cells_with_one"]["law2_s_per_K"]),
                           ("setpoint_10K", lf["s_setpoint_per_K"], lf["law2_s_setpoint_per_K"]),
                           ("drop_cell_without_T_cyc", lf["alt_drop_cells_without_T_cyc"]["s_per_K"], lf["alt_T_over_cells_with_one"]["law2_s_per_K"])):
        # compare against NCA and NMC under the same reading (setpoint: their setpoint s; otherwise their registered s)
        sN_, sM_ = (L["NCA_1C_25-35"]["s_setpoint_per_K"], L["NMC_1C_25-35"]["s_setpoint_per_K"]) if name == "setpoint_10K" else (sNCA, sNMC)
        ratio = sl / law2
        alt[name] = {"s_LFP": sl, "law2_s": law2, "ratio_to_law2": ratio, "s_NCA": sN_, "s_NMC": sM_,
                     "exceeds_both": bool(sl > sN_ and sl > sM_), "same_sign_as_law2": bool(np.sign(sl) == np.sign(law2)),
                     "within_factor_3": bool(1 / 3 <= ratio <= 3),
                     "verdict": "holds" if (sl > sN_ and sl > sM_ and np.sign(sl) == np.sign(law2) and 1 / 3 <= ratio <= 3) else "did not hold"}
    out["P4_alternatives"] = alt

    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps({"verdicts": v, "s_registered": s_reg, "P4_alt": alt}, indent=1, default=float))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
