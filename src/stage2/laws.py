#!/usr/bin/env python3
"""Stage 2: the two published laws over a temperature secant, evaluated with the Stage 1 code (src/thermal_lever.py,
imported unchanged).

Law 1 (Kucinskis et al. 2022): ln of the printed cycles to the SoH threshold, interpolated linearly in temperature
and in charge rate exactly as Stage 1 does (ln_L_kucinskis); the paper tested one discharge rate, so the discharge
rate does not enter. Law 2 (Schimpe et al. 2018): Stage 1's cycle coefficients for CC charge at the charge rate from
0 to 100 % SOC and CC discharge at the discharge rate (schimpe_cycle_coeffs), solved for the full cycles to the
capacity loss 1 - SoH threshold, a u + b u^2 = loss with u = sqrt(N), which is Stage 1's lifetime_law2 when the
threshold is 0.90. The secant is s = -(ln L(T2) - ln L(T1)) / (T2 - T1), per kelvin.
"""
import math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import thermal_lever as TL   # noqa: E402  (Stage 1, unchanged)


def ln_L_law1(T_C, soh, c_charge=0.5):
    return TL.ln_L_kucinskis(T_C, c_charge, soh)


def L_law2(T_C, c_discharge, soh, c_charge=0.5, terms=("cal", "high", "low", "hs")):
    a, b, _, _ = TL.schimpe_cycle_coeffs(T_C + 273.15, c_discharge, c_charge, None, terms)
    loss = 1.0 - soh
    u = loss / a if b <= 0 else (-a + math.sqrt(a * a + 4.0 * b * loss)) / (2.0 * b)
    return u * u


def ln_L_law2(T_C, c_discharge, soh, c_charge=0.5):
    return math.log(L_law2(T_C, c_discharge, soh, c_charge))


def secant(ln_fn, T1, T2):
    return -(ln_fn(T2) - ln_fn(T1)) / (T2 - T1)


def law_secants(T1, T2, c_discharge, soh, c_charge=0.5):
    return {"law_1_kucinskis2022": secant(lambda T: ln_L_law1(T, soh, c_charge), T1, T2),
            "law_2_schimpe2018": secant(lambda T: ln_L_law2(T, c_discharge, soh, c_charge), T1, T2)}


def law_curves(temps_C, c_discharge, soh, c_charge=0.5):
    return [{"T_C": T, "ln_L_law1": ln_L_law1(T, soh, c_charge), "ln_L_law2": ln_L_law2(T, c_discharge, soh, c_charge)}
            for T in temps_C]


def v3_against_stage1(lever_json):
    """V3: the Stage 2 secant over 34.5 to 35.5 C at Stage 1's conditions (C/2 charge, each Stage 1 rate, 90 % SoH)
    against Stage 1's registered s (results/lever.json), and Law 2's lifetime against Stage 1's lifetime_law2."""
    rows = []
    for law, entries in lever_json["levers"].items():
        for e in entries:
            c = e["crate"]
            if law == "law_1_kucinskis2022":
                s2 = secant(lambda T: ln_L_law1(T, 0.90, 0.5), 34.5, 35.5)
            else:
                s2 = secant(lambda T: ln_L_law2(T, c, 0.90, 0.5), 34.5, 35.5)
            rows.append({"law": law, "crate": c, "s_stage1": e["s_per_K"], "s_stage2": s2,
                         "rel_diff": abs(s2 - e["s_per_K"]) / abs(e["s_per_K"])})
    lt = [{"crate": c, "L_stage1": TL.lifetime_law2(308.15, c), "L_stage2": L_law2(35.0, c, 0.90)} for c in TL.CRATES]
    for r in lt:
        r["rel_diff"] = abs(r["L_stage2"] - r["L_stage1"]) / r["L_stage1"]
    return {"secants": rows, "law2_lifetimes_35C": lt,
            "max_rel_diff_secant": max(r["rel_diff"] for r in rows),
            "max_rel_diff_lifetime": max(r["rel_diff"] for r in lt)}
