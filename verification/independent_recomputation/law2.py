"""Law 2 (Schimpe et al. 2018, J. Electrochem. Soc. 165, A181, D4), written from the equation forms and Table IV
values transcribed in inputs/schimpe2018_parameters.json (the paper itself could not be fetched in this session).

  Q_loss = k_Cal(T,SOC) sqrt(t) + k_HighT(T) sqrt(Q_Tot) + k_LowT(T,I_Ch) sqrt(Q_Ch) + k_LowT,HighSOC(T,I_Ch,SOC) Q_Ch
  k_Cal   = k_Cal,Ref exp[-Ea_Cal/Rg (1/T - 1/Tref)] (exp[alpha F (U_a,Ref - U_a(SOC)) / (Rg Tref)] + k0)      (Eq. 9)
  k_HighT = k_HighT,Ref exp[-Ea_HighT/Rg (1/T - 1/Tref)]                                                       (Eq. 15)
  k_LowT  = k_LowT,Ref exp[+Ea_LowT/Rg (1/T - 1/Tref)] exp[beta_LowT (I_Ch - I_Ch,Ref)/C0]                     (Eq. 18)
  k_HS    = k_HS,Ref exp[+Ea_HS/Rg (1/T - 1/Tref)] exp[beta_HS (I_Ch - I_Ch,Ref)/C0] (sgn(SOC - SOC_Ref)+1)/2   (Eq. 21)
  U_a(x)  = Eq. A1;  x_a(SOC) = x0 + SOC (x100 - x0)  (Eq. A2, Table AI)

Protocol mapping (D4): CC charge at C/2 (1.5 A) from 0 to 100 % SOC, CC discharge at c, no rests. Per full cycle:
2 C0 total throughput, C0 charge throughput, (1 - SOC_Ref) C0 of it above SOC_Ref, cycle time 1/0.5 + 1/c h; the
calendar factor is its time average over the cycle (SOC uniform on [0, 1]). L = full cycles to Q_loss = 0.10.
"""
import math
import os

import numpy as np

from common import FLOORS, FLOOR_KEYS, RATES, RATE_NAMES, REPO, load_json, save

P = load_json(os.path.join(REPO, "inputs", "schimpe2018_parameters.json"))
TAB = P["table_IV_p_A190"]
APP = P["appendix_p_A192_A193"]
RG = P["other"]["R_g_J_per_mol_K"]["value"]
FARADAY = P["other"]["F_C_per_mol"]["value"]
C0 = P["other"]["C_0_Ah"]["value"]
SOC_REF = P["other"]["SOC_Ref"]["value"]


def params(text_values=False):
    p = {
        "kcal": TAB["k_Cal_Ref_h^-0.5"], "Ea_cal": TAB["E_a_Cal_J_per_mol"],
        "khi": TAB["k_Cyc_HighT_Ref_Ah^-0.5"], "Ea_hi": TAB["E_a_Cyc_HighT_J_per_mol"],
        "klo": TAB["k_Cyc_LowT_Ref_Ah^-0.5"], "Ea_lo": TAB["E_a_Cyc_LowT_J_per_mol"],
        "khs": TAB["k_Cyc_LowT_HighSOC_Ref_Ah^-1"], "Ea_hs": TAB["E_a_Cyc_LowT_HighSOC_J_per_mol"],
        "alpha": TAB["alpha"], "beta_lo": TAB["beta_LowT_h"], "beta_hs": TAB["beta_LowT_HighSOC_h"],
        "Tref": TAB["T_Ref_K"], "Iref": TAB["I_Ch_Ref_A"], "Uref": TAB["U_a_Ref_V"], "k0": TAB["k_0"],
    }
    if text_values:
        tv = P["text_values_for_comparison"]
        p.update(kcal=tv["k_Cal_Ref_h^-0.5"]["value"], Ea_cal=tv["E_a_Cal_J_per_mol"]["value"],
                 khi=tv["k_Cyc_HighT_Ref_Ah^-0.5"]["value"], klo=tv["k_Cyc_LowT_Ref_Ah^-0.5"]["value"],
                 Ea_hi=tv["E_a_Cyc_HighT_J_per_mol"]["value"], Ea_lo=tv["E_a_Cyc_LowT_J_per_mol"]["value"],
                 khs=tv["k_Cyc_LowT_HighSOC_Ref_Ah^-1"]["value"], Ea_hs=tv["E_a_Cyc_LowT_HighSOC_J_per_mol"]["value"],
                 beta_hs=tv["beta_LowT_HighSOC_h"]["value"])
    return p


def U_a(x):
    """Eq. A1 (graphite potential vs Li, Safari and Delacourt 2011), V."""
    x = np.asarray(x, dtype=float)
    return (0.6379 + 0.5416 * np.exp(-305.5309 * x) + 0.044 * np.tanh(-(x - 0.1958) / 0.1088)
            - 0.1978 * np.tanh((x - 1.0571) / 0.0854) - 0.6875 * np.tanh((x + 0.0117) / 0.0529)
            - 0.0175 * np.tanh((x - 0.5692) / 0.0875))


def x_a(soc):
    return APP["x_a_SOC0"] + soc * (APP["x_a_SOC100"] - APP["x_a_SOC0"])


# SOC average of the SOC part of Eq. 9, by composite Simpson on a fine uniform grid
_N = 400000
_SOC = np.linspace(0.0, 1.0, _N + 1)
_W = np.ones(_N + 1)
_W[1:-1:2] = 4.0
_W[2:-1:2] = 2.0
_W /= 3.0 * _N


def soc_factor_mean(p, rms=False):
    f = np.exp(p["alpha"] * FARADAY * (p["Uref"] - U_a(x_a(_SOC))) / (RG * p["Tref"])) + p["k0"]
    if rms:
        return math.sqrt(float(_W @ f ** 2))
    return float(_W @ f)


def arr(Ea, T, Tref, sign=-1.0):
    return math.exp(sign * Ea / RG * (1.0 / T - 1.0 / Tref))


def coefficients(T, c_dis, p, i_ch=1.5, rms=False):
    """Per-cycle coefficients of Q_loss(N) = a_cal sqrt(N) + a_hi sqrt(N) + a_lo sqrt(N) + b_hs N."""
    t_cyc = 1.0 / (i_ch / C0) + 1.0 / c_dis               # h per full cycle (charge at i_ch/C0 C-rate)
    kcal = p["kcal"] * arr(p["Ea_cal"], T, p["Tref"], -1) * soc_factor_mean(p, rms)
    khi = p["khi"] * arr(p["Ea_hi"], T, p["Tref"], -1)
    klo = p["klo"] * arr(p["Ea_lo"], T, p["Tref"], +1) * math.exp(p["beta_lo"] * (i_ch - p["Iref"]) / C0)
    khs = p["khs"] * arr(p["Ea_hs"], T, p["Tref"], +1) * math.exp(p["beta_hs"] * (i_ch - p["Iref"]) / C0)
    return {"cal": kcal * math.sqrt(t_cyc), "hi": khi * math.sqrt(2 * C0), "lo": klo * math.sqrt(C0),
            "hs": khs * (1.0 - SOC_REF) * C0}


def lifetime(T, c_dis, p, loss=0.10, rms=False, include_hs=True):
    co = coefficients(T, c_dis, p, rms=rms)
    a = co["cal"] + co["hi"] + co["lo"]
    b = co["hs"] if include_hs else 0.0
    if b == 0.0:
        x = loss / a
    else:
        x = (-a + math.sqrt(a * a + 4 * b * loss)) / (2 * b)     # x = sqrt(N)
    return x * x, co


def s_of(c_dis, p, T=308.15, d=0.25, **kw):
    return -(math.log(lifetime(T + d, c_dis, p, **kw)[0]) - math.log(lifetime(T - d, c_dis, p, **kw)[0])) / (2 * d)


def v1_loss(T, p, n_fec=2800.0, include_hs=True):
    """2800 full equivalent cycles at 1C CC (1 FEC = 6 Ah = 2 C0 throughput; cycle time 2 h; I_Ch = 3 A)."""
    t = n_fec * 2.0
    q_tot = n_fec * 2 * C0
    q_ch = n_fec * C0
    i_ch = 3.0
    kcal = p["kcal"] * arr(p["Ea_cal"], T, p["Tref"], -1) * soc_factor_mean(p)
    khi = p["khi"] * arr(p["Ea_hi"], T, p["Tref"], -1)
    klo = p["klo"] * arr(p["Ea_lo"], T, p["Tref"], +1) * math.exp(p["beta_lo"] * (i_ch - p["Iref"]) / C0)
    khs = p["khs"] * arr(p["Ea_hs"], T, p["Tref"], +1) * math.exp(p["beta_hs"] * (i_ch - p["Iref"]) / C0)
    parts = {"cal": kcal * math.sqrt(t), "high": khi * math.sqrt(q_tot), "low": klo * math.sqrt(q_ch),
             "high_soc": khs * (1 - SOC_REF) * q_ch}
    tot = sum(parts.values()) if include_hs else parts["cal"] + parts["high"] + parts["low"]
    return 100 * tot, {k: 100 * v for k, v in parts.items()}


def main():
    p = params()
    pt = params(text_values=True)
    T35 = 308.15
    out = {"soc_factor_mean": soc_factor_mean(p), "soc_factor_rms": soc_factor_mean(p, rms=True),
           "U_a_at_SOC0_V": float(U_a(x_a(0.0))), "U_a_at_SOC100_V": float(U_a(x_a(1.0))), "rates": {}}
    for c, name in zip(RATES, RATE_NAMES):
        L35, co = lifetime(T35, c, p)
        n = L35
        terms = {"cal": co["cal"] * math.sqrt(n), "high": co["hi"] * math.sqrt(n), "low": co["lo"] * math.sqrt(n),
                 "high_soc": co["hs"] * n}
        s_h = s_of(c, p, d=0.25)           # +-0.25 K (total width 0.5 K)
        s_q = s_of(c, p, d=0.5)            # +-0.5 K: primary, the reading that the published values follow
        s_fine = s_of(c, p, d=1e-3)        # essentially the exact derivative
        # term-only levers (each term alone sets L)
        s_terms = {}
        for key, label in (("cal", "calendar only"), ("hi", "high-temperature cycling only"),
                           ("lo", "low-temperature cycling only")):
            def lnL_only(T, key=key):
                cc = coefficients(T, c, p)
                return math.log((0.1 / cc[key]) ** 2)
            s_terms[label] = -(lnL_only(T35 + 0.5) - lnL_only(T35 - 0.5)) / 1.0

        def lnL_cyc(T):
            cc = coefficients(T, c, p)
            a = cc["hi"] + cc["lo"]
            b = cc["hs"]
            x = (-a + math.sqrt(a * a + 0.4 * b)) / (2 * b)
            return math.log(x * x)
        s_terms["cycling terms only"] = -(lnL_cyc(T35 + 0.5) - lnL_cyc(T35 - 0.5)) / 1.0
        lin = {k: FLOORS[k] / s_q for k in FLOOR_KEYS}
        exact = {k: 2 * math.atanh(FLOORS[k] / 2) / s_q for k in FLOOR_KEYS}
        sq3 = terms["cal"] + terms["high"] + terms["low"]
        out["rates"][name] = {
            "crate": c, "L_35C": L35, "s": s_q, "s_pm0.5K": s_q, "s_pm0.25K": s_h, "s_pm0.001K": s_fine,
            "dT_floor_K": lin, "dT_floor_exact_pair_K": exact,
            "share_of_loss_at_35C": {k: v / 0.1 for k, v in terms.items()},
            "share_over_three_sqrt_terms": {k: terms[k] / sq3 for k in ("cal", "high", "low")},
            "s_by_term": s_terms,
            "s_text_values": s_of(c, pt, d=0.5),
            "s_rms_calendar_variant": s_of(c, p, d=0.5, rms=True),
            "L_35C_rms_calendar_variant": lifetime(T35, c, p, rms=True)[0],
            "s_without_high_soc_term": s_of(c, p, d=0.5, include_hs=False),
        }
    # lifetime curve at a C/2 discharge, 15 to 55 C, with and without the fourth term
    out["curve_C2"] = [{"T_C": float(T), "L": lifetime(273.15 + T, 0.5, p)[0],
                        "L_without_high_soc": lifetime(273.15 + T, 0.5, p, include_hs=False)[0]}
                       for T in range(15, 56)]
    # V1
    v1 = {}
    for T_C in (25.0, 45.0):
        full, parts = v1_loss(273.15 + T_C, p, include_hs=True)
        three, _ = v1_loss(273.15 + T_C, p, include_hs=False)
        full_t, _ = v1_loss(273.15 + T_C, pt, include_hs=True)
        three_t, _ = v1_loss(273.15 + T_C, pt, include_hs=False)
        v1[str(T_C)] = {"loss_pct_full_eq23": full, "loss_pct_three_sqrt_terms": three, "parts_pct": parts,
                        "loss_pct_full_eq23_text_values": full_t, "loss_pct_three_terms_text_values": three_t}
    out["V1"] = v1
    save("law2.json", out)
    print("SOC factor mean", out["soc_factor_mean"], "rms", out["soc_factor_rms"])
    for name in RATE_NAMES:
        r = out["rates"][name]
        print(f"{name:5s} L35={r['L_35C']:.4f} s={r["s"]:.10f} (0.25K:{r["s_pm0.25K"]:.10f}, fine:{r['s_pm0.001K']:.10f})"
              f" dT8.7={r['dT_floor_K']['p95_8.7pct']:.6f} shares={ {k: round(v, 6) for k, v in r['share_of_loss_at_35C'].items()} }")
        print("       s_terms", {k: round(v, 8) for k, v in r["s_by_term"].items()}, "s_text", r["s_text_values"],
              "s_rms", r["s_rms_calendar_variant"], "s_no_hs", r["s_without_high_soc_term"])
    print("V1", v1)
    print("curve 15C, 25C:", out["curve_C2"][0], out["curve_C2"][10])


if __name__ == "__main__":
    main()
