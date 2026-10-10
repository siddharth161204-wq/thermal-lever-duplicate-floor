#!/usr/bin/env python3
"""
The thermal lever against the duplicate-cell floor: Step 1 (the lever) and Step 2(a) (the lumped module).

Step 1: for each ageing law, the lifetime L(T) to 90 % state of health as a function of steady cell temperature,
the sensitivity s = -d ln L / dT at 35 C, and the temperature difference dT_floor = f / s that produces a lifetime
difference equal to a floor f (the mean 2.99 %, the worst-of-8 shortfall 4.10 % and the 95th percentile 8.72 % of
doi.org/10.5281/zenodo.23205413, read from inputs/floor.json).

Step 2(a): an eight-node steady module network (two rows of four 21700 cells at 23.1 mm pitch): nickel strip and
air-gap conduction and grey radiation between neighbours, convection and radiation to 35 C air, I^2 R heat per cell;
solved by Newton iteration; reports the hottest-minus-coldest cell spread and where both cells are. A second,
secondary arrangement puts one row of eight on a cold plate at 35 C under the cell bases.

Every value comes from a named published source recorded in DECISIONS.md and inputs/SOURCES.md; the derived
resistance is labelled as derived. Nothing is fitted to the floor.

This file is the brief's thermal_lever.py, filled. The script as received is kept in
provenance/thermal_lever_as_received.py; every change is listed in DECISIONS.md (D1 to D16).

Usage:
  python3 src/thermal_lever.py --out results
"""
import argparse, csv, json, math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
INPUTS = os.path.join(ROOT, "inputs")

SIGMA = 5.670374419e-8
R_GAS = 8.314462618
T_AMB_C = 35.0
T_AMB = T_AMB_C + 273.15

# ---------------- the floor (inputs/floor.json, copied from the published results.json; D1) ----------------
with open(os.path.join(INPUTS, "floor.json")) as _f:
    FLOOR = json.load(_f)
FLOORS = dict(FLOOR["FLOORS"])          # fractions: mean_3.0pct, worst_of_8_4.1pct, p95_8.7pct

THRESHOLD_SOH = 0.90                    # lifetime defined at 90 % SOH, as in the floor
CRATES = [1 / 16, 1 / 10, 1 / 5, 1 / 2, 1.0, 2.0]   # the dataset's average rates, then what it would take
DATASET_CHARGE_CRATE = 0.5              # every protocol in the dataset charges CC-CV at C/2 (Geslin et al. 2025, Methods)


# =====================================================================================================
# Ageing law 1: Kucinskis et al. 2022 (NCA+LNO / Si-graphite 21700, 5 Ah), cycles to 90 % SoH, Fig. 7b (D3)
# =====================================================================================================
def _read_kucinskis():
    tab = {}
    with open(os.path.join(INPUTS, "kucinskis2022_fig7_cycles.csv")) as f:
        for r in csv.DictReader(f):
            if r["cycles_to_soh"] == "":
                continue
            key = (round(float(r["soh"]), 2), round(float(r["charge_c_rate"]), 1))
            tab.setdefault(key, {})[float(r["ambient_T_C"])] = float(r["cycles_to_soh"])
    return tab


KUC = _read_kucinskis()
LAW_1 = {"name": "Kucinskis et al. 2022, NCA+LNO / Si-graphite 21700 (5 Ah)",
         "source": "J. Power Sources 549 (2022) 232129, Fig. 7b, p. 8; doi:10.1016/j.jpowsour.2022.232129",
         "metric": "full cycles to 90 % SoH (charge CC-CV at the stated rate to 4.2 V, CV to 0.05C; discharge 0.5C CC to 2.5 V)",
         "temperatures_C": [15, 25, 35, 45, 55], "charge_rates_tabulated": [0.2, 0.4, 0.6]}


def ln_L_kucinskis(T_C, charge_crate=DATASET_CHARGE_CRATE, soh=THRESHOLD_SOH):
    """ln(cycles to soh): linear in T between the printed temperatures, linear in charge rate between printed rates."""
    rates = [0.2, 0.4, 0.6]

    def at_rate(c):
        row = KUC[(round(soh, 2), c)]
        ts = sorted(t for t in row if t >= 15.0)          # the high-temperature part of the table (D3)
        if T_C < ts[0] or T_C > ts[-1]:
            return float("nan")
        return float(np.interp(T_C, ts, [math.log(row[t]) for t in ts]))

    if charge_crate in rates:
        return at_rate(charge_crate)
    if not (rates[0] <= charge_crate <= rates[-1]):
        return float("nan")
    lo = max(c for c in rates if c <= charge_crate)
    hi = min(c for c in rates if c >= charge_crate)
    w = (charge_crate - lo) / (hi - lo)
    return (1 - w) * at_rate(lo) + w * at_rate(hi)


def lifetime_law1(T_K, discharge_crate, charge_crate=DATASET_CHARGE_CRATE):
    """Cycles to 90 % SoH. The paper tested one discharge rate (0.5C), so discharge_crate does not enter (D3)."""
    return math.exp(ln_L_kucinskis(T_K - 273.15, charge_crate))


# =====================================================================================================
# Ageing law 2: Schimpe et al. 2018 (LiFePO4/graphite, Sony US26650FTC1), the full model of Eq. 23 (D4)
# =====================================================================================================
with open(os.path.join(INPUTS, "schimpe2018_parameters.json")) as _f:
    SCH = json.load(_f)


def schimpe_params(which="table"):
    t = dict(SCH["table_IV_p_A190"])
    if which == "text":                                 # the more precise values printed in the text (D4)
        for k, v in SCH["text_values_for_comparison"].items():
            t[k] = v["value"]
    o = SCH["other"]
    a = SCH["appendix_p_A192_A193"]
    t.update({"C0": o["C_0_Ah"]["value"], "SOC_Ref": o["SOC_Ref"]["value"], "Rg": o["R_g_J_per_mol_K"]["value"],
              "F": o["F_C_per_mol"]["value"], "xa0": a["x_a_SOC0"], "xa100": a["x_a_SOC100"]})
    return t


def U_a(x):
    """Graphite open-circuit potential, Schimpe et al. 2018 Eq. A1 (from Safari and Delacourt 2011)."""
    return (0.6379 + 0.5416 * np.exp(-305.5309 * x) + 0.044 * np.tanh(-(x - 0.1958) / 0.1088)
            - 0.1978 * np.tanh((x - 1.0571) / 0.0854) - 0.6875 * np.tanh((x + 0.0117) / 0.0529)
            - 0.0175 * np.tanh((x - 0.5692) / 0.0875))


def schimpe_soc_factor_mean(p, soc_lo=0.0, soc_hi=1.0, n=4001):
    """Mean over a uniform SOC range of the calendar SOC factor (Eq. 9): exp[alpha F (U_ref - U_a)/(Rg T_ref)] + k0."""
    soc = np.linspace(soc_lo, soc_hi, n)
    xa = p["xa0"] + soc * (p["xa100"] - p["xa0"])
    g = np.exp(p["alpha"] * p["F"] * (p["U_a_Ref_V"] - U_a(xa)) / (p["Rg"] * p["T_Ref_K"])) + p["k_0"]
    return float(np.trapezoid(g, soc) / (soc_hi - soc_lo))


def schimpe_factors(T_K, p, i_ch_A):
    """Stress factors at temperature T for a charge current i_ch_A (Eqs. 9, 15, 18, 21)."""
    arr = lambda Ea, sign: math.exp(sign * Ea / p["Rg"] * (1.0 / T_K - 1.0 / p["T_Ref_K"]))
    k_cal_T = p["k_Cal_Ref_h^-0.5"] * arr(p["E_a_Cal_J_per_mol"], -1)                 # times the SOC factor
    k_high = p["k_Cyc_HighT_Ref_Ah^-0.5"] * arr(p["E_a_Cyc_HighT_J_per_mol"], -1)
    k_low = (p["k_Cyc_LowT_Ref_Ah^-0.5"] * arr(p["E_a_Cyc_LowT_J_per_mol"], +1)
             * math.exp(p["beta_LowT_h"] * (i_ch_A - p["I_Ch_Ref_A"]) / p["C0"]))
    k_hs = (p["k_Cyc_LowT_HighSOC_Ref_Ah^-1"] * arr(p["E_a_Cyc_LowT_HighSOC_J_per_mol"], +1)
            * math.exp(p["beta_LowT_HighSOC_h"] * (i_ch_A - p["I_Ch_Ref_A"]) / p["C0"]))
    return k_cal_T, k_high, k_low, k_hs


def schimpe_cycle_coeffs(T_K, discharge_crate, charge_crate=DATASET_CHARGE_CRATE, p=None, terms=("cal", "high", "low", "hs"),
                         soc_factor=None):
    """Loss after N full cycles = a sqrt(N) + b N for the dataset-style protocol (D4): CC charge at charge_crate from
    0 to 100 % SOC, CC discharge at discharge_crate back to 0 %, no rests. Per cycle: time 1/c_ch + 1/c_dis hours,
    total throughput 2 C0, charge throughput C0, of which (1 - SOC_Ref) C0 above SOC_Ref. The calendar factor is the
    time average over the cycle, in which SOC is uniform on [0, 1] in both ramps."""
    p = p or schimpe_params()
    C0 = p["C0"]
    k_cal_T, k_high, k_low, k_hs = schimpe_factors(T_K, p, charge_crate * C0)
    g = schimpe_soc_factor_mean(p) if soc_factor is None else soc_factor
    period_h = 1.0 / charge_crate + 1.0 / discharge_crate
    parts = {"cal": k_cal_T * g * math.sqrt(period_h), "high": k_high * math.sqrt(2.0 * C0), "low": k_low * math.sqrt(C0)}
    a = sum(v for k, v in parts.items() if k in terms)
    b = k_hs * (1.0 - p["SOC_Ref"]) * C0 if "hs" in terms else 0.0
    return a, b, parts, k_hs * (1.0 - p["SOC_Ref"]) * C0


def lifetime_law2(T_K, discharge_crate, charge_crate=DATASET_CHARGE_CRATE, p=None, terms=("cal", "high", "low", "hs")):
    """Full cycles (= equivalent full cycles) to a 10 % capacity loss: solve a u + b u^2 = 0.10 for u = sqrt(N)."""
    a, b, _, _ = schimpe_cycle_coeffs(T_K, discharge_crate, charge_crate, p, terms)
    loss = 1.0 - THRESHOLD_SOH
    u = loss / a if b <= 0 else (-a + math.sqrt(a * a + 4.0 * b * loss)) / (2.0 * b)
    return u * u


LAW_2 = {"name": "Schimpe et al. 2018, LiFePO4 / graphite 26650 (3 Ah), full model of Eq. 23",
         "source": "J. Electrochem. Soc. 165 (2018) A181-A193, Table IV p. A190, Eqs. 9, 15, 18, 21, 23, A1, A2; doi:10.1149/2.1181714jes",
         "metric": "equivalent full cycles to a 10 % capacity loss for CC charge at C/2 and CC discharge at the stated rate, 0 to 100 % SOC",
         "fitted_temperature_range_C": [0, 55]}


# =====================================================================================================
# The lever (generic for both laws)
# =====================================================================================================
def sensitivity(life_fn, crate, T_C=T_AMB_C, dT=0.5, **kw):
    """s = -d ln L / dT at T_C by a centred difference of +-dT kelvin (per kelvin)."""
    T = T_C + 273.15
    Lm, Lp = life_fn(T - dT, crate, **kw), life_fn(T + dT, crate, **kw)
    return -(math.log(Lp) - math.log(Lm)) / (2 * dT)


def dT_floor(f, s):
    """The temperature difference that gives a relative lifetime difference f: f / s (the brief's linear form).
    Negative where hotter cells live longer."""
    return f / s if s != 0 else float("inf")


def dT_floor_exact_pair(f, s):
    """Exact inversion for two cells whose ln L differs by s dT: 2 tanh(s dT / 2) = f."""
    return 2.0 * math.atanh(f / 2.0) / s if s != 0 else float("inf")


def lever(life_fn, crate, **kw):
    s = sensitivity(life_fn, crate, **kw)
    return {"crate": crate, "s_per_K": s,
            "dT_floor_K": {k: dT_floor(f, s) for k, f in FLOORS.items()},
            "dT_floor_exact_pair_K": {k: dT_floor_exact_pair(f, s) for k, f in FLOORS.items()},
            "L_at_35C": life_fn(T_AMB, crate, **kw)}


# =====================================================================================================
# Proxy cell and module (D5 to D12)
# =====================================================================================================
with open(os.path.join(INPUTS, "lgm50_resistance.json")) as _f:
    _RES = json.load(_f)
D_CELL, L_CELL = 0.0211, 0.0701          # m, LG M50 as in the two-cell study (doi:10.5281/zenodo.23253738)
CAP_AH = 5.0                             # Ah, rated capacity on 2.5 to 4.2 V (Chen et al. 2020, p. 9)
R_DC = None                              # set below from inputs/lgm50_resistance.json (D8)
PITCH = 0.0231                           # m, as in the two-cell study
GAP = PITCH - D_CELL                     # 2 mm
K_AIR = 0.03                             # W/(m K), as in the two-cell study
EPS = 0.8                                # can emissivity, as in the two-cell study
G_STRIP = 90.0 * 0.2e-3 * 8e-3 / PITCH   # W/K per nickel strip, the two-cell study's D24: 6.2 mW/K (two strips 12 mW/K)
H_LIST = [5.0, 10.0, 25.0]               # W/(m2 K)
H_BASE = 600.0                           # W/(m2 K), base-to-plate conductance of the published cell study (D12)
A_LAT = math.pi * D_CELL * L_CELL
A_END = math.pi * (D_CELL / 2) ** 2
A_OUTER = A_LAT + 2 * A_END
A_FACE = A_LAT * 0.25                    # facing-area fraction for the gap conduction, as the two-cell study used
X = 1.0 + GAP / D_CELL
F_VIEW = (math.sqrt(X * X - 1.0) + math.asin(1.0 / X) - X) / math.pi     # parallel equal cylinders, per whole lateral area
G_GAP = K_AIR * A_FACE / GAP

R_DC = float(_RES["registered"]["R_DC_ohm"])

# two rows of four: neighbour pairs (index: row*4 + column)
NEIGHBOURS = [(0, 1), (1, 2), (2, 3), (4, 5), (5, 6), (6, 7), (0, 4), (1, 5), (2, 6), (3, 7)]
POSITIONS = {i: ("corner" if i in (0, 3, 4, 7) else "inner") + f" (row {i // 4 + 1}, column {i % 4 + 1})" for i in range(8)}
# one row of eight on a cold plate (secondary arrangement, D12)
NEIGHBOURS_ROW = [(i, i + 1) for i in range(7)]
POSITIONS_ROW = {i: ("end" if i in (0, 7) else "inner") + f" (position {i + 1} of 8)" for i in range(8)}


def heat_per_cell(crate):
    return (crate * CAP_AH) ** 2 * R_DC


def _areas(n_nodes, pairs, base):
    n_nb = np.zeros(n_nodes)
    for a, b in pairs:
        n_nb[a] += 1
        n_nb[b] += 1
    base_area = A_END if base else 0.0
    a_conv = A_OUTER - base_area - n_nb * A_FACE                     # D6: no convection across a 2 mm gap
    a_rad = A_OUTER - base_area - n_nb * EPS * F_VIEW * A_LAT        # D5: first-order grey radiosity, room as black
    return a_conv, a_rad


def module_steady(q_cell, h, T_amb=T_AMB, g_strip=None, pairs=None, base=False, h_base=H_BASE, iters=100, tol=1e-11):
    """Newton solve of the steady node balance with conduction, radiation, convection and an optional cold plate.

    q_cell: heat per cell in W (scalar or one value per node). Returns T (K), the energy residual in W and the parts.
    """
    pairs = NEIGHBOURS if pairs is None else pairs
    g_strip = G_STRIP if g_strip is None else g_strip
    n = 1 + max(max(p) for p in pairs)
    q = np.full(n, float(q_cell)) if np.isscalar(q_cell) else np.asarray(q_cell, dtype=float)
    a_conv, a_rad = _areas(n, pairs, base)
    g_base = h_base * A_END if base else 0.0
    g_pair = g_strip + G_GAP
    r_pair = EPS ** 2 * SIGMA * F_VIEW * A_LAT                         # D5: exchange eps^2 sigma F A_lat (Ti^4 - Tj^4)
    T = np.full(n, T_amb + 1.0)
    for it in range(iters):
        F = q - h * a_conv * (T - T_amb) - EPS * SIGMA * a_rad * (T ** 4 - T_amb ** 4) - g_base * (T - T_amb)
        J = np.diag(-h * a_conv - 4 * EPS * SIGMA * a_rad * T ** 3 - g_base)
        for a, b in pairs:
            flow = g_pair * (T[a] - T[b]) + r_pair * (T[a] ** 4 - T[b] ** 4)
            F[a] -= flow
            F[b] += flow
            da = g_pair + 4 * r_pair * T[a] ** 3
            db = g_pair + 4 * r_pair * T[b] ** 3
            J[a, a] -= da; J[a, b] += db; J[b, a] += da; J[b, b] -= db
        dT = np.linalg.solve(J, -F)
        T = T + dT
        if np.abs(dT).max() < tol:
            break
    loss_conv = h * a_conv * (T - T_amb)
    loss_rad = EPS * SIGMA * a_rad * (T ** 4 - T_amb ** 4)
    loss_base = g_base * (T - T_amb)
    balance = float(q.sum() - (loss_conv + loss_rad + loss_base).sum())
    return T, balance, {"iterations": it + 1, "sources_W": float(q.sum()),
                        "losses_W": float((loss_conv + loss_rad + loss_base).sum())}


def module_row(crate, h, base=False):
    pairs = NEIGHBOURS_ROW if base else NEIGHBOURS
    pos = POSITIONS_ROW if base else POSITIONS
    q = heat_per_cell(crate)
    T, bal, info = module_steady(q, h, pairs=pairs, base=base)
    i_hot, i_cold = int(T.argmax()), int(T.argmin())
    return {"arrangement": "1x8 on cold plate" if base else "2x4 in air", "crate": crate, "q_cell_W": q, "h": h,
            "T_max_C": float(T.max() - 273.15), "T_min_C": float(T.min() - 273.15), "T_mean_C": float(T.mean() - 273.15),
            "dT_module_K": float(T.max() - T.min()), "dT_hot_minus_mean_K": float(T.max() - T.mean()),
            "i_hot": i_hot, "hot_position": pos[i_hot], "i_cold": i_cold, "cold_position": pos[i_cold],
            "balance_W": bal, "newton_iterations": info["iterations"], "T_cells_C": [float(x - 273.15) for x in T]}


def check_filled():
    missing = [k for k, v in (("R_DC", R_DC), ("CAP_AH", CAP_AH), ("G_STRIP", G_STRIP)) if v is None]
    if missing:
        raise SystemExit("FILL these from the named sources before running: " + ", ".join(missing))


LAWS = {"law_1_kucinskis2022": lifetime_law1, "law_2_schimpe2018": lifetime_law2}


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default="results"); a = ap.parse_args()
    check_filled(); os.makedirs(a.out, exist_ok=True)
    # Step 1
    levers = {name: [lever(fn, c) for c in CRATES] for name, fn in LAWS.items()}
    curves = {name: [{"T_C": t, "L": fn(t + 273.15, 0.5)} for t in np.arange(15.0, 55.01, 1.0)] for name, fn in LAWS.items()}
    # Step 2(a), registered arrangement, and the secondary cold-plate arrangement
    rows = [module_row(c, h) for c in CRATES for h in H_LIST]
    rows_plate = [module_row(c, h, base=True) for c in CRATES for h in H_LIST]
    keys = [k for k in rows[0] if k != "T_cells_C"] + [f"T_cell{i}_C" for i in range(8)]
    for name, rr in (("module_spread.csv", rows), ("module_spread_coldplate.csv", rows_plate)):
        with open(os.path.join(a.out, name), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
            for r in rr:
                d = {k: v for k, v in r.items() if k != "T_cells_C"}
                d.update({f"T_cell{i}_C": t for i, t in enumerate(r["T_cells_C"])})
                w.writerow(d)
    meta = {"floors": FLOORS, "floor_source": FLOOR["source"],
            "laws": {"law_1_kucinskis2022": LAW_1, "law_2_schimpe2018": LAW_2},
            "cell": {"D_m": D_CELL, "L_m": L_CELL, "capacity_Ah": CAP_AH, "R_DC_ohm": R_DC, "R_DC_definition": _RES["registered"]},
            "module": {"pitch_m": PITCH, "gap_m": GAP, "k_air": K_AIR, "eps": EPS, "G_strip_W_per_K": G_STRIP, "G_gap_W_per_K": G_GAP,
                       "F_view": F_VIEW, "A_lat_m2": A_LAT, "A_end_m2": A_END, "A_outer_m2": A_OUTER, "A_face_m2": A_FACE,
                       "h_list": H_LIST, "h_base": H_BASE, "T_amb_C": T_AMB_C}}
    with open(os.path.join(a.out, "lever.json"), "w") as f:
        json.dump({**meta, "levers": levers, "lifetime_curves_at_C2": curves, "module": rows, "module_coldplate": rows_plate}, f, indent=1)
    print("Step 1 (lever) and Step 2(a) (lumped module) written to", a.out)
