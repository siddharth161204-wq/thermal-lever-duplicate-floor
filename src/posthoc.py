#!/usr/bin/env python3
"""Post-hoc analyses, added after the registered results and NOT pre-registered (DECISIONS.md D24). Each one answers a
question the independent checks raised about the registered results; none changes a registered number.

1. The dataset's own discharge profiles. The brief's heat is the I^2 R heat of a constant current at the average rate
   (Step 2(a)), and D19 added the C/2 charge. The dataset's dynamic discharges carry more heat than a constant current
   of the same mean, by mean(I^2) / mean(I)^2. The authors publish, for every cell, the mean discharge current c and
   its "normalized current variance" NCV (data/protocol_features.csv of
   github.com/geslina/dynamic_cycling_Nature_Energy_2024 at commit 5b2f7f0, with the cell and protocol names in
   data/metadata.pkl). NCV is the variance over the mean, in C-rate units: the constant-current protocols with a
   storage rest are two-level profiles, for which c + Var/mean is the plateau current and Var/(mean * plateau) the rest
   fraction, and both identities hold in the data to round-off (checked below); for constant current without rest NCV
   is zero. So, per unit of discharged capacity, mean(I^2) = c^2 + c NCV over the discharge window (rests included),
   which lasts 1/c hours.
   The charge is the dataset's CC-CV charge at C/2 to 4.2 V with a C/20 cut-off, after a discharge to 3.1 V. Its
   duration t_ch and its integral of (I / I_1C)^2, e_ch, per unit of discharged capacity, are taken from PyBaMM's
   Doyle-Fuller-Newman model of the proxy cell (ORegan2022, isothermal 35 C), the model behind the registered
   resistance (D8). The cycle-averaged heat is then R I_1C^2 (c + NCV + e_ch) / (1/c + t_ch). D19's form, a 2 h
   constant-current charge (e_ch = 0.5, t_ch = 2), is kept alongside.
   The module (lumped, registered solver) is solved at that heat for every one of the 92 cells; the resolved cell's
   spread is scaled from its C/2 value (it is linear in the heat to 0.02 % between C/2 and 2C). The implied divergence
   uses the registered lever at 35 C; Law 2 evaluated at each cell's own temperature (D18) is reported alongside, as
   is the steady spread at the heat of the discharge window alone, which the cells approach while it lasts.
2. The heat at C/2 at which the module reaches Law 2's mean floor, as a multiple of the brief's C/2 heat.
3. The thermal time constants of a module cell, to say what the steady spread at the cycle-averaged heat means while
   the heat changes between the charge and the discharge.
4. V1 for Law 2 with all four terms of Eq. 23 (the registered V1 left out the high-SOC term).
5. The air conductivity, inherited from the two-cell study without a primary table: the C/2 module spread with it
   10 % lower and 10 % higher.

Usage: python3 src/posthoc.py --upstream PATH/TO/dynamic_cycling_Nature_Energy_2024 --results results
"""
import argparse, csv, hashlib, importlib, json, math, os, pickle, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import thermal_lever as TL

UPSTREAM_COMMIT = "5b2f7f04d05072fe1f9bd8af664f23eadbde6317"
EXPECTED_SHA256 = {
    "metadata.pkl": "7579f48cf3767552154fba3ad9acb048c3166088c2a46fadb06b2a221b3a1e1f",          # as in the floor study
    "protocol_features.csv": "ed6da2af4d8dc3267217573014a18a90d0d578d92a2a12c1272f2b33ea5c85da",
}
# The exact globals the released metadata pickle uses (the floor study's whitelist): data constructors only.
ALLOWED_GLOBALS = {
    ("builtins", "slice"), ("numpy", "dtype"), ("numpy", "ndarray"),
    ("numpy.core.multiarray", "_reconstruct"), ("numpy.core.multiarray", "scalar"), ("numpy.core.numeric", "_frombuffer"),
    ("pandas._libs.internals", "_unpickle_block"), ("pandas.core.frame", "DataFrame"),
    ("pandas.core.indexes.base", "Index"), ("pandas.core.indexes.base", "_new_Index"),
    ("pandas.core.indexes.multi", "MultiIndex"), ("pandas.core.indexes.range", "RangeIndex"),
    ("pandas.core.internals.managers", "BlockManager"),
    ("numpy._core.multiarray", "_reconstruct"), ("numpy._core.multiarray", "scalar"), ("numpy._core.numeric", "_frombuffer"),
    ("pandas", "DataFrame"), ("pandas", "Index"), ("pandas", "MultiIndex"), ("pandas", "StringDtype"),
    ("pandas.arrays", "StringArray"), ("pandas._libs.arrays", "__pyx_unpickle_NDArrayBacked"),
}
NOMINAL = {"Co16": 1 / 16, "Co10": 1 / 10, "Co5": 1 / 5, "Co2": 1 / 2}
LABEL = {"Co16": "C/16", "Co10": "C/10", "Co5": "C/5", "Co2": "C/2"}
L1, L2 = "law_1_kucinskis2022", "law_2_schimpe2018"


class _Whitelist(pickle.Unpickler):
    def find_class(self, module, name):
        if (module, name) not in ALLOWED_GLOBALS:
            raise pickle.UnpicklingError(f"refused global {module}.{name} in a data pickle")
        if module.startswith("numpy.core."):
            try:
                alt = "numpy._core." + module[len("numpy.core."):]
                importlib.import_module(alt)
                module = alt
            except ImportError:
                pass
        return super().find_class(module, name)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def load_dataset(upstream):
    d = os.path.join(upstream, "data")
    hashes = {name: sha256(os.path.join(d, name)) for name in EXPECTED_SHA256}
    for name, exp in EXPECTED_SHA256.items():
        if hashes[name] != exp:
            raise SystemExit(f"{name}: SHA-256 {hashes[name]} does not match the expected {exp}")
    with open(os.path.join(d, "metadata.pkl"), "rb") as f:
        meta = _Whitelist(f).load()
    with open(os.path.join(d, "protocol_features.csv")) as f:
        feats = list(csv.DictReader(f))
    if len(feats) != len(meta):
        raise SystemExit("protocol_features.csv and metadata.pkl differ in length")
    cells = []
    for i, ft in enumerate(feats):
        m = meta.iloc[i]
        c = float(ft["Average Current"])
        if abs(c - float(m["avg_crate_exp"])) > 1e-9:                       # the two files list the cells in one order
            raise SystemExit(f"row {i}: average current {c} differs from metadata {m['avg_crate_exp']}")
        cells.append({"cell": str(m["cell_name"]), "protocol": str(m["protocol_name"]), "type": str(m["protocol_type"]),
                      "nominal": str(m["avg_Crate"]), "c": c, "ncv": float(ft["Normalized Current Variance"]),
                      "i_max": float(ft["Maximum Discharge Current"]),
                      "rest_high_soc": float(ft["Rest Fraction at High SOC"])})
    return cells, hashes


def check_normalisation(cells):
    """Two-level constant-current protocols with a storage rest: plateau = c + Var/mean, which must equal the measured
    current of the same nominal rate without storage; with the rest at full charge, rest fraction = NCV / plateau."""
    plateau = {}
    for x in cells:
        if x["type"] == "CC" and "no_storage" in x["protocol"]:
            plateau.setdefault(x["nominal"], []).append(x["c"])
    plateau = {k: float(np.mean(v)) for k, v in plateau.items()}
    rows, dev_plateau, dev_rest = [], [], []
    for x in cells:
        if x["type"] == "CC" and "no_storage" not in x["protocol"]:
            p = x["c"] + x["ncv"]
            d1 = p - plateau[x["nominal"]]
            dev_plateau.append(abs(d1))
            row = {"cell": x["cell"], "protocol": x["protocol"], "c_plus_ncv": p, "plateau_without_storage": plateau[x["nominal"]],
                   "difference": d1}
            if "100%" in x["protocol"]:
                d2 = x["rest_high_soc"] - x["ncv"] / p
                dev_rest.append(abs(d2))
                row["rest_fraction_minus_ncv_over_plateau"] = d2
            rows.append(row)
    no_rest = [x["ncv"] for x in cells if x["type"] == "CC" and "no_storage" in x["protocol"]]
    return {"plateau_currents_without_storage": plateau, "max_abs_plateau_difference": max(dev_plateau),
            "max_abs_rest_fraction_difference": max(dev_rest), "max_ncv_without_rest": max(no_rest), "rows": rows}


def charge_phase(rates=(1 / 16, 1 / 10, 1 / 5, 1 / 2)):
    """The dataset's CC-CV charge (C/2 to 4.2 V, C/20 cut-off) after a discharge to 3.1 V, for the proxy cell in PyBaMM's
    DFN (ORegan2022, isothermal 35 C): per unit of discharged capacity, the charge time and the integral of (I/I_1C)^2."""
    import pybamm
    out = {}
    for c in rates:
        model = pybamm.lithium_ion.DFN()
        p = pybamm.ParameterValues("ORegan2022")
        p.update({"Ambient temperature [K]": TL.T_AMB, "Initial temperature [K]": TL.T_AMB})
        exp = pybamm.Experiment(["Charge at C/2 until 4.2 V", "Hold at 4.2 V until C/20", "Rest for 1 minute",
                                 f"Discharge at {c:.6f}C until 3.1 V", "Charge at C/2 until 4.2 V", "Hold at 4.2 V until C/20"])
        sol = pybamm.Simulation(model, parameter_values=p, experiment=exp).solve(initial_soc=0.5)
        steps = [cy.steps[0] for cy in sol.cycles]
        res = {}
        for k, name in ((3, "discharge"), (4, "cc"), (5, "cv")):
            t = steps[k]["Time [s]"].entries
            i = steps[k]["Current [A]"].entries
            res[name] = {"t_h": float((t[-1] - t[0]) / 3600.0), "Q_Ah": float(abs(np.trapezoid(i, t)) / 3600.0),
                         "int_I2_A2h": float(np.trapezoid(i ** 2, t) / 3600.0)}
        q_unit = res["discharge"]["Q_Ah"] / TL.CAP_AH
        res["per_unit_discharged_capacity"] = {
            "discharged_fraction_of_rated": q_unit,
            "t_ch_h": (res["cc"]["t_h"] + res["cv"]["t_h"]) / q_unit,
            "e_ch_Crate2_h": (res["cc"]["int_I2_A2h"] + res["cv"]["int_I2_A2h"]) / TL.CAP_AH ** 2 / q_unit,
            "cv_share_of_charge": res["cv"]["Q_Ah"] / (res["cc"]["Q_Ah"] + res["cv"]["Q_Ah"])}
        out[f"{c:.6f}"] = res
    return out


def law2_exact_pair(T_cells_K, c):
    """Law 2 at each cell's own temperature (D18), inside its fitted range (0 to 55 C): (L_cold - L_hot) / mean."""
    if not all(273.15 - 1e-9 <= t <= 328.15 + 1e-9 for t in T_cells_K):
        return float("nan")
    L = [TL.lifetime_law2(t, c) for t in T_cells_K]
    i_hot, i_cold = int(np.argmax(T_cells_K)), int(np.argmin(T_cells_K))
    return float((L[i_cold] - L[i_hot]) / (0.5 * (L[i_cold] + L[i_hot])))


def protocol_heat(cells, res_dir, charge):
    lev = json.load(open(os.path.join(res_dir, "lever.json")))
    cel = json.load(open(os.path.join(res_dir, "resolved_cell.json")))
    s = {law: {round(e["crate"], 6): e["s_per_K"] for e in lev["levers"][law]} for law in TL.LAWS}
    dtf = {law: {round(e["crate"], 6): e["dT_floor_K"] for e in lev["levers"][law]} for law in TL.LAWS}
    cell_c2 = {r["h"]: r for r in cel["rows"] if r["base_h"] is None and abs(r["crate"] - 0.5) < 1e-9}
    q_c2 = TL.heat_per_cell(0.5)
    k_heat = TL.R_DC * TL.CAP_AH ** 2
    k4 = "worst_of_8_4.1pct"
    rows = []
    for x in cells:
        c, v = x["c"], x["ncv"]
        cn = NOMINAL[x["nominal"]]
        ch = charge[f"{cn:.6f}"]["per_unit_discharged_capacity"]
        q = k_heat * (c + v + ch["e_ch_Crate2_h"]) / (1.0 / c + ch["t_ch_h"])          # with the CC-CV charge
        q_cc = k_heat * (c + v + 0.5) / (1.0 / c + 2.0)                                 # D19's 2 h constant-current charge
        q_dis = k_heat * (c * c + c * v)                                                # the discharge window alone
        row = {**{k: x[k] for k in ("cell", "protocol", "type")}, "nominal_rate": LABEL[x["nominal"]], "c_mean": c,
               "normalized_current_variance": v, "discharge_heat_factor": 1.0 + v / c, "q_cycle_avg_W": q,
               "q_cycle_avg_cc_charge_W": q_cc, "q_discharge_window_W": q_dis,
               "ratio_to_brief_heat_at_nominal_rate": q / TL.heat_per_cell(cn)}
        for h in TL.H_LIST:
            T, bal, _ = TL.module_steady(q, h)
            dm, dhm = float(T.max() - T.min()), float(T.max() - T.mean())
            row[f"dT_module_K_h{h:g}"] = dm
            row[f"dT_module_hot_minus_mean_K_h{h:g}"] = dhm
            row[f"dT_cell_K_h{h:g}"] = cell_c2[h]["dT_cell_K"] * q / q_c2
            row[f"T_mean_C_h{h:g}"] = float(T.mean() - 273.15)
            row[f"T_max_C_h{h:g}"] = float(T.max() - 273.15)
            Tcc, _, _ = TL.module_steady(q_cc, h)
            row[f"dT_module_cc_charge_K_h{h:g}"] = float(Tcc.max() - Tcc.min())
            Td, _, _ = TL.module_steady(q_dis, h)
            row[f"dT_module_discharge_window_K_h{h:g}"] = float(Td.max() - Td.min())
            for law in TL.LAWS:
                sv = abs(s[law][round(cn, 6)])
                row[f"implied_module_{law}_h{h:g}"] = sv * dm
                row[f"implied_module_cc_charge_{law}_h{h:g}"] = sv * row[f"dT_module_cc_charge_K_h{h:g}"]
                row[f"implied_cell_{law}_h{h:g}"] = sv * row[f"dT_cell_K_h{h:g}"]
                for kf in lev["floors"]:
                    row[f"module_exceeds_{kf}_{law}_h{h:g}"] = dm > abs(dtf[law][round(cn, 6)][kf])
                    row[f"module_cc_charge_exceeds_{kf}_{law}_h{h:g}"] = row[f"dT_module_cc_charge_K_h{h:g}"] > abs(dtf[law][round(cn, 6)][kf])
                row[f"module_hot_minus_mean_exceeds_{k4}_{law}_h{h:g}"] = dhm > abs(dtf[law][round(cn, 6)][k4])
            row[f"law2_exact_pair_h{h:g}"] = law2_exact_pair(T, cn)
        rows.append(row)
    return rows, lev["floors"]


def summarise(rows, floors):
    k4 = "worst_of_8_4.1pct"
    out = {}
    for nom in ("C/16", "C/10", "C/5", "C/2"):
        g = [r for r in rows if r["nominal_rate"] == nom]
        dyn = [r for r in g if r["type"] != "CC"]
        best = max(g, key=lambda r: r["q_cycle_avg_W"])
        d = {"cells": len(g), "dynamic_protocols": sorted({r["protocol"] for r in dyn}),
             "max_heat_ratio_to_brief": max(r["ratio_to_brief_heat_at_nominal_rate"] for r in g),
             "discharge_heat_factor_dynamic_range": [min(r["discharge_heat_factor"] for r in dyn), max(r["discharge_heat_factor"] for r in dyn)] if dyn else None,
             "hottest_protocol": best["protocol"], "hottest_protocol_c_mean": [r["c_mean"] for r in g if r["protocol"] == best["protocol"]],
             "max_q_cycle_avg_W": best["q_cycle_avg_W"]}
        for h in TL.H_LIST:
            d[f"max_dT_module_K_h{h:g}"] = max(r[f"dT_module_K_h{h:g}"] for r in g)
            d[f"max_dT_module_discharge_window_K_h{h:g}"] = max(r[f"dT_module_discharge_window_K_h{h:g}"] for r in g)
            d[f"max_dT_cell_K_h{h:g}"] = max(r[f"dT_cell_K_h{h:g}"] for r in g)
            d[f"T_mean_C_range_h{h:g}"] = [min(r[f"T_mean_C_h{h:g}"] for r in g), max(r[f"T_mean_C_h{h:g}"] for r in g)]
            for law in TL.LAWS:
                d[f"max_implied_module_{law}_h{h:g}"] = max(r[f"implied_module_{law}_h{h:g}"] for r in g)
                d[f"max_implied_module_cc_charge_{law}_h{h:g}"] = max(r[f"implied_module_cc_charge_{law}_h{h:g}"] for r in g)
                for kf in floors:
                    d[f"protocols_exceeding_{kf}_{law}_h{h:g}"] = sorted({r["protocol"] for r in g if r[f"module_exceeds_{kf}_{law}_h{h:g}"]})
                    d[f"protocols_exceeding_cc_charge_{kf}_{law}_h{h:g}"] = sorted({r["protocol"] for r in g if r[f"module_cc_charge_exceeds_{kf}_{law}_h{h:g}"]})
                d[f"protocols_exceeding_{k4}_hottest_minus_mean_{law}_h{h:g}"] = sorted(
                    {r["protocol"] for r in g if r[f"module_hot_minus_mean_exceeds_{k4}_{law}_h{h:g}"]})
            if dyn:
                d[f"implied_module_{L2}_dynamic_range_h{h:g}"] = [min(r[f"implied_module_{L2}_h{h:g}"] for r in dyn),
                                                                  max(r[f"implied_module_{L2}_h{h:g}"] for r in dyn)]
                ex = [r[f"law2_exact_pair_h{h:g}"] for r in dyn if math.isfinite(r[f"law2_exact_pair_h{h:g}"])]
                d[f"law2_exact_pair_dynamic_range_h{h:g}"] = [min(ex), max(ex)] if ex else None
        out[nom] = d
    return out


def heat_to_reach_mean_floor(lev):
    """Multiple of the brief's C/2 heat at which the module spread reaches Law 2's Delta T_floor(3.0 %), per h."""
    s2 = [e for e in lev["levers"][L2] if abs(e["crate"] - 0.5) < 1e-9][0]
    target = abs(s2["dT_floor_K"]["mean_3.0pct"])
    q0 = TL.heat_per_cell(0.5)
    out = {}
    for h in TL.H_LIST:
        lo, hi = 0.5, 5.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            T, _, _ = TL.module_steady(q0 * mid, h)
            if float(T.max() - T.min()) < target:
                lo = mid
            else:
                hi = mid
        out[f"h{h:g}"] = 0.5 * (lo + hi)
    return out


def time_constants():
    """Heat capacity of one cell at 35 C (jelly roll from inputs/oregan2022_properties.txt, interpolated between 25 and
    40 C; the two-cell study's steel can of about 9 g at 500 J/(kg K), values without a primary table) and the module's
    slowest time constant, C / G, with G the secant loss conductance per cell at C/2 for each h."""
    rows = []
    with open(os.path.join(TL.INPUTS, "oregan2022_properties.txt")) as f:
        for line in f:
            p = line.split()
            if len(p) == 6 and p[0].lstrip("-").isdigit():
                rows.append([float(v) for v in p])
    tab = {r[0]: r for r in rows}
    cp35 = tab[25.0][2] + (tab[40.0][2] - tab[25.0][2]) * (35.0 - 25.0) / 15.0
    m_roll = tab[25.0][5] * 1e-3
    c_cell = m_roll * cp35 + 0.009 * 500.0
    q = TL.heat_per_cell(0.5)
    out = {"cp_roll_35C_J_per_kgK": cp35, "roll_mass_kg": m_roll, "can_J_per_K": 0.009 * 500.0, "C_cell_J_per_K": c_cell}
    for h in TL.H_LIST:
        T, _, _ = TL.module_steady(q, h)
        g = q / float(T.mean() - TL.T_AMB)
        out[f"G_loss_secant_W_per_K_h{h:g}"] = g
        out[f"tau_h_h{h:g}"] = c_cell / g / 3600.0
    return out


def v1_law2_full():
    p = TL.schimpe_params("table")
    out = {}
    for T_C, stated in ((25.0, "less than 9%"), (45.0, "12%")):
        a, b, parts, _ = TL.schimpe_cycle_coeffs(T_C + 273.15, 1.0, 1.0, p, terms=("cal", "high", "low", "hs"))
        N = 2800.0
        out[f"{T_C:.0f}C"] = {"stated": stated, "three_terms_pct": 100 * a * math.sqrt(N),
                             "high_soc_term_pct": 100 * b * N, "four_terms_pct": 100 * (a * math.sqrt(N) + b * N)}
    return out


def k_air_sensitivity():
    q = TL.heat_per_cell(0.5)
    g0 = TL.G_GAP
    out = {}
    for fac in (0.9, 1.0, 1.1):
        TL.G_GAP = g0 * fac
        T, _, _ = TL.module_steady(q, 10.0)
        out[f"k_air_x{fac:.1f}"] = float(T.max() - T.min())
    TL.G_GAP = g0
    return out


def main(upstream, res_dir):
    cells, hashes = load_dataset(upstream)
    charge = charge_phase()
    rows, floors = protocol_heat(cells, res_dir, charge)
    with open(os.path.join(res_dir, "dataset_protocol_heat.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    lev = json.load(open(os.path.join(res_dir, "lever.json")))
    out = {"note": "Post hoc, not pre-registered (DECISIONS.md D24).",
           "upstream": {"repository": "https://github.com/geslina/dynamic_cycling_Nature_Energy_2024", "commit": UPSTREAM_COMMIT,
                        "sha256": hashes, "dataset": "Geslin et al., Nature Energy 10, 172-180 (2025); dataset doi.org/10.25740/td676xr4322 (CC BY 4.0)"},
           "normalisation_check": check_normalisation(cells),
           "charge_phase_proxy_DFN": charge,
           "dataset_protocol_heat_by_rate": summarise(rows, floors),
           "heat_multiple_of_brief_C2_to_reach_law2_mean_floor": heat_to_reach_mean_floor(lev),
           "time_constants": time_constants(),
           "v1_law2_four_terms": v1_law2_full(),
           "k_air_sensitivity_C2_h10_K": k_air_sensitivity()}
    with open(os.path.join(res_dir, "posthoc.json"), "w") as f:
        json.dump(out, f, indent=1)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--upstream", required=True)
    ap.add_argument("--results", default="results")
    a = ap.parse_args()
    o = main(a.upstream, a.results)
    print(json.dumps({k: v for k, v in o.items() if k not in ("normalisation_check", "charge_phase_proxy_DFN")}, indent=1)[:4000])
