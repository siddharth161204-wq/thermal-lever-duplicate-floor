#!/usr/bin/env python3
"""Post-hoc analyses, added after the registered results and NOT pre-registered (DECISIONS.md D24). Each one answers a
question the independent checks raised about the registered results; none changes a registered number.

1. The dataset's own discharge profiles. The brief's heat is the I^2 R heat of a constant current at the average rate
   (Step 2(a)), and D19 added the C/2 charge. The dataset's dynamic discharges carry more heat than a constant current
   of the same mean, by mean(I^2) / mean(I)^2. The authors publish, for every cell, the mean discharge current and its
   "normalized current variance" (data/protocol_features.csv of github.com/geslina/dynamic_cycling_Nature_Energy_2024
   at commit 5b2f7f0, with the cell and protocol names in data/metadata.pkl). The normalisation is variance over mean,
   in C-rate units: for the constant-current protocols with a storage rest, a two-level profile, Var/mean must equal
   I_max - mean, and it does (checked below); for constant current without rest it is zero. So, per unit capacity,
   mean(I^2) = c^2 + c * NCV over the discharge window (rests included), which lasts 1/c hours, and the cycle-averaged
   heat with the C/2 charge of D19 is R I_1C^2 (c + NCV + 0.5) / (1/c + 2). The module (lumped, registered solver) is
   solved at that heat for every cell's protocol; the resolved cell's spread is scaled from its C/2 value (it is
   linear in the heat to 0.02 % between C/2 and 2C). The CV part of the charge is ignored, as in D19.
2. The thermal time constants of a module cell, to say what the steady spread at the cycle-averaged heat means while
   the heat changes between the C/2 charge and the discharge.
3. V1 for Law 2 with all four terms of Eq. 23 (the registered V1 left out the high-SOC term).
4. The air conductivity, inherited from the two-cell study without a primary table: the C/2 module spread with it
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
        if exp is not None and hashes[name] != exp:
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
                      "regen_fraction": float(ft["Relative Charge Fraction"])})
    return cells, hashes


def check_normalisation(cells):
    """For the constant-current protocols with a storage rest (two current levels, I_max and 0), Var/mean = I_max - mean."""
    out = []
    for x in cells:
        if x["type"] == "CC" and "no_storage" not in x["protocol"]:
            out.append({"cell": x["cell"], "protocol": x["protocol"], "ncv": x["ncv"], "i_max_minus_mean": x["i_max"] - x["c"],
                        "ratio": x["ncv"] / (x["i_max"] - x["c"])})
        elif x["type"] == "CC":
            out.append({"cell": x["cell"], "protocol": x["protocol"], "ncv": x["ncv"], "i_max_minus_mean": 0.0, "ratio": None})
    return out


def protocol_heat(cells, res_dir):
    lev = json.load(open(os.path.join(res_dir, "lever.json")))
    cel = json.load(open(os.path.join(res_dir, "resolved_cell.json")))
    s = {law: {round(e["crate"], 6): e["s_per_K"] for e in lev["levers"][law]} for law in TL.LAWS}
    dtf = {law: {round(e["crate"], 6): e["dT_floor_K"] for e in lev["levers"][law]} for law in TL.LAWS}
    cell_c2 = {r["h"]: r for r in cel["rows"] if r["base_h"] is None and abs(r["crate"] - 0.5) < 1e-9}
    q_c2 = TL.heat_per_cell(0.5)
    i1 = TL.CAP_AH
    rows = []
    for x in cells:
        c, v = x["c"], x["ncv"]
        k_dis = 1.0 + v / c                                         # mean(I^2) / mean(I)^2 over the discharge window
        q = TL.R_DC * i1 ** 2 * (c + v + 0.5) / (1.0 / c + 2.0)       # cycle average with the C/2 charge (D19 form)
        q_cc = TL.R_DC * i1 ** 2 * (c + 0.5) / (1.0 / c + 2.0)        # the same with a constant-current discharge (D19)
        cn = NOMINAL[x["nominal"]]
        row = {**{k: x[k] for k in ("cell", "protocol", "type")}, "nominal_rate": LABEL[x["nominal"]], "c_mean": c,
               "normalized_current_variance": v, "discharge_heat_factor": k_dis, "q_cycle_avg_W": q,
               "ratio_to_D19_constant_current": q / q_cc, "ratio_to_brief_heat_at_nominal_rate": q / TL.heat_per_cell(cn)}
        for h in TL.H_LIST:
            T, bal, _ = TL.module_steady(q, h)
            dm = float(T.max() - T.min())
            dhm = float(T.max() - T.mean())
            dc = cell_c2[h]["dT_cell_K"] * q / q_c2
            row[f"dT_module_K_h{h:g}"] = dm
            row[f"dT_module_hot_minus_mean_K_h{h:g}"] = dhm
            row[f"dT_cell_K_h{h:g}"] = dc
            row[f"T_max_C_h{h:g}"] = float(T.max() - 273.15)
            for law in TL.LAWS:
                sv = abs(s[law][round(cn, 6)])
                row[f"implied_module_{law}_h{h:g}"] = sv * dm
                row[f"implied_cell_{law}_h{h:g}"] = sv * dc
                for kf in lev["floors"]:
                    row[f"module_exceeds_{kf}_{law}_h{h:g}"] = dm > abs(dtf[law][round(cn, 6)][kf])
        rows.append(row)
    return rows, lev["floors"]


def summarise(rows, floors):
    out = {}
    for nom in ("C/16", "C/10", "C/5", "C/2"):
        g = [r for r in rows if r["nominal_rate"] == nom]
        best = max(g, key=lambda r: r["q_cycle_avg_W"])
        d = {"cells": len(g), "max_heat_ratio_to_brief": max(r["ratio_to_brief_heat_at_nominal_rate"] for r in g),
             "max_heat_ratio_to_D19_constant_current": max(r["ratio_to_D19_constant_current"] for r in g),
             "max_discharge_heat_factor": max(r["discharge_heat_factor"] for r in g),
             "hottest_protocol": best["protocol"], "max_q_cycle_avg_W": best["q_cycle_avg_W"]}
        for h in TL.H_LIST:
            d[f"max_dT_module_K_h{h:g}"] = max(r[f"dT_module_K_h{h:g}"] for r in g)
            d[f"max_dT_cell_K_h{h:g}"] = max(r[f"dT_cell_K_h{h:g}"] for r in g)
            for law in TL.LAWS:
                d[f"max_implied_module_{law}_h{h:g}"] = max(r[f"implied_module_{law}_h{h:g}"] for r in g)
                for kf in floors:
                    exc = sorted({r["protocol"] for r in g if r[f"module_exceeds_{kf}_{law}_h{h:g}"]})
                    d[f"protocols_exceeding_{kf}_{law}_h{h:g}"] = exc
        out[nom] = d
    return out


def time_constants():
    """Heat capacity of one cell at 35 C (jelly roll from inputs/oregan2022_properties.txt, interpolated between 25 and
    40 C; the 9 g steel can at 500 J/(kg K) of the two-cell study) and the module's slowest time constant, C / G, with
    G the secant loss conductance per cell at C/2 for each h."""
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
    rows, floors = protocol_heat(cells, res_dir)
    with open(os.path.join(res_dir, "dataset_protocol_heat.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    out = {"note": "Post hoc, not pre-registered (DECISIONS.md D24).",
           "upstream": {"repository": "https://github.com/geslina/dynamic_cycling_Nature_Energy_2024", "commit": UPSTREAM_COMMIT,
                        "sha256": hashes, "dataset_licence": "CC BY 4.0 (doi.org/10.25740/td676xr4322)"},
           "normalisation_check_constant_current": check_normalisation(cells),
           "dataset_protocol_heat_by_rate": summarise(rows, floors),
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
    print(json.dumps({k: v for k, v in o.items() if k != "normalisation_check_constant_current"}, indent=1)[:6000])
