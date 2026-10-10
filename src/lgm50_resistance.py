#!/usr/bin/env python3
"""DC resistance and cycling heat of the proxy LG M50 cell from the PyBaMM "ORegan2022" parameter set.

Neither source named by the brief for the resistance reports a full-cell DC resistance that could be used:
Chen et al. 2020 (J. Electrochem. Soc. 167, 080534) gives only half-cell charge-transfer resistances (its Table VI),
and the full text of O'Regan et al. 2022 (Electrochim. Acta 425, 140700) was not accessible in this study
(DECISIONS.md D8). The resistance is therefore DERIVED: the Doyle-Fuller-Newman model of PyBaMM with the
ORegan2022 parameter set (which carries the temperature dependence O'Regan et al. measured), isothermal at the
dataset's 35 C, on the cell's rated 2.5 to 4.2 V window (Chen et al. 2020, p. 9: "The rated capacity of the LGM50
cell between 2.5 - 4.2 V is cited to be 5 Ah").

Reported (inputs/lgm50_resistance.json):
  R_DC pulse values: (V before the pulse - V at the end of the pulse) / I for a discharge pulse from rest at 50 % SOC,
        35 C, for pulse lengths 10 s and 30 s and currents C/2 and 1C (for comparison).
  R_eff: the time-averaged irreversible heat (total heating minus reversible heating) of a full constant-current
        discharge (4.2 to 2.5 V) and of a full constant-current charge (2.5 to 4.2 V, no CV), divided by I^2, at
        C/16, C/10, C/5, C/2, 1C and 2C. The REGISTERED resistance is the mean of the C/2 discharge and C/2 charge
        values: the resistance whose I^2 R reproduces the model's own heat of steady cycling at C/2 (D8).
The model has no contact or tab resistance (the parameter set's contact resistance is 0 ohm), so these values are
a lower estimate of the resistance of a real cell (D8).

Usage: python3 src/lgm50_resistance.py --out inputs/lgm50_resistance.json
"""
import argparse, json
import numpy as np
import pybamm

T_K = 35.0 + 273.15


def params():
    p = pybamm.ParameterValues("ORegan2022")
    p.update({"Ambient temperature [K]": T_K, "Initial temperature [K]": T_K,
              "Upper voltage cut-off [V]": 4.2, "Open-circuit voltage at 100% SOC [V]": 4.2,
              "Lower voltage cut-off [V]": 2.5, "Open-circuit voltage at 0% SOC [V]": 2.5})
    return p


def pulse_resistance(c_rate, seconds, soc=0.5):
    p = params()
    cap = p["Nominal cell capacity [A.h]"]
    model = pybamm.lithium_ion.DFN()
    exp = pybamm.Experiment([f"Rest for 600 seconds", f"Discharge at {c_rate}C for {seconds} seconds"],
                            period="0.1 seconds")
    sim = pybamm.Simulation(model, parameter_values=p, experiment=exp)
    sol = sim.solve(initial_soc=soc)
    v_rest = float(sol.cycles[0]["Voltage [V]"].entries[-1]) if sol.cycles else float(sol.sub_solutions[0]["Voltage [V]"].entries[-1])
    v_end = float(sol["Voltage [V]"].entries[-1])
    I = c_rate * cap
    return {"c_rate": c_rate, "pulse_s": seconds, "soc": soc, "T_C": T_K - 273.15, "current_A": I,
            "V_rest": v_rest, "V_end": v_end, "R_ohm": (v_rest - v_end) / I}


def cycle_heat(c_rate, direction):
    p = params()
    cap = p["Nominal cell capacity [A.h]"]
    model = pybamm.lithium_ion.DFN(options={"calculate heat source for isothermal models": "true"})
    if direction == "discharge":
        exp = pybamm.Experiment([f"Discharge at {c_rate}C until 2.5 V"], period="10 seconds")
        soc0 = 1.0
    else:
        exp = pybamm.Experiment([f"Charge at {c_rate}C until 4.2 V"], period="10 seconds")
        soc0 = 0.0
    sim = pybamm.Simulation(model, parameter_values=p, experiment=exp)
    sol = sim.solve(initial_soc=soc0)
    t = sol["Time [s]"].entries
    names = [n for n in ("Total heating [W]", "Ohmic heating [W]", "Irreversible electrochemical heating [W]",
                         "Reversible heating [W]") if n in sol.all_models[0].variables]
    out = {"c_rate": c_rate, "direction": direction, "T_C": T_K - 273.15, "current_A": c_rate * cap,
           "duration_s": float(t[-1] - t[0]), "charge_Ah": float(c_rate * cap * (t[-1] - t[0]) / 3600)}
    for n in names:
        q = sol[n].entries
        out["mean " + n] = float(np.trapezoid(q, t) / (t[-1] - t[0]))
    tot = out.get("mean Total heating [W]")
    rev = out.get("mean Reversible heating [W]", 0.0)
    if tot is not None:
        out["mean irreversible heat [W]"] = tot - rev
        out["R_eff_irreversible_ohm"] = (tot - rev) / (c_rate * cap) ** 2
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default="inputs/lgm50_resistance.json"); a = ap.parse_args()
    res = {"source": "PyBaMM DFN, parameter set ORegan2022 (O'Regan et al. 2022 with Chen et al. 2020), isothermal",
           "pybamm": pybamm.__version__, "window_V": [2.5, 4.2], "nominal_capacity_Ah": float(params()["Nominal cell capacity [A.h]"]),
           "contact_resistance_ohm": float(params()["Contact resistance [Ohm]"]),
           "pulses": [pulse_resistance(c, s) for c in (0.5, 1.0) for s in (10, 30)],
           "cycle_heat": [cycle_heat(round(c, 6), d) for c in (1 / 16, 1 / 10, 1 / 5, 1 / 2, 1.0, 2.0) for d in ("discharge", "charge")]}
    c2 = [x["R_eff_irreversible_ohm"] for x in res["cycle_heat"] if abs(x["c_rate"] - 0.5) < 1e-9]
    res["registered"] = {"R_DC_ohm": float(np.mean(c2)),
                         "definition": "mean of the C/2 discharge and C/2 charge values of (time-averaged irreversible heat) / I^2, "
                                       "DFN ORegan2022, isothermal 35 C, 2.5 to 4.2 V"}
    with open(a.out, "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))
