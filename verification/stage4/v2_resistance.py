#!/usr/bin/env python3
"""V2 of Stage 4 (DECISIONS.md D62): the registered resistance at 35, 25 and 10 C, rederived by separate code.

Written without src/lgm50_resistance.py or any Stage 4 function. The DFN of PyBaMM with the ORegan2022 set and the
isothermal heat-source option is built here; the initial stoichiometries for 100 % and 0 % of the 2.5 to 4.2 V
window are set explicitly from PyBaMM's electrode state-of-health solver; a constant current of 2.5 A (C/2 of the
set's 5 Ah) is applied through "Current function [A]" without PyBaMM's experiment interface, the voltage events end
the run, and the ohmic plus irreversible electrochemical heat is integrated on a 1 s grid. R = mean over the
discharge and the charge of (time-averaged heat) / I^2. Passes where it agrees with results/stage4_resistance.json
within 1 % at each ambient.

Usage: python3 verification/stage4/v2_resistance.py
"""
import json, os, sys

import numpy as np
import pybamm

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))


def one_direction(T_C, sign):
    p = pybamm.ParameterValues("ORegan2022")
    T = T_C + 273.15
    p.update({"Ambient temperature [K]": T, "Initial temperature [K]": T,
              "Upper voltage cut-off [V]": 4.2, "Lower voltage cut-off [V]": 2.5,
              "Open-circuit voltage at 100% SOC [V]": 4.2, "Open-circuit voltage at 0% SOC [V]": 2.5})
    cap = p["Nominal cell capacity [A.h]"]
    current = sign * 0.5 * cap                                 # positive discharges in PyBaMM
    soc0 = 1.0 if sign > 0 else 0.0
    x0, y0 = pybamm.lithium_ion.get_initial_stoichiometries(soc0, p)
    p.update({"Initial concentration in negative electrode [mol.m-3]": float(x0) * p["Maximum concentration in negative electrode [mol.m-3]"],
              "Initial concentration in positive electrode [mol.m-3]": float(y0) * p["Maximum concentration in positive electrode [mol.m-3]"],
              "Current function [A]": current})
    model = pybamm.lithium_ion.DFN(options={"calculate heat source for isothermal models": "true"})
    sim = pybamm.Simulation(model, parameter_values=p, solver=pybamm.IDAKLUSolver(rtol=1e-8, atol=1e-10))
    t_end = 4.0 * 3600.0
    sol = sim.solve(t_eval=[0.0, t_end], t_interp=np.arange(0.0, t_end + 1.0, 1.0))
    t = sol["Time [s]"].entries
    q = sol["Ohmic heating [W]"].entries + sol["Irreversible electrochemical heating [W]"].entries
    v = sol["Voltage [V]"].entries
    mean_q = float(np.trapezoid(q, t) / (t[-1] - t[0]))
    return {"T_C": T_C, "direction": "discharge" if sign > 0 else "charge", "current_A": current, "x0": float(x0), "y0": float(y0),
            "duration_s": float(t[-1] - t[0]), "samples": int(t.size), "V_start": float(v[0]), "V_end": float(v[-1]),
            "termination": sol.termination, "mean_heat_W": mean_q, "R_ohm": mean_q / current ** 2}


def main():
    ref = json.load(open(os.path.join(ROOT, "results", "stage4_resistance.json")))["R_used_ohm"]
    out = {"check": "V2 (D62): the resistance rederived by separate code, within 1 %", "pybamm": pybamm.__version__, "ambients": {}}
    for T in (35.0, 25.0, 10.0):
        dis, ch = one_direction(T, +1), one_direction(T, -1)
        R = 0.5 * (dis["R_ohm"] + ch["R_ohm"])
        r0 = ref[f"{T:g}"]
        out["ambients"][f"{T:g}"] = {"discharge": dis, "charge": ch, "R_ohm": R, "R_registered_ohm": r0,
                                     "rel_diff": (R - r0) / r0, "pass": abs(R - r0) / r0 <= 0.01}
        print(f"{T:g} C: separate {R:.9g} ohm, registered {r0:.9g} ohm, rel diff {(R - r0) / r0:+.2e}", flush=True)
    out["pass"] = all(v["pass"] for v in out["ambients"].values())
    json.dump(out, open(os.path.join(HERE, "v2.json"), "w"), indent=1)
    print("V2 pass:", out["pass"])


if __name__ == "__main__":
    main()
