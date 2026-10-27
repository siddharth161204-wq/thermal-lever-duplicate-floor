#!/usr/bin/env python3
"""Diagnostic: why the time-averaged irreversible heat of the C/2 charge at 10 C depends on the output grid.

Runs the C/2 charge at 10 C (and, for comparison, at 35 C) with the checker's constant-current set-up on output
grids of 10, 1 and 0.1 s, and prints the heat near the end of the charge, where the voltage reaches 4.2 V.
Also reruns the procedure's own experiment (period 10 s) with a 1 s period, to see whether the registered value
moves with the period.
Output: out/diag_charge_10C.json
"""
import json, os, sys

import numpy as np
import pybamm

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = os.path.join(HERE, "out")


def params(T_C):
    p = pybamm.ParameterValues("ORegan2022")
    T = T_C + 273.15
    p.update({"Ambient temperature [K]": T, "Initial temperature [K]": T, "Upper voltage cut-off [V]": 4.2,
              "Open-circuit voltage at 100% SOC [V]": 4.2, "Lower voltage cut-off [V]": 2.5,
              "Open-circuit voltage at 0% SOC [V]": 2.5})
    return p


def cc_run(T_C, direction, grid):
    p = params(T_C)
    p.update({"Current function [A]": 2.5 if direction == "discharge" else -2.5})
    model = pybamm.lithium_ion.DFN(options={"calculate heat source for isothermal models": "true"})
    sim = pybamm.Simulation(model, parameter_values=p)
    sol = sim.solve([0, 9360.0], initial_soc=1.0 if direction == "discharge" else 0.0, t_interp=np.arange(0, 9360.0, grid))
    t = sol["Time [s]"].entries
    q = sol["Ohmic heating [W]"].entries + sol["Irreversible electrochemical heating [W]"].entries
    return t, q, sol["Voltage [V]"].entries


def exp_run(T_C, direction, period):
    p = params(T_C)
    model = pybamm.lithium_ion.DFN(options={"calculate heat source for isothermal models": "true"})
    step = "Discharge at 0.5C until 2.5 V" if direction == "discharge" else "Charge at 0.5C until 4.2 V"
    exp = pybamm.Experiment([step], period=f"{period} seconds")
    sim = pybamm.Simulation(model, parameter_values=p, experiment=exp)
    sol = sim.solve(initial_soc=1.0 if direction == "discharge" else 0.0)
    t = sol["Time [s]"].entries
    q = sol["Total heating [W]"].entries - sol["Reversible heating [W]"].entries
    return float(np.trapezoid(q, t) / (t[-1] - t[0]) / 6.25), len(t), float(t[-1])


def main():
    res = {"cc": [], "experiment": []}
    for T_C in (10.0, 35.0):
        for d in ("charge", "discharge"):
            for g in (10.0, 1.0, 0.1):
                t, q, V = cc_run(T_C, d, g)
                cum = np.concatenate([[0.0], np.cumsum(0.5 * (q[1:] + q[:-1]) * np.diff(t))])
                rec = {"T_C": T_C, "direction": d, "grid_s": g, "points": len(t), "t_end_s": float(t[-1]),
                       "R_ohm": float(cum[-1] / (t[-1] - t[0]) / 6.25),
                       "share_last_30s": float((cum[-1] - np.interp(t[-1] - 30, t, cum)) / cum[-1]),
                       "share_last_2s": float((cum[-1] - np.interp(t[-1] - 2, t, cum)) / cum[-1]),
                       "q_last_W": float(q[-1]), "q_median_W": float(np.median(q)), "V_last": float(V[-1])}
                if g == 0.1:
                    tail = np.searchsorted(t, t[-1] - 5.0)
                    rec["tail_t_q"] = [[float(a), float(b)] for a, b in zip(t[tail::5], q[tail::5])]
                res["cc"].append(rec)
                print({k: v for k, v in rec.items() if k != "tail_t_q"}, flush=True)
            for per in (10, 1):
                R, npts, tend = exp_run(T_C, d, per)
                res["experiment"].append({"T_C": T_C, "direction": d, "period_s": per, "R_ohm": R, "points": npts, "t_end_s": tend})
                print(res["experiment"][-1], flush=True)
    with open(os.path.join(OUT, "diag_charge_10C.json"), "w") as f:
        json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()
