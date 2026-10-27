#!/usr/bin/env python3
"""Diagnostic, continued: the end of the C/2 charge at 10 C at fine output resolution.

Solves the checker's constant-current C/2 charge at 10 C with a coarse output grid up to 120 s before the end and
a fine one (0.1, 0.01 s) over the last 120 s; reports where the heat peaks, how large it is, the minimum electrolyte
concentration and the voltage there, and the average over the whole charge for each resolution. Also returns the
solver's own internal time points (no interpolation) when PyBaMM exposes them.
Output: out/diag_charge_10C_tail.json
"""
import json, os

import numpy as np
import pybamm

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
T_END = 3965.371110650174


def solve(t_interp, T_C=10.0):
    p = pybamm.ParameterValues("ORegan2022")
    T = T_C + 273.15
    p.update({"Ambient temperature [K]": T, "Initial temperature [K]": T, "Upper voltage cut-off [V]": 4.2,
              "Open-circuit voltage at 100% SOC [V]": 4.2, "Lower voltage cut-off [V]": 2.5,
              "Open-circuit voltage at 0% SOC [V]": 2.5, "Current function [A]": -2.5})
    model = pybamm.lithium_ion.DFN(options={"calculate heat source for isothermal models": "true"})
    sim = pybamm.Simulation(model, parameter_values=p)
    return sim.solve([0, 9360.0], initial_soc=0.0, t_interp=t_interp)


def main():
    res = {}
    for fine in (1.0, 0.1, 0.01):
        ti = np.concatenate([np.arange(0, T_END - 120, 10.0), np.arange(T_END - 120, T_END + 1, fine)])
        sol = solve(ti)
        t = sol["Time [s]"].entries
        q = sol["Ohmic heating [W]"].entries + sol["Irreversible electrochemical heating [W]"].entries
        V = sol["Voltage [V]"].entries
        ce = sol["Electrolyte concentration [mol.m-3]"].entries
        cmin = ce.min(axis=0)
        k = int(np.argmax(q))
        tail = t >= T_END - 120
        res[f"fine_{fine:g}s"] = {"points": int(len(t)), "t_end": float(t[-1]),
                                  "R_ohm_whole_charge": float(np.trapezoid(q, t) / (t[-1] - t[0]) / 6.25),
                                  "q_max_W": float(q[k]), "t_at_q_max": float(t[k]), "V_at_q_max": float(V[k]),
                                  "ce_min_at_q_max": float(cmin[k]), "ce_min_overall": float(cmin.min()),
                                  "t_at_ce_min": float(t[int(np.argmin(cmin))]),
                                  "heat_J_last_120s": float(np.trapezoid(q[tail], t[tail])),
                                  "heat_J_total": float(np.trapezoid(q, t))}
        if fine == 0.01:
            sel = np.where(tail)[0][::50]
            res["tail_samples_t_q_V_cemin"] = [[float(t[i]), float(q[i]), float(V[i]), float(cmin[i])] for i in sel]
        print(fine, res[f"fine_{fine:g}s"], flush=True)
    # the solver's own time points, if exposed
    try:
        sol = solve(None)
        ts = np.asarray(sol.t)
        res["solver_internal_points"] = int(len(ts))
        res["solver_last_20_t"] = [float(x) for x in ts[-20:]]
        q = sol["Ohmic heating [W]"].entries + sol["Irreversible electrochemical heating [W]"].entries
        res["solver_points_R_ohm"] = float(np.trapezoid(q, ts) / (ts[-1] - ts[0]) / 6.25)
        res["solver_points_q_last_20"] = [float(x) for x in q[-20:]]
        print("internal", res["solver_internal_points"], res["solver_points_R_ohm"], flush=True)
    except Exception as e:
        res["solver_internal_error"] = repr(e)[:300]
    with open(os.path.join(OUT, "diag_charge_10C_tail.json"), "w") as f:
        json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()
