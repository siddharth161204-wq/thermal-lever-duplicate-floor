#!/usr/bin/env python3
"""POST-FREEZE diagnostic (written after the freeze commit fd5d1e7, prompted by the pipeline's V2 code, which uses
IDAKLU rtol 1e-8, atol 1e-10). Is the sub-second heat spike of the C/2 charge at 10 C a solver-tolerance artifact?

The checker's constant-current C/2 charge at 10 C is solved with PyBaMM's default solver and with IDAKLU at
rtol 1e-8, atol 1e-10, each sampled at 0.01 s over the last 120 s; reported: the peak heat, the heat in the last 120 s
and the time-averaged irreversible heat over I^2 for the whole charge.
Output: out/postfreeze_diag_tolerance.json
"""
import json, os

import numpy as np
import pybamm

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
T_END = 3965.4


def run(solver):
    p = pybamm.ParameterValues("ORegan2022")
    T = 283.15
    p.update({"Ambient temperature [K]": T, "Initial temperature [K]": T, "Upper voltage cut-off [V]": 4.2,
              "Open-circuit voltage at 100% SOC [V]": 4.2, "Lower voltage cut-off [V]": 2.5,
              "Open-circuit voltage at 0% SOC [V]": 2.5, "Current function [A]": -2.5})
    model = pybamm.lithium_ion.DFN(options={"calculate heat source for isothermal models": "true"})
    sim = pybamm.Simulation(model, parameter_values=p, solver=solver)
    ti = np.concatenate([np.arange(0, T_END - 120, 1.0), np.arange(T_END - 120, T_END + 1, 0.01)])
    sol = sim.solve([0, 9360.0], initial_soc=0.0, t_interp=ti)
    t = sol["Time [s]"].entries
    q = sol["Ohmic heating [W]"].entries + sol["Irreversible electrochemical heating [W]"].entries
    tail = t >= t[-1] - 120
    k = int(np.argmax(q))
    return {"t_end": float(t[-1]), "q_max_W": float(q[k]), "t_at_q_max": float(t[k]),
            "heat_J_last_120s": float(np.trapezoid(q[tail], t[tail])),
            "R_ohm": float(np.trapezoid(q, t) / (t[-1] - t[0]) / 6.25)}


def main():
    res = {"default_solver": run(None), "idaklu_rtol1e-8_atol1e-10": run(pybamm.IDAKLUSolver(rtol=1e-8, atol=1e-10))}
    with open(os.path.join(OUT, "postfreeze_diag_tolerance.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
