#!/usr/bin/env python3
"""A separate derivation of the C/2 resistance, written by the checker without src/lgm50_resistance.py.

The DFN of PyBaMM with the ORegan2022 set and the same voltage window is built here, isothermal at the ambient,
and driven by a constant current through "Current function [A]" (no Experiment interface); the model's own
voltage cut-off events end the run. The irreversible heat is the model's ohmic plus irreversible electrochemical
heating, integrated by the trapezoidal rule on a 1 s grid; R = mean heat / I^2, and the value at each ambient is the
mean of the C/2 discharge (from 100 % SOC) and C/2 charge (from 0 % SOC). Also run: a sweep of the same quantity
from 10 to 90 C, used only for the checker's own sensitivity (the heat at each cell's own temperature), which is
not a registered quantity.

Usage: python3 resistance_check.py [--sweep]
Output: out/resistance_check.json
"""
import argparse, json, os, sys, time

import numpy as np
import pybamm

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")


def run(T_C, c_rate, direction):
    p = pybamm.ParameterValues("ORegan2022")
    T = T_C + 273.15
    p.update({"Ambient temperature [K]": T, "Initial temperature [K]": T,
              "Upper voltage cut-off [V]": 4.2, "Open-circuit voltage at 100% SOC [V]": 4.2,
              "Lower voltage cut-off [V]": 2.5, "Open-circuit voltage at 0% SOC [V]": 2.5})
    cap = p["Nominal cell capacity [A.h]"]
    I = c_rate * cap
    p.update({"Current function [A]": I if direction == "discharge" else -I})
    model = pybamm.lithium_ion.DFN(options={"calculate heat source for isothermal models": "true"})
    sim = pybamm.Simulation(model, parameter_values=p)
    t_end = 3600.0 / c_rate * 1.3
    sol = sim.solve([0, t_end], initial_soc=1.0 if direction == "discharge" else 0.0, t_interp=np.arange(0.0, t_end, 1.0))
    t = sol["Time [s]"].entries
    q = sol["Ohmic heating [W]"].entries + sol["Irreversible electrochemical heating [W]"].entries
    mean_q = float(np.trapezoid(q, t) / (t[-1] - t[0]))
    # spike screen: a sample more than 5 times the median of its 61-sample neighbourhood is replaced by that median
    # (the 10 C charge carries a sub-second heat spike of order 1e4 W with a smooth voltage, see diag_charge_10C_tail.py)
    w = 30
    med = np.array([np.median(q[max(0, i - w):i + w + 1]) for i in range(len(q))])
    bad = np.abs(q) > 5.0 * np.abs(med)
    qs = np.where(bad, med, q)
    mean_qs = float(np.trapezoid(qs, t) / (t[-1] - t[0]))
    return {"T_C": T_C, "c_rate": c_rate, "direction": direction, "duration_s": float(t[-1] - t[0]),
            "V_last": float(sol["Voltage [V]"].entries[-1]), "termination": str(sol.termination),
            "charge_Ah": float(I * (t[-1] - t[0]) / 3600.0), "mean_irreversible_W": mean_q, "R_ohm": mean_q / I ** 2,
            "spike_samples": int(bad.sum()), "spike_t_s": [float(x) for x in t[bad]][:20],
            "spike_q_W": [float(x) for x in q[bad]][:20], "R_ohm_spike_screened": mean_qs / I ** 2}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--sweep", action="store_true"); a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    temps = [35.0, 25.0, 10.0]
    if a.sweep:
        temps = [10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 70.0, 80.0, 90.0]
    rows = []
    for T in temps:
        rr = {}
        for d in ("discharge", "charge"):
            t0 = time.time()
            try:
                r = run(T, 0.5, d)
                r["status"] = "ok"
            except Exception as e:
                r = {"T_C": T, "c_rate": 0.5, "direction": d, "status": "failed", "error": repr(e)[:400]}
            r["seconds"] = time.time() - t0
            rr[d] = r
            print(r, flush=True)
        ok = all(x["status"] == "ok" for x in rr.values())
        R = 0.5 * (rr["discharge"]["R_ohm"] + rr["charge"]["R_ohm"]) if ok else None
        Rs = 0.5 * (rr["discharge"]["R_ohm_spike_screened"] + rr["charge"]["R_ohm_spike_screened"]) if ok else None
        rows.append({"T_C": T, "R_ohm": R, "R_ohm_spike_screened": Rs, "runs": rr})
    name = "resistance_sweep.json" if a.sweep else "resistance_check.json"
    with open(os.path.join(OUT, name), "w") as f:
        json.dump({"method": __doc__.strip().splitlines()[0], "rows": rows}, f, indent=1)
    for r in rows:
        print(r["T_C"], r["R_ohm"], r["R_ohm_spike_screened"],
              {d: r["runs"][d].get("spike_samples") for d in ("discharge", "charge")})


if __name__ == "__main__":
    main()
