#!/usr/bin/env python3
"""Stage 4, after the results: the heat traces behind the resistance procedure (D58), checked for isolated spikes.

Reruns, with Stage 1's parameter function unchanged and its temperature set to the ambient, the experiments of
lgm50_resistance.cycle_heat for the C/2 discharge and charge (which set the registered resistance) and the C/16
charge (whose values at 35 and 25 C stand apart) at 35, 25 and 10 C, and records per run: samples, the median and
largest irreversible heat (total minus reversible), where the largest falls, the time-averaged heat with and without
the samples above five times the median, the end voltage and the termination. Output:
verification/stage4/procedure_heat_check.json.

Usage: python3 verification/stage4/procedure_heat_check.py
"""
import json, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "src"))
import lgm50_resistance as LR           # noqa: E402  (Stage 1, unchanged)
import pybamm                           # noqa: E402

RUNS = [("C/2 discharge", 0.5, "Discharge at 0.5C until 2.5 V", 1.0), ("C/2 charge", 0.5, "Charge at 0.5C until 4.2 V", 0.0),
        ("C/16 charge", 0.0625, "Charge at 0.0625C until 4.2 V", 0.0)]


def main():
    out = {"check": "irreversible heat traces of the resistance procedure", "runs": []}
    for T in (35.0, 25.0, 10.0):
        LR.T_K = T + 273.15
        for name, c, step, soc in RUNS:
            p = LR.params()
            model = pybamm.lithium_ion.DFN(options={"calculate heat source for isothermal models": "true"})
            sim = pybamm.Simulation(model, parameter_values=p, experiment=pybamm.Experiment([step], period="10 seconds"))
            sol = sim.solve(initial_soc=soc)
            t = sol["Time [s]"].entries
            q = sol["Total heating [W]"].entries - sol["Reversible heating [W]"].entries
            med = float(np.median(q))
            i = int(np.argmax(q))
            mean = float(np.trapezoid(q, t) / (t[-1] - t[0]))
            mean_clip = float(np.trapezoid(np.minimum(q, 5 * med), t) / (t[-1] - t[0]))
            I = c * 5.0
            spikes = np.flatnonzero(q > 5 * med)
            out["runs"].append({"T_C": T, "run": name, "samples": int(t.size), "median_W": med, "max_W": float(q[i]),
                                "max_at_fraction_of_run": float((t[i] - t[0]) / (t[-1] - t[0])),
                                "samples_above_5x_median": int(spikes.size), "spike_times_s": [float(x) for x in t[spikes][:10]],
                                "mean_W": mean, "mean_without_spikes_W": mean_clip, "R_ohm": mean / I ** 2,
                                "R_without_spikes_ohm": mean_clip / I ** 2, "V_end": float(sol["Voltage [V]"].entries[-1]),
                                "termination": str(sol.termination)})
            r = out["runs"][-1]
            print(f"{T:g} C {name}: R {r['R_ohm']:.6g} ohm, without spikes {r['R_without_spikes_ohm']:.6g}, "
                  f"{r['samples_above_5x_median']} spikes, max {r['max_W']:.4g} W, V_end {r['V_end']:.4f}", flush=True)
    LR.T_K = 35.0 + 273.15
    json.dump(out, open(os.path.join(HERE, "procedure_heat_check.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
