#!/usr/bin/env python3
"""Stage 4: Stage 1's resistance procedure at 35, 25 and 10 C (DECISIONS.md D58).

The functions pulse_resistance and cycle_heat of src/lgm50_resistance.py are imported unchanged; that module's
temperature T_K is set to the ambient before each call, so the DFN with the ORegan2022 set is isothermal at the
ambient on 2.5 to 4.2 V. The registered resistance at each ambient is, as at 35 C, the mean of the C/2 discharge and
C/2 charge values of (time-averaged irreversible heat) / I^2. Every pulse and rate is run and recorded; a rate other
than C/2 that fails is recorded as failed and enters nothing. The 35 C cases of the module use the archived value.

Usage: python3 src/stage4/resistance.py --out results/stage4_resistance.json
"""
import argparse, json, os, sys, time, traceback

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as CM                     # noqa: E402
import lgm50_resistance as LR           # noqa: E402  (Stage 1's procedure, unchanged)
import pybamm                           # noqa: E402

RATES_PROC = (1 / 16, 1 / 10, 1 / 5, 1 / 2, 1.0, 2.0)   # as in Stage 1's main


def procedure(T_C):
    LR.T_K = T_C + 273.15
    p = LR.params()
    res = {"T_C": T_C, "source": "PyBaMM DFN, parameter set ORegan2022 (O'Regan et al. 2022 with Chen et al. 2020), isothermal",
           "pybamm": pybamm.__version__, "window_V": [2.5, 4.2],
           "nominal_capacity_Ah": float(p["Nominal cell capacity [A.h]"]),
           "contact_resistance_ohm": float(p["Contact resistance [Ohm]"]),
           "ambient_K_in_parameters": float(p["Ambient temperature [K]"]),
           "pulses": [], "cycle_heat": []}
    for c in (0.5, 1.0):
        for s in (10, 30):
            t0 = time.time()
            try:
                r = LR.pulse_resistance(c, s)
            except Exception as e:                                   # recorded, enters nothing
                r = {"c_rate": c, "pulse_s": s, "failed": repr(e), "trace": traceback.format_exc(limit=2)}
            r["wall_s"] = time.time() - t0
            res["pulses"].append(r)
    for c in RATES_PROC:
        for d in ("discharge", "charge"):
            t0 = time.time()
            try:
                r = LR.cycle_heat(round(c, 6), d)
            except Exception as e:
                r = {"c_rate": round(c, 6), "direction": d, "failed": repr(e), "trace": traceback.format_exc(limit=2)}
            r["wall_s"] = time.time() - t0
            res["cycle_heat"].append(r)
    c2 = [x["R_eff_irreversible_ohm"] for x in res["cycle_heat"]
          if abs(x["c_rate"] - 0.5) < 1e-9 and "failed" not in x]
    if len(c2) != 2:
        raise SystemExit(f"the C/2 heat did not solve at {T_C} C: the registered resistance is undefined")
    res["registered"] = {"R_DC_ohm": float(np.mean(c2)),
                         "definition": "mean of the C/2 discharge and C/2 charge values of (time-averaged irreversible heat) / I^2, "
                                       f"DFN ORegan2022, isothermal {T_C:g} C, 2.5 to 4.2 V"}
    return res


def main(out):
    r35 = CM.read_r35()
    runs = {}
    for T in CM.AMBIENTS_C:
        t0 = time.time()
        runs[f"{T:g}"] = procedure(T)
        print(f"{T:g} C: R = {runs[f'{T:g}']['registered']['R_DC_ohm']:.12g} ohm ({time.time() - t0:.0f} s)", flush=True)
    LR.T_K = 35.0 + 273.15                                            # restore Stage 1's value
    used = {"35": r35, "25": runs["25"]["registered"]["R_DC_ohm"], "10": runs["10"]["registered"]["R_DC_ohm"]}
    doc = {"procedure": "src/lgm50_resistance.py, pulse_resistance and cycle_heat unchanged, T_K set to the ambient (D58)",
           "archived_35C": {"R_DC_ohm": r35, "file": "inputs/lgm50_resistance.json"},
           "rerun_35C_minus_archived_ohm": runs["35"]["registered"]["R_DC_ohm"] - r35,
           "rerun_35C_identical_to_archived": runs["35"]["registered"]["R_DC_ohm"] == r35,
           "R_used_ohm": used,
           "ratio_to_35C": {k: v / r35 for k, v in used.items()},
           "runs": runs}
    with open(out, "w") as f:
        json.dump(doc, f, indent=1)
    print(json.dumps({k: doc[k] for k in ("R_used_ohm", "ratio_to_35C", "rerun_35C_minus_archived_ohm", "rerun_35C_identical_to_archived")}, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(CM.ROOT, "results", "stage4_resistance.json"))
    main(ap.parse_args().out)
