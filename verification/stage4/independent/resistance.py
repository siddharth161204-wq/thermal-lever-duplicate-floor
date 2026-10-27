#!/usr/bin/env python3
"""Independent check of Stage 4, part (a): the resistance at 35, 25 and 10 C by Stage 1's procedure (D58).

D58 says the procedure is run unchanged: the functions pulse_resistance and cycle_heat of src/lgm50_resistance.py
are imported and that module's T_K is set to the ambient before each call. That is what this script does. It adds
nothing to the computation; it only records, for every call, whether it raised, how long it took, and (through a
wrapper around pybamm.Simulation.solve that returns the solution untouched) the solver's termination message and the
last voltage, so that a run that stopped short of its cut-off can be seen. The registered value at each ambient is
the mean of the C/2 discharge and C/2 charge values of R_eff_irreversible_ohm (D8, D58).

Written by the independent checker from DECISIONS.md D56 to D63 and the brief; no Stage 4 pipeline file was read.

Usage: python3 verification/stage4/independent/resistance.py
Output: verification/stage4/independent/out/resistance.json
"""
import json, logging, math, os, sys, time, traceback, warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = os.path.join(HERE, "out")
sys.path.insert(0, os.path.join(ROOT, "src"))

import numpy as np                       # noqa: E402
import pybamm                            # noqa: E402
import lgm50_resistance as LR            # noqa: E402  (Stage 1's procedure, unchanged)

ARCHIVED_35C = 0.01871279570954531
AMBIENTS_C = [35.0, 25.0, 10.0]
PULSES = [(c, s) for c in (0.5, 1.0) for s in (10, 30)]
RATES = [round(c, 6) for c in (1 / 16, 1 / 10, 1 / 5, 1 / 2, 1.0, 2.0)]

# --- a transparent wrapper that keeps the last solution, so its termination can be recorded -------------------
_LAST = {}
_orig_solve = pybamm.Simulation.solve


def _solve(self, *args, **kwargs):
    sol = _orig_solve(self, *args, **kwargs)
    _LAST["sol"] = sol
    return sol


pybamm.Simulation.solve = _solve


class _Grab(logging.Handler):
    def __init__(self):
        super().__init__(level=logging.WARNING)
        self.msgs = []

    def emit(self, record):
        self.msgs.append(f"{record.levelname}: {record.getMessage()}"[:400])


def _diag(direction=None):
    sol = _LAST.get("sol")
    if sol is None:
        return {}
    d = {"termination": str(getattr(sol, "termination", "")),
         "V_last": float(sol["Voltage [V]"].entries[-1]),
         "t_last_s": float(sol["Time [s]"].entries[-1]),
         "n_steps_solved": len(getattr(sol, "steps", []) or [])}
    if direction == "discharge":
        d["reached_cutoff"] = bool(abs(d["V_last"] - 2.5) < 1e-3)
    elif direction == "charge":
        d["reached_cutoff"] = bool(abs(d["V_last"] - 4.2) < 1e-3)
    return d


def run_call(fn, args, direction=None):
    grab = _Grab()
    pybamm.logger.addHandler(grab)
    _LAST.clear()
    t0 = time.time()
    rec = {}
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            out = fn(*args)
        rec.update(out)
        rec["status"] = "ok"
        rec["warnings"] = [str(x.message)[:400] for x in w]
    except Exception as e:                                    # a failure is recorded and enters nothing (D58)
        rec["status"] = "failed"
        rec["error"] = repr(e)[:800]
        rec["traceback_tail"] = traceback.format_exc()[-1500:]
    rec["seconds"] = time.time() - t0
    rec["log_warnings"] = grab.msgs
    try:
        rec.update(_diag(direction))
    except Exception as e:
        rec["diag_error"] = repr(e)[:300]
    pybamm.logger.removeHandler(grab)
    return rec


def main():
    os.makedirs(OUT, exist_ok=True)
    res = {"procedure": "src/lgm50_resistance.py imported unchanged, T_K set to the ambient before each call (D58)",
           "pybamm": pybamm.__version__, "archived_35C_ohm": ARCHIVED_35C, "ambients": {}}
    for T_C in AMBIENTS_C:
        LR.T_K = T_C + 273.15
        amb = {"T_K": LR.T_K, "pulses": [], "cycle_heat": []}
        for c, s in PULSES:
            r = run_call(LR.pulse_resistance, (c, s))
            r.setdefault("c_rate", c); r.setdefault("pulse_s", s)
            amb["pulses"].append(r)
            print(f"{T_C:5.1f} C pulse {c}C {s}s: {r['status']} R={r.get('R_ohm')} ({r['seconds']:.1f} s)", flush=True)
        for c in RATES:
            for d in ("discharge", "charge"):
                r = run_call(LR.cycle_heat, (c, d), direction=d)
                r.setdefault("c_rate", c); r.setdefault("direction", d)
                amb["cycle_heat"].append(r)
                print(f"{T_C:5.1f} C {d:9s} {c:<7}: {r['status']} R_eff={r.get('R_eff_irreversible_ohm')} "
                      f"Ah={r.get('charge_Ah')} V_last={r.get('V_last')} term={r.get('termination')} ({r['seconds']:.1f} s)",
                      flush=True)
        c2 = [x for x in amb["cycle_heat"] if abs(x["c_rate"] - 0.5) < 1e-9]
        ok = all(x["status"] == "ok" and "R_eff_irreversible_ohm" in x for x in c2)
        amb["C2_values_ohm"] = {x["direction"]: x.get("R_eff_irreversible_ohm") for x in c2}
        amb["R_registered_ohm"] = float(np.mean([x["R_eff_irreversible_ohm"] for x in c2])) if ok else None
        amb["failed_calls"] = [f"{x.get('c_rate')} {x.get('direction', 'pulse ' + str(x.get('pulse_s')) + ' s')}"
                               for x in amb["pulses"] + amb["cycle_heat"] if x["status"] != "ok"]
        amb["short_of_cutoff"] = [f"{x['c_rate']} {x['direction']} (V_last {x.get('V_last')}, {x.get('charge_Ah')} Ah)"
                                  for x in amb["cycle_heat"] if x["status"] == "ok" and x.get("reached_cutoff") is False]
        res["ambients"][f"{T_C:g}"] = amb
        with open(os.path.join(OUT, "resistance.json"), "w") as f:
            json.dump(res, f, indent=1)
    R35, R25, R10 = (res["ambients"][k]["R_registered_ohm"] for k in ("35", "25", "10"))
    res["used"] = {"R_35C_ohm (archived, used by the 35 C cases)": ARCHIVED_35C,
                   "R_35C_rerun_ohm": R35, "rerun_minus_archived_ohm": None if R35 is None else R35 - ARCHIVED_35C,
                   "R_25C_ohm": R25, "R_10C_ohm": R10,
                   "ratio_25_over_35_archived": None if R25 is None else R25 / ARCHIVED_35C,
                   "ratio_10_over_35_archived": None if R10 is None else R10 / ARCHIVED_35C}
    with open(os.path.join(OUT, "resistance.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res["used"], indent=1))


if __name__ == "__main__":
    main()
