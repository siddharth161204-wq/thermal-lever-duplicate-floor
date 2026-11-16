#!/usr/bin/env python3
"""Stage 6, the scorecard of predictions_stage6.md (DECISIONS.md D96 (2)).

P1 and P6 were scored at the end of Part A0 (`a0_analysis.py`) and are carried over unchanged; P2 to P5, P7 and P8 are
scored here from the registered pooled results of Parts B and C.

Usage: python3 -I src/stage6/s6_score.py
Writes results/stage6_scorecard.json.
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s6_common as C                   # noqa: E402

P7_LO, P7_HI = C.LEVERS_D82[("secondary 2C", "25 to 35 C")][1:]
FLOORS = (C.SD_E_D82["Geslin"], C.SD_E_D82["KIT"])


def main():
    W = json.load(open(C.path("results/stage6_within.json")))["registered"]["pooled"]
    B = json.load(open(C.path("results/stage6_between.json")))["registered"]["pooled"]
    a0 = json.load(open(C.path("results/stage6_scorecard.json")))
    sw = W["B1"]["s_w_per_K"]
    sb = B["C1"]["s_b_per_K"]
    diff = B["C2"]["s_b_minus_s_w_per_K"]
    sc = {"scored_at": {"P1": "end of Part A0", "P6": "end of Part A0", "others": "after Parts B and C"},
          "P1": a0["P1"], "P6": a0["P6"]}
    sc["P2"] = {"quantity": "pooled s_w and its 95 % interval, per K", "value": sw["value"], "ci95": sw["ci95"],
                "holds": bool(sw["value"] > 0 and sw["ci95"][0] > 0)}
    b2 = W.get("B2", {})
    if "s_w_R_held_per_K" in b2:
        held, same = b2["s_w_R_held_per_K"]["value"], b2["s_w_same_cells_per_K"]["value"]
        keeps = (held > 0) == (same > 0) and held != 0 and abs(held) >= 0.5 * abs(same)
        sc["P3"] = {"quantity": "temperature's lever with ln R held, against s_w over the same cells, per K",
                    "value": held, "same_cells_s_w": same, "ratio": held / same if same else None,
                    "cells": b2["cells"], "holds": bool(keeps)}
    else:
        sc["P3"] = {"holds": None, "note": "not scorable: no cell with a valid resistance"}
    sc["P4"] = {"quantity": "pooled s_b and its 95 % interval, per K", "value": sb["value"], "ci95": sb["ci95"],
                "holds": bool(sb["value"] > 0 and sb["ci95"][0] > 0)}
    sc["P5"] = {"quantity": "s_b - s_w and its 95 % interval, per K", "value": diff["value"], "ci95": diff["ci95"],
                "holds": bool(diff["ci95"][0] > 0)}
    sc["P7"] = {"quantity": "pooled s_w against Sandia's 25 to 35 C 2C interval, per K", "value": sw["value"],
                "interval": [P7_LO, P7_HI], "holds": bool(P7_LO <= sw["value"] <= P7_HI)}
    s2 = W["B4"]["sqrt2_sd_cell"]
    sc["P8"] = {"quantity": "sqrt(2) times the within-set residual sd of ln L after B1 (pair scale)", "value": s2,
                "ci95": W["B4"]["sqrt2_sd_cell_ci95"], "floors": list(FLOORS),
                "holds": bool(FLOORS[0] <= s2 <= FLOORS[1])}
    sc["held"] = [k for k in ("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8") if sc[k].get("holds") is True]
    sc["failed"] = [k for k in ("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8") if sc[k].get("holds") is False]
    sc["not_scored"] = [k for k in ("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8") if sc[k].get("holds") is None]
    json.dump(sc, open(C.path("results/stage6_scorecard.json"), "w"), indent=1)
    print(json.dumps({k: sc[k].get("holds") for k in ("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8")}))


if __name__ == "__main__":
    main()
