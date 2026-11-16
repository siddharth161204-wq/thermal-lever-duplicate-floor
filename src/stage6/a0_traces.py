"""Stage 6, Part A0.4: each cell's temperature trace over cycles 2 to 100, and the contact rule of DECISIONS.md D88.
Temperatures only; no lifetime or capacity.

A cycle's swing is Tmax - Tmin. A thermocouple that reads the cell follows its heating in every charge and 4C
discharge; one that has lost contact, or never had it, reads mostly the chamber's air and swings far less than a cell
can on that policy. The rule flags a cell when its trace shows any of:
  (a) contact lost within the window: the median swing over cycles 91 to 100 below half that over cycles 2 to 11;
  (b) contact poor throughout: the median swing over cycles 2 to 100 below half the largest such median among the
      cells of its set (sets of two or more cells; every cell of a set runs the same policy);
  (c) a jump: a change of Tavg between two consecutive cycles larger than JUMP_K;
  (d) a flat line: FLAT_RUN or more consecutive cycles with an identical Tavg;
  (e) values below the chamber's: a cycle whose Tmin is more than BELOW_K below the chamber's set point (30 C), or
      a cycle of the window whose Tavg, Tmax or Tmin is missing or not finite;
  (f) the site's notes: a thermocouple found detached (2017-06-30, channels 7 and 21) or swapped with another cell's
      (channels 15 and 16 of 2017-05-12 and 2017-06-30, whose temperatures are exchanged before any measure, D87).
"""
import math

import numpy as np

CHAMBER_C = 30.0
HALF = 0.5
JUMP_K = 3.0
FLAT_RUN = 5
BELOW_K = 5.0
SITE_DETACHED = {"2017-06-30": {"7", "21"}}
SITE_SWAPPED = {"2017-05-12": ("15", "16"), "2017-06-30": ("15", "16")}


def per_cycle(rows):
    """Arrays of cycle, Tavg, Tmax and Tmin over cycles 2 to 100, and the number of incomplete cycles."""
    cy, a, x, m, bad = [], [], [], [], 0
    for r in rows:
        c = int(r["cycle"])
        if not 2 <= c <= 100:
            continue
        try:
            va, vx, vm = float(r["Tavg"]), float(r["Tmax"]), float(r["Tmin"])
        except (TypeError, ValueError):
            bad += 1
            continue
        if not (math.isfinite(va) and math.isfinite(vx) and math.isfinite(vm)):
            bad += 1
            continue
        cy.append(c); a.append(va); x.append(vx); m.append(vm)
    return np.array(cy), np.array(a), np.array(x), np.array(m), bad


def longest_flat(v):
    best = run = 1 if len(v) else 0
    for i in range(1, len(v)):
        run = run + 1 if v[i] == v[i - 1] else 1
        best = max(best, run)
    return best


def trace_stats(rows):
    cy, a, x, m, bad = per_cycle(rows)
    if len(cy) < 20:
        return {"n": int(len(cy)), "incomplete_cycles": bad}
    sw = x - m
    early, late = sw[cy <= 11], sw[cy >= 91]
    return {"n": int(len(cy)), "incomplete_cycles": bad, "swing_med": float(np.median(sw)),
            "swing_early": float(np.median(early)), "swing_late": float(np.median(late)),
            "floor_med": float(np.median(m)), "level_med": float(np.median(a)),
            "max_step": float(np.max(np.abs(np.diff(a)))), "flat_run": int(longest_flat(a)),
            "below_cycles": int(np.sum(m < CHAMBER_C - BELOW_K)), "min_Tmin": float(np.min(m))}


def flag(cells, stats):
    """Sets each entering cell's contact_flag, contact_reasons and site_note from the rules above."""
    best = {}
    for c in cells:
        if c["enters"] and c["set_size"] >= 2 and "swing_med" in stats[c["cell"]]:
            best[c["set"]] = max(best.get(c["set"], 0.0), stats[c["cell"]]["swing_med"])
    for c in cells:
        s = stats[c["cell"]]
        why = []
        if "swing_med" not in s:
            why.append("fewer than 20 complete cycles in 2 to 100")
        else:
            if s["swing_late"] < HALF * s["swing_early"]:
                why.append(f"contact lost: swing {s['swing_early']:.2f} K in cycles 2-11, {s['swing_late']:.2f} K in 91-100")
            if c["enters"] and c["set_size"] >= 2 and s["swing_med"] < HALF * best[c["set"]]:
                why.append(f"contact poor: swing {s['swing_med']:.2f} K against {best[c['set']]:.2f} K in its set")
            if s["max_step"] > JUMP_K:
                why.append(f"jump of {s['max_step']:.2f} K between cycles")
            if s["flat_run"] >= FLAT_RUN:
                why.append(f"flat line over {s['flat_run']} cycles")
            if s["below_cycles"] > 0:
                why.append(f"{s['below_cycles']} cycles with Tmin below {CHAMBER_C - BELOW_K:.0f} C")
        if s.get("incomplete_cycles", 0) > 0:
            why.append(f"{s['incomplete_cycles']} cycles with a missing temperature")
        ch = str(c["channel"]).strip()
        notes = []
        if ch in SITE_DETACHED.get(c["batch"], set()):
            notes.append("site: thermocouple found detached at unloading")
        if ch in SITE_SWAPPED.get(c["batch"], ()):
            notes.append("site: thermocouple swapped with channel " +
                         [x for x in SITE_SWAPPED[c["batch"]] if x != ch][0] + "; temperatures exchanged")
        c["contact_reasons"] = "; ".join(why + notes)
        c["contact_flag"] = bool(why or notes)
        c["site_note"] = "; ".join(notes)
        c.update({"swing_med": s.get("swing_med", float("nan")), "swing_early": s.get("swing_early", float("nan")),
                  "swing_late": s.get("swing_late", float("nan")), "floor_med": s.get("floor_med", float("nan"))})
