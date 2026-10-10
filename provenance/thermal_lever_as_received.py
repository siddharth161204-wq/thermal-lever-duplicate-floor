#!/usr/bin/env python3
"""
The thermal lever against the duplicate-cell floor: Step 1 (the lever) and Step 2(a) (the lumped module).

Step 1: for each ageing law, lifetime L(T) to a capacity threshold as a function of steady cell temperature,
the sensitivity s = -d ln L / dT at 35 C, and the temperature difference dT_floor = f / s that produces a
lifetime difference equal to a floor f (3.0 %, 4.1 %, 8.7 % from doi.org/10.5281/zenodo.23205413), with
Monte Carlo propagation of the law's stated parameter uncertainties.

Step 2(a): an eight-node steady module network (two rows of four 21700 cells at 23.1 mm pitch): busbar strip
conductance, air-gap conduction and radiation between neighbours, convection and radiation to 35 C air,
I^2 R heat per cell; solved by Newton iteration; reports the hottest-minus-coldest cell spread.

EVERY value marked FILL comes from a named published source recorded in DECISIONS.md. Never invent a value.
Do not fit anything to the floor.

Usage:
  python3 thermal_lever.py --out results
"""
import argparse, json, math, os, csv
import numpy as np

SIGMA = 5.670374419e-8
R_GAS = 8.314462618
T_AMB_C = 35.0

# ---------------- the floor (copy from the published re-analysis's results.json; cite its DOI) ----------------
FLOORS = {"mean_3.0pct": 0.030, "worst_of_8_4.1pct": 0.041, "p95_8.7pct": 0.087}   # fractions; confirm against results.json

# ---------------- ageing laws: FILL from the papers' tables ----------------
# Law A, Arrhenius-with-throughput form (the form of Wang et al. 2011, J. Power Sources 196, 3942-3948):
#   Q_loss (%) = B * exp(-(Ea - b * Crate) / (R T)) * Ah^z     (set b = 0 if the law has no rate term)
LAW_A = {"name": None, "source": None, "B": None, "Ea": None, "b": None, "z": None,
         "sigma_Ea": None,      # standard error of Ea if the paper gives one, else None
         "capacity_Ah": None}   # cell capacity to convert throughput to EFC
# Law B: a second law with its own expression. Fill `loss_B(T, efc, crate)` below from the paper, or set LAW_B to None.
LAW_B = {"name": None, "source": None}

THRESHOLD_LOSS_PCT = 10.0     # lifetime defined at 90 % SOH, as in the floor
CRATES = [1 / 16, 1 / 10, 1 / 5, 1 / 2, 1.0, 2.0]   # the dataset's average rates, then what it would take


def loss_A(T, efc, crate, p):
    Ah = efc * p["capacity_Ah"]
    return p["B"] * math.exp(-(p["Ea"] - p["b"] * crate) / (R_GAS * T)) * Ah ** p["z"]


def loss_B(T, efc, crate, p):
    raise NotImplementedError("FILL the second law's expression from its paper, or set LAW_B = None")


def lifetime(loss_fn, T, crate, p, lo=1.0, hi=1e6):
    """EFC at which the law's capacity loss reaches THRESHOLD_LOSS_PCT, by bisection."""
    f = lambda efc: loss_fn(T, efc, crate, p) - THRESHOLD_LOSS_PCT
    if f(hi) < 0:
        raise SystemExit(f"the law never reaches {THRESHOLD_LOSS_PCT} % loss within {hi:g} EFC at {T - 273.15:.1f} C; check the parameters and units")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0: hi = mid
        else: lo = mid
    return 0.5 * (lo + hi)


def sensitivity(loss_fn, p, crate, T_C=35.0, dT=0.5):
    T = T_C + 273.15
    Lm, Lp = lifetime(loss_fn, T - dT, crate, p), lifetime(loss_fn, T + dT, crate, p)
    return -(math.log(Lp) - math.log(Lm)) / (2 * dT)     # per kelvin


def lever(loss_fn, p, crate, n_mc=2000, seed=20261010):
    s = sensitivity(loss_fn, p, crate)
    out = {"crate": crate, "s_per_K": s, "dT_floor_K": {k: f / s for k, f in FLOORS.items()}}
    if p.get("sigma_Ea"):
        rng = np.random.default_rng(seed); ss = []
        for _ in range(n_mc):
            q = dict(p); q["Ea"] = rng.normal(p["Ea"], p["sigma_Ea"]); ss.append(sensitivity(loss_fn, q, crate))
        ss = np.array(ss)
        out["dT_floor_K_2.5_97.5"] = {k: [float(f / np.percentile(ss, 97.5)), float(f / np.percentile(ss, 2.5))] for k, f in FLOORS.items()}
    return out


# ---------------- proxy cell and module: FILL from the published sources ----------------
D_CELL, L_CELL = 0.0211, 0.0701          # LG M50, m
R_DC = None          # FILL ohm, DC resistance of the proxy cell at the recorded SOC and temperature (source in DECISIONS.md)
CAP_AH = None        # FILL Ah, capacity used to convert C-rate to current
PITCH = 0.0231       # m, as in the two-cell study
GAP = PITCH - D_CELL
K_AIR = 0.03         # W/(m K), as in the two-cell study
EPS = 0.8
G_STRIP = None       # FILL W/K per nickel strip, from the two-cell study's value (12 mW/K for two strips)
H_LIST = [5.0, 10.0, 25.0]   # W/(m2 K)
A_LAT = math.pi * D_CELL * L_CELL
A_END = math.pi * (D_CELL / 2) ** 2
A_FACE = A_LAT * 0.25        # facing-area fraction for the gap conduction, as the two-cell study used
X = 1.0 + GAP / D_CELL
F_VIEW = (math.sqrt(X * X - 1.0) + math.asin(1.0 / X) - X) / math.pi     # parallel equal cylinders

# two rows of four: neighbour pairs (index: row*4 + column)
NEIGHBOURS = [(0, 1), (1, 2), (2, 3), (4, 5), (5, 6), (6, 7), (0, 4), (1, 5), (2, 6), (3, 7)]


def exposed_area(i):
    """Lateral area not facing a neighbour, plus the two ends."""
    n_nb = sum(1 for a, b in NEIGHBOURS if i in (a, b))
    return A_LAT - n_nb * A_FACE + 2 * A_END


def module_steady(q_cell, h, T_amb=T_AMB_C + 273.15, iters=200):
    """Newton solve of the eight-node steady balance with conduction, radiation and convection."""
    T = np.full(8, T_amb + 5.0)
    g_gap = K_AIR * A_FACE / GAP
    for _ in range(iters):
        F = np.zeros(8); J = np.zeros((8, 8))
        for i in range(8):
            A = exposed_area(i)
            F[i] += q_cell - h * A * (T[i] - T_amb) - EPS * SIGMA * A * (T[i] ** 4 - T_amb ** 4)
            J[i, i] += -h * A - 4 * EPS * SIGMA * A * T[i] ** 3
        for a, b in NEIGHBOURS:
            cond = (G_STRIP + g_gap) * (T[a] - T[b])
            rad = EPS * SIGMA * F_VIEW * A_FACE * (T[a] ** 4 - T[b] ** 4)
            F[a] -= cond + rad; F[b] += cond + rad
            dca = (G_STRIP + g_gap) + 4 * EPS * SIGMA * F_VIEW * A_FACE * T[a] ** 3
            dcb = (G_STRIP + g_gap) + 4 * EPS * SIGMA * F_VIEW * A_FACE * T[b] ** 3
            J[a, a] -= dca; J[a, b] += dcb; J[b, a] += dca; J[b, b] -= dcb
        dT = np.linalg.solve(J, -F)
        T += dT
        if np.abs(dT).max() < 1e-9: break
    balance = 8 * q_cell - sum(h * exposed_area(i) * (T[i] - T_amb) + EPS * SIGMA * exposed_area(i) * (T[i] ** 4 - T_amb ** 4) for i in range(8))
    return T, balance


def check_filled():
    missing = [k for k, v in (("R_DC", R_DC), ("CAP_AH", CAP_AH), ("G_STRIP", G_STRIP)) if v is None]
    missing += [f"LAW_A.{k}" for k, v in LAW_A.items() if v is None and k not in ("sigma_Ea", "b")]
    if missing:
        raise SystemExit("FILL these from the named sources before running: " + ", ".join(missing))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default="results"); a = ap.parse_args()
    check_filled(); os.makedirs(a.out, exist_ok=True)
    if LAW_A["b"] is None: LAW_A["b"] = 0.0
    # Step 1
    levers = {"law_A": [lever(loss_A, LAW_A, c) for c in CRATES]}
    if LAW_B and LAW_B.get("name"):
        levers["law_B"] = [lever(loss_B, LAW_B, c) for c in CRATES]
    # Step 2(a)
    rows = []
    for c in CRATES:
        q = (c * CAP_AH) ** 2 * R_DC      # W per cell at the rate's current
        for h in H_LIST:
            T, bal = module_steady(q, h)
            rows.append({"crate": c, "q_cell_W": q, "h": h, "T_max_C": float(T.max() - 273.15), "T_min_C": float(T.min() - 273.15),
                         "dT_module_K": float(T.max() - T.min()), "i_hot": int(T.argmax()), "i_cold": int(T.argmin()), "balance_W": float(bal)})
    with open(os.path.join(a.out, "module_spread.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    with open(os.path.join(a.out, "lever.json"), "w") as f:
        json.dump({"floors": FLOORS, "laws": {"A": LAW_A, "B": LAW_B}, "levers": levers, "module": rows}, f, indent=1)
    print("Step 1 (lever) and Step 2(a) (lumped module) written to", a.out, "; Step 2(b) uses the two-cell study's resolved solver")
