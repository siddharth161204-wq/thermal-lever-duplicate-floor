#!/usr/bin/env python3
"""
Reference model: two LG M50 cells, cell 1 in thermal runaway, cell 2 the neighbour.
Both cells are LUMPED nodes here. This is the baseline that the RESOLVED cell-2 model
(src/resolved_cell2.py: anisotropic jelly roll + can shell, busbar strips entering at the
terminals, gap-facing surface receiving conduction + radiation) is compared against with
IDENTICAL inputs: same heat release in cell 1, same path conductances, same ambient losses,
same onset criteria. The resolved model imports every constant below from this file, so the
two models cannot drift apart. The only difference between them is whether cell 2 has a
temperature field or a single temperature.

Changes from the script as received on 7 Oct 2026 (provenance/lumped_two_cell_as_received.py),
each recorded in DECISIONS.md:
  D1, D2  band edges measured from Koenig, Zhao and Deng (2025) Fig. 7.
  D5      can-to-can radiation uses F over the WHOLE lateral area (the view factor of the
          parallel-cylinder formula is defined per unit of the whole lateral area), in the grey
          radiosity form with the room as a black third surface.
  D6      room convection excludes the quarter of the lateral area that faces the other cell.
  D13     T1(t) is written every 0.1 s from t = 0 on an integer step grid.
  Databank E_BODY_LIST filled from inputs/databank_m50_rows.csv (see inputs/extract_databank.py).

Paths from cell 1 to cell 2:
  busbar   : strip conductance in series with two weld joints  -> G_bus  [W/K]  (swept)
  gap      : air conduction across the 2 mm gap between the cans -> G_gap  [W/K]
  radiation: grey exchange between the facing cans, eps^2 sigma F A_lat (T1^4 - T2^4)
Losses from each cell to the room: convection h*A_conv + radiation eps*sigma*A_rad*(T^4 - Tamb^4),
with A_conv = A_outer - A_face and A_rad = A_outer - eps*F*A_lat (radiosity, first order in (1-eps)F).

Heat release in cell 1: a triangular pulse of total energy E_body (the Battery Failure
Databank's "cell body" share of the corrected total energy yield for the LG 21700-M50 at
100 % SOC) over a release duration tau_rel.

Onset criteria for cell 2 (Koenig, Zhao and Deng, Chem. Eng. J. 507 (2025) 160402,
Table A.1 main cluster and Fig. 7 two-sigma bands): T_initial (self-heating, 183.5 C) and
T_onset (internal short circuit, 278.9 C).

Usage:  python3 lumped_two_cell.py --out results     (runs the sweep, writes lumped_results.csv and t1_series/)
"""
import numpy as np, csv, math, os, argparse

# ---------------- cell (LG INR21700 M50), from O'Regan et al. 2022 via PyBaMM "ORegan2022" (D3) ----------------
#   jelly roll: density 2938.7 kg/m3; cp 842 J/kgK at 25 C, 999 at 60 C, 1110 at 100 C;
#   k radial (series) 1.23 to 1.38 W/mK; k axial and azimuthal (parallel) 24.9 to 25.1 W/mK;
#   jelly-roll mass 56.4 g (P2D volume 1.918e-5 m3) plus a steel can of about 9 g (cp 500 J/kgK).
#   Cell-average cp at 60 to 100 C is therefore about 950 J/kgK on 65 g. The component cp fits
#   are non-physical above about 150 C, so a constant value is used in both models (D3).
M_CELL   = 0.065      # kg        jelly roll 56 g + can 9 g (datasheet total 68 g)
CP_CELL  = 950.0      # J/(kg K)  cell average at 60 to 100 C
D_CELL   = 0.0211     # m
L_CELL   = 0.0701     # m
A_LAT    = math.pi * D_CELL * L_CELL                      # m^2, lateral area
A_CELL   = A_LAT + 2 * math.pi * (D_CELL / 2) ** 2       # m^2, total outer area
C_CELL   = M_CELL * CP_CELL   # J/K
M_CAN    = 0.009      # kg   steel can (resolved model: 0.25 mm shell carrying this heat capacity, D10)
CP_CAN   = 500.0      # J/(kg K)
T_CAN    = 0.25e-3    # m    can wall thickness (Finegan et al. 2024, Table 1: 250 um for the M50)
K_CAN    = 15.0       # W/(m K)
K_RADIAL = 1.3        # W/(m K)  jelly roll, through the layers (D11)
K_AXIAL  = 25.0       # W/(m K)  jelly roll, along the layers, axial and azimuthal (D11)

# ---------------- ambient ----------------
T_AMB    = 25.0 + 273.15
H_CONV   = 10.0       # W/(m^2 K) natural convection (brief)
EPS      = 0.8        # emissivity of the can surface (brief)
SIGMA    = 5.670e-8

# ---------------- gap between cans ----------------
GAP      = 0.002      # m (2 mm)
K_AIR    = 0.03       # W/(m K)
A_FACE   = 0.25 * A_LAT                                    # m^2, facing area: the 90 degree sector facing the neighbour
G_GAP    = K_AIR * A_FACE / GAP                            # W/K
X        = 1.0 + GAP / D_CELL
F_VIEW   = (math.sqrt(X * X - 1.0) + math.asin(1.0 / X) - X) / math.pi   # infinite parallel equal cylinders, per unit of the whole lateral area (D5)
A_CONV   = A_CELL - A_FACE                                 # D6: no room convection across the 2 mm gap
A_RAD    = A_CELL - EPS * F_VIEW * A_LAT                   # D5: room radiation, radiosity form

# ---------------- heat release in cell 1, from the Battery Failure Databank (Finegan et al. 2024) ----------------
E_BODY_LIST  = [8751.3, 21495.3, 33222.9]   # J, min / median / max of Energy-Fraction-Cell-Body-kJ over the 17 LG 21700-M50 (BV)
                                            # tests at 100 % SOC (pre-test OCV >= 4.1 V; 7 heater, 10 nail; no ISC-device tests
                                            # exist for this cell), Battery Failure Databank revision 2, inputs/extract_databank.py (D18, D19)
TAU_REL_LIST = [10.0, 30.0]         # s, release duration of the triangular pulse (brief; the Databank has energies, not heat-rate curves)

# ---------------- busbar path, swept over the design space ----------------
G_BUS_LIST = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]   # W/K  (wire bond < 0.005, nickel strip ~0.01, copper plate ~0.1 to 0.5)

# ---------------- onset criteria for cell 2 (Koenig, Zhao and Deng 2025), D1 and D2 ----------------
T_INITIAL_MEAN = 183.5 + 273.15
T_INITIAL_BAND = (183.5 * (1 - 0.158 / 2) + 273.15, 183.5 * (1 + 0.158 / 2) + 273.15)   # 169.0 to 198.0 C
T_ONSET_MEAN   = 278.9 + 273.15
T_ONSET_BAND   = (278.9 * (1 - 0.1665 / 2) + 273.15, 278.9 * (1 + 0.1665 / 2) + 273.15)  # 255.7 to 302.1 C
LEVELS = {"T_initial_low": T_INITIAL_BAND[0], "T_initial_mean": T_INITIAL_MEAN, "T_initial_high": T_INITIAL_BAND[1],
          "T_onset_low": T_ONSET_BAND[0], "T_onset_mean": T_ONSET_MEAN, "T_onset_high": T_ONSET_BAND[1]}

T_END     = 900.0  # s simulated
DT        = 0.01   # s, explicit Euler (as received)
DT_SERIES = 0.1    # s, T1(t) output interval for the resolved model (D13)


def q_release(t, e_body, tau):
    """Triangular pulse: rises to its peak at tau/2 and returns to zero at tau; integral = e_body."""
    if t < 0 or t > tau:
        return 0.0
    peak = 2.0 * e_body / tau
    return peak * (t / (tau / 2)) if t <= tau / 2 else peak * ((tau - t) / (tau / 2))


def loss_to_ambient(T):
    return H_CONV * A_CONV * (T - T_AMB) + EPS * SIGMA * A_RAD * (T ** 4 - T_AMB ** 4)


def coupling(T1, T2, g_bus):
    cond = (g_bus + G_GAP) * (T1 - T2)
    rad = EPS ** 2 * SIGMA * F_VIEW * A_LAT * (T1 ** 4 - T2 ** 4)
    return cond + rad


def run(e_body, tau, g_bus, series=None):
    """Explicit Euler at DT. If `series` is a list, (t, T1, T2) is appended every DT_SERIES from t = 0."""
    T1 = T2 = T_AMB
    n_steps = int(round(T_END / DT))
    every = int(round(DT_SERIES / DT))
    cross = {k: None for k in LEVELS}
    T1max = T2max = T_AMB
    if series is not None:
        series.append((0.0, T1 - 273.15, T2 - 273.15))
    for n in range(n_steps):
        t = n * DT
        q12 = coupling(T1, T2, g_bus)
        dT1 = (q_release(t, e_body, tau) - loss_to_ambient(T1) - q12) / C_CELL
        dT2 = (q12 - loss_to_ambient(T2)) / C_CELL
        T1 += dT1 * DT
        T2 += dT2 * DT
        t_new = (n + 1) * DT
        T1max = max(T1max, T1); T2max = max(T2max, T2)
        if series is not None and (n + 1) % every == 0:
            series.append((round(t_new, 3), T1 - 273.15, T2 - 273.15))
        for k, lv in LEVELS.items():
            if cross[k] is None and T2 >= lv:
                cross[k] = round(t_new, 2)
    return T1max, T2max, cross


def series_name(e, tau, g):
    return f"T1_E{int(round(e))}_tau{int(tau)}_G{g}.csv"


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results")
    args = ap.parse_args()
    if E_BODY_LIST is None:
        raise SystemExit("E_BODY_LIST is not set: fill it from the Battery Failure Databank before running.")
    os.makedirs(os.path.join(args.out, "t1_series"), exist_ok=True)
    rows = []
    for e in E_BODY_LIST:
        for tau in TAU_REL_LIST:
            for g in G_BUS_LIST:
                ser = []
                T1max, T2max, cross = run(e, tau, g, series=ser)
                # T1(t) is the boundary input for the resolved cell-2 model (one-way coupling, stated in the note)
                with open(os.path.join(args.out, "t1_series", series_name(e, tau, g)), "w", newline="") as fs:
                    ws = csv.writer(fs); ws.writerow(["t_s", "T1_C", "T2_lumped_C"])
                    ws.writerows([(f"{a:.1f}", f"{b:.4f}", f"{c:.4f}") for a, b, c in ser])
                rows.append({"E_body_J": round(e, 1), "tau_rel_s": tau, "G_bus_WK": g, "G_gap_WK": round(G_GAP, 5), "F_view": round(F_VIEW, 4),
                             "T1max_C": round(T1max - 273.15, 1), "T2max_C": round(T2max - 273.15, 1), **cross})
    with open(os.path.join(args.out, "lumped_results.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print("wrote lumped_results.csv with", len(rows), "cases; G_gap =", round(G_GAP, 4), "W/K; view factor =", round(F_VIEW, 4),
          "; A_conv =", round(A_CONV, 6), "m2; A_rad =", round(A_RAD, 6), "m2")
