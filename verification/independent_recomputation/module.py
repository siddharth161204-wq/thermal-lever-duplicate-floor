"""Lumped module networks (D9, D10, D11, D12, D19, D20), written from DECISIONS.md.

2x4 in air: two rows of four at 23.1 mm pitch, horizontal and vertical neighbours only (10 pairs). Per pair: one
nickel strip 90 W/mK x 0.2 mm x 8 mm / 23.1 mm, air-gap conduction 0.03 W/mK x (A_lat/4) / 2 mm, radiation
eps^2 sigma F A_lat (Ti^4 - Tj^4). Per cell: convection h (A_outer - n A_lat/4), radiation to the room
eps sigma (A_outer - n eps F A_lat)(T^4 - T_amb^4), n = number of neighbours; 35 C air; eps 0.8.
1x8 on a cold plate (D12): one row of eight, base to a 35 C plate through 600 W/(m2 K) x A_end; the rest of each
cell (outer area minus the base) in air with h and radiation as in D10.
Node index in the 2x4: i = 4 * row + column (row 0 = "row 1"), the indexing the published T_cells columns follow.
"""
import math

import numpy as np

from common import (A_END, A_LAT, A_OUTER, EPS, F_VIEW, G_AIRGAP, G_STRIP, H_BASE, H_LIST, RATES, RATE_NAMES,
                    SIGMA, T_AMB, q_cell, q_cycle_avg, save, R_DC)


def layout_2x4():
    pairs = []
    for r in range(2):
        for c in range(3):
            pairs.append((4 * r + c, 4 * r + c + 1))
    for c in range(4):
        pairs.append((c, 4 + c))
    return 8, pairs


def layout_1x8():
    return 8, [(i, i + 1) for i in range(7)]


def build(n, pairs, h, g_strip=G_STRIP, coldplate=False, as_received_radiation=False, f_view=F_VIEW):
    nb = np.zeros(n, dtype=int)
    for a, b in pairs:
        nb[a] += 1
        nb[b] += 1
    a_conv = A_OUTER - nb * A_LAT / 4
    if as_received_radiation == "eps2":
        a_rad = A_OUTER - nb * A_LAT / 4          # my first reading: eps^2, quarter area
        c_pair = EPS ** 2 * SIGMA * f_view * A_LAT / 4
    elif as_received_radiation:
        a_rad = A_OUTER - nb * A_LAT / 4          # received script (provenance/): exposed area for room radiation
        c_pair = EPS * SIGMA * f_view * A_LAT / 4  # and eps (not eps^2) times F times a quarter of the lateral area
    else:
        a_rad = A_OUTER - nb * EPS * f_view * A_LAT
        c_pair = EPS ** 2 * SIGMA * f_view * A_LAT
    g_base = np.zeros(n)
    if coldplate:
        a_conv = a_conv - A_END
        a_rad = a_rad - A_END
        g_base[:] = H_BASE * A_END
    return {"n": n, "pairs": pairs, "nb": nb, "g_conv": h * a_conv, "c_room": EPS * SIGMA * a_rad,
            "g_lin": g_strip + G_AIRGAP, "c_pair": c_pair, "g_base": g_base}


def residual_and_jacobian(T, net, q):
    n = net["n"]
    F = net["g_conv"] * (T - T_AMB) + net["c_room"] * (T ** 4 - T_AMB ** 4) + net["g_base"] * (T - T_AMB) - q
    J = np.diag(net["g_conv"] + 4 * net["c_room"] * T ** 3 + net["g_base"])
    for a, b in net["pairs"]:
        flow = net["g_lin"] * (T[a] - T[b]) + net["c_pair"] * (T[a] ** 4 - T[b] ** 4)
        F[a] += flow
        F[b] -= flow
        da = net["g_lin"] + 4 * net["c_pair"] * T[a] ** 3
        db = net["g_lin"] + 4 * net["c_pair"] * T[b] ** 3
        J[a, a] += da
        J[a, b] -= db
        J[b, b] += db
        J[b, a] -= da
    return F, J


def solve(net, q):
    q = np.broadcast_to(np.asarray(q, dtype=float), (net["n"],)).copy()
    T = np.full(net["n"], T_AMB + 1.0)
    for it in range(1, 60):
        F, J = residual_and_jacobian(T, net, q)
        dT = np.linalg.solve(J, -F)
        T += dT
        if np.max(np.abs(dT)) < 1e-12:
            break
    losses = (net["g_conv"] * (T - T_AMB) + net["c_room"] * (T ** 4 - T_AMB ** 4) + net["g_base"] * (T - T_AMB))
    balance = q.sum() - losses.sum()
    return T, it, balance


def pos_2x4(i):
    r, c = divmod(i, 4)
    kind = "corner" if c in (0, 3) else "inner"
    return f"{kind} (row {r + 1}, column {c + 1})"


def pos_1x8(i):
    kind = "end" if i in (0, 7) else "inner"
    return f"{kind} (position {i + 1} of 8)"


def summarise(T, it, balance, q, h, c, arrangement, posf):
    Tc = T - 273.15
    i_hot, i_cold = int(np.argmax(Tc)), int(np.argmin(Tc))
    return {"arrangement": arrangement, "crate": c, "q_cell_W": float(q), "h": h,
            "T_max_C": float(Tc.max()), "T_min_C": float(Tc.min()), "T_mean_C": float(Tc.mean()),
            "dT_module_K": float(Tc.max() - Tc.min()), "dT_hot_minus_mean_K": float(Tc.max() - Tc.mean()),
            "i_hot": i_hot, "hot_position": posf(i_hot), "i_cold": i_cold, "cold_position": posf(i_cold),
            "balance_W": float(balance), "newton_iterations": it, "T_cells_C": Tc.tolist()}


def two_node(h, q, linearised=False, T_lin=None):
    """Corner/inner reduction of the 2x4 (exact by symmetry when not linearised)."""
    n, pairs = layout_2x4()
    net = build(n, pairs, h)
    # corner: n = 2, inner: n = 3; corner-inner coupling one pair each; corner-corner and inner-inner pairs carry no flow
    gc_conv, gi_conv = net["g_conv"][0], net["g_conv"][1]
    cc_room, ci_room = net["c_room"][0], net["c_room"][1]
    if linearised:
        Tl = T_lin if T_lin is not None else T_AMB
        Gc = gc_conv + 4 * cc_room * Tl ** 3
        Gi = gi_conv + 4 * ci_room * Tl ** 3
        Gp = net["g_lin"] + 4 * net["c_pair"] * Tl ** 3
        d = (Gc - Gi) / (Gi + 2 * Gp)
        th_c = 2 * q / (Gi * (1 + d) + Gc)
        return {"G_c": Gc, "G_i": Gi, "G_p": Gp, "d": d, "theta_corner": th_c, "spread": d * th_c}
    T = np.array([T_AMB + 1.0, T_AMB + 1.0])
    for _ in range(60):
        Tc_, Ti_ = T
        flow = net["g_lin"] * (Ti_ - Tc_) + net["c_pair"] * (Ti_ ** 4 - Tc_ ** 4)
        F = np.array([gc_conv * (Tc_ - T_AMB) + cc_room * (Tc_ ** 4 - T_AMB ** 4) - flow - q,
                      gi_conv * (Ti_ - T_AMB) + ci_room * (Ti_ ** 4 - T_AMB ** 4) + flow - q])
        dflow_dTi = net["g_lin"] + 4 * net["c_pair"] * Ti_ ** 3
        dflow_dTc = -(net["g_lin"] + 4 * net["c_pair"] * Tc_ ** 3)
        J = np.array([[gc_conv + 4 * cc_room * Tc_ ** 3 - dflow_dTc, -dflow_dTi],
                      [dflow_dTc, gi_conv + 4 * ci_room * Ti_ ** 3 + dflow_dTi]])
        dT = np.linalg.solve(J, -F)
        T += dT
        if np.max(np.abs(dT)) < 1e-13:
            break
    return {"T_corner_C": T[0] - 273.15, "T_inner_C": T[1] - 273.15, "spread": T[1] - T[0]}


def main():
    out = {"constants": {"A_lat": A_LAT, "A_end": A_END, "A_outer": A_OUTER, "F_view": F_VIEW, "G_strip": G_STRIP,
                         "G_airgap": G_AIRGAP, "pair_rad_linearised_at_35C": 4 * EPS ** 2 * SIGMA * F_VIEW * A_LAT * T_AMB ** 3,
                         "R_DC": R_DC},
           "module": [], "module_coldplate": [], "two_node_C2": {}, "D19": [], "sensitivity_C2_h10": {}}
    n, pairs = layout_2x4()
    n8, pairs8 = layout_1x8()
    for c in RATES:
        for h in H_LIST:
            q = q_cell(c)
            T, it, bal = solve(build(n, pairs, h), q)
            out["module"].append(summarise(T, it, bal, q, h, c, "2x4 in air", pos_2x4))
            T, it, bal = solve(build(n8, pairs8, h, coldplate=True), q)
            out["module_coldplate"].append(summarise(T, it, bal, q, h, c, "1x8 on cold plate", pos_1x8))
    # two-node reduction at C/2
    for h in H_LIST:
        q = q_cell(0.5)
        exact = two_node(h, q)
        lin_amb = two_node(h, q, linearised=True, T_lin=T_AMB)
        T_mean_K = 273.15 + 0.5 * (exact["T_corner_C"] + exact["T_inner_C"])
        lin_mean = two_node(h, q, linearised=True, T_lin=T_mean_K)
        out["two_node_C2"][str(h)] = {"nonlinear_two_node": exact, "linearised_at_35C": lin_amb,
                                      "linearised_at_mean": lin_mean}
    # D19 cycle-averaged heat, dataset rates
    for c, name in zip(RATES[:4], RATE_NAMES[:4]):
        for h in H_LIST:
            q = q_cycle_avg(c)
            T, it, bal = solve(build(n, pairs, h), q)
            row = summarise(T, it, bal, q, h, c, "2x4 in air", pos_2x4)
            row.update(rate=name, ratio_to_brief_heat=q / q_cell(c))
            out["D19"].append(row)
    # sensitivities at C/2, h = 10
    q = q_cell(0.5)
    sens = {}
    for label, g in (("strip_0.00_mW_per_K", 0.0), ("strip_6.23_mW_per_K", G_STRIP), ("strip_12.47_mW_per_K", 2 * G_STRIP)):
        T, _, _ = solve(build(n, pairs, 10.0, g_strip=g), q)
        sens[label] = float(T.max() - T.min())
    T, _, _ = solve(build(n, pairs, 10.0, as_received_radiation=True), q)
    sens["radiation_as_received"] = float(T.max() - T.min())
    T, _, _ = solve(build(n, pairs, 10.0, as_received_radiation="eps2"), q)
    sens["radiation_quarter_area_eps2_variant"] = float(T.max() - T.min())
    for label, R in (("R_10s_pulse_15.47_mOhm", 0.015468155453132936), ("R_30s_pulse_20.37_mOhm", 0.020366730495318563)):
        T, _, _ = solve(build(n, pairs, 10.0), q * R / R_DC)
        sens[label] = float(T.max() - T.min())
    T, _, _ = solve(build(n, pairs, 10.0, f_view=0.160), q)
    sens["F_view_rounded_0.160"] = float(T.max() - T.min())
    out["sensitivity_C2_h10"] = sens
    # V2-style checks: strip at 1e3 W/K, 2C
    v2 = []
    for h in H_LIST:
        T, it, bal = solve(build(n, pairs, h, g_strip=1e3), q_cell(2.0))
        v2.append({"h": h, "max_minus_min_K": float(T.max() - T.min()), "balance_rel": float(bal / (8 * q_cell(2.0)))})
    out["V2_strip_1e3"] = v2
    save("module.json", out)
    for r in out["module"]:
        print(f"2x4 c={r['crate']:.4f} h={r['h']:4.0f} q={r['q_cell_W']:.6f} dT={r['dT_module_K']:.9f} Tmax={r['T_max_C']:.6f} "
              f"Tmin={r['T_min_C']:.6f} Tmean={r['T_mean_C']:.6f} hot={r['hot_position']} cold={r['cold_position']} bal={r['balance_W']:.1e}")
    for r in out["module_coldplate"]:
        print(f"1x8 c={r['crate']:.4f} h={r['h']:4.0f} dT={r['dT_module_K']:.9e} Tmax={r['T_max_C']:.6f} Tmean={r['T_mean_C']:.6f} "
              f"hot={r['hot_position']} cold={r['cold_position']} bal={r['balance_W']:.1e}")
    print("two-node", out["two_node_C2"])
    for r in out["D19"]:
        print("D19", r["rate"], r["h"], r["q_cell_W"], r["ratio_to_brief_heat"], r["dT_module_K"])
    print("sens", sens)
    print("V2", v2)
    print("constants", out["constants"])


if __name__ == "__main__":
    main()
