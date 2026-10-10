"""Resolved cell 2: three-dimensional transient conduction in scikit-fem, driven by cell 1's temperature history.

Every physical constant is imported from lumped_two_cell.py, so the resolved and lumped models share their
inputs by construction. See DECISIONS.md D5 to D16 for the choices made here.

Field equation (jelly roll, half cylinder):  rho_cp dT/dt = div(K grad T),  K = k_ax I + (k_rad - k_ax) e_r e_r^T
Can: 0.25 mm steel shell on every outer surface, in-plane conduction k_can t_can and heat capacity 4.5 J/K.
Boundary terms on the outer surface (per unit area, T1 from the lumped run, one-way coupling):
  positive strip, top disc r < 4 mm          : h_pos (T1 - T),  h_pos = G_pos / A_disc
  negative strip, bottom annulus r > 7.5 mm  : h_neg (T1 - T),  h_neg = G_neg / A_annulus
  gap-facing sector |theta| < 45 deg         : (k_air / gap) (T1 - T)
  lateral surface, radiation from cell 1     : eps^2 sigma F_loc(theta) (T1^4 - T_amb^4)
  room convection (outer minus sector)       : h (T_amb - T)
  room radiation (whole outer surface)       : eps sigma (T_amb^4 - T^4)
Symmetry plane y = 0: zero flux.
Time integration: backward Euler, row-sum (lumped) capacity, room radiation lagged one step (D13).
"""
import math
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from skfem import Basis, FacetBasis, ElementTetP1, BilinearForm, LinearForm, asm

import lumped_two_cell as LP
from mesh_cell2 import build_mesh, classify_facets, R_DISC, R_ANN_IN, SECTOR

R = LP.D_CELL / 2
LZ = LP.L_CELL
CENTRE_DIST = 2 * R + LP.GAP


def f_local(theta):
    """View factor from a strip of cell 2's lateral surface at angle theta (from the line of centres) to cell 1.

    Two-dimensional (infinitely long parallel cylinders): F = (sin b - sin a) / 2, where [a, b] is the range of
    directions, measured from the strip normal, in which cell 1 is seen, clipped to the strip's horizon.
    """
    th = np.asarray(theta, dtype=float)
    px, py = R * np.cos(th), R * np.sin(th)
    vx, vy = CENTRE_DIST - px, -py
    d = np.hypot(vx, vy)
    alpha = np.arcsin(np.clip(R / d, -1, 1))
    cos_psi = (np.cos(th) * vx + np.sin(th) * vy) / d
    psi = np.arccos(np.clip(cos_psi, -1, 1))
    lo = np.maximum(psi - alpha, -math.pi / 2)
    hi = np.minimum(psi + alpha, math.pi / 2)
    return np.where(hi > lo, 0.5 * (np.sin(hi) - np.sin(lo)), 0.0)


class _Solver:
    """Sparse symmetric positive definite solve: MKL PARDISO when available (METIS ordering), SuperLU otherwise."""

    def __init__(self, A):
        try:
            import pypardiso
            self.Au = sp.triu(A, format="csr")
            self.ps = pypardiso.PyPardisoSolver(mtype=2)
            self.ps.factorize(self.Au)
            self.kind = "pardiso"
        except ImportError:
            self.lu = spla.splu(A.tocsc(), permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0,
                                options=dict(SymmetricMode=True))
            self.kind = "superlu"

    def solve(self, b):
        if self.kind == "pardiso":
            return self.ps.solve(self.Au, b)
        return self.lu.solve(b)


class Cell2Model:
    def __init__(self, refine=1.0):
        self.mesh, self.mesh_info = build_mesh(R, LZ, refine)
        m = self.mesh
        e = ElementTetP1()
        self.vb = Basis(m, e)
        self.groups = classify_facets(m, R, LZ)
        fb = {k: FacetBasis(m, e, facets=v) for k, v in self.groups.items() if k != "symmetry" and len(v)}
        kr, kz = LP.K_RADIAL, LP.K_AXIAL

        @BilinearForm
        def conduction(u, v, w):
            x, y = w.x[0], w.x[1]
            r = np.sqrt(x ** 2 + y ** 2) + 1e-15
            ex, ey = x / r, y / r
            gu, gv = u.grad, v.grad
            ur = ex * gu[0] + ey * gu[1]
            vr = ex * gv[0] + ey * gv[1]
            return kz * (gu[0] * gv[0] + gu[1] * gv[1] + gu[2] * gv[2]) + (kr - kz) * ur * vr

        kt = LP.K_CAN * LP.T_CAN

        @BilinearForm
        def shell(u, v, w):
            n = w.n
            gu, gv = u.grad, v.grad
            un = gu[0] * n[0] + gu[1] * n[1] + gu[2] * n[2]
            vn = gv[0] * n[0] + gv[1] * n[1] + gv[2] * n[2]
            return kt * sum((gu[i] - un * n[i]) * (gv[i] - vn * n[i]) for i in range(3))

        @LinearForm
        def one(v, w):
            return v

        @LinearForm
        def view(v, w):
            th = np.arctan2(w.x[1], w.x[0])
            return f_local(th) * v

        self.K_roll = asm(conduction, self.vb).tocsr()
        self.K_shell = asm(shell, fb["outer"]).tocsr()
        self.K = (self.K_roll + self.K_shell).tocsr()
        w_vol = asm(one, self.vb)
        self.w_vol = w_vol
        self.w = {k: asm(one, b) for k, b in fb.items()}
        g_view = asm(view, fb["lateral"])

        # geometric fidelity of the mesh (half model, doubled to the full cell)
        full = lambda k: 2 * self.w[k].sum()
        self.geom = {
            "volume_mesh": 2 * w_vol.sum(), "volume_exact": math.pi * R ** 2 * LZ,
            "A_outer_mesh": full("outer"), "A_outer_exact": LP.A_CELL,
            "A_face_mesh": full("sector"), "A_face_exact": LP.A_FACE,
            "A_conv_mesh": full("conv"), "A_conv_exact": LP.A_CONV,
            "A_disc_mesh": full("top_disc"), "A_disc_exact": math.pi * R_DISC ** 2,
            "A_annulus_mesh": full("bottom_annulus"), "A_annulus_exact": math.pi * (R ** 2 - R_ANN_IN ** 2),
            "FA_lat_mesh": 2 * g_view.sum(), "FA_lat_exact": LP.F_VIEW * LP.A_LAT,
        }
        # capacities (D10): shell 4.5 J/K over the outer area, roll the remainder, total C_CELL
        c_shell_total = LP.M_CAN * LP.CP_CAN
        c_shell = c_shell_total * self.w["outer"] / full("outer")
        rho_cp_roll = (LP.C_CELL - c_shell_total) / self.geom["volume_mesh"]
        self.rho_cp_roll = rho_cp_roll
        self.C = rho_cp_roll * w_vol + c_shell              # nodal capacities, sum = C_CELL / 2
        self.C_total = self.C.sum()
        # boundary weights normalised so that a uniform field reproduces the lumped totals exactly (V0)
        self.d_conv_unit = LP.H_CONV * LP.A_CONV * self.w["conv"] / full("conv")
        self.d_gap_unit = LP.G_GAP * self.w["sector"] / full("sector")
        self.d_pos_unit = self.w["top_disc"] / full("top_disc")          # times G_pos
        self.d_neg_unit = self.w["bottom_annulus"] / full("bottom_annulus")  # times G_neg
        self.g_rad = LP.EPS ** 2 * LP.SIGMA * LP.F_VIEW * LP.A_LAT * g_view / (2 * g_view.sum())
        self.w_rad = LP.EPS * LP.SIGMA * LP.A_CELL * self.w["outer"] / full("outer")
        p = self.mesh.p
        self.r_node = np.hypot(p[0], p[1])
        self.th_node = np.arctan2(p[1], p[0])
        self.z_node = p[2]
        self._fact = {}

    # ---------------------------------------------------------------------------------------------
    def location(self, i):
        """Geometric description of node i (D14)."""
        r, th, z = self.r_node[i], self.th_node[i], self.z_node[i]
        facing = r >= R - 2e-3 and th <= SECTOR + math.radians(15)
        if z >= LZ - 3e-3 and r <= R_DISC + 3e-3:
            return "positive terminal end"
        if z <= 3e-3 and r >= R_ANN_IN - 3e-3:
            return "gap-facing rim of the negative end" if facing else "negative terminal end"
        if facing:
            return "gap-facing surface"
        return "elsewhere"

    def driver(self, i, T, T1, g_pos, g_neg, gap_on=True, rad_gain_on=True):
        """Which path heats node i directly (D14): the strips, or the air gap plus radiation from cell 1."""
        q_strip = (g_pos * self.d_pos_unit[i] + g_neg * self.d_neg_unit[i]) * (T1 - T[i])
        q_gap = (self.d_gap_unit[i] * (T1 - T[i]) if gap_on else 0.0) + \
                (self.g_rad[i] * (T1 ** 4 - LP.T_AMB ** 4) if rad_gain_on else 0.0)
        if q_strip <= 0 and q_gap <= 0:
            return "elsewhere"
        return "terminal end" if q_strip > q_gap else "gap-facing surface"

    def system(self, dt, g_pos, g_neg, gap_on=True, conv_on=True):
        key = (round(dt, 9), g_pos, g_neg, gap_on, conv_on)
        if key not in self._fact:
            d = self.C / dt
            if conv_on:
                d = d + self.d_conv_unit
            if gap_on:
                d = d + self.d_gap_unit
            d = d + g_pos * self.d_pos_unit + g_neg * self.d_neg_unit
            A = (self.K + sp.diags(d)).tocsr()
            self._fact = {}                                   # keep one factorisation in memory
            self._fact[key] = _Solver(A)
        return self._fact[key]

    def run(self, t1_of_t, t_end=900.0, dt=0.1, g_pos=0.0, g_neg=0.0, gap_on=True, rad_gain_on=True,
            conv_on=True, rad_loss_on=True, q_vol=None, q_shape="capacity", levels=None, record_every=0.5, T_init=None,
            volume_levels=None, cell1=None):
        """Integrate cell 2. t1_of_t(t) returns T1 in kelvin. q_vol(t) is an optional uniform internal source in W
        (full cell), used only by the verification runs. Returns a dict of crossings, peaks and a time series."""
        Ta = LP.T_AMB
        levels = levels if levels is not None else LP.LEVELS
        lu = self.system(dt, g_pos, g_neg, gap_on, conv_on)
        T = np.full(self.C.shape, Ta if T_init is None else T_init)
        n_steps = int(round(t_end / dt))
        every = max(1, int(round(record_every / dt)))
        Ctot = self.C_total
        cross_max = {k: None for k in levels}
        cross_mean = {k: None for k in levels}
        extra = {k: None for k in levels}
        Tmax_prev = Tmean_prev = T.max()
        peak = {"Tmax": T.max(), "t_Tmax": 0.0, "imax": None, "loc_peak": None, "geom_peak": None, "Tmean": T.max()}
        series = [(0.0, t1_of_t(0.0), T.max(), (self.C @ T) / Ctot, int(np.argmax(T)))]
        energy = {"in_bus": 0.0, "in_gap": 0.0, "in_rad": 0.0, "out_conv": 0.0, "out_rad": 0.0, "in_vol": 0.0}
        d_bus = g_pos * self.d_pos_unit + g_neg * self.d_neg_unit
        d_gap = self.d_gap_unit if gap_on else 0.0 * self.d_gap_unit
        d_conv = self.d_conv_unit if conv_on else 0.0 * self.d_conv_unit
        E0 = self.C @ T
        # two-way coupling (D25, exploratory): cell 1 becomes a lumped state driven by the resolved cell 2's heat draw
        if cell1 is not None:
            T1_state = Ta
            e1, tau1 = cell1["e_body"], cell1["tau"]

            def released(t):            # exact integral of the triangular pulse
                if t <= 0:
                    return 0.0
                if t >= tau1:
                    return e1
                if t <= tau1 / 2:
                    return e1 * 2 * (t / tau1) ** 2
                return e1 * (1 - 2 * ((tau1 - t) / tau1) ** 2)
            series[0] = (0.0, Ta, T.max(), (self.C @ T) / Ctot, int(np.argmax(T)))
        vol_total = self.w_vol.sum()
        vol_max = {lv: 0.0 for lv in (volume_levels or [])}
        vol_series = []
        for n in range(n_steps):
            t_new = (n + 1) * dt
            if cell1 is None:
                T1 = t1_of_t(t_new)
            else:
                q_out = 2 * float((d_bus + d_gap) @ (T1_state - T))                       # strips and air gap, to cell 2
                q_back = 2 * float(self.g_rad @ (T ** 4 - Ta ** 4))                       # radiation from cell 2, reciprocity
                q_room = LP.H_CONV * LP.A_CONV * (T1_state - Ta) + LP.EPS * LP.SIGMA * LP.A_CELL * (T1_state ** 4 - Ta ** 4)
                q_rel = (released(t_new) - released(t_new - dt)) / dt
                T1_state = T1_state + dt * (q_rel - q_out + q_back - q_room) / LP.C_CELL
                T1 = T1_state
            rhs = self.C / dt * T + d_conv * Ta + (d_gap + d_bus) * T1
            q_rad_gain = (T1 ** 4 - Ta ** 4) if rad_gain_on else 0.0
            rhs += self.g_rad * q_rad_gain
            rad_loss = self.w_rad * (T ** 4 - Ta ** 4) if rad_loss_on else 0.0 * T
            rhs -= rad_loss
            qv = 0.0
            if q_vol is not None:
                qv = q_vol(t_new) / 2.0                                     # half model
                shape = self.C / Ctot if q_shape == "capacity" else self.w_vol / self.w_vol.sum()
                rhs += qv * shape                                          # uniform per unit capacity or per unit roll volume
            Tn = lu.solve(rhs)
            # energy account (half model, doubled at the end)
            energy["in_bus"] += dt * float(d_bus @ (T1 - Tn))
            energy["in_gap"] += dt * float(d_gap @ (T1 - Tn))
            energy["in_rad"] += dt * float(self.g_rad.sum() * q_rad_gain)
            energy["out_conv"] += dt * float(d_conv @ (Tn - Ta))
            energy["out_rad"] += dt * float(np.sum(rad_loss))
            energy["in_vol"] += dt * qv
            T = Tn
            imax = int(np.argmax(T))
            Tmax = T[imax]
            Tmean = (self.C @ T) / Ctot
            for k, lv in levels.items():
                if cross_max[k] is None and Tmax >= lv:
                    tc = t_new - dt + dt * (lv - Tmax_prev) / (Tmax - Tmax_prev) if Tmax > Tmax_prev else t_new
                    cross_max[k] = tc
                    extra[k] = {"loc": self.driver(imax, T, T1, g_pos, g_neg, gap_on, rad_gain_on),
                                "geom": self.location(imax), "r_mm": 1e3 * self.r_node[imax],
                                "theta_deg": math.degrees(self.th_node[imax]), "z_mm": 1e3 * self.z_node[imax],
                                "Tmax_minus_Tmean_K": Tmax - Tmean}
                if cross_mean[k] is None and Tmean >= lv:
                    cross_mean[k] = t_new - dt + dt * (lv - Tmean_prev) / (Tmean - Tmean_prev) if Tmean > Tmean_prev else t_new
            if Tmax > peak["Tmax"]:
                peak.update(Tmax=Tmax, t_Tmax=t_new, imax=imax,
                            loc_peak=self.driver(imax, T, T1, g_pos, g_neg, gap_on, rad_gain_on),
                            geom_peak=self.location(imax))
            peak["Tmean"] = max(peak["Tmean"], Tmean)
            Tmax_prev, Tmean_prev = Tmax, Tmean
            if volume_levels:
                fr = {lv: float(self.w_vol[T >= lv].sum() / vol_total) for lv in volume_levels}
                for lv, v in fr.items():
                    vol_max[lv] = max(vol_max[lv], v)
            if (n + 1) % every == 0:
                series.append((t_new, T1, Tmax, Tmean, imax))
                if volume_levels:
                    vol_series.append([t_new] + [fr[lv] for lv in volume_levels])
        stored = (self.C @ T - E0)
        energy = {k: 2 * v for k, v in energy.items()}
        net = energy["in_bus"] + energy["in_gap"] + energy["in_rad"] + energy["in_vol"] - energy["out_conv"] - energy["out_rad"]
        return {"vol_fraction_max": vol_max, "vol_series": np.array(vol_series) if vol_series else None,
                "cross_max": cross_max, "cross_mean": cross_mean, "at_cross": extra, "peak": peak,
                "series": np.array(series), "energy": energy, "stored_J": 2 * stored, "net_in_J": net, "T_final": T}


def lumped_one_way(t1_of_t, t_end=900.0, dt=0.1, g_bus=0.0, gap_on=True, rad_gain_on=True, conv_on=True,
                   rad_loss_on=True, q_vol=None, levels=None, T_init=None):
    """The lumped cell-2 node integrated with the SAME scheme as the resolved model (backward Euler, room radiation
    lagged), driven by the same T1(t). Used for V1 and V3 and to check the script's explicit integration."""
    Ta = LP.T_AMB
    levels = levels if levels is not None else LP.LEVELS
    C = LP.C_CELL
    G = g_bus + (LP.G_GAP if gap_on else 0.0)
    hA = LP.H_CONV * LP.A_CONV if conv_on else 0.0
    T = Ta if T_init is None else T_init
    cross = {k: None for k in levels}
    Tprev = T
    out = [(0.0, T)]
    for n in range(int(round(t_end / dt))):
        t_new = (n + 1) * dt
        T1 = t1_of_t(t_new)
        q_in = (LP.EPS ** 2 * LP.SIGMA * LP.F_VIEW * LP.A_LAT * (T1 ** 4 - Ta ** 4)) if rad_gain_on else 0.0
        q_loss_rad = LP.EPS * LP.SIGMA * LP.A_CELL * (T ** 4 - Ta ** 4) if rad_loss_on else 0.0
        qv = q_vol(t_new) if q_vol is not None else 0.0
        T = (C / dt * T + G * T1 + hA * Ta + q_in - q_loss_rad + qv) / (C / dt + G + hA)
        for k, lv in levels.items():
            if cross[k] is None and T >= lv:
                cross[k] = t_new - dt + dt * (lv - Tprev) / (T - Tprev) if T > Tprev else t_new
        Tprev = T
        out.append((t_new, T))
    return cross, np.array(out)


def t1_interpolator(path):
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    t, T1 = data[:, 0], data[:, 1] + 273.15

    def f(tt):
        return float(np.interp(tt, t, T1))
    return f, data
