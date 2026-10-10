"""Half-cylinder tetrahedral mesh of cell 2 for the resolved model (D9, D12).

Cell 2 is centred on the z axis, 0 <= z <= L (z = 0 is the negative end, z = L the positive cap).
Cell 1 lies along +x, so the plane y = 0 through both axes is a plane of symmetry and only y >= 0 is meshed.
The cross-section is triangulated ring by ring (rings graded toward the can wall, with rings exactly at the
positive-disc radius 4 mm and the negative-annulus inner radius 7.5 mm, and every ring carrying a node at
45 degrees so the gap-facing sector edge lies on mesh edges). The triangles are extruded through graded
axial layers into prisms, and each prism is split into three tetrahedra with the diagonal rule that makes
neighbouring prisms conform.
"""
import math
import numpy as np
from skfem import MeshTet

R_DISC = 0.004      # m, positive strip contact: disc of 8 mm diameter on the cap (brief)
R_ANN_IN = 0.0075   # m, negative strip contact: annulus from 15 mm diameter to the can rim (brief, D7)
SECTOR = math.pi / 4  # half-angle of the gap-facing sector (brief)


def radial_rings(R, refine=1.0):
    """Ring radii in m (centre excluded). `refine` > 1 shrinks every spacing by that factor."""
    base_mm = [0.65, 1.3, 2.0, 2.65, 3.3, 4.0, 4.7, 5.4, 6.1, 6.8, 7.5,
               8.05, 8.55, 9.0, 9.4, 9.73, 10.0, 10.2, 10.35, 10.46, 10.55]
    rings = [0.0] + [x * 1e-3 for x in base_mm]
    rings[-1] = R
    if refine == 1.0:
        return np.array(rings[1:])
    out = []
    for a, b in zip(rings[:-1], rings[1:]):
        n = max(1, int(math.ceil(refine - 1e-9)))
        for k in range(1, n + 1):
            out.append(a + (b - a) * k / n)
    return np.array(out)


def n_segments(r, R, refine=1.0):
    """Segments over [0, pi] on a ring of radius r; always a multiple of 4 so 45 and 90 degrees are nodes."""
    ds = 0.5e-3 / refine if r >= 0.0089 else 0.55e-3 / refine
    n = 4 * int(math.ceil(math.pi * r / (4 * ds)))
    if r >= 0.0089:
        n = 4 * int(math.ceil(math.pi * R / (4 * 0.5e-3 / refine)))   # aligned outer rings
    return max(4, n)


def axial_layers(Lz, refine=1.0, first=0.2e-3, growth=1.3, dmax=2.5e-3):
    first /= refine
    dmax /= refine
    growth = growth ** (1.0 / refine)
    half = Lz / 2
    z = [0.0]
    d = first
    while z[-1] + d < half - 0.5 * min(d, dmax):
        z.append(z[-1] + d)
        d = min(d * growth, dmax)
    # stretch the graded half so it ends exactly at the mid-plane
    z = np.array(z + [half])
    lower = z
    upper = Lz - lower[::-1][1:]
    return np.concatenate([lower, upper])


def build_mesh(R, Lz, refine=1.0):
    rings = radial_rings(R, refine)
    pts2 = [(0.0, 0.0)]
    ring_nodes = []
    for r in rings:
        n = n_segments(r, R, refine)
        th = np.linspace(0.0, math.pi, n + 1)
        idx = list(range(len(pts2), len(pts2) + n + 1))
        pts2 += [(r * math.cos(t), r * math.sin(t)) for t in th]
        ring_nodes.append((idx, th))
    pts2 = np.array(pts2)
    tris = []
    # centre fan
    idx0, _ = ring_nodes[0]
    for j in range(len(idx0) - 1):
        tris.append((0, idx0[j], idx0[j + 1]))
    # zipper between consecutive rings
    for (ia, ta), (ib, tb) in zip(ring_nodes[:-1], ring_nodes[1:]):
        i = j = 0
        m, n = len(ia) - 1, len(ib) - 1
        while i < m or j < n:
            if i == m:
                tris.append((ia[i], ib[j], ib[j + 1])); j += 1
            elif j == n:
                tris.append((ia[i], ib[j], ia[i + 1])); i += 1
            elif ta[i + 1] < tb[j + 1] - 1e-12:
                tris.append((ia[i], ib[j], ia[i + 1])); i += 1
            else:
                tris.append((ia[i], ib[j], ib[j + 1])); j += 1
    tris = np.array(tris)
    z = axial_layers(Lz, refine)
    n2, nz = len(pts2), len(z)
    p = np.zeros((3, n2 * nz))
    for k in range(nz):
        p[0, k * n2:(k + 1) * n2] = pts2[:, 0]
        p[1, k * n2:(k + 1) * n2] = pts2[:, 1]
        p[2, k * n2:(k + 1) * n2] = z[k]
    tets = []
    for tri in tris:
        a, b, c = sorted(tri)
        for k in range(nz - 1):
            A0, B0, C0 = a + k * n2, b + k * n2, c + k * n2
            A1, B1, C1 = A0 + n2, B0 + n2, C0 + n2
            tets += [(A0, B0, C0, C1), (A0, B0, B1, C1), (A0, A1, B1, C1)]
    t = np.array(tets).T
    # orient positively
    P = p[:, t]
    vol = np.einsum("ij,ij->j", np.cross((P[:, 1] - P[:, 0]).T, (P[:, 2] - P[:, 0]).T).T, P[:, 3] - P[:, 0])
    neg = vol < 0
    t[[1, 2]] = np.where(neg, t[[2, 1]], t[[1, 2]])
    mesh = MeshTet(p, t)
    info = {"n_rings": len(rings), "n_layers": nz, "n_nodes_2d": n2, "n_tris_2d": len(tris),
            "min_abs_vol": float(np.abs(vol).min())}
    return mesh, info


def classify_facets(mesh, R, Lz):
    """Boundary facet groups of the half model (indices into mesh.facets)."""
    bf = mesh.boundary_facets()
    mid = mesh.p[:, mesh.facets[:, bf]].mean(axis=1)
    x, y, z = mid
    r = np.hypot(x, y)
    th = np.arctan2(y, x)
    tol = 1e-9
    top = np.abs(z - Lz) < tol
    bottom = np.abs(z) < tol
    sym = (np.abs(y) < tol) & ~top & ~bottom
    lateral = ~top & ~bottom & ~sym
    groups = {
        "top": bf[top], "bottom": bf[bottom], "symmetry": bf[sym], "lateral": bf[lateral],
        "top_disc": bf[top & (r < R_DISC)],
        "bottom_annulus": bf[bottom & (r > R_ANN_IN)],
        "sector": bf[lateral & (th < SECTOR)],
        "lateral_nonsector": bf[lateral & (th >= SECTOR)],
    }
    groups["outer"] = np.concatenate([groups["top"], groups["bottom"], groups["lateral"]])
    groups["conv"] = np.concatenate([groups["top"], groups["bottom"], groups["lateral_nonsector"]])
    return groups
