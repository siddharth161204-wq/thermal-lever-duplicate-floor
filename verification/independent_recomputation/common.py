"""Shared inputs for the independent check.

Everything here is read from the study's inputs/ folder or typed from DECISIONS.md (D6 to D12).
No code from the study's src/ scripts is used. The repository is opened read only.
"""
import csv
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.environ.get("THERMAL_LEVER_REPO", os.path.join(HERE, "..", "..")))   # the repository root
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)

# ----------------------------------------------------------------------------------------------
# floors (inputs/floor.json, D1)
# ----------------------------------------------------------------------------------------------
with open(os.path.join(REPO, "inputs", "floor.json")) as f:
    _floor = json.load(f)
FLOORS = dict(_floor["FLOORS"])          # keys: mean_3.0pct, worst_of_8_4.1pct, p95_8.7pct
FLOOR_KEYS = ["mean_3.0pct", "worst_of_8_4.1pct", "p95_8.7pct"]

# ----------------------------------------------------------------------------------------------
# rates and heat (D7, D8)
# ----------------------------------------------------------------------------------------------
RATES = [1 / 16, 1 / 10, 1 / 5, 1 / 2, 1.0, 2.0]
RATE_NAMES = ["C/16", "C/10", "C/5", "C/2", "1C", "2C"]
DATASET_RATES = RATES[:4]
H_LIST = [5.0, 10.0, 25.0]
CAPACITY_AH = 5.0
with open(os.path.join(REPO, "inputs", "lgm50_resistance.json")) as f:
    _res = json.load(f)
R_DC = float(_res["registered"]["R_DC_ohm"])


def q_cell(c):
    """Brief's heat per cell: (c x 5 A)^2 x R_DC [W]."""
    return (c * CAPACITY_AH) ** 2 * R_DC


def q_cycle_avg(c):
    """D19: CC charge at C/2 for 2 h, CC discharge at c for 1/c h, no rests, I^2 R heat averaged over the cycle."""
    i1 = CAPACITY_AH
    return R_DC * i1 ** 2 * (0.5 + c) / (2.0 + 1.0 / c)


# ----------------------------------------------------------------------------------------------
# cell geometry and module couplings (D6, D9, D10, D11, D12)
# ----------------------------------------------------------------------------------------------
D = 0.0211
L = 0.0701
RAD = D / 2
PITCH = 0.0231
GAP = PITCH - D                      # 2.0 mm
A_LAT = math.pi * D * L
A_END = math.pi * RAD ** 2
A_OUTER = A_LAT + 2 * A_END          # both ends included
T_AMB = 35.0 + 273.15
EPS = 0.8
SIGMA = 5.670374419e-8               # CODATA 2018 (the two-cell study used 5.670e-8; effect < 1e-4 relative)

K_NI, T_STRIP, W_STRIP = 90.0, 0.2e-3, 8e-3
G_STRIP = K_NI * T_STRIP * W_STRIP / PITCH          # one strip per pair (D9)
K_AIR = 0.03
A_FACE = A_LAT / 4
G_AIRGAP = K_AIR * A_FACE / GAP                     # per pair (D9)
_X = 1.0 + GAP / D
F_VIEW = (math.sqrt(_X * _X - 1.0) + math.asin(1.0 / _X) - _X) / math.pi   # parallel equal cylinders, per unit of the whole lateral area (D10)
H_BASE = 600.0                                      # D12

# resolved cell (D6, D13; geometry of the two-cell study's mesh: the whole cylinder is jelly roll,
# a zero-thickness steel shell with in-plane conduction k t on the lateral surface and both ends)
K_R = 1.3
K_Z = 25.0
K_CAN = 15.0
T_CAN = 0.25e-3


def rel(a, b):
    """Relative difference (mine - repo) / repo."""
    if b == 0:
        return float("nan") if a != 0 else 0.0
    return (a - b) / b


def load_json(path):
    with open(path) as f:
        return json.load(f)


def repo_json(name):
    return load_json(os.path.join(REPO, "results", name))


def repo_verif(name):
    return load_json(os.path.join(REPO, "verification", name))


def repo_csv(name):
    with open(os.path.join(REPO, "results", name)) as f:
        return list(csv.DictReader(f))


def save(name, obj):
    path = os.path.join(OUT, name)
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, default=float)
    return path
