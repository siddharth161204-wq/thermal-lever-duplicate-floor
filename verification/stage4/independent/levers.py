#!/usr/bin/env python3
"""The archived levers and floors of D56, read at full precision and checked against D56's printed values.

Bands (D57): warm = the twelve KIT 25 to 40 C levers + the Sandia NMC 25 to 35 C lever (13);
cold = the twelve KIT 10 to 25 C levers + the Sandia NMC, NCA and LFP 15 to 25 C levers (15).
s is per kelvin (a fraction), positive when the warmer cell reaches the threshold sooner.
"""
import csv, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))

# D56 (3) and (4), in % per K: point, interval low, interval high, at the rounding printed in D56
D56_KIT_WARM = {"0-100 0.33/0.33C": (+1.3029, +0.4122, +2.1937), "0-100 0.33/1.00C": (-0.6236, -1.5140, +0.2667),
                "0-100 1.00/1.00C": (-2.6850, -3.5876, -1.7825), "0-100 1.67/1.00C": (-3.4548, -4.4542, -2.4553),
                "10-100 0.33/0.33C": (+1.3598, +0.4592, +2.2605), "10-100 0.33/1.00C": (+0.6724, -0.2212, +1.5659),
                "10-100 1.00/1.00C": (-0.4222, -1.3097, +0.4653), "10-100 1.67/1.00C": (-1.9148, -2.8032, -1.0265),
                "10-90 0.33/0.33C": (+1.9097, +1.0043, +2.8152), "10-90 0.33/1.00C": (+1.0898, +0.1874, +1.9922),
                "10-90 1.00/1.00C": (-0.1935, -1.0698, +0.6828), "10-90 1.67/1.00C": (-0.9952, -1.8798, -0.1106)}
D56_KIT_COLD = {"0-100 0.33/0.33C": (-4.0180, -4.9813, -3.0547), "0-100 0.33/1.00C": (-3.1041, -4.0560, -2.1523),
                "0-100 1.00/1.00C": (-3.3410, -4.2788, -2.4033), "0-100 1.67/1.00C": (-3.5800, -4.6130, -2.5471),
                "10-100 0.33/0.33C": (-2.1876, -3.1272, -1.2481), "10-100 0.33/1.00C": (-3.6478, -4.5758, -2.7198),
                "10-100 1.00/1.00C": (-4.1691, -5.0866, -3.2517), "10-100 1.67/1.00C": (-11.5906, -12.5150, -10.6661),
                "10-90 0.33/0.33C": (-1.7564, -2.7006, -0.8121), "10-90 0.33/1.00C": (-2.2142, -3.1580, -1.2704),
                "10-90 1.00/1.00C": (-3.8750, -4.8438, -2.9062), "10-90 1.67/1.00C": (-5.2373, -6.1621, -4.3126)}
D56_SANDIA = {("NMC", "25 to 35 C"): (-2.4432, -4.3314, -0.4448), ("NMC", "15 to 25 C"): (-5.6227, -8.4354, -3.1004),
              ("NCA", "15 to 25 C"): (-1.3151, -5.2199, +0.7423), ("LFP", "15 to 25 C"): (+3.3718, +2.9491, +3.8270)}
D56_FLOORS = {"Geslin": (2.9910, 1.7457, 8.7176, 45), "KIT": (7.3323, 4.0902, 20.4751, 47)}


def _check(name, got_pct, printed, ndp=4):
    if round(got_pct, ndp) != printed and abs(got_pct - printed) > 0.5 * 10 ** -ndp + 1e-12:
        raise SystemExit(f"value differs from D56: {name}: {got_pct} vs {printed}")


def load():
    with open(os.path.join(ROOT, "results", "stage2_lever.json")) as f:
        kit = json.load(f)
    with open(os.path.join(ROOT, "results", "stage3_lever.csv")) as f:
        san = [r for r in csv.DictReader(f) if r["series"] == "registered 1C"]
    with open(os.path.join(ROOT, "inputs", "floor.json")) as f:
        gfl = json.load(f)
    warm, cold = [], []
    for band, key, ref in (("warm", "lever_25_40", D56_KIT_WARM), ("cold", "lever_10_25", D56_KIT_COLD)):
        rows = kit[key]
        assert len(rows) == 12 and all(r["determinable"] for r in rows)
        for r in rows:
            lv = {"band": band, "source": "KIT", "lever": f"KIT {r['protocol']}", "protocol": r["protocol"],
                  "secant": "25 to 40 C" if band == "warm" else "10 to 25 C",
                  "s": float(r["s_per_K"]), "lo": float(r["ci95"][0]), "hi": float(r["ci95"][1]),
                  "T_lo_C": float(r["T_lo_C"]), "T_hi_C": float(r["T_hi_C"]),
                  "file": "results/stage2_lever.json " + key}
            p = ref[r["protocol"]]
            for v, pv, nm in ((lv["s"], p[0], "s"), (lv["lo"], p[1], "lo"), (lv["hi"], p[2], "hi")):
                _check(f"{lv['lever']} {band} {nm}", 100 * v, pv)
            (warm if band == "warm" else cold).append(lv)
    for (chem, sec), p in D56_SANDIA.items():
        r = [x for x in san if x["chemistry"] == chem and x["secant"] == sec]
        assert len(r) == 1 and r[0]["determinable"] == "True"
        r = r[0]
        lv = {"band": "warm" if sec == "25 to 35 C" else "cold", "source": "Sandia", "lever": f"Sandia {chem} {sec}",
              "protocol": f"{chem} registered 1C", "secant": sec, "s": float(r["s_per_K"]), "lo": float(r["ci95_lo"]),
              "hi": float(r["ci95_hi"]), "T_lo_C": float(r["T_lo_C"]), "T_hi_C": float(r["T_hi_C"]),
              "file": "results/stage3_lever.csv registered 1C", "soh_threshold": float(r["soh_threshold"])}
        for v, pv, nm in ((lv["s"], p[0], "s"), (lv["lo"], p[1], "lo"), (lv["hi"], p[2], "hi")):
            _check(f"{lv['lever']} {nm}", 100 * v, pv)
        (warm if lv["band"] == "warm" else cold).append(lv)
    floors = {"Geslin": {"mean": float(gfl["mean"]), "median": float(gfl["median"]), "p95": float(gfl["p95"]),
                         "pairs": int(gfl["n_pairs"]), "file": "inputs/floor.json"},
              "KIT": {"mean": float(kit["floor"]["mean"]), "median": float(kit["floor"]["median"]),
                      "p95": float(kit["floor"]["p95"]), "pairs": int(kit["floor"]["pairs"]),
                      "file": "results/stage2_lever.json floor"}}
    for k, (m, md, q, n) in D56_FLOORS.items():
        _check(k + " mean", 100 * floors[k]["mean"], m)
        _check(k + " median", 100 * floors[k]["median"], md)
        _check(k + " p95", 100 * floors[k]["p95"], q)
        assert floors[k]["pairs"] == n
    assert len(warm) == 13 and len(cold) == 15
    return warm, cold, floors


if __name__ == "__main__":
    w, c, fl = load()
    for lv in w + c:
        print(f"{lv['band']:4s} {lv['lever']:28s} s={100*lv['s']:+.6f} [{100*lv['lo']:+.6f}, {100*lv['hi']:+.6f}]")
    print(fl)
