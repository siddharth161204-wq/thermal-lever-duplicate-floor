#!/usr/bin/env bash
# Independent Stage 5 check: every step in order, from the repository root.
# Writes verification/stage5/independent/out/ and record.json; touches nothing else.
set -euo pipefail
cd "$(dirname "$0")/../../.."
export PYTHONDONTWRITEBYTECODE=1
D=verification/stage5/independent
python3 $D/pairs.py                     # Part A: pairs, counts, hashes, shapes (out/pairs.csv, out/pairs.json)
python3 $D/partb.py                     # Part B (out/B.json, out/B.csv)
python3 $D/partc.py                     # Part C (out/C.json, out/C1.csv, out/C2.csv), about 1 minute
python3 -B $D/hand_v3.py                # V3 by hand in plain Python (out/v3_hand.json)
python3 $D/crosscheck_b.py              # Part B by a second method and its Monte Carlo spread (out/B_crosscheck.json)
python3 $D/crosscheck_stage4_power.py   # Part C calculation against Stage 4's power_t and n_required
python3 $D/mc_spread_c.py               # Part C Monte Carlo spread (out/C_mc_spread.json), about 5 minutes
python3 $D/influence.py                 # leave-one-set-out influence on beta (out/influence.json)
python3 $D/score.py                     # P1 to P8 (out/scorecard.json, out/scorecard.csv)
python3 $D/make_record.py               # record.json
