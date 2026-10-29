#!/bin/sh
# Stage 5 in order (DECISIONS.md D71 to D78): the synthetic tests, Part A, Part B, Part C with V4, the scoring,
# the checks V1 to V3 and the figure. Part C's simulations take a few minutes.
set -e
cd "$(dirname "$0")/../.."
python3 tests/stage5_synthetic.py
python3 src/stage5/s5_inputs.py
python3 src/stage5/s5_within.py
python3 src/stage5/s5_power.py
python3 src/stage5/s5_score.py
python3 verification/stage5/v1_inputs.py
python3 verification/stage5/v2_separate.py
python3 verification/stage5/v3_hand.py
python3 src/stage5/s5_figure.py
