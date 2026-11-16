#!/bin/sh
# Stage 6 in order (DECISIONS.md D83 to D96). Where the struct files are (the author's laptop), with Python 3.10,
# numpy and h5py, the readers were run as:
#   python3 src/stage6/a0_read.py <file.mat> <label> <directory>     (Part A0; labels b1, b2, b3, b9)
#   python3 src/stage6/l_read.py <file.mat> <label> <directory>      (the lifetimes, after Part A0 was pushed)
# and their outputs are committed in inputs/matr/. Everything below runs from those outputs.
set -e
cd "$(dirname "$0")/../.."
python3 -I tests/stage6_synthetic.py
python3 -I src/stage6/a0_cells.py
python3 -I src/stage6/a0_analysis.py
# after Part A0 was committed and pushed:
python3 -I src/stage6/l_table.py
python3 -I src/stage6/s6_parts.py
python3 -I src/stage6/s6_score.py
python3 -I verification/stage6/v1_files.py
python3 -I verification/stage6/v2_separate.py
python3 -I verification/stage6/v3_hand.py
python3 -I src/stage6/s6_paper.py
python3 -I src/stage6/s6_figure.py
