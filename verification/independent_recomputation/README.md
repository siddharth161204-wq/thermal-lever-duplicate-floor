# Independent numerical check of "The thermal lever against the duplicate-cell floor"

Everything here was written for this check, by a separate checker who had not seen the study being built, from the
study's DECISIONS.md (D1 to D22), predictions.md and the files in `inputs/`. No code from `src/thermal_lever.py`,
`resolved_cell_steady.py`, `compare.py` or `verify_pre.py` was read or used. `src/twocell/` was read only for the
resolved cell's geometry, and the received script in `provenance/` only for the form of its radiation terms (one
sensitivity row). The repository is opened read only.

The check was run on 10 October 2026 against commit 0e4ee5e, before the corrections of DECISIONS.md D23 (which renamed
some scorecard keys and left Law 2's exact forms at 2C blank), so `make_report.py` compares with the results of that
commit. To rerun it as it was run: `git worktree add ../tl-0e4ee5e 0e4ee5e` from the repository root, then set
`THERMAL_LEVER_REPO` to that folder. The other scripts read only `inputs/`, which D23 did not change. Its report is
`REPORT.md`.

Requirements: Python 3 with numpy, scipy and pyamg (pyamg is used only as a preconditioner on the finest
resolved-cell grid; SuperLU gives the same numbers to 10 digits).

## Run (in this folder, in this order)

    python3 law1.py            # Law 1 lever, printed-rate secants, exact pair form, assumed-noise Monte Carlo   -> out/law1.json
    python3 law2.py            # Law 2 (Eq. 23 with Eqs. 9, 15, 18, 21, A1, A2), lever, term shares, V1        -> out/law2.json
    python3 module.py          # 2x4 module in air, 1x8 on a cold plate, two-node reduction, D19, sensitivities -> out/module.json
    python3 resolved_cell.py   # axisymmetric finite-volume cell, five grids, base-cooled case                  -> out/resolved_cell.json
    python3 scoring.py         # implied divergence, Step 3 lowest rates, D19 implied, P1 to P5                 -> out/scoring.json
    python3 verify_fv_series.py  # analytic-series check of the cell solver, shell effect, k_r sensitivity      -> out/verify_fv_series.json
    python3 make_report.py     # compares every quantity with ../repo/results and ../repo/verification          -> out/comparison_all.json, out/tables.md, out/summary.md

`python3 resolved_cell.py quick` runs a three-grid version.

## Files

- `common.py`: inputs (floors, resistance, geometry, couplings) typed from DECISIONS.md or read from `inputs/`.
- `law1.py`, `law2.py`, `module.py`, `resolved_cell.py`, `scoring.py`: the computations.
- `verify_fv_series.py`: checks the cell solver against an analytic Bessel-series solution (no shell, linear Robin),
  the effect of the shell, and D6's k_r = 1.24 W/mK statement -> `out/verify_fv_series.json`.
- `make_report.py`: every comparison (repository value, this check, relative difference, verdict), written to
  `out/comparison_all.json` (machine-readable), `out/tables.md` (every compared quantity, by section),
  `out/summary.md` and `out/named_tables.json`; relative difference = (repository - this check) / this check.
  `python3 make_report.py --fill < template.md > REPORT.md` fills `{{name}}` placeholders in a report template
  with those tables.
- `REPORT.md`: the checker's report, filed here as returned; its tables are the ones in `out/named_tables.json`.
- `out/`: all outputs, plus `resolved_cell.log` from the last full run.
