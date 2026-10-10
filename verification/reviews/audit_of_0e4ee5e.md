Filed by the author with this header. The text below is the checker's report as it returned it; the checker could not write files itself, and the helper scripts it names stayed in its session and are not filed. It audited commit 0e4ee5e. Its three must-fix items, the sign of Law 1 (R1), the margin against the heat of dynamic profiles (R2) and the proxy cell in the note's limits (C7), and its should-fix items were taken up in DECISIONS.md D23 and D24, the note and the README, except these: C12 (the note keeps two headline items, each with its supporting numbers); C15 (the phrase it could not find is in the full brief, while the checker had a condensed copy); C17 (Zenodo's metadata file holds one licence; CC BY 4.0 for the note, figure and results is stated in the README and description); R16 (the commits carry the assistant's co-author lines).

# Independent audit: thermal-lever-duplicate-floor

Audited: the clone at /home/claude/independent_check/repo, branch main, HEAD 0e4ee5e (8 commits; the remote holds up to 8ede1ec). Rules: /home/claude/independent_check/BRIEF.md. Nothing in the repository was modified, and the original working copy was not touched. Helper scripts and saved evidence are in /home/claude/independent_check/audit/ (listed at the end).

Labels:
- **must fix**: wrong, misleading, or a brief rule broken.
- **should fix**: imprecise or incomplete, or a reviewer will ask about it.
- **fine**: checked and correct.

## Summary

- Every number printed in the note, README, .zenodo.json, CITATION.cff and the scored values of predictions.md was traced to its file and recomputed. None is wrong or mis-rounded.
- Separate implementations written for this audit reproduce results/lever.json: exactly or to 1e-12 for the module network and Law 1, and to 1e-9 for Law 2.
- Rerunning the study's own pipeline in a scratch clone:
  - Byte-identical: lever.json, module_spread*.csv, v1_v2.json, lgm50_resistance.json and the figure PNG.
  - Agreement to 1e-11: the resolved cell, V3, comparison.json and the scorecard.
- The floor matches its source byte for byte.
- The pre-registration is confirmed on GitHub.

Must fix (3):
1. The Law 1 sign claims ("so a hotter cell lives slightly longer", "the two laws disagree in sign") are contradicted by the registered Law 1 itself at every module temperature (R1).
2. "Below the 3.0 % mean floor under both laws at the dataset's conditions" holds for constant-current heat only (R2).
   - At C/2, 1.48 times that heat reaches the mean floor under Law 2.
   - The dataset's dynamic C/2 profiles plausibly give that much heat.
   - The note does not say so.
3. The note's Limits omit the proxy cell, which the brief lists first (C7).

Should fix: 26 items, listed per section below and collected at the end.

Not verified or not read:
- The Schimpe et al. Table IV values: no copy of the paper could be fetched.
- The 84 Kucinskis Fig. 7 counts: not present in the text extracted from the open PDF.
- https://api.crossref.org/works/10.2139/ssrn.7570666 could not be read. The WebFetch proxy returned HTTP 429 (rate limit) and asked that the page not be fetched again. That DOI, the source of the cold-plate base conductance, is unchecked.

---

## (1) Numbers

### N0. Every number in the note traced and recomputed (fine)

Sources: results/lever.json, comparison.json, implied_divergence.csv, resolved_cell.json, prediction_scorecard.json, verification/*.json, inputs/*. The recomputation is recompute.py, with output in recompute_out.json.

| Note line | Printed | Source and recomputed value | Status |
|---|---|---|---|
| 7 | 8.7 % | floor.json p95 0.087176 | fine |
| 7 | Law 2 needs 1.68 K at 35 C | dT_floor p95 1.6761 (C/16) to 1.6804 (C/2); 1.6821 at 2C | fine |
| 7 | 0.58 K for the 3.0 % mean | 0.5751 to 0.5766 | fine |
| 7 | Law 1 needs 8.8 K, opposite sign | -8.8137 K | fine (see R1 for the wording) |
| 7 | 102, 126 and 124 cycles at 25, 35, 45 C | sqrt(141 x 74) = 102.15; sqrt(169 x 94) = 126.04; sqrt(123 x 126) = 124.49 | fine |
| 8 | first offers 1.68 K at 2C, every h | module spread 1.377 to 1.519 K at 1C (below); 4.774 to 5.408 K at 2C | fine |
| 8 | 4.8 to 5.4 K | 4.774, 5.408, 5.246 | fine |
| 8 | cells 20 to 52 K above the air | at 2C: Tmin 55.18 C (h = 25), Tmax 86.73 C (h = 5) | fine |
| 8 | 1C gives 1.4 to 1.5 K | 1.407, 1.519, 1.377 | fine |
| 8 | at most 0.39 K; cell 0.09 K | 0.3913 (C/2, h = 10); 0.0921 (C/2, h = 5) | fine |
| 8 | at most 2.0 % (Law 2), 0.4 % (Law 1) | 0.051877 x 0.3913 = 2.030 %; 0.009891 x 0.3913 = 0.387 % | fine |
| 18 | 45 pairs; 2.99 %, 8.72 %, 4.10 % | n_pairs 45; 2.9910, 8.7176, 4.0975 | fine |
| 19 | 57 to 78 % calendar share | 0.5681 (2C) to 0.7792 (C/16) | fine |
| 19 | s = 0.051 against 0.052 | cycling terms only 0.05125; calendar only 0.05219; full model 0.0519 to 0.0520 | fine |
| 20 | 21.1 by 70.1 mm; 5 A; 18.71 mOhm | lever.json cell; D7; 0.0187128 ohm (also reproduced by rerunning PyBaMM) | fine |
| 21 | 6.2 mW/K, 17.4 mW/K, view factor 0.16 | recomputed 6.234, 17.425 mW/K, F = 0.1599; not in any results file (N1) | fine |
| 22 | 1.3 and 25 W/mK; 0.25 mm can | src/twocell/lumped_two_cell.py lines 56 to 59 | fine |
| 23 | 8.6 and 10.4 % | v1_v2.json 8.602, 10.353 | fine (see N2) |
| 23 | energy closes to 1e-12 | V2 max 4.9e-13; V3 max 7.5e-13 | fine |
| 23 | volume mean 1.8 to 5.2 % high | 1.773, 2.622, 5.167 % | fine |
| 23 | refinement 0.8 %; radial 0.13 % | 0.831 %; 0.127 % | fine (see R10) |
| 27 | 44.5 to 46.1 C; 1.5 K; 15 % | 44.542 and 46.061 C; 1.519 K; exact pair 14.78 % | fine (see R5) |
| 28 | 0.32 or 0.43 K | 0.3240, 0.4255 | fine |
| 28 | eight-fold at C/16; 0.05 K | ratio 8.000; 0.0437 to 0.0493 K | fine (see R9) |
| 29 | 0.65 K; 0.02 to 0.04 K; 0.29 K | 0.6545; 0.0156 to 0.0358; 0.2913 | fine |
| 35 to 39 | scorecard | identical to prediction_scorecard.json; "30 % high" = 0.12/0.092 = 1.30 | fine, except N3 |
| 43 | 0.6 K would suffice | 0.577 K | fine |

The PDF has one page (pdfinfo reports Pages 1, A4). Its text matches the markdown word for word (compare_note_pdf.py: 1071 words each). The only difference is the line-break hyphenation of "Si-graphite". Fine.

### N0b. README, .zenodo.json, CITATION.cff and predictions.md (fine)

README tables:
- Table 1 (lines 46 to 50): every s, L at 35 C and dT_floor is correctly rounded.
  - L at 35 C for Law 2: 1196.1, 1628.5, 2394.8, 3458.0.
  - The exact pair inversion differs by 0.063 %, under the 0.1 % stated on line 42.
- Table 2 (lines 56 to 61): all 48 entries are correctly rounded.
  - The heats are 1.83, 4.68, 18.7, 117, 468 and 1871 mW.
  - "Hottest always inner, coldest always corner" holds in all 18 cases.
- Table 3 (lines 69 to 74): all 24 entries are correctly rounded.
- Table 4 (lines 80 to 85): recomputed identically, and the answer is the same at every h.

Other README lines:
- Line 63: the means are 40.8 to 48.8 C at 1C and 57.8 to 84.3 C at 2C. At 2C every cell is above 55 C (minimum 55.18 C).
- Lines 89 to 95, secondary bullets: all correct.
  - Law 1 by charge rate: 0.036, 2.4 K, 0.007, 12.8 K, -0.027, -3.3 K.
  - Assumed-noise Monte Carlo: -0.0117 to -0.0081; 7.4 to 10.8 K; 0.0013 against a Monte Carlo SE of 0.0009 (my independent MC gives sd 0.00094).
  - Law 2 by term: 57 to 78 %; 0.052 and 0.051; text values within 0.1 % (largest difference 0.083 %).
  - Charge-heat correction: 8.0, 5.0, 2.5; 0.044 to 0.049 K; 0.26 %.
  - Cold plate: 0.016 to 0.036 K; 0.29 K; 1.16 K.
  - Sensitivities: 0.448, 0.347, 0.655, 0.324, 0.425 K; 15.47 and 20.37 mOhm.
  - Exploratory: 14.8 %; 31 cycles.
- Lines 99 to 101, verification bullets: all correct.

The other files:
- .zenodo.json description: all correct.
- predictions.md scored thresholds are the brief's. The outcomes recompute identically:
  - P1 failed, P2 held, P3 held.
  - P4 failed: the first clause held for both laws, the second failed for both.
  - P5 failed, with a ratio of 5.245 to 5.258.

### N1. lever.json loses the network constants (should fix)

- src/thermal_lever.py line 329 overwrites meta["module"] from line 325 (pitch, k_air, eps, G_strip, G_gap, F_view, areas) with "module": rows.
- So the note's 6.2 mW/K, 17.4 mW/K and 0.16 have no results-file source. The values themselves are correct.

### N2. V1 for Law 2 does not test the registered model (should fix)

- Printed in note line 23, README line 99 and predictions.md line 10: "8.6 and 10.4 %".
- src/verify_pre.py (lines 9 to 10, 59) runs V1 with terms=("cal","high","low"), leaving out the high-SOC term. The lever uses all four terms (D4).
- With all four terms (v1_terms_check.py): 8.91 % at 25 C and 10.35 % at 45 C. Both still pass.
- Fix: report these values, or say that V1 omits the term.

### N3. P5 central estimate printed as 5.2 (should fix)

- Printed: "(central 5.2)" (note line 39), "(5.2)" (README line 36), and "central_ratio": 5.2 (prediction_scorecard.json line 85, hard-coded at compare.py line 174).
- predictions.md line 31 registered "a factor of about 5". 5.2 was computed afterwards.

### N4. Repository creation time (should fix)

- README line 24 says "This repository was created and that commit pushed at 08:01:59 UTC".
- The GitHub API gives created_at 08:01:27Z. The push was at 08:01:59Z. See (5).

---

## (2) Compliance

### C1. No em or en dashes (fine)

- scan_dashes.py covered 42 tracked files, the text of both PDFs, all commit messages and the tags.
- The only matches were U+2212 minus signs, which are not dashes: note lines 19 and 35, and the minus signs in the figure's tick-label exponents.

### C2. A run-time estimate in the README (should fix)

- README line 133 says "(minutes)". This is a duration estimate, which the brief bans. It is also inaccurate: the step took 25 s in my rerun.
- No other schedule or timeline language was found.

### C3. No pilot labels (fine)

- There is no letter or numbered-pilot label in any file, commit message, title, PDF metadata, the GitHub description or the repository name.
- The word "pilot" appears only in the other studies' repository names and in the two-cell study's Zenodo title.

### C4. Title (fine)

- The exact title appears in note line 1, .zenodo.json line 2, CITATION.cff line 3 and the GitHub description.
- The README heading is shortened, with the subtitle on line 3. The brief allows this.
- The repository is named thermal-lever-duplicate-floor.

### C5. "Measurement of the thermal share" (fine)

Nothing calls the result a measurement of the thermal share; README line 18 says the opposite. See also R15.

### C6. The proxy statement is present (fine)

Note line 20 says the dataset's cells are "a different 21700 from the dataset's that sets magnitudes, not its values". README line 117 and D6 say the same. See F7 on whether the dataset's cells are 21700.

### C7. The note's Limits omit the proxy cell (must fix)

- Brief line 26 lists the limits as "proxy cell; steady cycling heat; laws fitted on other cells; lumped module, one geometry".
- The Limits section (lines 25 to 29) covers the other three. The proxy cell appears only in Method item 3.

### C8. The other three limits are present (fine)

### C9. Scorecards match the JSON and the scoring rule (fine)

The note and README scorecards match prediction_scorecard.json and the rule in predictions.md (lines 28 and 35). The one exception is N3.

### C10. The pre-registration commit is intact (fine)

- git diff 24ab8d3 HEAD -- predictions.md is empty.
- At 24ab8d3 there is no results/ or note/ content. The only output is verification/v1_v2.json, and the brief allows V1 and V2 before the predictions.
- No pushed commit contains results.

### C11. Commit order (fine)

The history is:
1. 0f9c22c: the brief's script as received.
2. 24ab8d3: the pre-registration.
3. 7791192: solver change (D21) and release metadata.
4. Push at 08:01:59Z.
5. 8ede1ec, pushed at 08:02:32Z.
6. e13a3f6: V3 check change (D22).
7. fbf5a18: results.
8. dfa3761: figure and note.
9. 0e4ee5e: README.

Since the push, the code changed only in the V3 check (e), the convergence flag and the figure layout. No registered calculation changed.

### C12. Note format against the brief (should fix)

- The brief asks for "two headline numbers". Lines 7 and 8 carry about 20 numbers.
- Method and Limits have six and three items, but each runs two to four lines.

### C13. L(T) saved at C/2 only (should fix)

- Brief Step 1 asks for L(T) from 25 to 55 C at C/16, C/10, C/5 and C/2.
- lever.json stores the curves only at C/2 (thermal_lever.py line 310). Law 2 differs by rate.

### C14. The V3 acceptance change is not listed as a change from the brief (should fix)

- The brief's V3 test (volume mean within 1 % of the lumped cell) fails, at 1.8 to 5.2 %.
- Acceptance moved to the surface mean (D15). This is disclosed under Verification, but not in "What changed from the brief" (README lines 103 to 108).

### C15. The README quotes words that are not in the brief (should fix)

README line 108 says the brief allowed "if it is the case". BRIEF.md contains no such phrase.

### C16. Three properties lack a primary published source (should fix)

- The values: k_air 0.03 W/(m K) (D9), the can conductivity 15 W/(m K) (two-cell lumped_two_cell.py line 57) and nickel k 90 W/(m K).
- All are inherited from the two-cell study without a primary source. That study's D24 calls the nickel value "illustrations only".
- The brief prescribes the strip value, so this is a recording gap, but the "never invent a physical parameter" rule asks for a paper, table and page.

### C17. Release metadata details (should fix)

- .zenodo.json sets "license": "MIT" only. CC BY 4.0 appears only in the description text. Make sure the Zenodo record lists both, as the floor study's record does.
- The version is "1.0.0" in .zenodo.json and "v1.0.0" in CITATION.cff.
- No DOI yet, which is expected before release.

### C18. Release metadata against the brief (fine)

- Upload type software; creator Satte, Siddharth, with ORCID 0009-0001-3584-1008.
- 10.5281/zenodo.23205413 and 10.5281/zenodo.23253738 are marked isDerivedFrom.
- The law papers' DOIs are marked references, plus Geslin, O'Regan, Chen and PyBaMM.
- CITATION.cff is valid YAML, lists MIT and CC-BY-4.0, and carries the ORCID.

### C19. Deliverables (fine)

- Present: the results CSVs (module_spread, module_spread_coldplate, resolved_cell, implied_divergence), the results JSONs, the figure as PNG and PDF, the one-page note in markdown and PDF, and the release metadata.
- BRIEF.md does not describe the figure's panels. The three panels map onto Step 1, Step 2(a)/(b) and Step 2(c)/Step 3.
- The figure matches its caption and regenerates pixel for pixel from the committed results.

---

## (3) Floor and inputs

### F1. The floor matches its source (fine)

- I cloned geslin-floor-pilot. Tag v1.0.1 is commit 103f7e7f9f19d267cfe0c2e5f4a4a702d80e410b.
- results/results.json there has SHA-256 be63577f328a2c49850b9ed6d9242459c89569e7b1a065ccce7880d3d825b8ac, the value recorded in the study.
- The mean, median, p95 and worst_of_8 values are bit-identical.
- The bootstrap intervals equal the headline_percentile_bootstrap_95 keys, which floor.json does not name.
- D16's sigma (2.65 %) equals A2_cell_noise_sigma.soh0.900.sigma_from_mean, 0.026507.

### F2. The two-cell study files match their source (fine)

- src/twocell/*.py and oregan2022_properties.txt match SHA256SUMS and the files at m50-two-cell-pilot tag v1.0.0 (2b35ac7aa1a2cc822ed59f61fc73128b220bdd1c).
- The strip value is that study's D24. The radiation and convection forms are its D5 and D6.

### F3. The derived resistance reproduces (fine)

- Rerunning lgm50_resistance.py with PyBaMM 26.9.0.0 gives a byte-identical file: R = 18.713 mOhm.
- The effective resistance ranges from 18.1 to 20.2 mOhm between C/10 and 2C. The pulse values are 15.47 and 20.37 mOhm.

### F4. Schimpe Table IV could not be verified

| Attempted source | Result |
|---|---|
| WebFetch iopscience.iop.org article and /pdf | ROBOTS_DISALLOWED |
| opus4.kobv.de (Kempten open copy) | Anubis bot-check page |
| mediatum.ub.tum.de/node?id=1446747 | Anubis: access denied |
| OSTI search and API | ROBOTS_DISALLOWED |
| Semantic Scholar API | client error |
| Internet Archive Scholar | rate limited |
| Wayback Machine | site blocked |
| CORE | ROBOTS_DISALLOWED |
| Shell (curl) | the shell cannot reach these hosts |

Only internal consistency could be checked:
- Every Table IV value is the rounding of the text value listed beside it: 3.69e-4/3.694e-4, 2.06e4/20592, 1.46e-4/1.456e-4, 4.01e-4/4.009e-4, 3.27e4/32699, 5.55e4/55546, 2.03e-6/2.031e-6, 2.33e5/2.3e5, 7.84/7.8.
- The text values change s by less than 0.1 %.

**Status: unverified.** These still need checking against the PDF: Table IV, SOC_Ref, I_Ch,Ref, U_a,Ref, alpha, k_0, Table AI and Eq. 24.

### F5. Kucinskis Fig. 7 counts could not be checked

- The dspace.lu.lv item and its bitstream 3fb0395b were read with WebFetch. The Fig. 7 numbers are not in the extracted text, and a shell download was refused (proxy CONNECT 403).
- What could be checked matches the inputs:
  - The values quoted on p. 7 (662, 115, 387) equal the CSV.
  - "Approximately half ... at 15 °C" holds: 332/662 = 0.50.
  - The row-maximum statements match.
  - The cell is a 5 Ah 21,700-type NCA LiNi0.90Co0.05Al0.05O2 with about 6 % LNO and about 3 % Si.
  - The 1 kHz resistance is 12.9 mOhm.
  - Two or more cells per condition, but one cell at 55 C.

### F6. Law 2 calendar factor: arithmetic mean or RMS (should fix: check against Eq. 24)

- D4 uses the arithmetic mean of the calendar SOC factor (1.0989).
- If Eq. 24 integrates the square-root law in state form (dQ/dt = k^2/(2Q)), a cycling SOC gives the RMS instead (1.2145).
- Effect:
  - s changes by under 0.05 %, so no dT_floor moves.
  - L at 35 C would be 12 to 15 % lower (3045 against 3458 cycles at C/2).
  - The calendar share would be 59 to 80 % rather than 57 to 78 %.
  - V1 at 25 C would be 9.2 %.
- This is a question to check, not a finding of error.

### F7. The dataset cells' 21700 format is not in Geslin et al.'s main text (should fix: give a source)

- Printed: README line 7, note line 20, D2, D6.
- The Methods (p. 177) say: "Ninety-two commercial EV energy lithium-ion cells (silicon oxide graphite / nickel cobalt aluminium) were cycled ... with four-point contact cylindrical cell fixtures".
- The format is not given on the Stanford Digital Repository (SDR) landing page or in the authors' code README either. The Supplementary Information was robots-blocked.

---

## (4) Citations

Crossref was read with WebFetch. The Geslin pages come from the nature.com page, because Crossref has no volume or pages for it. The two Zenodo records were read through DataCite, because zenodo.org/api/records is robots-disallowed for WebFetch.

| Reference | Authors | Title | Journal, volume, pages, year | Status |
|---|---|---|---|---|
| Geslin et al. 10.1038/s41560-024-01675-8 | 6, in order | match | Nature Energy 10, 172 to 180 (2025) per nature.com; Crossref has the online date only | fine |
| Kucinskis et al. 10.1016/j.jpowsour.2022.232129 | 6, match | match, but the publisher's subtitle separator is an en dash and README/SOURCES use a colon | J. Power Sources 549, 232129 (2022) | fine (forced by the no-dash rule) |
| Schimpe et al. 10.1149/2.1181714jes | 6, match | match (sentence case) | JES 165(2), A181 to A193 (2018) | fine |
| O'Regan et al. 10.1016/j.electacta.2022.140700 | 4, match | match | Electrochim. Acta 425, 140700 (2022) | fine |
| Chen et al. 10.1149/1945-7111/ab9050 | 6, match | match | JES 167(8), 080534 (2020) | fine |
| Sulzer et al. 10.5334/jors.309 | 5, match | match | JORS 9(1), 14 (2021) | fine |
| Satte 10.5281/zenodo.23205413 | Satte, Siddharth | "The cell-to-cell lifetime floor in the Geslin et al. dynamic-cycling dataset" | version 1.0.1, Software, 2026 | fine |
| Satte 10.5281/zenodo.23253738 | Satte, Siddharth (ORCID) | "Where and when a neighbour cell first reaches thermal-runaway onset temperatures: resolved against lumped, a pre-registered two-cell LG M50 pilot" | version v1.0.0, Software, 2026 | fine |

- Each DOI string is identical everywhere it appears. Law 1's and Law 2's DOIs agree between SOURCES.md (lines 11, 13) and DECISIONS.md (lines 11, 13).
- 10.5281/zenodo.23053419 resolves to "Data and code for: Current collection and event duration set the error of uniform heat-generation models ..." v1.0 by Satte. Fine.
- Unchecked:
  - 10.1016/j.jpowsour.2014.03.112 (Waldmann): Crossref returned a client error.
  - 10.2139/ssrn.7570666: HTTP 429.
  - Neither appears in the README references or .zenodo.json.
- Outside this repository: the floor study's DataCite record has no creator ORCID.

---

## (5) Public record

WebFetch to api.github.com failed (CLIENT_ERROR). The same endpoints were read with curl and saved in audit/evidence/.

| Endpoint | Result |
|---|---|
| repository | public, visibility public, created_at 2026-10-10T08:01:27Z, pushed_at 08:02:32Z |
| commits | only 8ede1ec, 7791192, 24ab8d3, 0f9c22c are on GitHub |
| activity | branch_creation at 08:01:59Z (from 0000000 to 7791192, which contains 24ab8d3); push at 08:02:32Z (7791192 to 8ede1ec) |
| events | an empty array |

- User-level events are blocked by the session proxy (403). I can confirm the repository is public now and that 24ab8d3 was pushed at 08:01:59Z, but not that it was public from creation. Nothing suggests a private period.
- README line 24's push time is right, but the creation time should read 08:01:27 (N4). The note's "public at 08:02 UTC" is the rounded push time, which is fine.
- Results commits: e13a3f6 and fbf5a18 at 08:26:20; dfa3761 and 0e4ee5e at 08:26:21. They are local and record the order of work, not run times.
  - PDF creation stamps support the order: figure D:20261010081215Z, note D:20261010082237+00'00'.
  - Once pushed, the GitHub push time becomes the public proof.
- 8ede1ec's commit message omits that it also added note.css and render_note.sh (styling files only). This is fine and cannot be changed without rewriting pushed history.

---

## (6) Reviewer flags

### R1. Law 1's sign is unresolved, and the module case contradicts the stated sign (must fix)

Printed:
- Note line 7: "so a hotter cell lives slightly longer".
- Note line 43: "where the two laws disagree in sign".
- README lines 16 and 18 say the same.

Evidence (law1_sign_check.py):
- The registered s = -0.0099/K is the centred secant from 25 to 45 C. The one-sided slopes at a C/2 charge are -0.0210/K (25 to 35 C) and +0.0012/K (35 to 45 C). On the upper segment a hotter cell lives shorter, as under Law 2, and dT_floor(8.7 %) would be 70 K.
- Every module cell is above 35 C. Evaluated with the registered piecewise law at each cell's own temperature (D18), the hotter cell lives shorter in all 12 cases from C/16 to C/2. At C/2, h = 10: -0.048 % against +0.387 % from s x dT.
- implied_divergence.csv stores that pair difference unsigned, which hides the sign.
- The change from 35 to 45 C (-1.2 %) is smaller than D16's assumed per-value noise (1.9 %).
- README line 89 already says the sign is unknown.

Fix: say the sign is undetermined. The thesis conclusion still holds.

### R2. The headline margin is a factor of 1.5 in heat (must fix)

Printed in note lines 8 and 43 and README lines 14 and 17: "below the 3.0 % mean floor".

Evidence (margin.py):
- The module reaches Law 2's dT_floor(3.0 %) = 0.577 K at 0.61C (h = 10). At h = 5 and h = 25 this happens at 0.63C and 0.64C.
- At C/2 this needs 1.48 times the brief's heat, which a dynamic discharge supplies at I_rms/I_mean of 1.40 or more.
- Geslin et al. (p. 175) report discharge current peaks up to 1,800 % of the average.
- D19 names the mean-square factor but does not quantify it, and the note is silent. Its limit 2 reassures only about the C/16 charge correction.
- The missing contact resistance and the 30 s pulse value (+9 % spread) push in the same direction.

Fix: state the margin and the dynamic-profile caveat, or compute the mean-square factor from the dataset's public current traces.

### R3. The 2C comparison lies outside both laws (should fix)

- Note lines 8 and 43 compare a spread at 55 to 87 C with the lever taken at 35 C.
- At 2C all cells are above 55 C, outside both laws' fitted ranges. README line 63 says so; the note does not.

### R4. The worst-of-8 floor is compared with the wrong spread (should fix)

- README Table 4 (lines 81 and 84) and figure panel (c) compare the worst-of-8 floor with hottest minus coldest.
- D17's matching quantity is hottest minus mean. With it:
  - Law 2's module row moves from 1C to 2C.
  - Law 1's moves from 2C to "none up to 2C".

### R5. The Law 1 exploratory 15 % is fragile (should fix)

- The note line 27 pair straddles 45 C: the coldest cell is at 44.54 C.
- The result depends on h: 19.5 % (h = 5), 14.8 % (h = 10), 0.17 % (h = 25).
- The steep segment rests on single-cell 55 C points, where 95 % and 90 % SoH tie at 31 cycles.
- At the cells' own temperatures, Law 1 exceeds 8.7 % at 1C for h = 5 and 10. That conflicts with "and under Law 2 alone" (line 43).

### R6. "Fixtures of unreported temperature" (should fix)

- Note line 43 and README lines 18 and 108 use this phrase.
- Geslin et al. write: "The batteries were held at a constant temperature of 35 °C inside a CSZ temperature chamber". Only the cells' own temperatures are unreported.

### R7. "It shows that ... a module test could not separate it from cell-to-cell chance" (should fix)

- README line 18.
- This is a model estimate under two laws fitted on other cells and a proxy module.
- A systematic inner-against-corner difference of about 2 % could be separated statistically by repeating the module test.

### R8. "The closest published match" (should fix)

- Note line 7, README lines 16 and 113, and .zenodo.json. D3 says "found", which these drop.
- "(... charged at C/2)" reads as if the paper tested a C/2 charge. It tested 0.2C, 0.4C and 0.6C; C/2 is interpolated.

### R9. Low-rate spreads are time averages (should fix)

- Note limit 2 (line 28) and README line 92 quote time-averaged spreads.
- Transient check (charge_phase_spread.py, with 61.75 J/K per cell):
  - The network's time constants are 0.09 to 0.39 h.
  - During each 2 h C/2 charge the spread averages 0.28 to 0.30 K, at every rate from C/16 to C/5.
  - The time-weighted average at C/16 is 0.049 K.
- Charge-driven ageing therefore sees about six times the quoted spread. It is still below every dT_floor.

### R10. The V3(e) metric (should fix)

- Note line 23 and README line 101 say the radial profile is "within 0.13 %". That figure is the error over the maximum rise.
- The centre-to-surface difference, which is what dT_cell measures, is 0.1005 K numerically against 0.1021 K analytically: 1.6 % low.

### R11. The figure shows magnitudes without saying so (should fix)

- Panels (a) and (c) plot absolute values.
- The caption (line 14) says "s × ΔT" and never says that Law 1's values are negative.

### R12. Expected failures are unmarked in the note's scorecard (should fix)

- predictions.md (lines 17, 27, 31) expected P1, P5 and P4's second clause to fail.
- README line 38 says so. The note's scorecard (lines 35 to 39) does not.

### R13. Stale cross-references (should fix)

- DECISIONS.md line 3 says "D1 to D19". D20 was already present at 24ab8d3.
- thermal_lever.py line 19 says "(D1 to D16)".
- lgm50_resistance.py line 7 cites D7; it should be D8.
- twocell/resolved_cell2.py line 4 cites "D5 to D16", meaning the two-cell study's decisions. That file is copied unchanged, so note this in the README Files section instead.

### R14. Informal name (should fix)

DECISIONS.md line 29 says "which Sid decided". It should name the author.

### R15. Computed values labelled "measured" (should fix)

prediction_scorecard.json uses measured_K, measured and measured_ratio (lines 4, 30, 41, 51, 79) for computed values.

### R16. AI co-author trailers (fine, author's call)

All 8 commits carry a Claude co-author trailer and a session URL. The README and note mention no AI assistance.

### R17. Law 2 at C/16 and C/10 (fine)

Law 2 is applied at C/16 and C/10, below the paper's 0.25C to 1C tests. D4 states this, and the calendar term dominates at those rates.

### R18. Voltage windows (fine)

Kucinskis cycles 4.2 to 2.5 V; the dataset cycles 4.2 to 3.1 V (Geslin p. 173). This is covered by note limit 1.

### R19. "C/2 cycling heat" (fine)

The note means irreversible heat. Including the net reversible heat, the average is 0.124 W against 0.117 W, 6 % more.

---

## All should-fix items

N1, N2, N3, N4, C2, C12, C13, C14, C15, C16, C17, F6, F7, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15.

## Files in /home/claude/independent_check/audit/

- recompute.py and recompute_out.json: the independent Law 1, Law 2 and module network; Tables 1 to 4; the secondary results; Step 3 with hottest minus mean.
- scan_dashes.py: the dash scan.
- check_floor.py: floor values and hash against the source.
- compare_note_pdf.py: note markdown against the PDF.
- law1_sign_check.py: Law 1's signed pair at each module cell temperature (R1).
- margin.py: the rate and heat at which Law 2 reaches the mean floor (R2).
- v1_terms_check.py: V1 with and without the high-SOC term, and mean against RMS.
- charge_phase_spread.py: the transient check (R9).
- rerun_pipeline.sh: the reproduction run in a scratch clone.
- evidence/: the GitHub API responses (repo, commits, activity, events), the upstream floor results.json at v1.0.1, and the rerun logs.

REPORT.md itself is not in the folder, because the subagent harness refused to write it.
