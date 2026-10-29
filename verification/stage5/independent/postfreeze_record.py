#!/usr/bin/env python3
"""Stage 5 independent check: the record of the comparison phase (after the freeze), built from out_compare/*.json so
that it cannot differ from them. The frozen record.json is left as it was at ea3a5b4.

Output: postfreeze_record.json.
Usage: PYTHONDONTWRITEBYTECODE=1 python3 verification/stage5/independent/postfreeze_record.py
"""
import hashlib, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ind_common as C   # noqa: E402

OUTC = os.path.join(C.HERE, "out_compare")


def main():
    cmp_ = json.load(open(os.path.join(OUTC, "compare.json")))
    txt = json.load(open(os.path.join(OUTC, "text_review.json")))
    rec = {
        "check": "Stage 5 independent check, comparison phase",
        "freeze": {"commit": "ea3a5b46a891e74d504c818d30b35d34b4351ff5", "committed_utc": "2026-10-29T04:40:07Z",
                   "pushed_utc": "2026-10-29T04:40:14Z", "branch": "origin/stage5-independent"},
        "opened_after_the_freeze": {
            "fetch": "git fetch origin +refs/heads/stage5:refs/remotes/origin/stage5, after the freeze push",
            "read_with_git_show": {"101ca54a199def64ab7227a41277e2408b7f5e91": "registered results, V1 to V4 and the figure (pushed 2026-10-29T04:22:03Z)",
                                   "157827f1a77ef43e1f4545fc1c6fb24c2fb42270": "D79, D80, results/stage5_posthoc.json, results/stage5_paper_numbers.json",
                                   "67685ae92dd66fd10359456d21fde8e66a87bc1c": "README, note, CITATION.cff, .zenodo.json, D80 rounding",
                                   "ee4a403eb7bc23e41e1fdf35facd289995a46cb0": "the code change D79 (1) describes (diff against 101ca54)"},
            "branch_never_switched": True},
        "deterministic_all_agree": cmp_["deterministic_all_agree"],
        "table": cmp_["table"],
        "verdicts": {"P1_to_P8_mismatches": cmp_["verdicts"]["P1_to_P8"]["mismatches"],
                     "floor_verdict_mismatches": cmp_["verdicts"]["floor_verdicts"]["mismatches"],
                     "v4_pass_differences": cmp_["verdicts"]["v4_differences"]},
        "monte_carlo_summary": cmp_["monte_carlo_summary"],
        "text_review": {"numbers_checked": txt["numbers_checked"],
                        "numbers_not_matching": [(c["where"], c["quoted"], c["pipeline"]) for c in txt["numbers_not_matching"]],
                        "findings": [{k: f[k] for k in ("id", "kind", "where", "suggest")} for f in txt["findings"]]},
        "post_freeze_files": {}}
    for fn in sorted(os.listdir(OUTC)):
        rec["post_freeze_files"]["out_compare/" + fn] = hashlib.sha256(open(os.path.join(OUTC, fn), "rb").read()).hexdigest()
    for fn in sorted(os.listdir(C.HERE)):
        if fn.startswith(("compare", "postfreeze")) and fn.endswith(".py"):
            rec["post_freeze_files"][fn] = hashlib.sha256(open(os.path.join(C.HERE, fn), "rb").read()).hexdigest()
    with open(os.path.join(C.HERE, "postfreeze_record.json"), "w") as f:
        json.dump(rec, f, indent=1)
    print(json.dumps({k: rec[k] for k in ("deterministic_all_agree", "verdicts")}, indent=1))


if __name__ == "__main__":
    main()
