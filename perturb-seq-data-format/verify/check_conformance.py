"""Runs the Python validator over conformance/cases.json and writes what it
reports for each data file to conformance/results-python.json: the rules, and
whether the file conforms. A file that does not conform and must still be read
(S4.4) is then read, and must give the same counts as the conforming file it was
made from.

    python verify/check_conformance.py
"""

import json
import sys
from pathlib import Path

import mudata as md

sys.path.insert(0, "format")
import perturb_seq_data as sf


def same_counts(a, b) -> bool:
    x, y = md.read_h5mu(a), md.read_h5mu(b)
    return all((x[m].X != y[m].X).nnz == 0 for m in y.mod)


cases = json.loads(Path("conformance/cases.json").read_text())
results, wrong = {}, 0
for case in cases:
    found = sf.validate_screen(case["path"])
    result = {"rules": sorted(set(found["rule"])), "conforms": len(found) == 0}
    if "reads_as" in case:
        result["reads_correctly"] = same_counts(case["path"], case["reads_as"])
    if "conforming_with_index_limit" in case:
        # the 64-bit branch of S4: no test file reaches the real limit, so the limit is lowered.
        # Only the modality that passes the lowered limit may stay wide; the others then break S4.
        lowered = sf.validate_screen(case["path"], index_limit=case["conforming_with_index_limit"])
        result["wide_branch"] = not any("mod/rna/X" in w for w in lowered["where"])
    if case.get("reader_refuses"):
        try:
            sf.check_screen(md.read_h5mu(case["path"]))
            result["reader_refused"] = False
        except ValueError:
            result["reader_refused"] = True
    results[case["name"]] = result
    failed = [k for k in ("reads_correctly", "reader_refused", "wide_branch") if result.get(k) is False]
    if result["rules"] != case["expect"] or result["conforms"] != case["conforms"] or failed:
        wrong += 1
        print(f"MISMATCH {case['name']}: expected {case['expect']}, reported {result}")
        print(found.to_string(index=False))
Path("conformance/results-python.json").write_text(json.dumps(results, indent=2) + "\n")
print(f"Python: {len(cases) - wrong} of {len(cases)} cases as expected")
sys.exit(1 if wrong else 0)
