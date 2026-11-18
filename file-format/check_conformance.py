"""Runs the Python validators over conformance/cases.json and writes what they
report for each case to conformance/results-python.json: the rules, and whether
the file is conforming, readable or not readable. A "readable" file is then read,
and must give the same counts as the conforming file it was made from.

    python check_conformance.py
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, "format")
import sceptre_format as sf


def same_counts(a, b) -> bool:
    x, y = sf.read_screen(a), sf.read_screen(b)
    return all((x[m].X != y[m].X).nnz == 0 for m in y.mod)


cases = json.loads(Path("conformance/cases.json").read_text())
results, wrong = {}, 0
for case in cases:
    if case["kind"] == "screen":
        found = sf.validate_screen(case["path"])
    else:
        found = sf.validate_analysis(case["path"], case.get("screen"))
    result = {"rules": sorted(set(found["rule"])), "level": sf.conformance_level(found)}
    if "reads_as" in case:
        result["reads_correctly"] = same_counts(case["path"], case["reads_as"])
    if case.get("reader_refuses"):
        try:
            sf.read_screen(case["path"])
            result["reader_refused"] = False
        except ValueError:
            result["reader_refused"] = True
    results[case["name"]] = result
    if result["rules"] != case["expect"] or result["level"] != case["level"] or result.get("reads_correctly") is False or result.get("reader_refused") is False:
        wrong += 1
        print(f"MISMATCH {case['name']}: expected {case['expect']} ({case['level']}), reported {result}")
        print(found.to_string(index=False))
Path("conformance/results-python.json").write_text(json.dumps(results, indent=2) + "\n")
print(f"Python: {len(cases) - wrong} of {len(cases)} cases as expected")
sys.exit(1 if wrong else 0)
