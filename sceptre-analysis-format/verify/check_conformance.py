"""Runs the Python validator over conformance/cases.json and writes what it
reports for each analysis folder to conformance/results-python.json: the rules,
and whether the folder conforms. Each case names the data file the analysis is
checked against, in the data file format's own conformance set.

    python verify/check_conformance.py
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, "format")
import sceptre_analysis as sf

cases = json.loads(Path("conformance/cases.json").read_text())
results, wrong = {}, 0
for case in cases:
    found = sf.validate_analysis(case["path"], case.get("screen"))
    result = {"rules": sorted(set(found["rule"])), "conforms": len(found) == 0}
    if "fits_used" in case:
        # whether a reader takes fits from the cache of this folder (T6.5)
        result["fits_as_expected"] = (sf.FITS in sf.read_analysis(case["path"])["tables"]) == case["fits_used"]
    if "fits_for" in case:
        fits = sf.read_fits(case["path"], case["fits_for"]["identity"])
        result["fits_as_expected"] = fits is not None and len(fits["table"]) == case["fits_for"]["rows"]
    if case.get("reader_refuses"):
        try:
            sf.read_analysis(case["path"])
            result["reader_refused"] = False
        except ValueError:
            result["reader_refused"] = True
    results[case["name"]] = result
    failed = [k for k in ("reader_refused", "fits_as_expected") if result.get(k) is False]
    if result["rules"] != case["expect"] or result["conforms"] != case["conforms"] or failed:
        wrong += 1
        print(f"MISMATCH {case['name']}: expected {case['expect']}, reported {result}")
        print(found.to_string(index=False))
Path("conformance/results-python.json").write_text(json.dumps(results, indent=2) + "\n")
print(f"Python: {len(cases) - wrong} of {len(cases)} cases as expected")
sys.exit(1 if wrong else 0)
