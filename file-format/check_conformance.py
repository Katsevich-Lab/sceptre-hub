"""Runs the Python validators over conformance/cases.json and writes what they
report for each case to conformance/results-python.json: the rules, and whether
the file is conforming, readable or not readable. A "readable" file is then read,
and must give the same counts as the conforming file it was made from.

    python check_conformance.py
"""

import json
import shutil
import sys
import tempfile
from pathlib import Path

import h5py

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
    if "conforming_with_index_limit" in case:
        # the 64-bit branch of S4: no test file reaches the real limit, so the limit is lowered.
        # Only the modality that passes the lowered limit may stay wide; the others then break S4.
        lowered = sf.validate_screen(case["path"], index_limit=case["conforming_with_index_limit"])
        result["wide_branch"] = not any("mod/rna/X" in w for w in lowered["where"])
    if case.get("reader_refuses"):
        try:
            sf.read_screen(case["path"]) if case["kind"] == "screen" else sf.read_analysis(case["path"])
            result["reader_refused"] = False
        except ValueError:
            result["reader_refused"] = True
    if case.get("copy_keeps_modalities"):
        with tempfile.TemporaryDirectory() as folder:
            copy = Path(folder) / "copy.h5mu"
            sf.write_screen(sf.read_screen(case["path"], case["copy_keeps_modalities"]), copy)
            with h5py.File(case["path"]) as a, h5py.File(copy) as b:
                result["copy_kept_modalities"] = sorted(a["mod"]) == sorted(b["mod"]) and len(sf.validate_screen(copy)) == 0
    results[case["name"]] = result
    failed = [k for k in ("reads_correctly", "reader_refused", "copy_kept_modalities", "wide_branch") if result.get(k) is False]
    if result["rules"] != case["expect"] or result["level"] != case["level"] or failed:
        wrong += 1
        print(f"MISMATCH {case['name']}: expected {case['expect']} ({case['level']}), reported {result}")
        print(found.to_string(index=False))
# A writer's own behaviour: a stage that did not finish leaves a folder without a
# record; writing the stage again must leave nothing of the first attempt (A2.3).
with tempfile.TemporaryDirectory() as folder:
    retry = Path(folder) / "analysis"
    shutil.copytree("conformance/valid/analyses/left_sided_written_by_python", retry)
    saved = sf.read_analysis(retry)
    sf.clear_stage(retry, "assignment")
    (retry / "assignment").mkdir()
    shutil.copy("conformance/valid/analyses/left_sided_written_by_python/qc/cells.parquet", retry / "assignment/flagged_cells.parquet")
    record = saved["stages"]["assignment"]
    sf.write_stage(retry, "assignment", {"settings": record["settings"]}, {"assignment/assignments": saved["tables"]["assignment/assignments"]}, record["provenance"])
    clean = not (retry / "assignment/flagged_cells.parquet").exists() and len(sf.validate_analysis(retry)) == 0
    results["(writer) a retried stage leaves nothing of the first attempt"] = {"rules": [], "level": "conforming", "holds": clean}
    if not clean:
        wrong += 1
        print("MISMATCH: a retried stage left files of its first attempt")
# The writer's 64-bit branch, with the limit lowered: both index arrays of every
# modality come out 64-bit, the file conforms under that limit, and it reads back equal.
with tempfile.TemporaryDirectory() as folder:
    wide = Path(folder) / "wide.h5mu"
    source = "conformance/valid/types_written_by_python.h5mu"
    sf.write_screen(sf.read_screen(source), wide, index_limit=1)
    with h5py.File(wide) as f:
        types = {str(f[f"mod/{m}/X/{a}"].dtype) for m in f["mod"] for a in ("indices", "indptr")}
    holds = types == {"int64"} and len(sf.validate_screen(wide, index_limit=1)) == 0 and same_counts(wide, source)
    results["(writer) past the index limit both index arrays are 64-bit"] = {"rules": [], "level": "conforming", "holds": holds}
    if not holds:
        wrong += 1
        print("MISMATCH: the writer's 64-bit branch")
Path("conformance/results-python.json").write_text(json.dumps(results, indent=2) + "\n")
print(f"Python: {len(cases) - wrong} of {len(cases)} cases as expected")
sys.exit(1 if wrong else 0)
