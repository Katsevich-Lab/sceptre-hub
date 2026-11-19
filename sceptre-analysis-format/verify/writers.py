"""What a writer does beyond producing a conforming file, for the Python writer
in format/. verify/writers.R is the same for the R writer.

    python verify/writers.py
"""

import json
import shutil
import sys
import tempfile
from pathlib import Path

import anndata as ad
import h5py
import mudata as md

sys.path.insert(0, "format")
sys.path.insert(0, "verify")
sys.path.insert(0, "../perturb-seq-data-format/verify")
import sceptre_analysis as sf
from common import Results

results = Results("writers")
DATA = "../perturb-seq-data-format"   # the data file format's folder: its conformance files are used here


def same_counts(a, b) -> bool:
    x, y = md.read_h5mu(a), md.read_h5mu(b)
    return all((x[m].X != y[m].X).nnz == 0 for m in y.mod)


# A copy of a data file keeps its `id` and every modality (S10.4, S2).
for case in json.loads(Path(f"{DATA}/conformance/cases.json").read_text()):
    if case.get("copy_keeps_modalities"):
        with tempfile.TemporaryDirectory() as folder:
            copy = Path(folder) / "copy.h5mu"
            sf.write_screen(md.read_h5mu(f"{DATA}/{case['path']}"), copy)
            with h5py.File(f"{DATA}/{case['path']}") as a, h5py.File(copy) as b:
                results.add(f"{case['name']}: a copy keeps the id", ad.io.read_elem(a["uns/dataset_id"]) == ad.io.read_elem(b["uns/dataset_id"]))
                results.add(f"{case['name']}: a copy keeps every modality", sorted(a["mod"]) == sorted(b["mod"]))
            results.add(f"{case['name']}: the copy conforms", len(sf.validate_screen(copy)) == 0)

# A stage that did not finish leaves a folder without a record; writing the stage
# again leaves nothing of the first attempt (A2.3).
with tempfile.TemporaryDirectory() as folder:
    retry = Path(folder) / "analysis"
    shutil.copytree("conformance/valid/analyses/left_sided_written_by_python", retry)
    saved = sf.read_analysis(retry)
    sf.clear_stage(retry, "assignment")
    (retry / "assignment").mkdir()
    shutil.copy("conformance/valid/analyses/left_sided_written_by_python/qc/cells.parquet", retry / "assignment/flagged_cells.parquet")
    record = saved["stages"]["assignment"]
    sf.write_stage(retry, "assignment", {"settings": record["settings"]}, {"assignment/assignments": saved["tables"]["assignment/assignments"]}, record["provenance"])
    results.add("a retried stage leaves nothing of the first attempt", not (retry / "assignment/flagged_cells.parquet").exists())
    results.add("the analysis conforms after the retry", len(sf.validate_analysis(retry)) == 0)

# Past the index limit (lowered here) both index arrays of every modality are
# 64-bit, the file conforms under that limit, and it reads back equal (S4).
with tempfile.TemporaryDirectory() as folder:
    wide = Path(folder) / "wide.h5mu"
    source = f"{DATA}/conformance/valid/types_written_by_python.h5mu"
    sf.write_screen(md.read_h5mu(source), wide, index_limit=1)
    with h5py.File(wide) as f:
        types = {str(f[f"mod/{m}/X/{a}"].dtype) for m in f["mod"] for a in ("indices", "indptr")}
    results.add("past the index limit both index arrays are 64-bit", types == {"int64"}, types)
    results.add("the wide file conforms under the lowered limit", len(sf.validate_screen(wide, index_limit=1)) == 0)
    results.add("the wide file reads back equal", same_counts(wide, source))

results.save()
