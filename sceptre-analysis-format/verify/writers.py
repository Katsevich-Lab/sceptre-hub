"""What a writer does beyond producing a conforming analysis, for the Python
writer in format/. verify/writers.R is the same for the R writer.

    python verify/writers.py
"""

import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, "format")
sys.path.insert(0, "../perturb-seq-data-format/verify")
import sceptre_analysis as sf
from common import Results

results = Results("writers")

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

# The parts of the fit cache that carry one identity are merged into one part:
# the same fits, the same identity, and the old parts are gone (T6.4).
with tempfile.TemporaryDirectory() as folder:
    merged = Path(folder) / "analysis"
    shutil.copytree("conformance/valid/analyses/fits_from_two_jobs", merged)
    before = sf.read_fits(merged)
    written = sf.merge_fits(merged)
    after = sf.read_fits(merged)
    results.add("merging the fit cache leaves one part", len(written) == 1 and len(sf.fit_parts(merged)) == 1)
    results.add("the merged part holds the same fits, under the same identity",
                sf.same_identity(before["identity"], after["identity"]) and sf.same_table(sf.rows_sorted(before["table"]), sf.rows_sorted(after["table"])))
    results.add("the analysis conforms after the merge", len(sf.validate_analysis(merged)) == 0)

# A table with a list per row is written in row groups of at most 10,000 rows (T5.7).
with tempfile.TemporaryDirectory() as folder:
    import pandas as pd
    import pyarrow.parquet as pq
    rows = 25_000
    table = pd.DataFrame({"grna_target": [f"t{i // 100:04d}" for i in range(rows)], "response_id": [f"g{i % 100:03d}" for i in range(rows)],
                          "z_null": [[0.5, -0.5]] * rows})
    sf.write_table(table, Path(folder) / "resamples.parquet", ["grna_target", "response_id"], table_name="discovery/resamples")
    file = pq.ParquetFile(Path(folder) / "resamples.parquet")
    largest = max(file.metadata.row_group(i).num_rows for i in range(file.metadata.num_row_groups))
    results.add("resampled statistics are written in row groups of at most 10,000 rows", largest <= 10_000 and file.metadata.num_rows == rows, largest)

results.save()
