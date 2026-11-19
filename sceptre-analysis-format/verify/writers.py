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

results.save()
