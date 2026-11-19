# What a writer does beyond producing a conforming analysis, for the R writer in
# format/. verify/writers.py is the same for the Python writer.
#     Rscript verify/writers.R
source("format/setup.R")
source("../perturb-seq-data-format/verify/common.R")
# A stage that did not finish leaves a folder without a record; writing the stage
# again leaves nothing of the first attempt (A2.3).
retry <- file.path(tempfile(), "analysis")
invisible(dir.create(dirname(retry)))
source_analysis <- "conformance/valid/analyses/left_sided_written_by_r"
invisible(file.copy(source_analysis, dirname(retry), recursive = TRUE))
invisible(file.rename(file.path(dirname(retry), basename(source_analysis)), retry))
saved <- fmt_read_analysis(retry)
fmt_clear_stage(retry, "assignment")
invisible(dir.create(file.path(retry, "assignment")))
invisible(file.copy(file.path(source_analysis, "qc/cells.parquet"), file.path(retry, "assignment/flagged_cells.parquet")))
record <- saved$stages$assignment
fmt_write_stage(retry, "assignment", list(settings = record$settings), list("assignment/assignments" = saved$tables[["assignment/assignments"]]), record$provenance)
check("a retried stage leaves nothing of the first attempt", !file.exists(file.path(retry, "assignment/flagged_cells.parquet")))
check("the analysis conforms after the retry", nrow(validate_analysis(retry)) == 0L)
save_results("writers")
