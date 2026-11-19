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
# The parts of the fit cache that carry one identity are merged into one part:
# the same fits, the same identity, and the old parts are gone (A13.4).
merged <- file.path(tempfile(), "analysis")
invisible(dir.create(dirname(merged)))
invisible(file.copy("conformance/valid/analyses/fits_from_two_jobs", dirname(merged), recursive = TRUE))
invisible(file.rename(file.path(dirname(merged), "fits_from_two_jobs"), merged))
before <- fmt_read_fits(merged)
written <- fmt_merge_fits(merged)
after <- fmt_read_fits(merged)
check("merging the fit cache leaves one part", length(written) == 1L && length(fit_parts(merged)) == 1L)
check("the merged part holds the same fits, under the same identity",
    same_identity(before$identity, after$identity) && same_table(rows_sorted(before$table), rows_sorted(after$table)))
check("the analysis conforms after the merge", nrow(validate_analysis(merged)) == 0L)

# A table with a list per row is written in row groups of at most 10,000 rows (A12.7).
rows <- 25000L
table <- data.frame(grna_target = sprintf("t%04d", (seq_len(rows) - 1L) %/% 100L), response_id = sprintf("g%03d", (seq_len(rows) - 1L) %% 100L))
table$z_null <- rep(list(c(0.5, -0.5)), rows)
path <- tempfile(fileext = ".parquet")
write_table(table, path, c("grna_target", "response_id"), name = "discovery/resamples")
groups <- nanoparquet::read_parquet_metadata(path)$row_groups$num_rows
check("resampled statistics are written in row groups of at most 10,000 rows", max(groups) <= 10000L && sum(groups) == rows, max(groups))

save_results("writers")
