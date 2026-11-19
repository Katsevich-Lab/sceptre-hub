# What a writer does beyond producing a conforming file, for the R writer in
# format/. verify/writers.py is the same for the Python writer.
#     Rscript verify/writers.R
source("format/setup.R")
source("verify/as_sce.R")
source("verify/common.R")
same_counts <- function(a, b) {
    x <- read_as_sce(a)
    y <- read_as_sce(b)
    rhdf5::h5closeAll()   # MuData's reader leaves a handle open
    identical(counts(x), counts(y)) && identical(counts(altExp(x, "gRNA")), counts(altExp(y, "gRNA")))
}
# A copy of a data file keeps its `id` and every modality (S10.4, S2).
for (case in jsonlite::fromJSON("conformance/cases.json", simplifyVector = FALSE)) {
    if (is.null(case$copy_keeps_modalities)) next
    copy <- tempfile(fileext = ".h5mu")
    fmt_write_screen(read_as_sce(case$path, case$copy_keeps_modalities), copy)
    modalities <- function(file) sort(subset(rhdf5::h5ls(file, recursive = 2L), group == "/mod")$name)
    check(paste0(case$name, ": a copy keeps the id"), identical(as.vector(rhdf5::h5read(case$path, "uns/sceptre/id")), as.vector(rhdf5::h5read(copy, "uns/sceptre/id"))))
    check(paste0(case$name, ": a copy keeps every modality"), identical(modalities(case$path), modalities(copy)))
    check(paste0(case$name, ": the copy conforms"), nrow(validate_screen(copy)) == 0L)
    rhdf5::h5closeAll()
    unlink(copy)
}
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
# Past the index limit (lowered here) both index arrays of every modality are
# 64-bit, the file conforms under that limit, and it reads back equal (S4).
wide <- tempfile(fileext = ".h5mu")
wide_source <- "conformance/valid/types_written_by_r.h5mu"
fmt_write_screen(read_as_sce(wide_source), wide, index_limit = 1)
listing <- rhdf5::h5ls(wide, all = TRUE)
types <- listing$dtype[grepl("/X$", listing$group) & listing$name %in% c("indices", "indptr")]
check("past the index limit both index arrays are 64-bit", length(types) == 4L && all(grepl("^H5T_STD_I64", types)), paste(types, collapse = " "))
check("the wide file conforms under the lowered limit", nrow(validate_screen(wide, index_limit = 1)) == 0L)
check("the wide file reads back equal", same_counts(wide, wide_source))
rhdf5::h5closeAll()
save_results("writers")
