# Runs the R validators over conformance/cases.json and writes what they report
# for each case to conformance/results-r.json: the rules, and whether the file
# is conforming, readable or not readable. A "readable" file is then read, and
# must give the same counts as the conforming file it was made from.
#     Rscript check_conformance.R
source("format/setup.R")
same_counts <- function(a, b) {
    x <- fmt_read_screen(a)
    y <- fmt_read_screen(b)
    rhdf5::h5closeAll()   # MuData's reader leaves a handle open
    identical(counts(x), counts(y)) && identical(counts(altExp(x, "gRNA")), counts(altExp(y, "gRNA")))
}
cases <- jsonlite::fromJSON("conformance/cases.json", simplifyVector = FALSE)
wrong <- 0L
results <- list()
for (case in cases) {
    found <- if (case$kind == "screen") validate_screen(case$path) else validate_analysis(case$path, case$screen)
    result <- list(rules = I(sort(unique(found$rule), method = "radix")), level = conformance_level(found))
    if (!is.null(case$reads_as)) result$reads_correctly <- same_counts(case$path, case$reads_as)
    if (!is.null(case$conforming_with_index_limit)) {
        # the 64-bit branch of S4: no test file reaches the real limit, so the limit is lowered.
        # Only the modality that passes the lowered limit may stay wide; the others then break S4.
        lowered <- validate_screen(case$path, index_limit = case$conforming_with_index_limit)
        result$wide_branch <- !any(grepl("mod/rna/X", lowered$where, fixed = TRUE))
    }
    if (isTRUE(case$reader_refuses)) {
        result$reader_refused <- tryCatch({
            if (case$kind == "screen") fmt_read_screen(case$path) else fmt_read_analysis(case$path)
            FALSE
        }, error = function(e) TRUE)
    }
    if (!is.null(case$copy_keeps_modalities)) {
        copy <- tempfile(fileext = ".h5mu")
        fmt_write_screen(fmt_read_screen(case$path, case$copy_keeps_modalities), copy)
        modalities <- function(file) sort(subset(rhdf5::h5ls(file, recursive = 2L), group == "/mod")$name)
        result$copy_kept_modalities <- identical(modalities(case$path), modalities(copy)) && nrow(validate_screen(copy)) == 0L
        rhdf5::h5closeAll()
        unlink(copy)
    }
    results[[case$name]] <- result
    expected <- sort(as.character(unlist(case$expect)), method = "radix")
    if (!identical(as.character(result$rules), expected) || result$level != case$level || isFALSE(result$reads_correctly) ||
        isFALSE(result$reader_refused) || isFALSE(result$copy_kept_modalities) || isFALSE(result$wide_branch)) {
        wrong <- wrong + 1L
        cat("MISMATCH ", case$name, ": expected ", paste(expected, collapse = ", "), " (", case$level, "); reported ",
            paste(result$rules, collapse = ", "), " (", result$level, ")\n", sep = "")
        print(found, row.names = FALSE)
    }
}
# A writer's own behaviour: a stage that did not finish leaves a folder without a
# record; writing the stage again must leave nothing of the first attempt (A2.3).
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
clean <- !file.exists(file.path(retry, "assignment/flagged_cells.parquet")) && nrow(validate_analysis(retry)) == 0L
results[["(writer) a retried stage leaves nothing of the first attempt"]] <- list(rules = I(character()), level = "conforming", holds = clean)
if (!clean) {
    wrong <- wrong + 1L
    cat("MISMATCH: a retried stage left files of its first attempt\n")
}
# The writer's 64-bit branch, with the limit lowered: both index arrays of every
# modality come out 64-bit, the file conforms under that limit, and it reads back equal.
wide <- tempfile(fileext = ".h5mu")
wide_source <- "conformance/valid/types_written_by_r.h5mu"
fmt_write_screen(fmt_read_screen(wide_source), wide, index_limit = 1)
listing <- rhdf5::h5ls(wide, all = TRUE)
types <- listing$dtype[grepl("/X$", listing$group) & listing$name %in% c("indices", "indptr")]
holds <- length(types) == 4L && all(grepl("^H5T_STD_I64", types)) && nrow(validate_screen(wide, index_limit = 1)) == 0L && same_counts(wide, wide_source)
rhdf5::h5closeAll()
results[["(writer) past the index limit both index arrays are 64-bit"]] <- list(rules = I(character()), level = "conforming", holds = holds)
if (!holds) {
    wrong <- wrong + 1L
    cat("MISMATCH: the writer's 64-bit branch\n")
}
writeLines(jsonlite::toJSON(results, pretty = TRUE, auto_unbox = TRUE), "conformance/results-r.json")
cat("R: ", length(cases) - wrong, " of ", length(cases), " cases as expected\n", sep = "")
if (wrong > 0L) quit(status = 1L)
