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
        isFALSE(result$reader_refused) || isFALSE(result$copy_kept_modalities)) {
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
writeLines(jsonlite::toJSON(results, pretty = TRUE, auto_unbox = TRUE), "conformance/results-r.json")
cat("R: ", length(cases) - wrong, " of ", length(cases), " cases as expected\n", sep = "")
if (wrong > 0L) quit(status = 1L)
