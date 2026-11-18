# Runs the R validators over conformance/cases.json and writes what they report
# for each case to conformance/results-r.json: the rules, and whether the file
# conforms. A file that does not conform and must still be read (S4.4) is then read, and
# must give the same counts as the conforming file it was made from.
#     Rscript verify/check_conformance.R
source("format/setup.R")
source("verify/as_sce.R")
same_counts <- function(a, b) {
    x <- read_as_sce(a)
    y <- read_as_sce(b)
    rhdf5::h5closeAll()   # MuData's reader leaves a handle open
    identical(counts(x), counts(y)) && identical(counts(altExp(x, "gRNA")), counts(altExp(y, "gRNA")))
}
cases <- jsonlite::fromJSON("conformance/cases.json", simplifyVector = FALSE)
wrong <- 0L
results <- list()
for (case in cases) {
    found <- if (case$kind == "screen") validate_screen(case$path) else validate_analysis(case$path, case$screen)
    result <- list(rules = I(sort(unique(found$rule), method = "radix")), conforms = nrow(found) == 0L)
    if (!is.null(case$reads_as)) result$reads_correctly <- same_counts(case$path, case$reads_as)
    if (!is.null(case$conforming_with_index_limit)) {
        # the 64-bit branch of S4: no test file reaches the real limit, so the limit is lowered.
        # Only the modality that passes the lowered limit may stay wide; the others then break S4.
        lowered <- validate_screen(case$path, index_limit = case$conforming_with_index_limit)
        result$wide_branch <- !any(grepl("mod/rna/X", lowered$where, fixed = TRUE))
    }
    if (!is.null(case$fits_used)) {
        # whether a reader takes fits from the cache of this folder (T6.5)
        result$fits_as_expected <- identical(!is.null(fmt_read_analysis(case$path)$tables[[FITS]]), case$fits_used)
    }
    if (!is.null(case$fits_for)) {
        fits <- fmt_read_fits(case$path, case$fits_for$identity)
        result$fits_as_expected <- !is.null(fits) && nrow(fits$table) == case$fits_for$rows
    }
    if (isTRUE(case$reader_refuses)) {
        result$reader_refused <- tryCatch({
            if (case$kind == "screen") check_screen(suppressWarnings(MuData::readH5MU(case$path))) else fmt_read_analysis(case$path)
            FALSE
        }, error = function(e) TRUE)
        rhdf5::h5closeAll()
    }
    if (!is.null(case$copy_keeps_modalities)) {
        copy <- tempfile(fileext = ".h5mu")
        fmt_write_screen(read_as_sce(case$path, case$copy_keeps_modalities), copy)
        modalities <- function(file) sort(subset(rhdf5::h5ls(file, recursive = 2L), group == "/mod")$name)
        same_id <- identical(as.vector(rhdf5::h5read(case$path, "uns/sceptre/id")), as.vector(rhdf5::h5read(copy, "uns/sceptre/id")))
        result$copy_kept_modalities <- identical(modalities(case$path), modalities(copy)) && same_id && nrow(validate_screen(copy)) == 0L
        rhdf5::h5closeAll()
        unlink(copy)
    }
    results[[case$name]] <- result
    expected <- sort(as.character(unlist(case$expect)), method = "radix")
    if (!identical(as.character(result$rules), expected) || !identical(result$conforms, case$conforms) || isFALSE(result$reads_correctly) ||
        isFALSE(result$reader_refused) || isFALSE(result$copy_kept_modalities) || isFALSE(result$wide_branch) || isFALSE(result$fits_as_expected)) {
        wrong <- wrong + 1L
        cat("MISMATCH ", case$name, ": expected ", paste(expected, collapse = ", "), "; reported ",
            paste(result$rules, collapse = ", "), "\n", sep = "")
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
results[["(writer) a retried stage leaves nothing of the first attempt"]] <- list(rules = I(character()), conforms = TRUE, holds = clean)
if (!clean) {
    wrong <- wrong + 1L
    cat("MISMATCH: a retried stage left files of its first attempt\n")
}
# The writer's 64-bit branch, with the limit lowered: both index arrays of every
# modality come out 64-bit, the file conforms under that limit, and it reads back equal.
wide <- tempfile(fileext = ".h5mu")
wide_source <- "conformance/valid/types_written_by_r.h5mu"
fmt_write_screen(read_as_sce(wide_source), wide, index_limit = 1)
listing <- rhdf5::h5ls(wide, all = TRUE)
types <- listing$dtype[grepl("/X$", listing$group) & listing$name %in% c("indices", "indptr")]
holds <- length(types) == 4L && all(grepl("^H5T_STD_I64", types)) && nrow(validate_screen(wide, index_limit = 1)) == 0L && same_counts(wide, wide_source)
rhdf5::h5closeAll()
results[["(writer) past the index limit both index arrays are 64-bit"]] <- list(rules = I(character()), conforms = TRUE, holds = holds)
if (!holds) {
    wrong <- wrong + 1L
    cat("MISMATCH: the writer's 64-bit branch\n")
}
writeLines(jsonlite::toJSON(results, pretty = TRUE, auto_unbox = TRUE), "conformance/results-r.json")
cat("R: ", length(cases) - wrong, " of ", length(cases), " cases as expected\n", sep = "")
if (wrong > 0L) quit(status = 1L)
