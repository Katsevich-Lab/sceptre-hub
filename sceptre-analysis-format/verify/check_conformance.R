# Runs the R validator over conformance/cases.json and writes what it reports
# for each analysis folder to conformance/results-r.json: the rules, and whether
# the folder conforms. Each case names the data file the analysis is checked
# against, in the data file format's own conformance set.
#     Rscript verify/check_conformance.R
source("format/setup.R")
cases <- jsonlite::fromJSON("conformance/cases.json", simplifyVector = FALSE)
wrong <- 0L
results <- list()
for (case in cases) {
    found <- validate_analysis(case$path, case$screen)
    result <- list(rules = I(sort(unique(found$rule), method = "radix")), conforms = nrow(found) == 0L)
    if (!is.null(case$fits_used)) {
        # whether a reader takes fits from the cache of this folder (A13.5)
        result$fits_as_expected <- identical(!is.null(fmt_read_analysis(case$path)$tables[[FITS]]), case$fits_used)
    }
    if (!is.null(case$fits_for)) {
        fits <- fmt_read_fits(case$path, case$fits_for$identity)
        result$fits_as_expected <- !is.null(fits) && nrow(fits$table) == case$fits_for$rows
    }
    if (isTRUE(case$reader_refuses)) {
        result$reader_refused <- tryCatch({
            fmt_read_analysis(case$path)
            FALSE
        }, error = function(e) TRUE)
        rhdf5::h5closeAll()
    }
    results[[case$name]] <- result
    expected <- sort(as.character(unlist(case$expect)), method = "radix")
    if (!identical(as.character(result$rules), expected) || !identical(result$conforms, case$conforms) ||
        isFALSE(result$reader_refused) || isFALSE(result$fits_as_expected)) {
        wrong <- wrong + 1L
        cat("MISMATCH ", case$name, ": expected ", paste(expected, collapse = ", "), "; reported ",
            paste(result$rules, collapse = ", "), "\n", sep = "")
        print(found, row.names = FALSE)
    }
}
writeLines(jsonlite::toJSON(results, pretty = TRUE, auto_unbox = TRUE), "conformance/results-r.json")
cat("R: ", length(cases) - wrong, " of ", length(cases), " cases as expected\n", sep = "")
if (wrong > 0L) quit(status = 1L)
