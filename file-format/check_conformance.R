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
    results[[case$name]] <- result
    expected <- sort(as.character(unlist(case$expect)), method = "radix")
    if (!identical(as.character(result$rules), expected) || result$level != case$level || isFALSE(result$reads_correctly)) {
        wrong <- wrong + 1L
        cat("MISMATCH ", case$name, ": expected ", paste(expected, collapse = ", "), " (", case$level, "); reported ",
            paste(result$rules, collapse = ", "), " (", result$level, ")\n", sep = "")
        print(found, row.names = FALSE)
    }
}
writeLines(jsonlite::toJSON(results, pretty = TRUE, auto_unbox = TRUE), "conformance/results-r.json")
cat("R: ", length(cases) - wrong, " of ", length(cases), " cases as expected\n", sep = "")
if (wrong > 0L) quit(status = 1L)
