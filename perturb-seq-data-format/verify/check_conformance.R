# Runs the R validator over conformance/cases.json and writes what it reports
# for each data file to conformance/results-r.json: the rules, and whether the
# file conforms. A file that does not conform and must still be read (D4.4) is
# then read, and must give the same counts as the conforming file it was made from.
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
    found <- validate_screen(case$path)
    result <- list(rules = I(sort(unique(found$rule), method = "radix")), conforms = nrow(found) == 0L)
    if (!is.null(case$reads_as)) result$reads_correctly <- same_counts(case$path, case$reads_as)
    if (!is.null(case$conforming_with_index_limit)) {
        # the 64-bit branch of D4: no test file reaches the real limit, so the limit is lowered.
        # Only the modality that passes the lowered limit may stay wide; the others then break D4.
        lowered <- validate_screen(case$path, index_limit = case$conforming_with_index_limit)
        result$wide_branch <- !any(grepl("mod/rna/X", lowered$where, fixed = TRUE))
    }
    if (isTRUE(case$reader_refuses)) {
        result$reader_refused <- tryCatch({
            check_screen(suppressWarnings(MuData::readH5MU(case$path)))
            FALSE
        }, error = function(e) TRUE)
        rhdf5::h5closeAll()
    }
    results[[case$name]] <- result
    expected <- sort(as.character(unlist(case$expect)), method = "radix")
    if (!identical(as.character(result$rules), expected) || !identical(result$conforms, case$conforms) || isFALSE(result$reads_correctly) ||
        isFALSE(result$reader_refused) || isFALSE(result$wide_branch)) {
        wrong <- wrong + 1L
        cat("MISMATCH ", case$name, ": expected ", paste(expected, collapse = ", "), "; reported ",
            paste(result$rules, collapse = ", "), "\n", sep = "")
        print(found, row.names = FALSE)
    }
}
writeLines(jsonlite::toJSON(results, pretty = TRUE, auto_unbox = TRUE), "conformance/results-r.json")
cat("R: ", length(cases) - wrong, " of ", length(cases), " cases as expected\n", sep = "")
if (wrong > 0L) quit(status = 1L)
