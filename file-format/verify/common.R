# Shared by the package-free checks: recording results, and making sure that
# the sceptre package was not loaded.
dir.create("verify/out/results", recursive = TRUE, showWarnings = FALSE)
results <- list()
# `ok` is TRUE (pass), FALSE (fail) or text (an observation, not a verdict)
check <- function(name, ok, detail = "") {
    result <- if (isTRUE(ok)) "pass" else if (identical(ok, FALSE)) "fail" else as.character(ok)
    results[[length(results) + 1L]] <<- list(check = name, result = result, detail = as.character(detail))
    invisible(result)
}
save_results <- function(name) {
    check("neither sceptre nor pysceptre was loaded", !"sceptre" %in% loadedNamespaces())
    writeLines(jsonlite::toJSON(results, auto_unbox = TRUE, pretty = TRUE), file.path("verify/out/results", paste0(name, "-r.json")))
    failed <- Filter(function(r) r$result == "fail", results)
    cat(name, " (R): ", length(results) - length(failed), " of ", length(results), " checks pass\n", sep = "")
    for (r in failed) cat("  FAIL ", r$check, " - ", r$detail, "\n", sep = "")
}
