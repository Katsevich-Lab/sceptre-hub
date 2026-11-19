# Encoding matrix for the tables, R side. Writes one Parquet file holding each
# kind of column with nanoparquet (not the prototype's writer), and reads the
# file Python wrote with pyarrow, comparing with what was put in. See
# verify/encodings.py.
#     Rscript verify/encodings.R
if (dir.exists("../perturb-seq-data-format/rlib")) .libPaths(c("../perturb-seq-data-format/rlib", .libPaths()))
source("../perturb-seq-data-format/verify/common.R")
folder <- "verify/out/encodings"
dir.create(file.path(folder, "r"), recursive = TRUE, showWarnings = FALSE)
table <- data.frame(text = c("a", "b", "c"), text_dictionary = factor(c("a", "b", "a")), boolean_with_missing = c(TRUE, NA, FALSE),
    integer_32 = 1:3, double = c(0.5, NA, Inf), stringsAsFactors = FALSE)
table$list_of_doubles <- list(c(1.5, 2.5), numeric(), 3.5)
nanoparquet::write_parquet(table, file.path(folder, "r", "table.parquet"))

# read what Python wrote
expected <- jsonlite::fromJSON("verify/expected_encodings.json", simplifyVector = FALSE)
record <- function(name, allowed, outcome) {
    if (allowed) {
        check(name, startsWith(outcome, "correct"), outcome)
    } else if (!is.null(expected[[name]])) {
        check(name, if (identical(expected[[name]], outcome)) TRUE else paste("CHANGED:", outcome), paste("forbidden kind; as expected:", expected[[name]]))
    } else {
        check(name, outcome, "forbidden kind; no expectation recorded")
    }
}
put <- list(text = c("a", "b", "c"), text_dictionary = c("a", "b", "a"), boolean_with_missing = c(TRUE, NA, FALSE), integer_32 = 1:3,
    integer_64 = c(1, 2, 2^40), double = c(0.5, NA, Inf))
readers <- list(
    nanoparquet = function(p) as.data.frame(nanoparquet::read_parquet(p)),
    arrow = function(p) as.data.frame(arrow::read_parquet(p))
)
for (reader in names(readers)) {
    frame <- tryCatch(readers[[reader]](file.path(folder, "python", "table.parquet")), error = function(e) NULL)
    if (is.null(frame)) {
        check(paste0("table: Python (pyarrow) writes, ", reader, " reads"), FALSE, "cannot read")
        next
    }
    for (name in names(put)) {
        got <- frame[[name]]
        plain <- if (is.factor(got)) as.character(got) else if (inherits(got, "integer64")) as.numeric(got) else got
        ok <- length(plain) == length(put[[name]]) && all((is.na(plain) & is.na(put[[name]])) | (!is.na(plain) & !is.na(put[[name]]) & plain == put[[name]]))
        if (name == "integer_64") record(paste0("table, integer_64: Python (pyarrow) writes, ", reader, " reads"), FALSE, paste0(if (ok) "correct" else "values differ", " (read as ", class(got)[1], ")"))
        else check(paste0("table, ", name, ": Python (pyarrow) writes, ", reader, " reads"), isTRUE(ok), paste("read as", class(got)[1]))
    }
    lists <- lapply(frame$list_of_doubles, as.numeric)
    check(paste0("table, list_of_doubles: Python (pyarrow) writes, ", reader, " reads"), identical(lists, list(c(1.5, 2.5), numeric(), 3.5)), paste("read as", class(frame$list_of_doubles)[1]))
}
save_results("encodings")
