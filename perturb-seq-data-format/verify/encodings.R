# Encoding matrix, R side. Writes each kind of value with the standard R
# libraries (MuData, nanoparquet; not the prototype's writer), and reads the
# files Python wrote with its standard libraries, comparing with what was put in.
# See verify/encodings.py.
#     Rscript verify/encodings.R
if (dir.exists("rlib")) .libPaths(c("rlib", .libPaths()))
suppressPackageStartupMessages({ library(SingleCellExperiment); library(Matrix) })
source("verify/common.R")
folder <- "verify/out/encodings"
dir.create(file.path(folder, "r"), recursive = TRUE, showWarnings = FALSE)
cells <- paste0("cell_", 1:6)
kinds <- list(
    double = list(allowed = TRUE, values = c(0.5, NaN, 1.5, Inf, -Inf, 3.5)),
    integer_32 = list(allowed = TRUE, values = 1:6),
    factor_levels_not_alphabetical = list(allowed = TRUE, values = factor(c("b", "a", NA, "a", "b", "a"), levels = c("b", "a"))),
    text_in_feature_table = list(allowed = TRUE, values = c("x", "y")),
    boolean = list(allowed = FALSE, values = c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE)),
    integer_64 = list(allowed = FALSE, values = c(1, 2, 3, 4, 5, 2^40)),
    integer_with_missing = list(allowed = FALSE, values = c(1L, NA, 3L, 4L, NA, 6L)),
    boolean_with_missing = list(allowed = FALSE, values = c(TRUE, NA, TRUE, TRUE, NA, FALSE)),
    text_with_missing_in_feature_table = list(allowed = FALSE, values = c("x", NA)),
    text_covariate = list(allowed = FALSE, values = c("L2", "L1", "L2", "L1", "L2", "L1"))
)
counts <- as(Matrix::sparseMatrix(i = c(1, 1, 2, 2, 1, 2, 1, 2), j = c(1, 2, 3, 4, 4, 6, 6, 3), x = c(3, 1, 2, 1, 5, 4, 2, 2), dims = c(2, 6)), "RsparseMatrix")
data_file <- function(values, path) {
    experiment <- function(ids, row_values = NULL, col_values = NULL) {
        SingleCellExperiment(
            assays = list(counts = counts),
            colData = if (is.null(col_values)) S4Vectors::DataFrame(row.names = cells) else S4Vectors::DataFrame(value = col_values, row.names = cells),
            rowData = if (is.null(row_values)) S4Vectors::DataFrame(row.names = ids) else S4Vectors::DataFrame(value = row_values, row.names = ids)
        )
    }
    feature <- length(values) == 2L
    rna <- experiment(c("gene_1", "gene_2"), if (feature) values, if (!feature) values)
    unlink(path)
    MuData::writeH5MU(MultiAssayExperiment::MultiAssayExperiment(list(rna = rna, grna = experiment(c("grna_1", "grna_2")))), path)
}
for (kind in setdiff(names(kinds), "integer_64")) data_file(kinds[[kind]]$values, file.path(folder, "r", paste0(kind, ".h5mu")))
table <- data.frame(text = c("a", "b", "c"), text_dictionary = factor(c("a", "b", "a")), boolean_with_missing = c(TRUE, NA, FALSE),
    integer_32 = 1:3, double = c(0.5, NA, Inf), stringsAsFactors = FALSE)
table$list_of_doubles <- list(c(1.5, 2.5), numeric(), 3.5)
nanoparquet::write_parquet(table, file.path(folder, "r", "table.parquet"))

# read what Python wrote
# "correct": the same values in the same kind of vector. "same values, read as
# X": the values survive and the kind does not.
same_values <- function(read, put) {
    read <- if (is.factor(read)) as.character(read) else as.vector(read)
    put <- if (is.factor(put)) as.character(put) else put
    if (length(read) != length(put) || is.character(read) != is.character(put)) return(FALSE)
    both_na <- is.na(read) & is.na(put)
    all(both_na | (!is.na(read) & !is.na(put) & read == put)) && (is.character(read) || identical(is.nan(read), is.nan(put)))
}
describe <- function(read, put) {
    kind <- paste0(class(read)[1], if (is.factor(read)) paste0(", levels ", paste(levels(read), collapse = ",")) else "")
    same_kind <- identical(class(read)[1], class(put)[1]) && (!is.factor(put) || identical(levels(read), levels(put)))
    if (!isTRUE(same_values(read, put))) return(paste0("values differ (read as ", kind, ")"))
    if (same_kind) paste0("correct (read as ", kind, ")") else paste0("same values, read as ", kind)
}
expected_file <- "verify/expected_encodings.json"
expected <- if (file.exists(expected_file)) jsonlite::fromJSON(expected_file, simplifyVector = FALSE) else list()
observed <- list()
record <- function(name, allowed, outcome) {
    if (allowed) {
        check(name, startsWith(outcome, "correct"), outcome)
    } else if (!is.null(expected[[name]])) {
        check(name, if (identical(expected[[name]], outcome)) TRUE else paste("CHANGED:", outcome), paste("forbidden kind; as expected:", expected[[name]]))
    } else {
        check(name, outcome, "forbidden kind; no expectation recorded")
    }
    observed[[name]] <<- outcome
}
read_value <- function(path, feature) {
    mae <- suppressWarnings(MuData::readH5MU(path))
    frame <- if (feature) SummarizedExperiment::rowData(mae[["rna"]]) else SummarizedExperiment::colData(mae[["rna"]])
    frame$value
}
for (kind in c(names(kinds), "default_string_settings")) {
    put <- if (kind == "default_string_settings") kinds$double else kinds[[kind]]
    path <- file.path(folder, "python", paste0(kind, ".h5mu"))
    outcome <- tryCatch(describe(read_value(path, length(put$values) == 2L), put$values), error = function(e) "cannot open")
    rhdf5::h5closeAll()
    record(paste0("data file, ", kind, ": Python writes, R reads"), put$allowed && kind != "default_string_settings", outcome)
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
if (length(expected) == 0L) writeLines(jsonlite::toJSON(observed, auto_unbox = TRUE, pretty = TRUE), file.path(folder, "observed-r.json"))
save_results("encodings")
