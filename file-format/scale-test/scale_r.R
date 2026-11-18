# Scale test, R side: validate the Python-written data file, read it in memory
# and feature by feature from disk, write it again from R, and read and rewrite
# the assignment table.
suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(Matrix)
})
Sys.setenv(SCEPTRE_FORMAT_DIR = getwd())
source("format/sceptre_format.R")
results <- list()
step <- function(name, expr) {
    invisible(gc(reset = TRUE))
    start <- Sys.time()
    value <- force(expr)
    results[[name]] <<- list(seconds = round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1), peak_gb = round(sum(gc()[, 6]) / 1024, 1))
    cat(name, ":", results[[name]]$seconds, "s,", results[[name]]$peak_gb, "GB\n")
    jsonlite::write_json(results, "out/results-r.json", auto_unbox = TRUE, pretty = TRUE)
    value
}
file <- "out/screen-python.h5mu"
python <- jsonlite::fromJSON("out/results-python.json", simplifyVector = FALSE)$dataset

problems <- step("validate the data file", validate_screen(file))
results$validation <- list(level = conformance_level(problems), problems = problems)

picked <- unlist(jsonlite::fromJSON("out/picked-features.json")) + 1L
by_feature <- step("read 200 features from disk, one at a time (rhdf5)", {
    pointers <- as.vector(rhdf5::h5read(file, "mod/rna/X/indptr", bit64conversion = "double"))
    lapply(picked, function(j) {
        if (pointers[j + 1] == pointers[j]) return(list(cells = integer(), counts = integer()))
        rows <- list(seq(pointers[j] + 1, pointers[j + 1]))
        list(cells = as.vector(rhdf5::h5read(file, "mod/rna/X/indices", index = rows)), counts = as.vector(rhdf5::h5read(file, "mod/rna/X/data", index = rows)))
    })
})
rhdf5::h5closeAll()
backed_warnings <- character()
backed <- tryCatch(withCallingHandlers(step("open disk-backed (HDF5Array) and read 200 features", {
    matrix <- HDF5Array::H5SparseMatrix(file, "mod/rna/X")
    # the file is cells x features; HDF5Array may present it either way round
    by_column <- ncol(matrix) == python$responses
    part <- as(if (by_column) matrix[, picked] else matrix[picked, ], "dgCMatrix")
    list(dim = dim(matrix), class = class(matrix), sums = if (by_column) Matrix::colSums(part) else Matrix::rowSums(part))
}), warning = function(w) { backed_warnings <<- c(backed_warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
error = function(e) list(error = conditionMessage(e)))
results$`disk-backed` <- list(dim = backed$dim, class = backed$class, error = backed$error, warnings = unique(backed_warnings))
rhdf5::h5closeAll()

sce <- step("read the data file (in memory)", fmt_read_screen(file))
rhdf5::h5closeAll()
counts <- counts(sce)
grna <- altExp(sce, "gRNA")
results$`read back` <- list(
    `same shape and nonzeros` = identical(as.numeric(c(ncol(sce), nrow(sce), nrow(grna), length(counts@x), length(counts(grna)@x))),
        as.numeric(c(python$cells, python$responses, python$grnas, python$rna_nonzeros, python$grna_nonzeros))),
    `same count sums` = sum(counts@x) == python$rna_count_sum && sum(counts(grna)@x) == python$grna_count_sum,
    `200 features from disk equal the in-memory ones` = all(vapply(seq_along(picked), function(i) {
        column <- counts[picked[i], ]
        cells <- which(column != 0)
        identical(as.numeric(by_feature[[i]]$cells) + 1, as.numeric(cells)) && identical(as.numeric(by_feature[[i]]$counts), as.numeric(column[cells]))
    }, logical(1))),
    `disk-backed sums equal the in-memory ones` = if (is.null(backed$error)) isTRUE(all.equal(unname(backed$sums), unname(Matrix::rowSums(counts[picked, ])))) else backed$error,
    `disk-backed sums, first five` = if (is.null(backed$error)) unname(backed$sums)[1:5],
    `in-memory sums, first five` = unname(Matrix::rowSums(counts[picked[1:5], ])),
    covariates = paste(names(colData(sce)), collapse = ", ")
)
rm(counts); invisible(gc())

unlink("out/screen-r.h5mu")
step("write the data file", fmt_write_screen(sce, "out/screen-r.h5mu"))
rhdf5::h5closeAll()
results$`R-written file` <- list(gb = round(file.size("out/screen-r.h5mu") / 1e9, 2))
rm(sce); invisible(gc())

table <- step("read the assignment table", read_table("out/assignments-python.parquet"))
results$assignments <- list(rows = nrow(table))
step("write the assignment table", write_table(table, "out/assignments-r.parquet", c("grna_id", "cell_id")))
jsonlite::write_json(results, "out/results-r.json", auto_unbox = TRUE, pretty = TRUE)
cat("R side done\n")
