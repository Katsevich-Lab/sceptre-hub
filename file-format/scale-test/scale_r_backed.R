# The disk-backed path the sceptre3 branch uses (R/backed_counts.R): MuData's
# backed reader gives a features x cells DelayedMatrix; one response is read at
# a time with DelayedArray::extract_array(). Checked against direct HDF5 reads.
suppressPackageStartupMessages({ library(SingleCellExperiment); library(Matrix) })
file <- "out/screen-python.h5mu"
out <- list()
timed <- function(name, expr) {
    invisible(gc(reset = TRUE)); start <- Sys.time(); value <- force(expr)
    out[[name]] <<- list(seconds = round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1), peak_gb = round(sum(gc()[, 6]) / 1024, 1))
    cat(name, ":", out[[name]]$seconds, "s,", out[[name]]$peak_gb, "GB\n"); value
}
picked <- unlist(jsonlite::fromJSON("out/picked-features.json")) + 1L
pointers <- as.vector(rhdf5::h5read(file, "mod/rna/X/indptr", bit64conversion = "double"))
direct <- vapply(picked, function(j) {
    if (pointers[j + 1] == pointers[j]) return(0)
    sum(as.vector(rhdf5::h5read(file, "mod/rna/X/data", index = list(seq(pointers[j] + 1, pointers[j + 1])))))
}, numeric(1))
rhdf5::h5closeAll()
warnings_seen <- character()
result <- tryCatch(withCallingHandlers({
    mae <- timed("open with MuData's backed reader", MuData::readH5MU(file, backed = TRUE))
    counts <- SummarizedExperiment::assay(mae[["rna"]], 1L)
    out$matrix <- list(class = class(counts)[1], seed = class(DelayedArray::seed(counts))[1], dim = dim(counts))
    sums <- timed("read 200 responses, one at a time", vapply(picked, function(i) sum(as.numeric(DelayedArray::extract_array(counts, list(i, NULL)))), numeric(1)))
    out$`sums equal direct HDF5 reads` <- identical(as.numeric(sums), as.numeric(direct))
    out$`first five, backed` <- sums[1:5]; out$`first five, direct` <- direct[1:5]
    block <- timed("read the first 500 responses as one in-memory block", as(counts[1:500, ], "dgCMatrix"))
    out$`block nonzeros equal the pointers` <- length(block@x) == pointers[501]
    "ok"
}, warning = function(w) { warnings_seen <<- c(warnings_seen, conditionMessage(w)); invokeRestart("muffleWarning") }),
error = function(e) conditionMessage(e))
out$outcome <- result; out$warnings <- unique(warnings_seen)
jsonlite::write_json(out, "out/results-r-backed.json", auto_unbox = TRUE, pretty = TRUE)
print(jsonlite::toJSON(out, auto_unbox = TRUE, pretty = TRUE))
