# Readers we did not write, R side. The tables of the example analyses are read
# with arrow, beside nanoparquet.
#     Rscript verify/independent_readers.R
if (dir.exists("../perturb-seq-data-format/rlib")) .libPaths(c("../perturb-seq-data-format/rlib", .libPaths()))
source("../perturb-seq-data-format/verify/common.R")
plain <- function(frame) {
    frame <- as.data.frame(frame, stringsAsFactors = FALSE)
    class(frame) <- "data.frame"
    frame[] <- lapply(frame, function(column) if (is.factor(column)) as.character(column) else if (is.list(column)) lapply(column, as.numeric) else as.vector(column))
    frame
}
for (folder in c("example/highmoi/analyses/left_sided", "example/highmoi/analyses/singleton", "example/lowmoi/analyses/default",
    "conformance/valid/analyses/left_sided_written_by_python")) {
    for (path in sort(list.files(folder, pattern = "\\.parquet$", recursive = TRUE, full.names = TRUE))) {
        if (grepl("/fits/", path)) next
        check(paste0("arrow reads ", sub(paste0(dirname(folder), "/"), "", path), " to the same table as nanoparquet"),
            isTRUE(all.equal(plain(arrow::read_parquet(path)), plain(nanoparquet::read_parquet(path)), check.attributes = FALSE)))
    }
    # the fit cache is a folder of part files: arrow reads them as one dataset
    parts <- list.files(file.path(folder, "fits/response"), pattern = "[.]parquet$", full.names = TRUE)
    if (length(parts) > 0L) {
        rows <- sum(vapply(parts, function(p) nanoparquet::read_parquet_info(p)$num_rows, numeric(1)))
        check(paste0("arrow reads the part files of the fit cache of ", folder, " as one dataset"), arrow::open_dataset(parts, format = "parquet")$num_rows == rows)
        # the folder holds part files and nothing else (A13.1), so a dataset reader opens the folder itself
        as_folder <- tryCatch(arrow::open_dataset(file.path(folder, "fits/response"), format = "parquet")$num_rows == rows, error = function(e) FALSE)
        check(paste0("arrow opens the fit cache folder of ", folder, " itself"), as_folder)
        # the identity of the fits is in each part's own metadata, which arrow shows too
        seen <- arrow::read_parquet(parts[1], as_data_frame = FALSE)$metadata[["sceptre_fit_identity"]]
        check(paste0("arrow shows the identity inside a part of the fit cache of ", folder), !is.null(seen) && !is.null(jsonlite::fromJSON(seen)$cells_in_use_sha256))
    }
}
save_results("independent-readers")
