# Readers we did not write, R side. The example files are opened with tools
# outside the prototype: HDF5Array for the counts, read from disk, and arrow for
# the tables, beside nanoparquet.
#     Rscript verify/independent_readers.R
if (dir.exists("rlib")) .libPaths(c("rlib", .libPaths()))
suppressPackageStartupMessages(library(Matrix))
source("verify/common.R")
files <- c("example/highmoi/screen.h5mu" = "written by R", "example/lowmoi/screen.h5mu" = "written by R",
    "conformance/valid/highmoi_written_by_python.h5mu" = "written by Python")
for (file in names(files)) {
    for (modality in c("rna", "grna")) {
        # HDF5Array presents the stored cells x features matrix as features x cells, read from disk
        counts <- HDF5Array::H5SparseMatrix(file, paste0("mod/", modality, "/X"))
        shape <- rev(as.integer(rhdf5::h5readAttributes(file, paste0("mod/", modality, "/X"))$shape))
        check(paste0("HDF5Array opens ", modality, " of ", file, " (", files[[file]], ") from disk"), identical(dim(counts), shape), paste(dim(counts), collapse = " x "))
        in_memory <- methods::as(counts, "dgCMatrix")
        stored <- function(name) as.numeric(rhdf5::h5read(file, paste0("mod/", modality, "/obs/", name)))
        check(paste0("HDF5Array's per-cell totals equal the stored n_umis, ", modality, " of ", file), isTRUE(all.equal(unname(Matrix::colSums(in_memory)), stored("n_umis"))))
        check(paste0("HDF5Array's per-cell detected features equal the stored n_nonzero, ", modality, " of ", file), isTRUE(all.equal(unname(Matrix::colSums(in_memory != 0)), stored("n_nonzero"))))
        rhdf5::h5closeAll()
    }
}
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
        # pointed at the folder itself, a dataset reader also meets identity.json
        as_folder <- tryCatch({ arrow::open_dataset(file.path(folder, "fits/response"), format = "parquet")$num_rows; "works" },
            error = function(e) "FAILS (the folder also holds identity.json, which is not Parquet)")
        check(paste0("arrow pointed at the fit cache folder of ", folder, " itself"), as_folder, "an observation: the record of how the fits were made is inside the folder of parts")
    }
}
save_results("independent-readers")
