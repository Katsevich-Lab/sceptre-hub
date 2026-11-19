# Readers we did not write, R side. The example data files are opened with
# HDF5Array, which reads the counts from disk with none of this folder's code.
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
save_results("independent-readers")
