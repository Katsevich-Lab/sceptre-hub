# The 64-bit branch of S4, R side. verify/wide_index.py wrote the highmoi example
# with both index arrays 64-bit (the writer's limit lowered). R's readers open
# it: the prototype's reader in memory, MuData's backed reader (the path the
# sceptre3 branch uses), and HDF5Array. R's own writer is run the same way. This
# shows the types are handled; it cannot show that values above 2,147,483,647
# are, which needs a file of that size.
#     Rscript verify/wide_index.R
source("format/setup.R")
source("verify/common.R")
wide <- "verify/out/wide/screen-wide.h5mu"
narrow <- "example/highmoi/screen.h5mu"
types_of <- function(file) {
    listing <- rhdf5::h5ls(file, all = TRUE)
    listing$dtype[grepl("/X$", listing$group) & listing$name %in% c("indices", "indptr")]
}
check("the file from Python has 64-bit index arrays", all(grepl("^H5T_STD_I64", types_of(wide))), paste(types_of(wide), collapse = ", "))
caught <- character()
outcome <- function(expr) tryCatch(withCallingHandlers(expr, warning = function(w) {
    caught <<- c(caught, conditionMessage(w))
    invokeRestart("muffleWarning")
}), error = function(e) paste0("FAILS (", conditionMessage(e), ")"))

# 1. the prototype's reader, in memory
a <- outcome(fmt_read_screen(wide))
b <- fmt_read_screen(narrow)
rhdf5::h5closeAll()
ok <- !is.character(a) && identical(SingleCellExperiment::counts(a), SingleCellExperiment::counts(b)) &&
    identical(SingleCellExperiment::counts(SingleCellExperiment::altExp(a, "gRNA")), SingleCellExperiment::counts(SingleCellExperiment::altExp(b, "gRNA")))
check("the prototype's reader returns the same counts as from the 32-bit file", ok, if (is.character(a)) a else "")

# 2. MuData's backed reader: one response at a time, and a block
for (modality in c("rna", "grna")) {
    result <- outcome({
        x <- SummarizedExperiment::assay(MuData::readH5MU(wide, backed = TRUE)[[modality]], 1L)
        y <- SummarizedExperiment::assay(MuData::readH5MU(narrow, backed = TRUE)[[modality]], 1L)
        rows <- seq_len(min(20L, nrow(x)))
        one <- vapply(rows, function(i) identical(as.numeric(DelayedArray::extract_array(x, list(i, NULL))), as.numeric(DelayedArray::extract_array(y, list(i, NULL)))), logical(1))
        all(one) && identical(dim(x), dim(y)) && identical(as.matrix(x), as.matrix(y))
    })
    rhdf5::h5closeAll()
    check(paste0("MuData's backed reader gives the same ", modality, " counts, one feature at a time and whole"), isTRUE(result), if (is.character(result)) result else "")
    # 3. HDF5Array, from disk
    result <- outcome({
        x <- HDF5Array::H5SparseMatrix(wide, paste0("mod/", modality, "/X"))
        y <- HDF5Array::H5SparseMatrix(narrow, paste0("mod/", modality, "/X"))
        identical(as(x, "dgCMatrix"), as(y, "dgCMatrix"))
    })
    check(paste0("HDF5Array gives the same ", modality, " counts"), isTRUE(result), if (is.character(result)) result else "")
}
# warnings the same readers give on the 32-bit file are not about the index type
baseline <- character()
invisible(withCallingHandlers({
    fmt_read_screen(narrow)
    MuData::readH5MU(narrow, backed = TRUE)
}, warning = function(w) {
    baseline <<- c(baseline, conditionMessage(w))
    invokeRestart("muffleWarning")
}))
rhdf5::h5closeAll()
extra <- setdiff(unique(caught), unique(baseline))
check("no reader gave a warning it does not give on the 32-bit file", length(extra) == 0L, paste(extra, collapse = " | "))

# 4. R's writer on the same branch, and its file against Python's
mine <- "verify/out/wide/screen-wide-r.h5mu"
unlink(mine)
fmt_write_screen(fmt_read_screen(narrow), mine, index_limit = 1)
check("R's writer, with the limit lowered, writes 64-bit index arrays", all(grepl("^H5T_STD_I64", types_of(mine))), paste(types_of(mine), collapse = ", "))
check("R's validator passes R's file under the lowered limit", nrow(validate_screen(mine, index_limit = 1)) == 0L)
check("R's validator passes Python's file under the lowered limit", nrow(validate_screen(wide, index_limit = 1)) == 0L)
found <- validate_screen(mine)
check("under the real limit R's file does not conform, breaking S4.2 and S4.3 only", identical(sort(unique(found$rule)), c("S4.2", "S4.3")))
same <- all(vapply(c("mod/rna/X/data", "mod/rna/X/indices", "mod/rna/X/indptr", "mod/grna/X/data", "mod/grna/X/indices", "mod/grna/X/indptr"), function(path) {
    identical(as.vector(rhdf5::h5read(mine, path, bit64conversion = "double")), as.vector(rhdf5::h5read(wide, path, bit64conversion = "double")))
}, logical(1)))
check("R's and Python's wide files hold the same arrays", same)
rhdf5::h5closeAll()
save_results("wide-index")
