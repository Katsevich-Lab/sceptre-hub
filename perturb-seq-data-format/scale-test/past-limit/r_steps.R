# R steps of the past-the-limit test.
#   Rscript r_steps.R {validate|hdf5array|mudata|write} ROOT
args <- commandArgs(trailingOnly = TRUE)
step <- args[1]
root <- args[2]
Sys.setenv(PERTURB_SEQ_DATA_FORMAT_DIR = file.path(root, "code"))
source(file.path(root, "code", "format", "setup.R"))
file <- file.path(root, "data", "screen.h5mu")
out <- file.path(root, "results", paste0("r-", step, ".json"))
record <- list(step = step, file = file, checks = list(), timings = list(), started = format(Sys.time()),
               versions = lapply(c(MuData = "MuData", rhdf5 = "rhdf5", HDF5Array = "HDF5Array", Matrix = "Matrix", R = "base"), function(p) as.character(packageVersion(p))))
stamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
say <- function(...) cat("[", stamp(), "] ", ..., "\n", sep = "")
save <- function() {
    record$finished <<- format(Sys.time())
    record$gc_max_mb <<- sum(gc()[, 6])
    jsonlite::write_json(record, out, auto_unbox = TRUE, pretty = TRUE, digits = NA)
}
check <- function(name, ok, detail = "") {
    record$checks[[length(record$checks) + 1L]] <<- list(check = name, ok = isTRUE(ok), detail = paste(detail, collapse = " "))
    say(if (isTRUE(ok)) "PASS " else "FAIL ", name, " ", paste(detail, collapse = " "))
    save()
}
# runs expr, records time, warnings and the exact error; returns the value or NULL
timed <- function(name, expr) {
    say("start ", name)
    warns <- character()
    t <- proc.time()[["elapsed"]]
    value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
    }), error = function(e) structure(list(message = conditionMessage(e), call = paste(deparse(conditionCall(e)), collapse = " ")), class = "step_error"))
    dt <- proc.time()[["elapsed"]] - t
    entry <- list(seconds = round(dt, 1), ok = !inherits(value, "step_error"), warnings = unique(warns))
    if (inherits(value, "step_error")) {
        entry$error <- value$message
        entry$error_call <- value$call
        say("ERROR in ", name, " after ", round(dt, 1), "s: ", value$message, " [call: ", value$call, "]")
        value <- NULL
    } else say("done ", name, " in ", round(dt, 1), "s")
    if (length(warns)) say("warnings: ", paste(unique(warns), collapse = " | "))
    record$timings[[name]] <<- entry
    save()
    value
}
truth <- utils::read.csv(file.path(root, "data", "truth_genes.csv"))
samples_file <- file.path(root, "data", "truth_samples.h5")
picks <- as.integer(rhdf5::h5read(samples_file, "genes", bit64conversion = "double"))   # 0-based
sample_of <- function(g) list(pos = as.numeric(rhdf5::h5read(samples_file, paste0("pos/", g), bit64conversion = "double")),
                              val = as.integer(rhdf5::h5read(samples_file, paste0("val/", g))))

if (step == "validate") {
    found <- timed("validate_screen (real limit)", validate_screen(file))
    if (!is.null(found)) {
        record$problems <- found
        check("R validate_screen reports nothing", nrow(found) == 0L, paste(utils::capture.output(print(found)), collapse = "\n"))
    }
}

if (step == "hdf5array") {
    x <- timed("HDF5Array::H5SparseMatrix open", HDF5Array::H5SparseMatrix(file, "mod/rna/X"))
    if (!is.null(x)) {
        record$h5sparse <- list(class = class(x)[1], dim = dim(x), seed_class = class(DelayedArray::seed(x))[1])
        say("dim ", paste(dim(x), collapse = " x "), " class ", class(x)[1])
        by_col <- ncol(x) == nrow(truth)   # features are columns
        per <- list()
        for (g in picks) {
            t <- proc.time()[["elapsed"]]
            v <- timed(paste0("read gene ", g), {
                if (by_col) as.vector(DelayedArray::extract_array(x, list(NULL, g + 1L))) else as.vector(DelayedArray::extract_array(x, list(g + 1L, NULL)))
            })
            per[[as.character(g)]] <- round(proc.time()[["elapsed"]] - t, 2)
            if (!is.null(v)) {
                s <- sample_of(g)
                nz <- which(v != 0)
                ok <- identical(as.numeric(nz - 1L), s$pos) && identical(as.integer(v[nz]), s$val)
                check(paste0("HDF5Array gene ", g, " equals regenerated truth exactly"), ok,
                      c("nnz", length(nz), "vs", truth$nnz[g + 1L], "sum", sum(as.numeric(v)), "vs", truth$sum[g + 1L]))
            }
        }
        record$per_gene_seconds <- per
        # a block of genes as a sparse object, to see what the R side gives for a slice
        blk <- timed("HDF5Array: 50 genes from the end as SparseArray", {
            j <- (nrow(truth) - 49L):nrow(truth)
            b <- if (by_col) DelayedArray::extract_sparse_array(x, list(NULL, j)) else DelayedArray::extract_sparse_array(x, list(j, NULL))
            list(class = class(b)[1], nnz = SparseArray::nzcount(b), sums = if (by_col) colSums(b) else rowSums(b))
        })
        if (!is.null(blk)) {
            j <- (nrow(truth) - 49L):nrow(truth)
            check("HDF5Array last 50 genes: column sums equal truth", isTRUE(all.equal(unname(as.numeric(blk$sums)), as.numeric(truth$sum[j]), tolerance = 0)),
                  c(blk$class, "nnz", blk$nnz, "truth", sum(as.numeric(truth$nnz[j]))))
        }
        # what happens if one asks HDF5Array for the whole matrix in memory
        whole <- timed("HDF5Array: as(x, 'dgCMatrix') whole matrix", {
            m <- as(x, "dgCMatrix")
            list(class = class(m)[1], nnz = length(m@x))
        })
        if (!is.null(whole)) record$whole_dgCMatrix <- whole
    }
}

if (step == "mudata") {
    mae <- timed("MuData::readH5MU in memory", MuData::readH5MU(file))
    if (!is.null(mae)) {
        rna <- mae[["rna"]]
        a <- SummarizedExperiment::assay(rna, 1L)
        record$mudata <- list(class = class(a)[1], dim = dim(a), nnz = if (methods::is(a, "sparseMatrix")) length(a@x) else NA)
        say("assay class ", class(a)[1], " dim ", paste(dim(a), collapse = " x "))
        timed("compare MuData in-memory to truth", {
            by_row <- nrow(a) == nrow(truth)
            sums <- if (by_row) Matrix::rowSums(a) else Matrix::colSums(a)
            nnz <- if (methods::is(a, "dgCMatrix")) { if (by_row) tabulate(a@i + 1L, nrow(a)) else diff(a@p) } else if (methods::is(a, "dgRMatrix")) { if (by_row) diff(a@p) else tabulate(a@j + 1L, ncol(a)) } else NA
            check("MuData in memory: per-gene sums equal truth", isTRUE(all(as.numeric(sums) == as.numeric(truth$sum))), c("sum total", sum(as.numeric(sums)), "truth", sum(as.numeric(truth$sum))))
            check("MuData in memory: per-gene nnz equal truth", isTRUE(all(as.numeric(nnz) == as.numeric(truth$nnz))), c("total", sum(as.numeric(nnz)), "truth", sum(as.numeric(truth$nnz))))
        })
    }
    rhdf5::h5closeAll()
    # the backed reader, for comparison: one gene at a time
    b <- timed("MuData::readH5MU backed", MuData::readH5MU(file, backed = TRUE))
    if (!is.null(b)) {
        a <- SummarizedExperiment::assay(b[["rna"]], 1L)
        record$mudata_backed <- list(class = class(a)[1], dim = dim(a))
        by_row <- nrow(a) == nrow(truth)
        for (g in picks[c(1L, length(picks) %/% 2L, length(picks))]) {
            v <- timed(paste0("backed gene ", g), as.vector(if (by_row) DelayedArray::extract_array(a, list(g + 1L, NULL)) else DelayedArray::extract_array(a, list(NULL, g + 1L))))
            if (!is.null(v)) {
                s <- sample_of(g)
                nz <- which(v != 0)
                check(paste0("MuData backed gene ", g, " equals truth exactly"), identical(as.numeric(nz - 1L), s$pos) && identical(as.integer(v[nz]), s$val), c("nnz", length(nz)))
            }
        }
    }
}

if (step == "write") {
    # Get the data into R without an in-memory sparse matrix: a disk-backed
    # H5SparseMatrix (features x cells, via t()) as the counts assay, and the
    # small gRNA modality in memory. n_umis / n_nonzero are given so the writer
    # does not need to sum the disk-backed matrix.
    sce <- timed("build SingleCellExperiment with disk-backed counts", {
        x <- HDF5Array::H5SparseMatrix(file, "mod/rna/X")
        counts <- if (nrow(x) == nrow(truth)) x else t(x)    # features x cells
        g <- HDF5Array::H5SparseMatrix(file, "mod/grna/X")
        gcounts <- as(if (ncol(g) == length(as.character(rhdf5::h5read(file, "mod/grna/obs/_index")))) g else t(g), "dgCMatrix")
        cells <- as.character(rhdf5::h5read(file, "mod/rna/obs/_index"))
        genes <- as.character(rhdf5::h5read(file, "mod/rna/var/_index"))
        grnas <- as.character(rhdf5::h5read(file, "mod/grna/var/_index"))
        dimnames(counts) <- list(genes, cells)
        dimnames(gcounts) <- list(grnas, cells)
        gvar <- S4Vectors::DataFrame(grna_target = as.character(rhdf5::h5read(file, "mod/grna/var/grna_target")), grna_type = as.character(rhdf5::h5read(file, "mod/grna/var/grna_type")), row.names = grnas)
        rna_obs <- rhdf5::h5read(file, "mod/rna/obs")
        grna_obs <- rhdf5::h5read(file, "mod/grna/obs")
        cd <- S4Vectors::DataFrame(rna_n_umis = as.numeric(rna_obs$n_umis), rna_n_nonzero = as.integer(rna_obs$n_nonzero),
                                   grna_n_umis = as.numeric(grna_obs$n_umis), grna_n_nonzero = as.integer(grna_obs$n_nonzero), row.names = cells)
        s <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts), colData = cd,
                                                        rowData = S4Vectors::DataFrame(gene_name = paste0("GENE", seq_len(nrow(counts)) - 1L), row.names = genes))
        SingleCellExperiment::altExp(s, "gRNA") <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = gcounts), rowData = gvar)
        S4Vectors::metadata(s) <- list(title = "Synthetic screen past 2^31-1 nonzero counts", perturbation_type = "CRISPRi", organism = "Homo sapiens")
        s
    })
    if (!is.null(sce)) {
        record$assay_class <- class(SummarizedExperiment::assay(sce, "counts"))[1]
        mine <- file.path(root, "data", "screen-r.h5mu")
        unlink(mine)
        w <- timed("fmt_write_screen (real limit), disk-backed counts", fmt_write_screen(sce, mine))
        record$r_file_exists_after <- file.exists(mine)
        if (file.exists(mine)) record$r_file_size <- file.size(mine)
        # where exactly it fails: the coercion fmt_write_screen makes (perturb_seq_data.R line 210)
        timed("the writer's coercion alone: as(drop0(as(counts, 'generalMatrix')), 'RsparseMatrix')",
              methods::as(Matrix::drop0(methods::as(SummarizedExperiment::assay(sce, "counts"), "generalMatrix")), "RsparseMatrix"))
        # and whether any in-memory Matrix sparse class can hold the counts
        timed("as(counts, 'dgCMatrix') (Matrix's in-memory sparse class)", {
            m <- as(SummarizedExperiment::assay(sce, "counts"), "dgCMatrix")
            length(m@x)
        })
    }
}
save()
say("STEP-END ", step)
