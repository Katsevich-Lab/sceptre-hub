# Two writers, one file: the R half. Reads the random inputs that
# verify/random_inputs.py made, writes them with the R library, checks both
# languages' files with the R validator, and reads Python's files back to
# compare them with the source. Also a plain round trip of an analysis folder.
#     Rscript verify/two_writers.R
source("format/setup.R")
source("verify/as_sce.R")
source("verify/common.R")
formula <- "~log(rna_n_umis) + batch"
settings <- list(formula = formula, side = "both", resampling_mechanism = "permutations", resampling_approximation = "spa",
    control_group = "complement", grna_integration_strategy = "union", B1 = 499L, B2 = 4999L, B3 = 0L,
    multiple_testing_method = "BH", multiple_testing_alpha = 0.1, seed = 4L, response_fit_method = "batched", grna_fit_method = "pysceptre")
missing_to <- function(x, na) vapply(x, function(v) if (is.null(v)) na else v, na)

for (folder in sort(Sys.glob("verify/out/random/seed-*"))) {
    seed <- basename(folder)
    inputs <- jsonlite::fromJSON(file.path(folder, "input.json"), simplifyVector = FALSE)
    cells <- unlist(inputs$cells); genes <- unlist(inputs$genes); grnas <- unlist(inputs$grnas)
    count_matrix <- function(triplets, ids) {
        t <- matrix(unlist(triplets), ncol = 3L, byrow = TRUE)
        Matrix::sparseMatrix(i = t[, 2] + 1L, j = t[, 1] + 1L, x = t[, 3], dims = c(length(ids), length(cells)), dimnames = list(ids, cells))
    }
    sce <- SingleCellExperiment(assays = list(counts = count_matrix(inputs$rna, genes)), colData = S4Vectors::DataFrame(
        batch = factor(unlist(inputs$batch$values), levels = unlist(inputs$batch$levels)),
        quality = missing_to(inputs$quality, NaN), row.names = cells
    ))
    altExp(sce, "gRNA") <- SingleCellExperiment(assays = list(counts = count_matrix(inputs$grna, grnas)), rowData = S4Vectors::DataFrame(
        grna_target = unlist(inputs$grna_target), target_chr = missing_to(inputs$target_chr, NA_character_),
        target_start = missing_to(inputs$target_start, NA_integer_), target_end = missing_to(inputs$target_end, NA_integer_), row.names = grnas
    ))
    metadata(sce) <- list(sceptre = list(moi = inputs$moi, id = inputs$id, genome_build = inputs$genome_build))
    screen <- file.path(folder, "screen-r.h5mu")
    unlink(screen)
    fmt_write_screen(sce, screen)
    rhdf5::h5closeAll()

    analysis <- file.path(folder, "analysis-r")
    unlink(analysis, recursive = TRUE)
    provenance <- list(package = "verify", version = "0", language = "R", written_at = "2026-01-01T00:00:00Z")
    back <- read_as_sce(screen)
    rhdf5::h5closeAll()
    fmt_init_analysis(analysis, screen_link(back), "rna")
    pairs <- data.frame(grna_target = vapply(inputs$pairs, `[[`, "", 1L), response_id = vapply(inputs$pairs, `[[`, "", 2L))
    results_table <- data.frame(
        grna_target = vapply(inputs$results, function(r) r$grna_target, ""), response_id = vapply(inputs$results, function(r) r$response_id, ""),
        pass_qc = vapply(inputs$results, function(r) r$pass_qc, logical(1)), p_value = missing_to(lapply(inputs$results, `[[`, "p_value"), NA_real_),
        fold_change = missing_to(lapply(inputs$results, `[[`, "fold_change"), NA_real_), se_fold_change = missing_to(lapply(inputs$results, `[[`, "se_fold_change"), NA_real_),
        log_2_fold_change = missing_to(lapply(inputs$results, `[[`, "log_2_fold_change"), NA_real_)
    )
    fmt_write_stage(analysis, "parameters", list(settings = settings, design = list(terms = I(c("(Intercept)", "log(rna_n_umis)", "batchb1")))),
        list("parameters/discovery_pairs" = pairs, "parameters/positive_control_pairs" = pairs[0, ]), provenance)
    fmt_write_stage(analysis, "assignment", list(settings = list(method = "thresholding", hyperparameters = list(threshold = 5L))),
        list("assignment/assignments" = data.frame(grna_id = vapply(inputs$assignments, `[[`, "", 1L), cell_id = vapply(inputs$assignments, `[[`, "", 2L))), provenance)
    fmt_write_stage(analysis, "qc", list(settings = list(response_n_umis_range = c(0.01, 0.99), response_n_nonzero_range = c(0.01, 0.99), p_mito_threshold = 0.2, pairwise_qc = TRUE, B3 = 0L)),
        list("qc/cells" = data.frame(cell_id = cells, in_use = unlist(inputs$in_use)),
            "qc/discovery_pair_qc" = results_table[c("grna_target", "response_id", "pass_qc")],
            "qc/positive_control_pair_qc" = results_table[0, c("grna_target", "response_id", "pass_qc")]), provenance)
    called <- results_table$p_value[!is.na(results_table$p_value) & results_table$p_value < 0.2]
    fmt_write_stage(analysis, "discovery", list(settings = stats::setNames(list(), character()), significance_threshold = if (length(called)) max(called)),
        list("discovery/results" = results_table), provenance)
    fits <- data.frame(response_id = vapply(inputs$fits, function(f) f$response_id, ""), theta = vapply(inputs$fits, function(f) f$theta, numeric(1)))
    fits$coefficients <- lapply(inputs$fits, function(f) as.numeric(unlist(f$coefficients)))
    fmt_write_fits(analysis, fits[c("response_id", "coefficients", "theta")], list(screen_id = inputs$id, formula = formula,
        cells_in_use_sha256 = sha256_ids(cells[unlist(inputs$in_use)]), response_fit_method = "batched", control_group = "complement"))

    for (language in c("r", "python")) {
        file <- file.path(folder, paste0("screen-", language, ".h5mu"))
        problems <- validate_screen(file)
        check(paste0(seed, ": R validator passes the file ", language, " wrote"), nrow(problems) == 0L, paste(problems$rule, problems$message, collapse = "; "))
        problems <- validate_analysis(file.path(folder, paste0("analysis-", language)), file)
        check(paste0(seed, ": R validator passes the analysis ", language, " wrote"), nrow(problems) == 0L, paste(problems$rule, problems$message, collapse = "; "))
    }

    # Python's file read back and compared with the source
    python <- read_as_sce(file.path(folder, "screen-python.h5mu"))
    rhdf5::h5closeAll()
    triplets <- function(matrix) {
        t <- methods::as(methods::as(matrix, "generalMatrix"), "TsparseMatrix")
        out <- cbind(t@j, t@i, t@x)
        out[order(out[, 1], out[, 2]), , drop = FALSE]
    }
    source_triplets <- function(x) {
        t <- matrix(unlist(x), ncol = 3L, byrow = TRUE)
        t[order(t[, 1], t[, 2]), , drop = FALSE]
    }
    check(paste0(seed, ": rna counts read from Python's file equal the source"), isTRUE(all.equal(unname(triplets(counts(python))), unname(source_triplets(inputs$rna)))))
    check(paste0(seed, ": grna counts read from Python's file equal the source"), isTRUE(all.equal(unname(triplets(counts(altExp(python, "gRNA")))), unname(source_triplets(inputs$grna)))))
    check(paste0(seed, ": identifiers read from Python's file equal the source"), identical(colnames(python), cells) && identical(rownames(python), genes) && identical(rownames(altExp(python, "gRNA")), grnas))
    check(paste0(seed, ": factor levels and values read from Python's file equal the source"),
        identical(levels(python$batch), unlist(inputs$batch$levels)) && identical(as.character(python$batch), unlist(inputs$batch$values)))
    quality <- missing_to(inputs$quality, NA_real_)
    check(paste0(seed, ": doubles with missing values read from Python's file equal the source"), identical(is.na(python$quality), is.na(quality)) && identical(python$quality[!is.na(quality)], quality[!is.na(quality)]))
    table <- rowData(altExp(python, "gRNA"))
    check(paste0(seed, ": gRNA table read from Python's file equals the source"),
        identical(as.character(table$grna_target), unlist(inputs$grna_target)) && isTRUE(all.equal(as.numeric(table$target_start), as.numeric(missing_to(inputs$target_start, NA_integer_)))) &&
            identical(as.character(table$target_chr), missing_to(inputs$target_chr, NA_character_)))
    saved <- fmt_read_analysis(file.path(folder, "analysis-python"))
    got <- saved$tables[["discovery/results"]]
    want <- results_table[byte_order(results_table, c("grna_target", "response_id")), ]
    equal <- all(vapply(names(want), function(column) identical(is.na(got[[column]]), is.na(want[[column]])) && identical(got[[column]][!is.na(want[[column]])], want[[column]][!is.na(want[[column]])]), logical(1)))
    check(paste0(seed, ": results read from Python's analysis equal the source, to the last bit"), equal)
}

# A plain round trip of an analysis folder in R: read it, write every stage and
# the fits again to a new folder, and compare.
source_folder <- "example/highmoi/analyses/left_sided"
saved <- fmt_read_analysis(source_folder)
copy <- file.path("verify/out", "analysis-rewritten-by-r")
unlink(copy, recursive = TRUE)
fmt_init_analysis(copy, saved$screen, saved$response_modality)
for (stage in names(saved$stages)) {
    record <- saved$stages[[stage]]
    names_listed <- paste0(stage, "/", sub("\\.parquet$", "", vapply(record$files, function(f) f$path, character(1))))
    fmt_write_stage(copy, stage, record[setdiff(names(record), c("stage", "provenance", "files"))], saved$tables[names_listed], record$provenance)
}
fmt_write_fits(copy, saved$tables[[FITS]], saved$fits_identity)
again <- fmt_read_analysis(copy)
check("R: an analysis folder read and written again has identical tables", identical(saved$tables[sort(names(saved$tables))], again$tables[sort(names(again$tables))]))
check("R: an analysis folder read and written again has equal stage records", isTRUE(all.equal(saved$stages, again$stages)))
check("R: the rewritten analysis folder conforms", nrow(validate_analysis(copy, "example/highmoi/screen.h5mu")) == 0L)
save_results("two-writers")
