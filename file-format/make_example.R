# Builds example/: two datasets and their analyses, written by R, and
# example-objects.rds, the sceptre_objects they were written from.
# Run from this folder: Rscript make_example.R
source("format/setup.R")
quiet <- function(object, ..., through = "discovery") {
    object <- object |> assign_grnas(parallel = FALSE, ...) |> run_qc()
    if (through == "qc") return(object)
    object |>
        run_calibration_check(parallel = FALSE, print_progress = FALSE) |>
        run_power_check(parallel = FALSE, print_progress = FALSE) |>
        run_discovery_analysis(parallel = FALSE, print_progress = FALSE)
}

# ---- high MOI: 500 cells, 100 responses, 60 gRNAs
data(highmoi_example_data)
data(grna_target_data_frame_highmoi)
high <- import_data(
    response_matrix = highmoi_example_data$response_matrix,
    grna_matrix = highmoi_example_data$grna_matrix,
    grna_target_data_frame = grna_target_data_frame_highmoi,
    moi = "high",
    extra_covariates = highmoi_example_data$extra_covariates,
    response_names = highmoi_example_data$gene_names
)
positive_control_pairs <- construct_positive_control_pairs(high)
discovery_pairs <- construct_cis_pairs(high, positive_control_pairs = positive_control_pairs, distance_threshold = 5e6)
parameters <- function(...) {
    set_analysis_parameters(high, discovery_pairs = discovery_pairs, positive_control_pairs = positive_control_pairs, ...)
}
highmoi <- list(
    imported = high,
    left_sided = parameters(side = "left") |> quiet(method = "thresholding"),
    two_sided = parameters(side = "both") |> quiet(method = "thresholding"),
    # each gRNA tested on its own instead of pooled by target
    singleton = parameters(side = "left", grna_integration_strategy = "singleton") |> quiet(method = "thresholding"),
    # each gRNA tested on its own, then combined per target
    bonferroni = parameters(side = "left", grna_integration_strategy = "bonferroni") |> quiet(method = "thresholding"),
    stopped_after_qc = parameters(side = "left") |> quiet(method = "thresholding", through = "qc"),
    # discovery only, on 40 pairs, keeping every resampled statistic (sceptre's
    # largest output); with all 1,406 pairs this build of sceptre fails
    with_resamples = set_analysis_parameters(
        high, discovery_pairs = discovery_pairs[1:40, ], positive_control_pairs = positive_control_pairs,
        side = "left", resampling_approximation = "skew_normal"
    ) |>
        quiet(method = "thresholding", through = "qc") |>
        run_discovery_analysis(parallel = FALSE, print_progress = FALSE, output_amount = 3)
)

# ---- low MOI: 1,000 cells, 100 responses, 50 gRNAs
data(lowmoi_example_data)
low <- import_data(
    response_matrix = lowmoi_example_data$response_matrix,
    grna_matrix = lowmoi_example_data$grna_matrix,
    grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
    moi = "low",
    extra_covariates = lowmoi_example_data$extra_covariates
)
low_positive <- construct_positive_control_pairs(low)
low_discovery <- construct_trans_pairs(low, positive_control_pairs = low_positive, pairs_to_exclude = "pc_pairs")
lowmoi <- list(
    imported = low,
    # sceptre's low-MOI defaults: the "maximum" assignment method and the
    # non-targeting cells as the control group
    default = set_analysis_parameters(low, discovery_pairs = low_discovery, positive_control_pairs = low_positive) |> quiet()
)

# ---- a tiny synthetic dataset with every kind of column the format allows
cells <- paste0("cell_", 1:6)
types_counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 1, 2, 3, 1), j = c(1, 1, 2, 3, 4, 5, 6), x = c(3, 1, 2, 5, 1, 4, 2),
    dims = c(3, 6), dimnames = list(paste0("gene_", 1:3), cells)
)
types_grna <- Matrix::sparseMatrix(
    i = c(1, 2, 1, 2, 1, 2), j = 1:6, x = c(7, 9, 4, 6, 8, 5), dims = c(2, 6),
    dimnames = list(c("grna_a", "grna_nt"), cells)
)
types <- SingleCellExperiment(assays = list(counts = types_counts), colData = S4Vectors::DataFrame(
    rna_n_umis = as.double(Matrix::colSums(types_counts)),           # double
    rna_n_nonzero = as.integer(Matrix::colSums(types_counts > 0)),   # integer
    grna_n_umis = as.double(Matrix::colSums(types_grna)),
    grna_n_nonzero = as.integer(Matrix::colSums(types_grna > 0)),
    batch = factor(c("b", "a", "b", "a", "b", "a"), levels = c("b", "a")), # factor, levels not alphabetical
    lane = c("L2", "L1", "L2", "L1", "L2", "L1"),                     # text
    treated = c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE),               # true/false
    score = c(0.5, NA, 1.5, 2.5, NA, 3.5),                            # double with missing values
    S.Score = c(0.5, -1.25, 1.5, 2.5, 0.125, 3.5),                    # a name with a dot
    row.names = cells
))
altExp(types, "gRNA") <- SingleCellExperiment(assays = list(counts = types_grna), rowData = S4Vectors::DataFrame(
    grna_target = c("site_1", "non-targeting"),
    target_chr = c("chr1", NA),                                       # text with a missing value
    target_start = c(100L, NA),                                       # integer with a missing value
    row.names = rownames(types_grna)
))
metadata(types) <- list(sceptre = list(moi = "low"))

objects <- list(highmoi = highmoi, lowmoi = lowmoi)
saveRDS(objects, "example-objects.rds")
unlink("example", recursive = TRUE)
write_sceptre_screen(high, "example/highmoi/screen.h5mu", response_names = highmoi_example_data$gene_names)
write_sceptre_screen(low, "example/lowmoi/screen.h5mu")
fmt_write_screen(types, "example/types/screen.h5mu")
for (dataset in names(objects)) {
    screen <- file.path("example", dataset, "screen.h5mu")
    for (name in setdiff(names(objects[[dataset]]), "imported")) {
        write_sceptre_analysis(objects[[dataset]][[name]], screen, file.path("example", dataset, "analyses", name))
    }
}

# The same objects in the October layout of the sceptre3 branch, which the
# Python notebook uses as a reference for pysceptre's own reader.
unlink("reference-october-layout", recursive = TRUE)
for (dataset in names(objects)) for (name in names(objects[[dataset]])) {
    write_sceptre_object(objects[[dataset]][[name]], file.path("reference-october-layout", dataset, name))
}
print(file_listing("example"))
