# The example analyses, run by sceptre on its own example data: a list of
# sceptre_objects per dataset. Sourced by make_example.R and check.R.
run_stages <- function(object, ..., through = "discovery", output_amount = 1) {
    object <- object |> assign_grnas(parallel = FALSE, ...) |> run_qc()
    if (through == "qc") return(object)
    object |>
        run_calibration_check(parallel = FALSE, print_progress = FALSE, output_amount = output_amount) |>
        run_power_check(parallel = FALSE, print_progress = FALSE, output_amount = output_amount) |>
        run_discovery_analysis(parallel = FALSE, print_progress = FALSE, output_amount = output_amount)
}
example_objects <- function() {
    # high MOI: 500 cells, 100 responses, 60 gRNAs
    data(highmoi_example_data, package = "sceptre", envir = environment())
    data(grna_target_data_frame_highmoi, package = "sceptre", envir = environment())
    high <- import_data(
        response_matrix = highmoi_example_data$response_matrix, grna_matrix = highmoi_example_data$grna_matrix,
        grna_target_data_frame = grna_target_data_frame_highmoi, moi = "high",
        extra_covariates = highmoi_example_data$extra_covariates, response_names = highmoi_example_data$gene_names
    )
    positive <- construct_positive_control_pairs(high)
    discovery <- construct_cis_pairs(high, positive_control_pairs = positive, distance_threshold = 5e6)
    parameters <- function(..., pairs = discovery) set_analysis_parameters(high, discovery_pairs = pairs, positive_control_pairs = positive, ...)
    # low MOI: 1,000 cells, 100 responses, 50 gRNAs
    data(lowmoi_example_data, package = "sceptre", envir = environment())
    low <- import_data(
        response_matrix = lowmoi_example_data$response_matrix, grna_matrix = lowmoi_example_data$grna_matrix,
        grna_target_data_frame = lowmoi_example_data$grna_target_data_frame, moi = "low",
        extra_covariates = lowmoi_example_data$extra_covariates
    )
    low_positive <- construct_positive_control_pairs(low)
    low_discovery <- construct_trans_pairs(low, positive_control_pairs = low_positive, pairs_to_exclude = "pc_pairs")
    list(
        highmoi = list(
            left_sided = parameters(side = "left") |> run_stages(method = "thresholding"),
            two_sided = parameters(side = "both") |> run_stages(method = "thresholding"),
            # each gRNA tested on its own instead of pooled by target
            singleton = parameters(side = "left", grna_integration_strategy = "singleton") |> run_stages(method = "thresholding"),
            # each gRNA tested on its own, then combined per target
            bonferroni = parameters(side = "left", grna_integration_strategy = "bonferroni") |> run_stages(method = "thresholding"),
            stopped_after_qc = parameters(side = "left") |> run_stages(method = "thresholding", through = "qc"),
            # sceptre's default assignment method at high MOI
            mixture_assignment = parameters(side = "left") |> run_stages(),
            # 40 pairs, keeping every resampled statistic (sceptre's largest output)
            with_resamples = parameters(side = "left", pairs = discovery[1:40, ]) |> run_stages(method = "thresholding", output_amount = 3)
        ),
        # sceptre's low-MOI defaults: the "maximum" assignment method and the
        # non-targeting cells as the control group
        lowmoi = list(default = set_analysis_parameters(low, discovery_pairs = low_discovery, positive_control_pairs = low_positive) |> run_stages())
    )
}
screen_of <- function(dataset) file.path("..", "perturb-seq-data-format", "example", dataset, "screen.h5mu")

# What sceptre held, for the package-free checks (verify/derived.R, verify/derived.py).
write_expected <- function(object, cells, folder) {
    ids_text <- function(lists) paste(paste0(names(lists), "\t", vapply(lists, function(idx) paste(cells[idx], collapse = ","), character(1))), collapse = "\n")
    sha256_text <- function(text) unname(tools::sha256sum(bytes = charToRaw(enc2utf8(text))))
    done <- called(object)
    expected <- list(captured_from = paste("sceptre", as.character(utils::packageVersion("sceptre"))), stages = I(done))
    if ("assignment" %in% done) {
        expected$cells_by_grna_sha256 <- sha256_text(ids_text(slot_of(object, "initial_grna_assignment_list")))
        raw <- slot_of(object, "grna_assignments_raw")
        expected$cells_by_group_sha256 <- sha256_text(ids_text(raw$grna_group_idxs))
        expected$cells_by_control_grna_sha256 <- sha256_text(ids_text(raw$indiv_nt_grna_idxs))
        per_cell <- slot_of(object, "grnas_per_cell")
        if (length(per_cell) > 0L) expected$grnas_per_cell_sha256 <- sha256_text(paste(per_cell, collapse = ","))
    }
    if ("qc" %in% done) {
        expected$cells_in_use_sha256 <- sha256_ids(cells[slot_of(object, "cells_in_use")])
        expected$n_ok_discovery_pairs <- slot_of(object, "n_ok_discovery_pairs")
        expected$n_ok_positive_control_pairs <- slot_of(object, "n_ok_positive_control_pairs")
    }
    if ("calibration" %in% done) expected$control_groups <- I(sort(unique(as.character(slot_of(object, "negative_control_pairs")$grna_group)), method = "radix"))
    for (stage in intersect(c("calibration", "power", "discovery"), done)) {
        result <- as.data.frame(slot_of(object, result_slot[[stage]]), stringsAsFactors = FALSE)
        keep <- intersect(c("grna_target", "grna_id", "response_id", "p_value", "log_2_fold_change", "significant"), names(result))
        write_table(result[keep], file.path(folder, paste0(stage, ".parquet")), unit_key(result))
    }
    write_json(expected, file.path(folder, "expected.json"))
}
write_examples <- function(objects, root) {
    for (dataset in names(objects)) {
        cells <- colnames(read_as_sce(screen_of(dataset), "rna"))
        for (name in names(objects[[dataset]])) {
            suppressWarnings(write_sceptre_object(objects[[dataset]][[name]], screen_of(dataset), file.path(root, dataset, "analyses", name)))
            write_expected(objects[[dataset]][[name]], cells, file.path(root, dataset, "expected", name))
        }
    }
}
