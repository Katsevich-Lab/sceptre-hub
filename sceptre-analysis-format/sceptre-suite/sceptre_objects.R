# A sceptre_object of sceptre 0.99.1 to and from the files of the two formats.
# It reads and sets the object's slots, which are internal to sceptre, so it
# works with that version only. It is test code, not part of the format: a
# package that adopts the format has a writer and a reader of its own.
slot_of <- function(object, name) methods::slot(object, name)
swap_prefix <- function(text, from, to) gsub(paste0("(?<![A-Za-z0-9_.])", from, "_"), paste0(to, "_"), text, perl = TRUE)
to_file_names <- function(text, modality) swap_prefix(text, "response", modality)
to_package_names <- function(text, modality) swap_prefix(text, modality, "response")
SIDES <- c(left = -1L, both = 0L, right = 1L)
formula_text <- function(x) {
    if (inherits(x, "formula")) return(paste(deparse(x, width.cutoff = 500L), collapse = ""))
    if (is.list(x)) return(lapply(x, formula_text))
    x
}
same <- function(a, b) isTRUE(all.equal(a, b, check.attributes = FALSE))
must_follow <- function(derived, held, what) {
    if (!same(derived, held)) stop("The object's `", what, "` does not follow from what the files hold.", call. = FALSE)
}
called <- function(object) {
    done <- slot_of(object, "functs_called")
    names(STAGE_FUNCTIONS)[STAGE_FUNCTIONS %in% names(done)[done]]
}
unit_columns <- function(grna_group, grna_table, strategy) {
    if (identical(strategy, "union")) return(data.frame(grna_target = as.character(grna_group), stringsAsFactors = FALSE))
    data.frame(grna_target = grna_table$grna_target[match(grna_group, grna_table$grna_id)], grna_id = as.character(grna_group), stringsAsFactors = FALSE)
}
unit_group <- function(table, strategy) if (identical(strategy, "union")) table$grna_target else table$grna_id
result_slot <- c(calibration = "calibration_result", power = "power_result", discovery = "discovery_result")

# ------------------------------------------------------------------ writing
write_sceptre_object <- function(object, screen_file, directory, modality = "rna") {
    done <- called(object)
    sce <- read_as_sce(screen_file, modality)
    cells <- colnames(sce)
    grna_ids <- rownames(SingleCellExperiment::altExp(sce, "gRNA"))
    must_follow(colnames(sceptre::get_response_matrix(object)), cells, "cells")
    must_follow(rownames(sceptre::get_response_matrix(object)), rownames(sce), "responses")
    # the covariates sceptre computed are the ones the data file stores
    held <- slot_of(object, "covariate_data_frame")
    stored <- as.data.frame(SummarizedExperiment::colData(sce), optional = TRUE)
    names(stored) <- to_package_names(names(stored), modality)
    for (name in names(held)) must_follow(as.vector(stored[[name]]), as.vector(held[[name]]), paste0("covariate ", name))
    fmt_init_analysis(directory, screen_link(sce), modality)
    provenance <- list(package = "sceptre", version = as.character(utils::packageVersion("sceptre")), language = "R", written_at = "2026-11-19T00:00:00Z")
    write <- function(stage, record, tables = list()) fmt_write_stage(directory, stage, record, tables[!vapply(tables, is.null, logical(1))], provenance)
    low_moi <- isTRUE(slot_of(object, "low_moi"))
    strategy <- slot_of(object, "grna_integration_strategy")
    grna_table <- as.data.frame(slot_of(object, "grna_target_data_frame"), stringsAsFactors = FALSE)

    if ("parameters" %in% done) {
        must_follow(grna_groups(grna_table, strategy), grna_table$grna_group, "grna_group")
        formula <- to_file_names(formula_text(slot_of(object, "formula_object")), modality)
        if (!formula_in_subset(formula)) stop("The formula ", formula, " is outside the additive subset.")
        design <- slot_of(object, "covariate_matrix")
        pairs <- lapply(c(discovery = "discovery_pairs", positive_control = "positive_control_pairs"), function(name) {
            held <- as.data.frame(slot_of(object, name), stringsAsFactors = FALSE)
            table <- unique(data.frame(grna_target = as.character(held$grna_target), response_id = as.character(held$response_id), stringsAsFactors = FALSE))
            if (nrow(table) > 0L) must_follow(rows_sorted(expand_pairs(table, grna_table)), rows_sorted(held[c("grna_target", "response_id", "grna_group")]), name)
            table
        })
        write("parameters", list(
            settings = list(
                moi = if (low_moi) "low" else "high",
                formula = formula,
                side = names(SIDES)[match(slot_of(object, "side_code"), SIDES)],
                grna_integration_strategy = strategy,
                control_group = if (slot_of(object, "control_group_complement")) "complement" else "nt_cells",
                resampling_mechanism = if (slot_of(object, "run_permutations")) "permutations" else "crt",
                resampling_approximation = slot_of(object, "resampling_approximation"),
                multiple_testing_method = slot_of(object, "multiple_testing_method"),
                multiple_testing_alpha = slot_of(object, "multiple_testing_alpha"),
                B1 = slot_of(object, "B1"), B2 = slot_of(object, "B2"), B3 = slot_of(object, "B3")
            ),
            design = list(terms = I(to_file_names(colnames(design), modality)))
        ), list("parameters/discovery_pairs" = pairs$discovery, "parameters/positive_control_pairs" = pairs$positive_control))
    }

    if ("assignment" %in% done) {
        initial <- slot_of(object, "initial_grna_assignment_list")
        must_follow(names(initial), grna_ids, "the order of initial_grna_assignment_list")
        must_follow(group_cells(initial, grna_table), slot_of(object, "grna_assignments_raw"), "grna_assignments_raw")
        must_follow(mean(lengths(initial)), slot_of(object, "mean_cells_per_grna"), "mean_cells_per_grna")
        method <- slot_of(object, "grna_assignment_method")
        if (!identical(method, "maximum")) must_follow(grnas_per_cell(initial, length(cells)), slot_of(object, "grnas_per_cell"), "grnas_per_cell")
        flagged <- setdiff(slot_of(object, "cells_w_zero_or_twoplus_grnas"), if (low_moi) not_one_grna(initial, length(cells)))
        hyperparameters <- formula_text(slot_of(object, "grna_assignment_hyperparameters"))
        if (!is.null(hyperparameters$formula_object)) hyperparameters$formula_object <- to_file_names(hyperparameters$formula_object, modality)
        hyperparameters <- lapply(hyperparameters, function(x) if (length(x) > 1L) I(x) else x)
        write("assignment", list(settings = list(method = method, hyperparameters = hyperparameters)), list(
            "assignment/assignments" = data.frame(grna_id = rep(names(initial), lengths(initial)), cell_id = cells[unlist(initial, use.names = FALSE)], stringsAsFactors = FALSE),
            "assignment/flagged_cells" = if (length(flagged) > 0L) data.frame(cell_id = cells[flagged], stringsAsFactors = FALSE)
        ))
    }

    if ("qc" %in% done) {
        in_use <- slot_of(object, "cells_in_use")
        must_follow(qc_assignments(slot_of(object, "grna_assignments_raw"), in_use, !slot_of(object, "control_group_complement")), slot_of(object, "grna_assignments"), "grna_assignments")
        pair_qc <- lapply(c(discovery = "discovery", positive_control = "positive_control"), function(type) {
            info <- as.data.frame(slot_of(object, paste0(type, "_pairs_with_info")), stringsAsFactors = FALSE)
            if (nrow(info) == 0L) info <- data.frame(response_id = character(), grna_group = character(), n_nonzero_trt = integer(), n_nonzero_cntrl = integer(), pass_qc = logical(), stringsAsFactors = FALSE)
            table <- cbind(unit_columns(info$grna_group, grna_table, strategy), response_id = as.character(info$response_id),
                n_nonzero_trt = as.integer(info$n_nonzero_trt), n_nonzero_cntrl = as.integer(info$n_nonzero_cntrl), pass_qc = info$pass_qc, stringsAsFactors = FALSE)
            if (nrow(table) > 0L) must_follow(n_pairs_passing(table, strategy), slot_of(object, paste0("n_ok_", type, "_pairs")), paste0("n_ok_", type, "_pairs"))
            table
        })
        removed <- slot_of(object, "cell_removal_metrics")
        thresholds <- slot_of(object, "cellwise_qc_thresholds")
        write("qc", list(
            settings = list(
                n_nonzero_trt_thresh = slot_of(object, "n_nonzero_trt_thresh"), n_nonzero_cntrl_thresh = slot_of(object, "n_nonzero_cntrl_thresh"),
                response_n_umis_range = I(thresholds$response_n_umis_range), response_n_nonzero_range = I(thresholds$response_n_nonzero_range),
                p_mito_threshold = thresholds$p_mito_threshold, B3 = slot_of(object, "B3")
            ),
            summaries = list(cells_removed = as.list(stats::setNames(removed, sub("^n_cells_rm_", "", names(removed)))))
        ), list(
            "qc/cells" = data.frame(cell_id = cells, in_use = seq_along(cells) %in% in_use, stringsAsFactors = FALSE),
            "qc/discovery_pair_qc" = pair_qc$discovery, "qc/positive_control_pair_qc" = pair_qc$positive_control
        ))
    }

    result_of <- function(stage) as.data.frame(slot_of(object, result_slot[[stage]]), stringsAsFactors = FALSE)
    result_table <- function(stage) {
        result <- result_of(stage)
        if (!is.null(result$significant)) {
            must_follow(add_significant(result, significance_threshold(result))$significant, result$significant, "significant")
            result$significant <- NULL
        }
        if (!is.null(result$fold_change)) must_follow(log_2_of(result$fold_change), result$log_2_fold_change, "log_2_fold_change")
        result[setdiff(names(result), resample_columns(result))]
    }
    resamples <- function(stage) resamples_table(result_of(stage))
    if ("calibration" %in% done) {
        result <- result_of("calibration")
        negative <- as.data.frame(slot_of(object, "negative_control_pairs"), stringsAsFactors = FALSE)
        must_follow(
            rows_sorted(data.frame(response_id = result$response_id, grna_group = unit_group(result, strategy), n_nonzero_trt = result$n_nonzero_trt, n_nonzero_cntrl = result$n_nonzero_cntrl, pass_qc = result$pass_qc, stringsAsFactors = FALSE)),
            rows_sorted(negative[c("response_id", "grna_group", "n_nonzero_trt", "n_nonzero_cntrl", "pass_qc")]), "negative_control_pairs"
        )
        groups <- unique(as.character(negative$grna_group))
        members <- strsplit(groups, "&", fixed = TRUE)
        write("calibration", list(
            settings = list(calibration_group_size = slot_of(object, "calibration_group_size"), n_calibration_pairs = slot_of(object, "n_calibration_pairs"), B3 = slot_of(object, "B3")),
            significance_threshold = significance_threshold(result)
        ), list(
            "calibration/negative_control_groups" = data.frame(grna_group = rep(groups, lengths(members)), grna_id = unlist(members), stringsAsFactors = FALSE),
            "calibration/results" = result_table("calibration"), "calibration/resamples" = resamples("calibration")
        ))
    }
    if ("power" %in% done) write("power", list(settings = stats::setNames(list(), character())), list("power/results" = result_table("power"), "power/resamples" = resamples("power")))
    if ("discovery" %in% done) {
        write("discovery", list(settings = stats::setNames(list(), character()), significance_threshold = significance_threshold(result_of("discovery"))),
            list("discovery/results" = result_table("discovery"), "discovery/resamples" = resamples("discovery")))
    }
    fits <- slot_of(object, "response_precomputations")
    if (length(fits) > 0L && "qc" %in% done) {
        terms <- colnames(slot_of(object, "covariate_matrix"))
        table <- data.frame(response_id = names(fits), theta = vapply(fits, function(f) f$theta, numeric(1)), stringsAsFactors = FALSE)
        table$coefficients <- lapply(fits, function(f) { must_follow(names(f$fitted_coefs), terms, "the terms of a saved fit"); unname(f$fitted_coefs) })
        fmt_write_fits(directory, table[c("response_id", "coefficients", "theta")], fit_identity(object, sce, modality))
    }
    invisible(directory)
}
fit_identity <- function(object, sce, modality) {
    list(
        screen_id = S4Vectors::metadata(sce)$dataset_id,
        formula = to_file_names(formula_text(slot_of(object, "formula_object")), modality),
        cells_in_use_sha256 = sha256_ids(colnames(sce)[slot_of(object, "cells_in_use")]),
        control_group = if (slot_of(object, "control_group_complement")) "complement" else "nt_cells"
    )
}

# ------------------------------------------------------------------ reading
# The data of a data file as import_data() takes it.
import_from_file <- function(screen_file, moi, modality = "rna") {
    sce <- read_as_sce(screen_file, modality)
    grna <- SingleCellExperiment::altExp(sce, "gRNA")
    table <- as.data.frame(SummarizedExperiment::rowData(grna), optional = TRUE)
    grna_table <- data.frame(grna_id = rownames(grna), grna_target = as.character(table$grna_target), stringsAsFactors = FALSE)
    if (!is.null(table$target_chr)) {
        grna_table$chr <- as.character(table$target_chr)
        grna_table$start <- as.integer(table$target_start)
        grna_table$end <- as.integer(table$target_end)
    }
    covariates <- as.data.frame(SummarizedExperiment::colData(sce), optional = TRUE)
    extra <- covariates[!grepl(paste0("^(", modality, "|grna)_"), names(covariates))]
    object <- sceptre::import_data(
        response_matrix = SummarizedExperiment::assay(sce, "counts"), grna_matrix = SummarizedExperiment::assay(grna, "counts"),
        grna_target_data_frame = grna_table, moi = moi, extra_covariates = if (ncol(extra) > 0L) extra else data.frame(),
        response_names = if (is.null(SummarizedExperiment::rowData(sce)$response_name)) NA_character_ else as.character(SummarizedExperiment::rowData(sce)$response_name)
    )
    list(object = object, sce = sce)
}

# A data file and one analysis of it as a sceptre_object. The parameters are
# set through sceptre's own function, from the settings in the record; every
# later stage's outputs are put into the object's slots from the tables.
read_sceptre_object <- function(screen_file, directory) {
    analysis <- fmt_read_analysis(directory)
    modality <- analysis$response_modality
    records <- analysis$stages
    tables <- analysis$tables
    done <- names(records)
    settings <- records$parameters$settings
    imported <- import_from_file(screen_file, settings$moi, modality)
    object <- imported$object
    sce <- imported$sce
    check_link(analysis$screen, sce, directory)
    cells <- colnames(sce)
    response_ids <- rownames(sce)
    grna_ids <- rownames(SingleCellExperiment::altExp(sce, "gRNA"))
    set <- function(name, value) methods::slot(object, name, check = FALSE) <<- value
    strategy <- settings$grna_integration_strategy
    object <- sceptre::set_analysis_parameters(object,
        discovery_pairs = tables[["parameters/discovery_pairs"]], positive_control_pairs = tables[["parameters/positive_control_pairs"]],
        side = settings$side, grna_integration_strategy = strategy,
        formula_object = stats::as.formula(to_package_names(settings$formula, modality), env = globalenv()),
        resampling_approximation = settings$resampling_approximation, control_group = settings$control_group,
        resampling_mechanism = settings$resampling_mechanism, multiple_testing_method = settings$multiple_testing_method,
        multiple_testing_alpha = settings$multiple_testing_alpha)
    if (!identical(to_file_names(colnames(slot_of(object, "covariate_matrix")), modality), as.character(unlist(records$parameters$design$terms)))) {
        stop("The design matrix rebuilt from the formula has other columns than the analysis recorded.")
    }
    for (name in c("B1", "B2", "B3")) must_follow(slot_of(object, name), as.integer(settings[[name]]), name)
    grna_table <- as.data.frame(slot_of(object, "grna_target_data_frame"), stringsAsFactors = FALSE)
    low_moi <- identical(settings$moi, "low")

    if ("assignment" %in% done) {
        a <- records$assignment$settings
        hyper <- lapply(as.list(a$hyperparameters), function(x) if (is.list(x)) unlist(x) else x)
        if (!is.null(hyper$formula_object)) hyper$formula_object <- stats::as.formula(to_package_names(hyper$formula_object, modality), env = globalenv())
        set("grna_assignment_method", a$method)
        set("grna_assignment_hyperparameters", hyper)
        initial <- cells_by_grna(tables[["assignment/assignments"]], grna_ids, cells)
        raw <- group_cells(initial, grna_table)
        set("initial_grna_assignment_list", initial)
        set("grna_assignments_raw", raw)
        set("mean_cells_per_grna", mean(lengths(initial)))
        if (!identical(a$method, "maximum")) set("grnas_per_cell", grnas_per_cell(initial, length(cells)))
        flagged <- match(tables[["assignment/flagged_cells"]]$cell_id, cells)
        if (low_moi) set("cells_w_zero_or_twoplus_grnas", if (identical(a$method, "maximum")) flagged else union(which(grnas_per_cell(initial, length(cells)) >= 2L), which(grnas_per_cell(initial, length(cells)) == 0L)))
    }
    if ("qc" %in% done) {
        q <- records$qc$settings
        set("cellwise_qc_thresholds", list(response_n_umis_range = as.numeric(unlist(q$response_n_umis_range)), response_n_nonzero_range = as.numeric(unlist(q$response_n_nonzero_range)), p_mito_threshold = q$p_mito_threshold))
        set("n_nonzero_trt_thresh", as.integer(q$n_nonzero_trt_thresh))
        set("n_nonzero_cntrl_thresh", as.integer(q$n_nonzero_cntrl_thresh))
        set("B3", as.integer(q$B3))
        removed <- unlist(records$qc$summaries$cells_removed)
        set("cell_removal_metrics", stats::setNames(as.integer(removed), paste0("n_cells_rm_", names(removed))))
        table_cells <- tables[["qc/cells"]][match(cells, tables[["qc/cells"]]$cell_id), , drop = FALSE]
        in_use <- which(table_cells$in_use)
        set("cells_in_use", in_use)
        set("grna_assignments", qc_assignments(raw, in_use, identical(settings$control_group, "nt_cells")))
        for (type in c("discovery", "positive_control")) {
            pair_qc <- tables[[paste0("qc/", type, "_pair_qc")]]
            if (nrow(pair_qc) == 0L) next
            info <- data.frame(response_id = pair_qc$response_id, grna_group = unit_group(pair_qc, strategy), n_nonzero_trt = pair_qc$n_nonzero_trt,
                n_nonzero_cntrl = pair_qc$n_nonzero_cntrl, pass_qc = pair_qc$pass_qc, stringsAsFactors = FALSE)
            set(paste0(type, "_pairs_with_info"), info)
            set(paste0("n_ok_", type, "_pairs"), as.integer(n_pairs_passing(pair_qc, strategy)))
        }
    }
    result <- function(stage) {
        table <- add_resamples(tables[[paste0(stage, "/results")]], tables[[paste0(stage, "/resamples")]])
        if (stage != "power") table <- add_significant(table, records[[stage]]$significance_threshold)
        table
    }
    if ("calibration" %in% done) {
        s <- records$calibration$settings
        set("calibration_group_size", as.integer(s$calibration_group_size))
        set("n_calibration_pairs", as.integer(s$n_calibration_pairs))
        set("B3", as.integer(s$B3))
        calibration <- tables[["calibration/results"]]
        set("negative_control_pairs", data.frame(response_id = calibration$response_id, grna_group = unit_group(calibration, strategy),
            n_nonzero_trt = calibration$n_nonzero_trt, n_nonzero_cntrl = calibration$n_nonzero_cntrl, pass_qc = calibration$pass_qc, stringsAsFactors = FALSE))
    }
    for (stage in intersect(c("calibration", "power", "discovery"), done)) set(result_slot[[stage]], data.table::as.data.table(result(stage)))
    if ("qc" %in% done) {
        fits <- fmt_read_fits(directory, fit_identity(object, sce, modality))
        if (!is.null(fits)) {
            terms <- colnames(slot_of(object, "covariate_matrix"))
            set("response_precomputations", lapply(stats::setNames(seq_len(nrow(fits$table)), fits$table$response_id), function(i) {
                list(fitted_coefs = stats::setNames(fits$table$coefficients[[i]], terms), theta = fits$table$theta[i])
            }))
        }
    }
    status <- stats::setNames(c(TRUE, STAGE_FUNCTIONS %in% STAGE_FUNCTIONS[done]), c("import_data", STAGE_FUNCTIONS))
    set("functs_called", status[names(slot_of(object, "functs_called"))])
    set("last_function_called", unname(STAGE_FUNCTIONS[[utils::tail(done, 1L)]]))
    object
}

# One slot of two objects: identical, equal (to rounding, or up to row or name
# order), or different.
compare_slot <- function(x, y) {
    plain <- function(x) {
        if (inherits(x, "formula")) return(paste(deparse(x, width.cutoff = 500L), collapse = ""))
        if (is.data.frame(x)) { x <- as.data.frame(x, stringsAsFactors = FALSE); return(if (nrow(x) > 0L) rows_sorted(x) else x[0, 0]) }
        if (is.list(x) && !is.null(names(x))) return(lapply(x[sort(names(x))], plain))
        x
    }
    if (identical(x, y)) return("identical")
    x <- plain(x); y <- plain(y)
    if (identical(x, y)) return("identical up to row or name order")
    difference <- tryCatch(all.equal(x, y, check.attributes = FALSE), error = function(e) conditionMessage(e))
    if (isTRUE(difference)) "equal" else paste("DIFFERENT:", paste(utils::head(difference, 2), collapse = "; "))
}
