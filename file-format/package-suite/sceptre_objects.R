# Files of the sceptre format <-> sceptre_object.
#
# This is the package-dependent half of the prototype. It needs the sceptre3
# build of the sceptre package, an unreviewed branch build, and reaches into its
# internals. It is here to test that the files hold everything an analysis
# needs, and is meant to move into the package when the package adopts the
# format. Nothing in format/ or verify/ uses it.
#
#   write_sceptre_screen(), write_sceptre_analysis(), read_sceptre_object()
#
# The sceptre package calls the response modality `response` whatever was
# measured, and the gRNA table's target coordinates `chr`, `start`, `end`. The
# functions here translate between those names and the file's.

sceptre_ns <- function() asNamespace("sceptre")
field <- function(object, name) sceptre_ns()$object_field(object, name)
`field<-` <- function(object, name, value) {
    setter <- get("object_field<-", envir = sceptre_ns())
    setter(object, name, value)
}

# The package's names <-> the file's names. The package prefixes the response
# modality's covariates with `response_`; the file uses the modality's name.
swap_prefix <- function(text, from, to) {
    gsub(paste0("(?<![A-Za-z0-9_.])", from, "_"), paste0(to, "_"), text, perl = TRUE)
}
to_file_names <- function(text, modality) swap_prefix(text, "response", modality)
to_package_names <- function(text, modality) swap_prefix(text, modality, "response")
rename_columns <- function(frame, map) {
    hit <- match(names(frame), names(map))
    names(frame)[!is.na(hit)] <- unname(map[hit[!is.na(hit)]])
    frame
}

SIDES <- c(left = -1L, both = 0L, right = 1L)
# The user-facing name "spa" is the default saddlepoint route of each mechanism.
SPA_DEFAULT <- c(crt = "crt_spa_fast", permutations = "rpt_spa_fast")

formula_text <- function(x) {
    if (inherits(x, "formula")) return(paste(deparse(x, width.cutoff = 500L), collapse = ""))
    if (is.list(x)) return(lapply(x, formula_text))
    x
}
text_formula <- function(x) {
    if (is.list(x)) {
        is_formula <- names(x) %in% "formula_object"
        x[is_formula] <- lapply(x[is_formula], text_formula)
        return(x)
    }
    stats::as.formula(x, env = globalenv())
}

same <- function(a, b) isTRUE(all.equal(a, b, check.attributes = FALSE))
must_follow <- function(derived, held, what) {
    if (!same(derived, held)) fmt_stop("The object's `", what, "` does not follow from what the files hold.")
}

package_provenance <- function() {
    list(package = "sceptre", version = as.character(utils::packageVersion("sceptre")), language = "R")
}

# The data of a sceptre_object as the SingleCellExperiment fmt_write_screen takes.
object_screen <- function(object, response_names = NULL, id = NULL, response_modality = DEFAULT_RESPONSE_MODALITY,
                          metadata = attr(object, "screen_metadata")) {
    ns <- sceptre_ns()
    screen <- ns$file_screen(object)
    meta <- S4Vectors::metadata(screen)
    grna_table <- as.data.frame(meta$grna_target_data_frame, stringsAsFactors = FALSE)
    # a gRNA under two targets has no place in the gRNA table (S8)
    if (anyDuplicated(grna_table$grna_id) > 0L) fmt_stop("A gRNA listed under two targets cannot be written in this format.")
    grna <- SingleCellExperiment::altExp(screen, "gRNA")
    grna_table <- grna_table[match(rownames(grna), grna_table$grna_id), , drop = FALSE]
    grna_table$grna_group <- NULL
    held <- as.list(SummarizedExperiment::colData(grna))
    if (length(held) > 0L) {
        must_follow(top_grna(SummarizedExperiment::assay(grna, "counts")), held[c("max_grna", "max_grna_frac_umis")], "top gRNA per cell")
    }
    SummarizedExperiment::colData(grna) <- NULL
    SummarizedExperiment::rowData(grna) <- S4Vectors::DataFrame(
        rename_columns(grna_table[setdiff(names(grna_table), "grna_id")], GRNA_COORDINATES),
        row.names = grna_table$grna_id, check.names = FALSE
    )
    SingleCellExperiment::altExps(screen) <- list()
    if (!is.null(response_names)) SummarizedExperiment::rowData(screen)$response_name <- as.character(response_names)
    SingleCellExperiment::altExp(screen, "gRNA") <- grna
    names(SummarizedExperiment::colData(screen)) <- to_file_names(names(SummarizedExperiment::colData(screen)), response_modality)
    # `metadata`: the optional entries (name, source, perturbation_type, organism, genome_build)
    S4Vectors::metadata(screen) <- list(sceptre = c(
        list(moi = if (isTRUE(meta$low_moi)) "low" else "high", id = id, response_modality = response_modality),
        as.list(metadata)[intersect(OPTIONAL_METADATA, names(metadata))]
    ))
    screen
}

write_sceptre_screen <- function(object, file, response_names = NULL, id = attr(object, "screen_id"),
                                 response_modality = DEFAULT_RESPONSE_MODALITY, metadata = attr(object, "screen_metadata")) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    fmt_write_screen(object_screen(object, response_names, id, response_modality, metadata), file)
}

called <- function(object) {
    done <- field(object, "functs_called")
    names(STAGE_FUNCTIONS)[STAGE_FUNCTIONS %in% names(done)[done]]
}

unit_columns <- function(grna_group, grna_table, strategy) {
    if (identical(strategy, "union")) return(data.frame(grna_target = as.character(grna_group)))
    data.frame(
        grna_target = grna_table$grna_target[match(grna_group, grna_table$grna_id)],
        grna_id = as.character(grna_group)
    )
}
unit_group <- function(table, strategy) if (identical(strategy, "union")) table$grna_target else table$grna_id

# Write the analysis of a sceptre_object into a new folder, as the stages
# would have written it one after another.
write_sceptre_analysis <- function(object, screen_file, directory, provenance = attr(object, "format_provenance"),
                                   response_modality = NULL) {
    ns <- sceptre_ns()
    done <- called(object)
    sce <- fmt_read_screen(screen_file, response_modality)
    modality <- S4Vectors::metadata(sce)$sceptre$response_modality
    cells <- colnames(sce)
    response_ids <- rownames(sce)
    grna_ids <- rownames(SingleCellExperiment::altExp(sce, "gRNA"))
    must_follow(colnames(ns$file_screen(object)), cells, "cells")
    must_follow(rownames(ns$file_screen(object)), response_ids, "responses")
    fmt_init_analysis(directory, screen_link(sce, screen_file, directory), modality)
    stage_provenance <- function(stage) if (is.null(provenance[[stage]])) package_provenance() else provenance[[stage]]
    write <- function(stage, record, tables = list()) {
        fmt_write_stage(directory, stage, record, tables[!vapply(tables, is.null, logical(1))], stage_provenance(stage))
    }
    low_moi <- isTRUE(field(object, "low_moi"))

    if ("parameters" %in% done) {
        strategy <- field(object, "grna_integration_strategy")
        grna_table <- as.data.frame(field(object, "grna_target_data_frame"), stringsAsFactors = FALSE)
        must_follow(grna_groups(grna_table, strategy), grna_table$grna_group, "grna_group")
        mechanism <- if (field(object, "run_permutations")) "permutations" else "crt"
        approximation <- field(object, "resampling_approximation")
        if (identical(approximation, SPA_DEFAULT[[mechanism]])) approximation <- "spa"
        # in the file, a covariate carries its modality's name
        formula <- to_file_names(formula_text(field(object, "formula_object")), modality)
        if (!formula_in_subset(formula)) {
            fmt_stop("The formula ", formula, " cannot be written: the format takes additive terms only, each a covariate name, log(x) or log(x + 1).")
        }
        design <- field(object, "covariate_matrix")
        rebuilt <- ns$convert_covariate_df_to_design_matrix(field(object, "covariate_data_frame"), field(object, "formula_object"))
        if (!identical(unname(design), unname(rebuilt))) fmt_stop("The object's design matrix is not the one its formula and covariates give.")
        pairs <- lapply(c(discovery = "discovery_pairs", positive_control = "positive_control_pairs"), function(name) {
            held <- as.data.frame(field(object, name), stringsAsFactors = FALSE)
            table <- unique(data.frame(grna_target = as.character(held$grna_target), response_id = as.character(held$response_id)))
            must_follow(rows_sorted(expand_pairs(table, grna_table)), rows_sorted(held[c("grna_target", "response_id", "grna_group")]), name)
            table
        })
        write("parameters", list(
            settings = list(
                formula = formula,
                side = names(SIDES)[match(field(object, "side_code"), SIDES)],
                resampling_mechanism = mechanism,
                resampling_approximation = approximation,
                control_group = if (field(object, "control_group_complement")) "complement" else "nt_cells",
                grna_integration_strategy = strategy,
                B1 = field(object, "B1"), B2 = field(object, "B2"), B3 = field(object, "B3"),
                multiple_testing_method = field(object, "multiple_testing_method"),
                multiple_testing_alpha = field(object, "multiple_testing_alpha"),
                seed = field(object, "seed"),
                response_fit_method = field(object, "response_fit_method"),
                grna_fit_method = field(object, "grna_fit_method")
            ),
            design = list(terms = I(to_file_names(colnames(design), modality)), sha256 = design_sha256(design))
        ), list("parameters/discovery_pairs" = pairs$discovery, "parameters/positive_control_pairs" = pairs$positive_control))
    }

    if ("assignment" %in% done) {
        initial <- field(object, "initial_grna_assignment_list")
        must_follow(names(initial), grna_ids, "the order of initial_grna_assignment_list")
        must_follow(lapply(initial, sort), initial, "the order of cells in initial_grna_assignment_list")
        must_follow(group_cells(initial, grna_table), field(object, "grna_assignments_raw"), "grna_assignments_raw")
        must_follow(mean(lengths(initial)), field(object, "mean_cells_per_grna"), "mean_cells_per_grna")
        held <- field(object, "grnas_per_cell")
        if (length(held) > 0L) must_follow(grnas_per_cell(initial, length(cells)), held, "grnas_per_cell")
        # cells the method itself flags as not carrying exactly one gRNA, beyond
        # what the assignment table shows (the "maximum" method does this)
        flagged <- setdiff(field(object, "cells_w_zero_or_twoplus_grnas"), if (low_moi) not_one_grna(initial, length(cells)))
        if (!low_moi && length(flagged) > 0L) fmt_stop("Cells are flagged in a high-MOI analysis.")
        hyperparameters <- formula_text(field(object, "grna_assignment_hyperparameters"))
        write("assignment", list(settings = list(
            method = field(object, "grna_assignment_method"),
            hyperparameters = if (length(hyperparameters) == 0L) stats::setNames(list(), character()) else hyperparameters
        )), list(
            "assignment/assignments" = data.frame(
                grna_id = rep(names(initial), lengths(initial)), cell_id = cells[unlist(initial, use.names = FALSE)]
            ),
            "assignment/flagged_cells" = if (length(flagged) > 0L) data.frame(cell_id = cells[flagged])
        ))
    }

    if ("qc" %in% done) {
        in_use <- field(object, "cells_in_use")
        nt_cells_control <- !field(object, "control_group_complement")
        must_follow(
            qc_assignments(field(object, "grna_assignments_raw"), in_use, nt_cells_control),
            field(object, "grna_assignments"), "grna_assignments"
        )
        pair_qc <- lapply(c(discovery = "discovery", positive_control = "positive_control"), function(type) {
            info <- as.data.frame(field(object, paste0(type, "_pairs_with_info")), stringsAsFactors = FALSE)
            table <- cbind(unit_columns(info$grna_group, grna_table, strategy), response_id = as.character(info$response_id), pass_qc = info$pass_qc)
            if (nrow(table) > 0L) must_follow(n_pairs_passing(table, strategy), field(object, paste0("n_ok_", type, "_pairs")), paste0("n_ok_", type, "_pairs"))
            table
        })
        removed <- field(object, "cell_removal_metrics")
        write("qc", list(
            settings = c(field(object, "cellwise_qc_thresholds"), list(pairwise_qc = field(object, "pairwise_qc"), B3 = field(object, "B3"))),
            summaries = list(cells_removed = as.list(stats::setNames(removed, sub("^n_cells_rm_", "", names(removed)))))
        ), list(
            "qc/cells" = data.frame(cell_id = cells, in_use = seq_along(cells) %in% in_use),
            "qc/discovery_pair_qc" = pair_qc$discovery, "qc/positive_control_pair_qc" = pair_qc$positive_control
        ))
    }

    fits <- ns$response_fits(object)
    fit_table <- function() {
        terms <- colnames(field(object, "covariate_matrix"))
        table <- data.frame(response_id = names(fits), theta = vapply(fits, function(f) f$theta, numeric(1)))
        table$coefficients <- lapply(fits, function(f) {
            must_follow(names(f$fitted_coefs), terms, "the terms of a saved fit")
            unname(f$fitted_coefs)
        })
        table[c("response_id", "coefficients", "theta")]
    }
    result_table <- function(stage) {
        result <- as.data.frame(ns$analysis_result(object, stage), stringsAsFactors = FALSE)
        if (!is.null(result$significant)) {
            must_follow(add_significant(result, significance_threshold(result))$significant, result$significant, "significant")
            result$significant <- NULL
        }
        # stored only where it does not follow from `fold_change` (T5)
        if (log_2_follows(result)) result$log_2_fold_change <- NULL
        result[setdiff(names(result), resample_columns(result))]
    }
    resamples <- function(stage) resamples_table(as.data.frame(ns$analysis_result(object, stage), stringsAsFactors = FALSE))

    if ("calibration" %in% done) {
        result <- as.data.frame(ns$analysis_result(object, "calibration"), stringsAsFactors = FALSE)
        negative <- as.data.frame(field(object, "negative_control_pairs"), stringsAsFactors = FALSE)
        # a control group is named in the column that names the tested unit:
        # `grna_target` under union, `grna_id` otherwise
        must_follow(
            rows_sorted(data.frame(response_id = result$response_id, grna_group = unit_group(result, strategy), pass_qc = result$pass_qc)),
            rows_sorted(data.frame(response_id = negative$response_id, grna_group = negative$grna_group, pass_qc = negative$pass_qc)),
            "negative_control_pairs"
        )
        # today a control group's members exist only in its name (T4)
        groups <- unique(as.character(negative$grna_group))
        members <- strsplit(groups, "&", fixed = TRUE)
        write("calibration", list(
            settings = list(
                calibration_group_size = field(object, "calibration_group_size"),
                n_calibration_pairs = field(object, "n_calibration_pairs"), B3 = field(object, "B3")
            ),
            significance_threshold = significance_threshold(result)
        ), list(
            "calibration/negative_control_groups" = data.frame(grna_group = rep(groups, lengths(members)), grna_id = unlist(members)),
            "calibration/results" = result_table("calibration"), "calibration/resamples" = resamples("calibration")
        ))
    }
    if ("power" %in% done) {
        write("power", list(settings = stats::setNames(list(), character())),
            list("power/results" = result_table("power"), "power/resamples" = resamples("power")))
    }
    if ("discovery" %in% done) {
        write("discovery", list(
            settings = stats::setNames(list(), character()),
            significance_threshold = significance_threshold(as.data.frame(ns$analysis_result(object, "discovery")))
        ), list("discovery/results" = result_table("discovery"), "discovery/resamples" = resamples("discovery")))
    }
    # the fits the object holds go to the cache, which no stage owns (T6)
    if (length(fits) > 0L && "qc" %in% done) fmt_write_fits(directory, fit_table(), fit_identity(object, sce))
    invisible(directory)
}

# What makes a saved fit valid (T6).
fit_identity <- function(object, sce) {
    list(
        screen_id = S4Vectors::metadata(sce)$sceptre$id,
        formula = to_file_names(formula_text(field(object, "formula_object")), S4Vectors::metadata(sce)$sceptre$response_modality),
        cells_in_use_sha256 = sha256_ids(colnames(sce)[field(object, "cells_in_use")]),
        response_fit_method = field(object, "response_fit_method"),
        control_group = if (field(object, "control_group_complement")) "complement" else "nt_cells"
    )
}

# Read a data file, and optionally one analysis of it, as a sceptre_object.
# Without `directory` the object is as after import_data(). The response
# modality is the one the analysis names; without an analysis, `response_modality`
# or the file's only modality besides `grna`.
read_sceptre_object <- function(screen_file, directory = NULL, response_modality = NULL) {
    ns <- sceptre_ns()
    if (!is.null(directory)) response_modality <- fmt_response_modality(directory)
    sce <- fmt_read_screen(screen_file, response_modality)
    modality <- S4Vectors::metadata(sce)$sceptre$response_modality
    grna <- SingleCellExperiment::altExp(sce, "gRNA")
    grna_table <- data.frame(
        grna_id = rownames(grna),
        rename_columns(
            as.data.frame(SummarizedExperiment::rowData(grna), optional = TRUE),
            stats::setNames(names(GRNA_COORDINATES), GRNA_COORDINATES)
        ),
        stringsAsFactors = FALSE
    )
    rownames(grna_table) <- NULL
    cells <- colnames(sce)
    response_ids <- rownames(sce)
    low_moi <- identical(S4Vectors::metadata(sce)$sceptre$moi, "low")
    covariates <- SummarizedExperiment::colData(sce)
    names(covariates) <- to_package_names(names(covariates), modality)
    file_screen <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = ns$as_row_counts(SummarizedExperiment::assay(sce, "counts"))),
        colData = covariates,
        altExps = list(gRNA = SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = ns$as_row_counts(SummarizedExperiment::assay(grna, "counts"))),
            colData = S4Vectors::DataFrame(top_grna(SummarizedExperiment::assay(grna, "counts")), row.names = cells)
        )),
        metadata = list(low_moi = low_moi, grna_target_data_frame = grna_table, covariate_row_names = cells)
    )
    object <- ns$object_from_file_screen(file_screen)
    # what import_data() sets beyond the data
    field(object, "covariate_names") <- sort(names(field(object, "covariate_data_frame")))
    field(object, "nuclear") <- FALSE
    field(object, "nf_pipeline") <- FALSE
    field(object, "elements_to_analyze") <- NA_character_
    done <- character()
    provenance <- list()
    if (!is.null(directory)) {
        analysis <- fmt_read_analysis(directory)
        check_link(analysis$screen, sce, directory)
        done <- names(analysis$stages)
        tables <- analysis$tables
        records <- analysis$stages
        provenance <- lapply(records, function(r) r$provenance)
    }
    set <- function(name, value) field(object, name) <<- value

    if ("parameters" %in% done) {
        settings <- records$parameters$settings
        strategy <- settings$grna_integration_strategy
        grna_table$grna_group <- grna_groups(grna_table, strategy)
        approximation <- settings$resampling_approximation
        if (identical(approximation, "spa")) approximation <- SPA_DEFAULT[[settings$resampling_mechanism]]
        set("grna_target_data_frame", grna_table)
        set("formula_object", text_formula(to_package_names(settings$formula, modality)))
        set("side_code", SIDES[[settings$side]])
        set("run_permutations", identical(settings$resampling_mechanism, "permutations"))
        set("resampling_approximation", approximation)
        set("control_group_complement", identical(settings$control_group, "complement"))
        set("grna_integration_strategy", strategy)
        for (name in c("B1", "B2", "B3", "seed")) set(name, as.integer(settings[[name]]))
        set("multiple_testing_method", settings$multiple_testing_method)
        set("multiple_testing_alpha", as.double(settings$multiple_testing_alpha))
        set("response_fit_method", settings$response_fit_method)
        set("grna_fit_method", settings$grna_fit_method)
        design <- ns$convert_covariate_df_to_design_matrix(field(object, "covariate_data_frame"), field(object, "formula_object"))
        if (!identical(to_file_names(colnames(design), modality), as.character(records$parameters$design$terms)) ||
            !identical(design_sha256(design), records$parameters$design$sha256)) {
            fmt_stop("The design matrix rebuilt from the formula is not the one the analysis was run with.")
        }
        set("covariate_matrix", design)
        for (type in c("discovery", "positive_control")) {
            pairs <- expand_pairs(tables[[paste0("parameters/", type, "_pairs")]], grna_table)
            rownames(pairs) <- NULL
            set(paste0(type, "_pairs"), pairs)
            # counted as the user gave them, before any expansion by gRNA
            set(paste0("n_", type, "_pairs"), nrow(tables[[paste0("parameters/", type, "_pairs")]]))
        }
    }

    if ("assignment" %in% done) {
        settings <- records$assignment$settings
        set("grna_assignment_method", settings$method)
        set("grna_assignment_hyperparameters", text_formula(as.list(settings$hyperparameters)))
        initial <- cells_by_grna(tables[["assignment/assignments"]], rownames(grna), cells)
        raw <- group_cells(initial, grna_table)
        set("initial_grna_assignment_list", initial)
        set("grna_assignments_raw", raw)
        set("mean_cells_per_grna", mean(lengths(initial)))
        # the "maximum" method does not count gRNAs per cell
        if (!identical(settings$method, "maximum")) set("grnas_per_cell", grnas_per_cell(initial, length(cells)))
        if (low_moi) {
            flagged <- match(tables[["assignment/flagged_cells"]]$cell_id, cells)
            set("cells_w_zero_or_twoplus_grnas", sort(union(not_one_grna(initial, length(cells)), flagged)))
        }
    }

    if ("qc" %in% done) {
        settings <- records$qc$settings
        set("cellwise_qc_thresholds", settings[c("response_n_umis_range", "response_n_nonzero_range", "p_mito_threshold")])
        set("pairwise_qc", isTRUE(settings$pairwise_qc))
        set("B3", as.integer(settings$B3))
        removed <- unlist(records$qc$summaries$cells_removed)
        set("cell_removal_metrics", stats::setNames(as.integer(removed), paste0("n_cells_rm_", names(removed))))
        table_cells <- tables[["qc/cells"]][match(cells, tables[["qc/cells"]]$cell_id), , drop = FALSE]
        in_use <- which(table_cells$in_use)
        set("cells_in_use", in_use)
        set("grna_assignments", qc_assignments(raw, in_use, identical(records$parameters$settings$control_group, "nt_cells")))
        for (type in c("discovery", "positive_control")) {
            pair_qc <- tables[[paste0("qc/", type, "_pair_qc")]]
            info <- data.frame(response_id = pair_qc$response_id, grna_group = unit_group(pair_qc, strategy), pass_qc = pair_qc$pass_qc)
            # the order a package tests in: by response, in the data file's order
            info <- info[order(match(info$response_id, response_ids)), , drop = FALSE]
            rownames(info) <- NULL
            set(paste0(type, "_pairs_with_info"), data.table::as.data.table(info))
            if (nrow(info) > 0L) set(paste0("n_ok_", type, "_pairs"), as.integer(n_pairs_passing(pair_qc, strategy)))
        }
        # per response, the cells of the control pool that express it: a pass
        # over the counts, done with the package's own routine
        object <- restore_pool_nonzero(object)
    }

    result <- function(stage) {
        table <- add_log_2_fold_change(tables[[paste0(stage, "/results")]])
        table <- add_resamples(table, tables[[paste0(stage, "/resamples")]])
        if (stage != "power") table <- add_significant(table, records[[stage]]$significance_threshold)
        # as sceptre returns results: by p-value, then response
        table <- table[do.call(order, c(list(table$p_value, table$response_id), unname(as.list(table[unit_key(table)])), list(method = "radix"))), , drop = FALSE]
        rownames(table) <- NULL
        data.table::as.data.table(table)
    }
    set_result <- function(stage) {
        setter <- get("analysis_result<-", envir = ns)
        object <<- setter(object, stage, value = result(stage))
    }
    if ("calibration" %in% done) {
        settings <- records$calibration$settings
        set("calibration_group_size", as.integer(settings$calibration_group_size))
        set("n_calibration_pairs", as.integer(settings$n_calibration_pairs))
        set("B3", as.integer(settings$B3))
        calibration <- tables[["calibration/results"]]
        set("negative_control_pairs", data.table::as.data.table(data.frame(
            response_id = calibration$response_id, grna_group = unit_group(calibration, strategy), pass_qc = calibration$pass_qc
        )))
        set_result("calibration")
    }
    if ("power" %in% done) set_result("power")
    if ("discovery" %in% done) set_result("discovery")

    # cached fits are reused only if they were made as this analysis would make them (T6)
    saved <- if (is.null(directory)) NULL else analysis$fits_identity
    if (!is.null(saved) && "qc" %in% done && !is.null(tables[[FITS]])) {
        reason <- fit_mismatch(saved, fit_identity(object, sce))
        if (!is.null(reason)) {
            message("The fits saved in ", FITS, " are not reused: ", reason, ".")
        } else {
            table <- tables[[FITS]]
            terms <- colnames(field(object, "covariate_matrix"))
            fits <- lapply(stats::setNames(seq_len(nrow(table)), table$response_id), function(i) {
                list(fitted_coefs = stats::setNames(table$coefficients[[i]], terms), theta = table$theta[i])
            })
            setter <- get("response_fits<-", envir = ns)
            object <- setter(object, value = c(ns$response_fits(object), fits))
        }
    }

    status <- stats::setNames(c(TRUE, STAGE_FUNCTIONS %in% STAGE_FUNCTIONS[done]), c("import_data", STAGE_FUNCTIONS))
    field(object, "functs_called") <- status
    field(object, "last_function_called") <- if (length(done) == 0L) "import_data" else STAGE_FUNCTIONS[[utils::tail(done, 1L)]]
    attr(object, "format_provenance") <- provenance
    attr(object, "screen_id") <- S4Vectors::metadata(sce)$sceptre$id
    stamp <- S4Vectors::metadata(sce)$sceptre
    attr(object, "screen_metadata") <- stamp[intersect(OPTIONAL_METADATA, names(stamp))]
    object
}

# Why cached fits cannot be reused, or NULL. A method name means the same
# estimator in every language that lists it; a method this language does not
# share is never reused.
fit_mismatch <- function(saved, current) {
    for (name in names(current)) {
        if (!identical(identity_value(saved, name), identity_value(current, name))) return(paste0("their `", name, "` differs"))
    }
    method <- current$response_fit_method
    if (!"R" %in% FIT_METHODS[[method]]) return(paste0("the method '", method, "' is not one R shares"))
    NULL
}

# Per response, how many cells of the control pool express it. sceptre computes
# this during pairwise QC; it is recomputed here with sceptre's own routine.
restore_pool_nonzero <- function(object) {
    ns <- sceptre_ns()
    info <- rbind(
        as.data.frame(field(object, "discovery_pairs_with_info")),
        as.data.frame(field(object, "positive_control_pairs_with_info"))
    )
    if (nrow(info) == 0L) return(object)
    response_ids <- rownames(ns$get_response_matrix(object))
    groups <- names(field(object, "grna_assignments")$grna_group_idxs)
    # as sceptre calls it: pairs ordered by response
    info <- info[order(match(info$response_id, response_ids)), , drop = FALSE]
    margins <- ns$pair_margins(object, match(info$response_id, response_ids), match(info$grna_group, groups))
    field(object, "n_nonzero_tot_vector") <- margins$pool_nonzero
    object
}

# Compare two sceptre_objects field by field. Tables are compared with their
# rows sorted and the fit cache by response, since the files keep neither order.
compare_objects <- function(a, b) {
    ns <- sceptre_ns()
    fields <- setdiff(ns$object_field_names, "response_precomputations")
    plain <- function(x, name) {
        if (inherits(x, "formula")) return(paste(deparse(x, width.cutoff = 500L), collapse = ""))
        if (is.data.frame(x)) return(rows_sorted(x))
        if (methods::is(x, "sparseMatrix")) return(methods::as(methods::as(x, "generalMatrix"), "CsparseMatrix"))
        if (name == "cells_w_zero_or_twoplus_grnas") return(sort(x))
        x
    }
    verdict <- vapply(fields, function(name) {
        x <- plain(field(a, name), name)
        y <- plain(field(b, name), name)
        if (identical(x, y)) return("identical")
        difference <- all.equal(x, y, check.attributes = FALSE)
        if (isTRUE(difference)) "equal" else paste(difference, collapse = "; ")
    }, character(1))
    by_name <- function(fits) fits[sort(names(fits))]
    fits <- if (identical(by_name(ns$response_fits(a)), by_name(ns$response_fits(b)))) {
        "identical"
    } else {
        difference <- all.equal(by_name(ns$response_fits(a)), by_name(ns$response_fits(b)), check.attributes = FALSE)
        if (isTRUE(difference)) "equal" else paste(difference, collapse = "; ")
    }
    c(verdict, response_fits = fits)
}

# A one-line summary of compare_objects().
comparison_summary <- function(verdict) {
    different <- names(verdict)[!verdict %in% c("identical", "equal")]
    paste0(
        sum(verdict == "identical"), " fields identical, ", sum(verdict == "equal"), " equal, ", length(different), " different",
        if (length(different) > 0L) paste0(": ", paste(different, collapse = ", ")) else ""
    )
}
