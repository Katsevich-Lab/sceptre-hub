# The checks that need sceptre: that the format holds what a sceptre analysis
# holds. sceptre runs the example analyses again, and then
#   1. the files written from its objects conform and equal the committed examples;
#   2. a sceptre_object rebuilt from the data file and the folder alone equals
#      the object the folder was written from, slot by slot;
#   3. sceptre, run again from only the settings in the records, gives the
#      same results.
# Run from the folder that holds format/:  Rscript sceptre-suite/check.R
source("sceptre-suite/setup.R")
failures <- 0L
n_checks <- 0L
check <- function(name, ok, detail = "") {
    n_checks <<- n_checks + 1L
    if (!isTRUE(ok)) {
        failures <<- failures + 1L
        cat("FAIL ", name, if (nzchar(detail)) paste0(" - ", detail), "\n", sep = "")
    }
}
objects <- suppressMessages(example_objects())
out <- "verify/out/sceptre"
unlink(out, recursive = TRUE)
write_examples(objects, out)

# Slots sceptre computes again from the data, which the format does not store.
NOT_STORED <- c("M_matrix", "n_nonzero_tot_vector")
# A set of cells, held in no particular order.
UNORDERED <- "cells_w_zero_or_twoplus_grnas"
# Numbers are compared to a relative tolerance: a result computed on another
# machine may differ in its last bits.
close_enough <- function(a, b) isTRUE(all.equal(a, b, check.attributes = FALSE, tolerance = 1e-8))
same_files <- function(a, b) {
    # a part of the fit cache has a name of its writer's choosing; each example has one part
    files <- lapply(list(a, b), function(folder) sort(list.files(folder, recursive = TRUE)))
    if (!identical(sub("part-[0-9a-f]+", "part", files[[1]]), sub("part-[0-9a-f]+", "part", files[[2]]))) return("the folders hold different files")
    for (i in seq_along(files[[1]])) {
        paths <- c(file.path(a, files[[1]][i]), file.path(b, files[[2]][i]))
        same <- if (grepl("\\.parquet$", paths[1])) {
            close_enough(as.data.frame(nanoparquet::read_parquet(paths[1])), as.data.frame(nanoparquet::read_parquet(paths[2])))
        } else {
            close_enough(read_json(paths[1]), read_json(paths[2]))
        }
        if (!same) return(paste(files[[1]][i], "differs"))
    }
    TRUE
}
replay <- function(screen, directory) {
    analysis <- fmt_read_analysis(directory)
    records <- analysis$stages
    modality <- analysis$response_modality
    p <- records$parameters$settings
    object <- import_from_file(screen, p$moi, modality)$object |> set_analysis_parameters(
        discovery_pairs = analysis$tables[["parameters/discovery_pairs"]], positive_control_pairs = analysis$tables[["parameters/positive_control_pairs"]],
        side = p$side, grna_integration_strategy = p$grna_integration_strategy,
        formula_object = stats::as.formula(to_package_names(p$formula, modality), env = globalenv()),
        resampling_approximation = p$resampling_approximation, control_group = p$control_group, resampling_mechanism = p$resampling_mechanism,
        multiple_testing_method = p$multiple_testing_method, multiple_testing_alpha = p$multiple_testing_alpha
    )
    a <- records$assignment$settings
    hyper <- lapply(a$hyperparameters, function(x) if (is.list(x)) unlist(x) else x)
    if (!is.null(hyper$formula_object)) hyper$formula_object <- stats::as.formula(to_package_names(hyper$formula_object, modality), env = globalenv())
    object <- do.call(assign_grnas, c(list(object, method = a$method, parallel = FALSE), hyper))
    q <- records$qc$settings
    object <- run_qc(object, n_nonzero_trt_thresh = q$n_nonzero_trt_thresh, n_nonzero_cntrl_thresh = q$n_nonzero_cntrl_thresh,
        response_n_umis_range = unlist(q$response_n_umis_range), response_n_nonzero_range = unlist(q$response_n_nonzero_range), p_mito_threshold = q$p_mito_threshold)
    if (is.null(records$calibration)) return(object)
    k <- records$calibration$settings
    object |>
        run_calibration_check(parallel = FALSE, print_progress = FALSE, n_calibration_pairs = k$n_calibration_pairs, calibration_group_size = k$calibration_group_size) |>
        run_power_check(parallel = FALSE, print_progress = FALSE) |>
        run_discovery_analysis(parallel = FALSE, print_progress = FALSE)
}
for (dataset in names(objects)) for (name in names(objects[[dataset]])) {
    label <- paste0(dataset, "/", name)
    object <- objects[[dataset]][[name]]
    directory <- file.path(out, dataset, "analyses", name)
    problems <- validate_analysis(directory, screen_of(dataset))
    check(paste0(label, ": the folder written from sceptre's object conforms"), nrow(problems) == 0L, paste(problems$rule, collapse = ", "))
    for (part in c("analyses", "expected")) {
        verdict <- same_files(file.path(out, dataset, part, name), file.path("example", dataset, part, name))
        check(paste0(label, ": ", part, " equal to the committed example"), verdict, as.character(verdict))
    }
    back <- suppressWarnings(suppressMessages(read_main_object(screen_of(dataset), directory)))
    for (slot in setdiff(methods::slotNames(object), NOT_STORED)) {
        held <- get(object, slot)
        rebuilt <- get(back, slot)
        if (slot %in% UNORDERED) { held <- sort(held); rebuilt <- sort(rebuilt) }
        verdict <- compare_slot(held, rebuilt)
        check(paste0(label, ": slot `", slot, "` of the rebuilt object"), !startsWith(verdict, "DIFFERENT"), verdict)
    }
    again <- suppressWarnings(suppressMessages(replay(screen_of(dataset), directory)))
    check(paste0(label, ": cells in use, run again from the recorded settings"), identical(get(object, "cells_in_use"), get(again, "cells_in_use")))
    for (stage in intersect(c("calibration", "power", "discovery"), called(object))) {
        first <- as.data.frame(get(object, result_slot[[stage]]))
        second <- as.data.frame(get(again, result_slot[[stage]]))
        first <- first[names(second)]   # the rerun keeps sceptre's default columns
        check(paste0(label, ": ", stage, " results, run again from the recorded settings"), close_enough(rows_sorted(first), rows_sorted(second)))
    }
}
cat("sceptre suite: ", n_checks - failures, " of ", n_checks, " checks pass\n", sep = "")
if (failures > 0L) quit(status = 1L)
