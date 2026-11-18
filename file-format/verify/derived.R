# Derived values against frozen expectations, R side. See verify/derived.py.
# Loads neither sceptre nor anything built on it.
#     Rscript verify/derived.R
source("format/setup.R")
source("verify/as_sce.R")
source("verify/common.R")
ids_sha <- function(lists, cells) {
    text <- paste(paste0(names(lists), "\t", vapply(lists, function(idx) paste(cells[idx], collapse = ","), character(1))), collapse = "\n")
    unname(tools::sha256sum(bytes = charToRaw(enc2utf8(text))))
}
for (expected_file in sort(Sys.glob("example/*/expected/*/expected.json"))) {
    parts <- strsplit(expected_file, "/", fixed = TRUE)[[1]]
    dataset <- parts[2]; name <- parts[4]
    label <- paste0(dataset, "/", name)
    expected <- read_json(expected_file)
    sce <- read_as_sce(file.path("example", dataset, "screen.h5mu"))
    rhdf5::h5closeAll()
    saved <- fmt_read_analysis(file.path("example", dataset, "analyses", name))
    tables <- saved$tables; records <- saved$stages
    cells <- colnames(sce)
    grna <- altExp(sce, "gRNA")
    strategy <- records$parameters$settings$grna_integration_strategy
    grna_table <- data.frame(grna_id = rownames(grna), as.data.frame(rowData(grna)), stringsAsFactors = FALSE)
    grna_table$grna_group <- grna_groups(grna_table, strategy)
    if (!is.null(expected$cells_by_grna_sha256)) {
        by_grna <- cells_by_grna(tables[["assignment/assignments"]], rownames(grna), cells)
        pooled <- group_cells(by_grna, grna_table)
        check(paste0(label, ": cells of each gRNA"), identical(ids_sha(by_grna, cells), expected$cells_by_grna_sha256))
        check(paste0(label, ": cells of each target (pooled over its gRNAs)"), identical(ids_sha(pooled$grna_group_idxs, cells), expected$cells_by_group_sha256))
        check(paste0(label, ": cells of each non-targeting gRNA"), identical(ids_sha(pooled$indiv_nt_grna_idxs, cells), expected$cells_by_control_grna_sha256))
        if (!is.null(expected$grnas_per_cell_sha256)) {
            per_cell <- paste(grnas_per_cell(by_grna, length(cells)), collapse = ",")
            check(paste0(label, ": gRNAs per cell"), identical(unname(tools::sha256sum(bytes = charToRaw(per_cell))), expected$grnas_per_cell_sha256))
        }
    }
    if (!is.null(expected$cells_in_use_sha256)) {
        kept <- tables[["qc/cells"]]
        kept <- kept$in_use[match(cells, kept$cell_id)]
        check(paste0(label, ": cells in use"), identical(sha256_ids(cells[kept]), expected$cells_in_use_sha256))
        for (type in c("discovery", "positive_control")) {
            want <- expected[[paste0("n_ok_", type, "_pairs")]]
            if (is.numeric(want) && length(want) == 1L) {
                check(paste0(label, ": ", type, " pairs passing QC"), n_pairs_passing(tables[[paste0("qc/", type, "_pair_qc")]], strategy) == want)
            }
        }
    }
    if (!is.null(expected$control_groups)) {
        groups <- split(tables[["calibration/negative_control_groups"]]$grna_id, tables[["calibration/negative_control_groups"]]$grna_group)
        check(paste0(label, ": the control groups sceptre formed"), identical(sort(names(groups), method = "radix"), as.character(unlist(expected$control_groups))))
        non_targeting <- grna_table$grna_id[grna_table$grna_target == "non-targeting"]
        check(paste0(label, ": every control group's members are non-targeting gRNAs, and are the ones its name lists"),
            all(mapply(function(group, members) all(members %in% non_targeting) && setequal(members, strsplit(group, "&", fixed = TRUE)[[1]]), names(groups), groups)))
    }
    for (stage in c("calibration", "power", "discovery")) {
        path <- file.path(dirname(expected_file), paste0(stage, ".parquet"))
        if (!file.exists(path)) next
        want <- read_table(path)
        got <- add_log_2_fold_change(tables[[paste0(stage, "/results")]])
        if (stage != "power") got <- add_significant(got, records[[stage]]$significance_threshold)
        key <- unit_key(want)
        both <- merge(want, got, by = key, suffixes = c("_expected", ""))
        check(paste0(label, ": ", stage, " results cover the same pairs"), nrow(both) == nrow(want) && nrow(want) == nrow(got))
        for (column in c("p_value", "log_2_fold_change")) {
            check(paste0(label, ": ", stage, " ", column, " identical to sceptre's, to the last bit"), identical(both[[column]], both[[paste0(column, "_expected")]]))
        }
        if (!is.null(both$significant_expected)) {
            check(paste0(label, ": ", stage, " `significant`, from the stored threshold, identical to sceptre's"), identical(both$significant, both$significant_expected))
        }
    }
}
save_results("derived")
